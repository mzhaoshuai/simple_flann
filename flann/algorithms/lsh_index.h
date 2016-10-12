/***********************************************************************
 * Software License Agreement (BSD License)
 /***********************************************************************
 *	Author: Vincent Rabaud
 *	redefine zhaoshuai
 *************************************************************************/

#ifndef FLANN_LSH_INDEX_H_
#define FLANN_LSH_INDEX_H_

#include <algorithm>
#include <cassert>
#include <cstring>
#include <map>
#include <vector>


#include "../general.h"
//#include "../algorithms/nn_index.h"
#include "../util/matrix.h"
#include "../util/result_set.h"
#include "../util/heap.h"
#include "../util/lsh_table.h"
#include "../util/allocator.h"
#include "../util/random.h"
#include "../util/saving.h"


#ifndef min
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#endif

namespace flann
{

	struct LshIndexParams //: public IndexParams
	{
	public:
		/*default constructor*/
		LshIndexParams(unsigned int table_number_, unsigned int key_size_, float gap_w_, unsigned int table_size_)
			:table_number(table_number_), key_size(key_size_), gap_w(gap_w_), table_size(table_size_) {}
		unsigned int table_number;
		unsigned int key_size;
		float gap_w;
		unsigned int table_size;
	};

	/**
	 * Randomized kd-tree index
	 *
	 * Contains the k-d trees and other information for indexing a set of points
	 * for nearest-neighbor matching.
	 */
	template<typename Distance>
	class LshIndex //: public NNIndex<Distance>
	{
	public:
		typedef typename Distance::ElementType ElementType;
		typedef typename Distance::ResultType DistanceType;
		//typedef float ElementType;
		//typedef float DistanceType;
	private:
		/** The different hash tables */
		std::vector<lsh::LshTable<ElementType> > tables_;
		/*特征向量数据，points_中按顺序存储每个特征向量的首地址，大小应为size_*/
		/*debug:10.06.15.54 points_°值正确*/
		std::vector<ElementType*> points_;
		/*数据集中特征向量的个数*/
		size_t size_;
		/*数据集中每个特征向量的维数*/
		size_t veclen_;
		/*hash表的个数*/
		unsigned int table_number_;
		/** key size */
		unsigned int key_size_;
		/*投影时的间隔*/
		float m_gap_w;
		/*哈希表的大小，一般取特征向量总数*/
		unsigned int m_table_size;
		/*距离函子*/
		Distance m_distance;
		/*lsh索引参数*/
		LshIndexParams index_params_;
		/*上一次索引建立时数据集中特征向量的个数，用以重建数据集*/
		size_t size_at_build_;
		/*距离计算次数*/
		size_t m_distance_cnt;
		/*距离因子*/
		//Distance distance_;
		/*寻求的近邻个数*/
		//size_t m_knn;
		/** How far should we look for neighbors in multi-probe LSH */
		//unsigned int multi_probe_level_;
		/** The XOR masks to apply to a key to get the neighboring buckets */
		//std::vector<lsh::BucketKey> xor_masks_;
		//typedef NNIndex<Distance> BaseClass;

	public:

		/** Constructor
		 * @param params parameters passed to the LSH algorithm
		 * @param d the distance used
		 */
		LshIndex(const LshIndexParams& params = LshIndexParams(), Distance d = Distance())
			:m_distance(d), size_(0), size_at_build_(0), veclen_(0), index_params_(params)
		{
			table_number_ = params.table_number;
			key_size_ = params.key_size;
			m_gap_w = params.gap_w;
			m_table_size = params.table_size;
			m_distance_cnt = 0;
			std::cout << "LshIndex第一种构造函数" << std::endl;
		}

		/*
		** Constructor
		** @param input_data dataset with the input features
		** @param params parameters passed to the LSH algorithm
		** @param d the distance used
		*/
		LshIndex(Matrix<ElementType>& input_data, const LshIndexParams& params = LshIndexParams(), Distance d = Distance())
			: m_distance(d), size_(0), size_at_build_(0), veclen_(0), index_params_(params)
		{
			table_number_ = params.table_number;
			key_size_ = params.key_size;
			m_gap_w = params.gap_w;
			setDataset(input_data);
			m_distance_cnt = 0;
			m_table_size = size_;
			//std::cout << "LshIndex第二种构造函数"<< std::endl;
		}

		/*
		**@brief	Constructor,由已有的Index初始化
		*/
		LshIndex(const LshIndex& other) :
			m_distance(other.m_distance),
			m_gap_w(other.m_gap_w),
			size_(other.size_),
			size_at_build_(other.size_at_build_),
			table_number_(other.table_number_),
			key_size_(other.key_size_),
			veclen_(other.veclen_),
			index_params_(other.index_params_),
			//points_(other.points_),
			m_table_size(other.m_table_size)
		{
		}

		LshIndex& operator=(LshIndex other)
		{
			this->swap(other);
			return *this;
		}

		 ~LshIndex() {	freeIndex();	}


		LshIndex* clone() const
		{
			return new LshIndex(*this);
		}

		/*释放掉所有内存*/
		void freeIndex()
		{
			tables_.clear();
		
		}

		/*
		**@brief 返回距离计算次数
		*/
		size_t get_m_distance_cnt()
		{
			return(m_distance_cnt);
		}

		/*
		**@brief	Builds the index
		*/
		void buildIndex( )
		{
			//freeIndex();
			//cleanRemovedPoints();
			// building index
			buildIndexImpl( );
			size_at_build_ = size_;
		}

		/*
		**@brief	用于重建索引
		*/
		void buildIndex(Matrix<ElementType>& input_data)
		{
			setDataset(input_data);
			m_distance_cnt = 0;
			m_table_size = size_;
			buildIndexImpl();
			size_at_build_ = size_;
			m_distance_cnt = 0;
		}


	protected:

		/*
		**@brief	Builds the index
		*/
		void buildIndexImpl( )
		{
			tables_.clear();
			tables_.resize(table_number_);									//根据哈希表的数目分配内存
			//std::vector<std::pair<size_t, ElementType*> > features;
			//features.reserve(points_.size());
			//features.resize(size_);
			//features.reserve(size_);										//size_	-特征向量个数
			unsigned int i = 0;
			/**
			while (i < points_.size() - 9)									//依次为所有特征向量编号
			{
				features.push_back(std::make_pair(i, points_[i]));
				features.push_back(std::make_pair(i + 1, points_[i + 1]));
				features.push_back(std::make_pair(i + 2, points_[i + 2]));
				features.push_back(std::make_pair(i + 3, points_[i + 3]));
				features.push_back(std::make_pair(i + 4, points_[i + 4]));
				features.push_back(std::make_pair(i + 5, points_[i + 5]));
				features.push_back(std::make_pair(i + 6, points_[i + 6]));
				features.push_back(std::make_pair(i + 7, points_[i + 7]));
				features.push_back(std::make_pair(i + 8, points_[i + 8]));
				features.push_back(std::make_pair(i + 9, points_[i + 9]));
				i += 5;
			}
			while (i < points_.size())
			{
				features.push_back(std::make_pair(i, points_[i]));
				i++;
			}
			**/
			for (i = 0; i < table_number_; ++i)								//根据hash表的数目，为每个表压入哈希后的关键字以及特征向量编号
			{
				//lsh::LshTable<ElementType>& table = tables_[i];
				//table = lsh::LshTable<ElementType>(veclen_, key_size_, m_gap_w, m_table_size);
				tables_[i]=lsh::LshTable<ElementType>(veclen_, key_size_, m_gap_w, m_table_size);
				// Add the features to the table
				//tables_[i].add(features);
				tables_[i].add(points_);
			}
		}

		

		/*
		**@brief	由输入的数据dataset修改索引的信息，从nn_index复制
		**
		**@param	dataset	-输入的数据集
		*/
		void setDataset(Matrix<ElementType>& dataset)
		{
			size_ = dataset.rows;
			veclen_ = dataset.cols;
			//last_id_ = 0;

			//ids_.clear();
			//removed_points_.clear();
			//removed_ = false;
			//removed_count_ = 0;
			points_.clear();
			points_.resize(size_);
			size_t i = 0;
			while( i < size_ - 4 )											//一次算5个
			{
				points_[i] = dataset[i];									//dataset[]重载后返回行向量的首地址
				points_[i + 1] = dataset[i + 1];
				points_[i + 2] = dataset[i + 2];
				points_[i + 3] = dataset[i + 3];
				points_[i + 4] = dataset[i + 4];
				i += 5;
			}
			while (i < size_)
			{
				points_[i] = dataset[i];
				i++;
			}
			/*nothing to do here*/
		}

		/*
		**@brief	扩展数据集
		*/
		void extendDataset(const Matrix<ElementType>& new_points)
		{
			//size_t new_size = size_ + new_points.rows;
			//if (removed_) 
			//{
			//	removed_points_.resize(new_size);
			//	ids_.resize(new_size);
			//}
			//points_.resize(new_size);
			//for (size_t i = size_; i < new_size; ++i)
			//{
			//	points_[i] = new_points[i - size_];
				//if (removed_) 
				//{
				//	ids_[i] = last_id_++;
				//	removed_points_.reset(i);
				//}
			//}
			//size_ = new_size;
		}

	public:

		/*
		**@brief 添加数据集
		*/
		void addPoints(const Matrix<ElementType>& points, float rebuild_threshold = 2)
		{
			/**
			assert(points.cols == veclen_);
			size_t old_size = size_;

			extendDataset(points);											//其中重定义points_的空间大小

			if (rebuild_threshold > 1 && size_at_build_*rebuild_threshold < size_)
			{
				buildIndex();
			}
			else
			{
				for (unsigned int i = 0; i < table_number_; ++i)
				{
					lsh::LshTable<ElementType>& table = tables_[i];
					for (size_t i = old_size; i < size_; ++i)
					{
						table.add(i, points_[i]);
					}
				}
			}
			**/
		}


		//flann_algorithm_t getType() const
		//{
		//	return FLANN_INDEX_LSH;
		//}


		/**
		template<typename Archive>
		void serialize(Archive& ar)
		{
		ar.setObject(this);

		ar & *static_cast<NNIndex<Distance>*>(this);

		ar & table_number_;
		ar & key_size_;
		ar & multi_probe_level_;

		ar & xor_masks_;
		ar & tables_;

		if (Archive::is_loading::value) {
		index_params_["algorithm"] = getType();
		index_params_["table_number"] = table_number_;
		index_params_["key_size"] = key_size_;
		index_params_["multi_probe_level"] = multi_probe_level_;
		}
		}
		**/

		void saveIndex(FILE* stream)
		{
			//	serialization::SaveArchive sa(stream);
			//	sa & *this;
		}

		void loadIndex(FILE* stream)
		{
			//	serialization::LoadArchive la(stream);
			//	la & *this;
		}

		/**
		 * Computes the index memory usage
		 * Returns: memory used by the index
		 */
		int usedMemory() const
		{
			//return size_ * sizeof(int);
			int size_kb = size_*table_number_*sizeof(ElementType);
			size_kb = size_kb >> 10;
			return(size_kb);
		}

		/**
		 * \brief Perform k-nearest neighbor search
		 * \param[in] queries The query points for which to find the nearest neighbors
		 * \param[out] indices The indices of the nearest neighbors found
		 * \param[out] dists Distances to the nearest neighbors found
		 * \param[in] knn Number of nearest neighbors to return
		 * \param[in] params Search parameters
		 */
		int knnSearch(
			const Matrix<ElementType>& queries,
			Matrix<size_t>& indices,
			Matrix<DistanceType>& dists,
			size_t knn)//const
			//const SearchParams& params) const
		{
			assert(queries.cols == veclen_);
			assert(indices.rows >= queries.rows);
			assert(dists.rows >= queries.rows);
			assert(indices.cols >= knn);
			assert(dists.cols >= knn);
			size_t dis_cnt0 = 0, dis_cnt1 = 0, dis_cnt2 = 0, dis_cnt3 = 0, dis_cnt4 = 0;
			//m_knn = knn;
			//size_t n = 0;
			//int count = 0;
			//if (params.use_heap == FLANN_True)
			//{
//#pragma omp parallel num_threads(params.cores)
			//	{
					//KNNUniqueResultSet<DistanceType> resultSet(knn);
					int i = 0;
//#pragma omp for schedule(static) reduction(+:count)
					while (i < (int)(queries.rows - 4))
					{
						dis_cnt0 = getNeighbors(indices[i], dists[i], queries[i], knn);

						dis_cnt1 = getNeighbors(indices[i + 1], dists[i + 1], queries[i + 1], knn);

						dis_cnt2 = getNeighbors(indices[i + 2], dists[i + 2], queries[i + 2], knn);
						
						dis_cnt3 = getNeighbors(indices[i + 3], dists[i + 3], queries[i + 3], knn);

						dis_cnt4 = getNeighbors(indices[i + 4], dists[i + 4], queries[i + 4], knn);

						i += 5;
						m_distance_cnt += (dis_cnt0 + dis_cnt1 + dis_cnt2 + dis_cnt3 + dis_cnt4);
					}
					while (i < (int)queries.rows)
					{
						dis_cnt0 = getNeighbors(indices[i], dists[i], queries[i], knn);
						m_distance_cnt += dis_cnt0;
						i++;
					}
			//	}
			//}
			//else
			//{
//#pragma omp parallel num_threads(params.cores)
			//	{
			//		KNNResultSet<DistanceType> resultSet(knn);
//#pragma omp for schedule(static) reduction(+:count)
			//		for (int i = 0; i < (int)queries.rows; i++)
			//		{
			//			resultSet.clear();
			//			//findNeighbors(resultSet, queries[i], params);
			//			getNeighbors(queries[i], resultSet);
			//			n = (std::min)(resultSet.size(), knn);
			//			resultSet.copy(indices[i], dists[i], n, params.sorted);
			//			//indices_to_ids(indices[i], indices[i], n);
			//			count += n;
			//		}
			//	}
			//}

			return 0;
		}

	private:

		/*
		**	Performs the approximate nearest-neighbor search.
		**	This is a slower version than the above as it uses the ResultSet
		**	@param vec -the feature to analyze
		*/
		size_t getNeighbors(
			size_t* indices,
			DistanceType*  dists,
			ElementType* vec, 
			size_t knn) const
		{
			typename std::vector<lsh::LshTable<ElementType> >::const_iterator table = tables_.begin();
			typename std::vector<lsh::LshTable<ElementType> >::const_iterator table_end = tables_.end();	//迭代器，表头表尾
			//unsigned int worst_score = std::numeric_limits<unsigned int>::max();
			size_t distance_cnt = 0;
			lsh::BucketKey bucket_key;
			std::vector<size_t>indices_one;
			std::vector<DistanceType>dists_one;
			//size_t* indices_one = new size_t[knn];
			//DistanceType* dists_one = new DistanceType[knn];
			size_t i = 0;
			unsigned char first_flag = 0;
			for (; table != table_end; table++)																//每个表都要进行计算并存储数据
			{
				i = 0;
				bucket_key = std::make_pair(0, 0);
				bucket_key = table->getKey(vec);

				const lsh::Bucket* bucket_p = table->getBucketFromKey(bucket_key);
				if (bucket_p == 0)
					continue;
				// Go over each descriptor index
				std::vector<lsh::FeatureIndex>::const_iterator training_index = bucket_p->begin();				//
				std::vector<lsh::FeatureIndex>::const_iterator last_training_index = bucket_p->end();
				DistanceType euclidean_distance = 0;
				// Process the rest of the candidates
				for (; training_index < last_training_index; training_index++,i++)								//算出查询点到所有点额距离并存储
				{
					euclidean_distance = 0;
					if (*training_index >= size_)
						continue;
					// Compute the euclidean distance
					euclidean_distance = m_distance(vec, points_[*training_index], veclen_);
					distance_cnt++;
					dists_one.push_back(euclidean_distance);
					indices_one.push_back(*training_index);
					//if (i < knn)
					//{
					//	*(dists_one + i) = euclidean_distance;
					//	*(indices_one + i) = *training_index;														//将距离压入结果数据集,所有表的数据压到一起
					//}
					//else
					//{
					//	if (i == knn)
					//		bubble_swap(indices_one, dists_one, knn);
					//	if ( euclidean_distance> *(indices_one+knn-1) )
					//		continue;
					//	else
					//	{
					//		*(dists_one + knn - 1) = euclidean_distance;
					//		*(indices_one + knn - 1) = *training_index;
					//		bubble_swap(indices_one, dists_one, knn);
					//	}
					//}
				}
				
			}
			for (i = dists_one.size(); i < knn; i++)
			{
				dists_one.push_back(E2LSH_BigPrimeNum);
				indices_one.push_back(E2LSH_BigPrimeNum);
			}
			bubble_swap(indices_one, dists_one, knn);
			for (i = 0; i < knn; i++)
			{
				*(indices + i) = indices_one[i];
				*(dists + i) = dists_one[i];
			}
			dists_one.clear();
			indices_one.clear();
			return(distance_cnt);
		}


		/*
		**	@brief	对dist_one中的数据进行冒泡排序，将较小数冒到前面,同时交换indices_one中的数据
		**	
		**	@param	
		*/
		void bubble_swap(std::vector<size_t>&indices_one, std::vector<DistanceType>&dists_one, size_t knn)const
		{
			size_t i = 0, j = 0, i_temp = 0;
			float d_temp = 0;
			int total_num = 0;
			total_num = indices_one.size();
			for (i = total_num - 1; i > total_num - 1 - knn; i--)
			{
				for (j = i; j>0; j--)
				{
					if (dists_one[j] < dists_one[j-1])
					{
						/*距离交换*/
						d_temp = dists_one[j];
						dists_one[j] = dists_one[j-1];
						dists_one[j-1] = d_temp;
						/*索引交换*/
						i_temp = indices_one[j];
						indices_one[j] = indices_one[j - 1];
						indices_one[j - 1] = i_temp;
					}
				}
			}
		}

		void swap(LshIndex& other)
		{
			//BaseClass::swap(other);
			std::swap(tables_, other.tables_);
			std::swap(points_, other.points_);
			std::swap(size_, other.size_);
			std::swap(veclen_, other.veclen_);
			std::swap(table_number_, other.table_number_);
			std::swap(key_size_, other.key_size_);
			std::swap(m_gap_w, other.m_gap_w);
			std::swap(m_table_size, other.m_table_size);
			std::swap(m_distance, other.m_distance);
			std::swap(index_params_, other.index_params_);
			std::swap(size_at_build_, other.size_at_build_);
			//std::swap(xor_masks_, other.xor_masks_);
		}
		//USING_BASECLASS_SYMBOLS
	};
}

#endif //FLANN_LSH_INDEX_H_
