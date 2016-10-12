/*
**@brief	LSH和KDTrees的主要测试文件
**
**@author	zhaoshuai
*/



#include "flann/flann.hpp"
//#include "flann/io/hdf5.h"											//去除hdf5文件的读写库
#include "flann/io/save_result_txt.h"
#include "flann/io/dataset_read.h"
#include "flann/algorithms/lsh_index.h"
#include "flann/algorithms/mean_precision.h"

#include <stdio.h>
#include <cstdlib>
#include <iostream>
#include <windows.h>

/*为0选择kdtree方式，为1选择lsh方式*/
#define LSH_INDEX_ENABLE (1)		
/*为0选择SIFT数据集，为1选择GIST数据集*/
//#define SIFT_GIST_CHOOSE	(0)

/*SIFT数据集的定义*/
#define ANN_SIFT1M_BASE_NUM (1000000)
#define ANN_SIFT1M_QUERY_NUM (1000)
#define ANN_SIFT1M_GROUND_TRUTH_NUM (1000)

/*GIST数据集的定义*/
#define ANN_GIST1M_BASE_NUM_1 (350000)									//内存有限，分成三段读取
#define ANN_GIST1M_BASE_NUM_2 (700000)
#define ANN_GIST1M_BASE_NUM_3 (1000000)
#define ANN_GIST1M_QUERY_NUM (1000)
#define ANN_GIST1M_GROUND_TRUTH_NUM (1000)

using namespace flann;




/*
**@brief	kd-tree在SIFT1M数据集下的测试函数
**
**@params	knn			-搜索的近邻个数
**			tree_num	-要搜索时要使用的树的个数
*/
void test_kdtree_sift(
	size_t knn,
	int tree_num)
{

	Matrix<float> dataset;
	Matrix<float> query;
	Matrix<size_t> ground_truth;

	DWORD start_time, end_time, spend_time;
	float mean_precision = 0;

	/*读取数据信息*/
	flann::fvecs_ivecs_read(dataset, "dataset/sift/sift_base.fvecs", 1, ANN_SIFT1M_BASE_NUM);
	flann::fvecs_ivecs_read(query, "dataset/sift/sift_query.fvecs", 1, ANN_SIFT1M_QUERY_NUM);
	flann::fvecs_ivecs_read(ground_truth, "dataset/sift/sift_groundtruth.ivecs", 1, ANN_SIFT1M_GROUND_TRUTH_NUM);


	Matrix<size_t> indices(new size_t[query.rows*knn], query.rows, knn);
	Matrix<float> dists(new float[query.rows*knn], query.rows, knn);


	start_time = GetTickCount();
	std::cout << std::endl;
	std::cout << "KD-TREE方式检索,采用" << tree_num << "棵树" << std::endl;
	Index<L2<float>> index(dataset, flann::KDTreeIndexParams(tree_num));					//FLANN_INDEX_KDTREE
	index.buildIndex();
	end_time = GetTickCount();
	std::cout << "建立索引完成" << std::endl;
	std::cout << "建立索引所用时间: " << end_time - start_time << "ms" << std::endl;
	// do a knn search, using 128 checks
	start_time = GetTickCount();
	index.knnSearch(query, indices, dists, knn, flann::SearchParams(128));
	end_time = GetTickCount();
	std::cout << "查询完成" << std::endl;
	std::cout << "查询所用时间: " << end_time - start_time << "ms" << std::endl;
	std::cout << "查询结果的平均精度：";
	mean_precision = flann::mean_precision(indices, ground_truth);
	std::cout << std::setprecision(4) << mean_precision * 100 << "%" << std::endl;
	std::cout << "内存开销:" << index.usedMemory()<<"Bytes"<<std::endl;
	std::cout << std::endl;
	/*删除new的空间*/
	delete[] indices.ptr();
	delete[] dists.ptr();
	delete[] query.ptr();
	delete[] ground_truth.ptr();
	delete[] dataset.ptr();
	return;
}


/*
**@brief	LSH在SIFT1M数据集下的测试函数
**
**@params	knn				-搜索的近邻个数
**			table_number	-用的哈希表个数
**			key_size		-关键字的个数。第一次哈希后得到的向量维数。
**			gap_w			-分割间隔
*/
void test_lsh_sift(
	size_t knn,
	unsigned int table_number,
	unsigned int key_size,
	float gap_w)
{

	Matrix<float> dataset;
	Matrix<float> query;
	Matrix<size_t> ground_truth;

	DWORD start_time, end_time, spend_time;
	float mean_precision = 0;
	/*读取数据信息*/
	std::cout << std::endl;
	flann::fvecs_ivecs_read(dataset, "dataset/sift/sift_base.fvecs", 1, ANN_SIFT1M_BASE_NUM);
	flann::fvecs_ivecs_read(query, "dataset/sift/sift_query.fvecs", 1, ANN_SIFT1M_QUERY_NUM);
	flann::fvecs_ivecs_read(ground_truth, "dataset/sift/sift_groundtruth.ivecs", 1, ANN_SIFT1M_GROUND_TRUTH_NUM);


	Matrix<size_t> indices(new size_t[query.rows*knn], query.rows, knn);
	Matrix<float> dists(new float[query.rows*knn], query.rows, knn);

	std::cout << std::endl;
	std::cout << "LSH方式检索" << std::endl;
	std::cout << "哈希表数:" << table_number << "\t";
	std::cout << "关键字个数:" << key_size << "\t";
	std::cout << std::setprecision(4) << "分割间隔:" << gap_w << "\t";

	start_time = GetTickCount();
	/*参数*/
	unsigned int table_size = dataset.rows;

	LshIndexParams lsh_param(table_number, key_size, gap_w, table_size);
	LshIndex<L2<float>> lsh_index(dataset, lsh_param, flann::L2<float>());
	lsh_index.buildIndex();

	end_time = GetTickCount();
	std::cout << "建立索引完成" << std::endl;
	std::cout << "建立索引所用时间: " << end_time - start_time << "ms" << std::endl;
	start_time = GetTickCount();
	lsh_index.knnSearch(query, indices, dists, knn);
	end_time = GetTickCount();
	std::cout << "查询完成" << std::endl;
	std::cout << "查询所用时间: " << end_time - start_time << "ms" << std::endl;
	std::cout << "查询结果的平均精度：";
	mean_precision = flann::mean_precision(indices, ground_truth);
	std::cout << std::setprecision(4) << mean_precision * 100 << "%" << std::endl;
	std::cout << "内存开销:" << lsh_index.usedMemory() << "KB" << std::endl;
	//std::cout << std::endl;
	std::cout << "距离计算次数:" << lsh_index.get_m_distance_cnt() << "次" << std::endl;
	std::cout << std::endl;

	delete[] indices.ptr();
	delete[] dists.ptr();
	delete[] query.ptr();
	delete[] ground_truth.ptr();
	delete[] dataset.ptr();
	return;

}

/*
**@brief 用于gist数据集最后的合并工作
*/
void gist_combine(
	Matrix<size_t>& indices,
	Matrix<float>& dists,
	Matrix<size_t>& indices_1,
	Matrix<float>& dists_1,
	Matrix<size_t>& indices_2,
	Matrix<float>& dists_2
	)
{
	int cols_1 = indices.cols;
	int cols_2 = 2 * indices.cols;
	int cols_3 = 3 * indices.cols;
	int i = 0, j = 0, k = 0;
	size_t* res_indice = new size_t[cols_3];						//索引备用空间
	float* res_dist = new float[cols_3];							//距离备用空间
	size_t i_temp = 0;
	float d_temp = 0;

	for (i = 0; i < indices.rows; i++)
	{
		for (j = 0; j < cols_1; j++)								//复制到一个空间中
		{
			*(res_indice + j) = *(indices[i] + j);
			*(res_dist + j) = *(dists[i] + j);
		}
		for (j = cols_1; j < cols_2; j++)
		{
			*(res_indice + j) = *(indices_1[i] + j - cols_1);
			*(res_dist + j) = *(dists_1[i] + j - cols_1);
		}
		for (j = cols_2; j < cols_3; j++)
		{
			*(res_indice + j) = *(indices_2[i] + j - cols_2);
			*(res_dist + j) = *(dists_2[i] + j - cols_2);
		}

		for (j = cols_3 - 1; j>cols_2 - 1; j--)						//冒泡，将较小的cols_1个数冒到前面来，同时交换标志
		{
			for (k = j; k > 0; k--)
			{
				if (*(res_dist + k) < *(res_dist + k - 1))
				{
					//交换距离
					d_temp = *(res_dist + k);
					*(res_dist + k) = *(res_dist + k - 1);
					*(res_dist + k - 1) = d_temp;
					//交换索引
					i_temp = *(res_indice + k);
					*(res_indice + k) = *(res_indice + k - 1);
					*(res_indice + k - 1) = i_temp;
				}
			}
		}
		for (j = 0; j < cols_1; j++)						//冒泡的的结果复制到原来存储空间中
		{
			*(indices[i] + j) = *(res_indice + j);
			*(dists[i] + j) = *(res_dist + j);
		}

	}
	delete[]res_indice;
	delete[]res_dist;
	return;
}


/*
**@brief	kd-tree在GIST1M数据集下的测试函数
**
**@params	knn			-搜索的近邻个数
**			tree_num	-要搜索时要使用的树的个数
*/
void test_kdtree_gist(size_t knn,int tree_num)
{
	Matrix<float> dataset;
	Matrix<float> dataset_1;
	Matrix<float> dataset_2;
	Matrix<float> query;
	Matrix<size_t> ground_truth;
	
	std::cout << std::endl;
	flann::fvecs_ivecs_read(query, "dataset/gist/gist_query.fvecs", 1, ANN_GIST1M_QUERY_NUM);
	flann::fvecs_ivecs_read(ground_truth, "dataset/gist/gist_groundtruth.ivecs", 1, ANN_GIST1M_GROUND_TRUTH_NUM);


	Matrix<size_t> indices(new size_t[query.rows*knn], query.rows, knn);
	Matrix<float> dists(new float[query.rows*knn], query.rows, knn);
	Matrix<size_t> indices_1(new size_t[query.rows*knn], query.rows, knn);
	Matrix<float> dists_1(new float[query.rows*knn], query.rows, knn);
	Matrix<size_t> indices_2(new size_t[query.rows*knn], query.rows, knn);
	Matrix<float> dists_2(new float[query.rows*knn], query.rows, knn);

	DWORD start_time, end_time, spend_time;
	DWORD buildtime, buildtime_1, buildtime_2;
	DWORD searchtime, searchtime_1, searchtime_2;
	int mem = 0, mem_1, mem_2;
	float mean_precision = 0;

	
	std::cout << std::endl;
	std::cout << "KD-TREE方式检索,采用" << tree_num << "棵树" << std::endl;

	
	/*部分搜索1*/
	flann::fvecs_ivecs_read(dataset, "dataset/gist/gist_base.fvecs", 1, ANN_GIST1M_BASE_NUM_1);
	start_time = GetTickCount();
	Index<L2<float>> index(dataset, flann::KDTreeIndexParams(tree_num));					//FLANN_INDEX_KDTREE
	index.buildIndex();
	end_time = GetTickCount();
	buildtime = end_time - start_time;
	start_time = GetTickCount();
	index.knnSearch(query, indices, dists, knn, flann::SearchParams(128));
	end_time = GetTickCount();
	searchtime = end_time - start_time;
	mem = index.usedMemory();
	delete[] dataset.ptr();
	/*部分搜索2*/
	flann::fvecs_ivecs_read(dataset_1, "dataset/gist/gist_base.fvecs", ANN_GIST1M_BASE_NUM_1, ANN_GIST1M_BASE_NUM_2);
	start_time = GetTickCount();
	//Index<L2<float>> index_1(dataset_1, flann::KDTreeIndexParams(tree_num));					//FLANN_INDEX_KDTREE
	index.buildIndex(dataset_1);
	end_time = GetTickCount();
	buildtime_1 = end_time - start_time;

	start_time = GetTickCount();
	index.knnSearch(query, indices_1, dists_1, knn, flann::SearchParams(128));
	end_time = GetTickCount();
	searchtime_1 = end_time - start_time;

	mem_1 = index.usedMemory(); 
	delete[] dataset_1.ptr();
	/*部分搜索3*/
	flann::fvecs_ivecs_read(dataset_2, "dataset/gist/gist_base.fvecs", ANN_GIST1M_BASE_NUM_2, ANN_GIST1M_BASE_NUM_3);
	start_time = GetTickCount();
	//Index<L2<float>> index_2(dataset_2, flann::KDTreeIndexParams(tree_num));					//FLANN_INDEX_KDTREE
	index.buildIndex(dataset_2);
	end_time = GetTickCount();
	buildtime_2 = end_time - start_time;

	start_time = GetTickCount();
	index.knnSearch(query, indices_2, dists_2, knn, flann::SearchParams(128));
	end_time = GetTickCount();
	searchtime_2 = end_time - start_time;

	mem_2 = index.usedMemory();
	delete[] dataset_2.ptr();

	std::cout << "建立索引完成" << std::endl;
	std::cout << "建立索引所用时间: " << (buildtime + buildtime_1 + buildtime_2) << "ms" << std::endl;	
	std::cout << "查询完成" << std::endl;
	std::cout << "查询所用时间: " << (searchtime + searchtime_1 + searchtime_2) << "ms" << std::endl;
	std::cout << "查询结果的平均精度：";
	gist_combine(indices, dists, indices_1, dists_1, indices_2, dists_2);
	mean_precision = flann::mean_precision(indices, ground_truth);
	std::cout << std::setprecision(4) << mean_precision * 100 << "%" << std::endl;
	std::cout << "内存开销:" << (mem + mem_1 + mem_2) << "Bytes" << std::endl;
	std::cout << std::endl;

	delete[] indices.ptr();
	delete[] indices_1.ptr();
	delete[] indices_2.ptr();
	delete[] dists.ptr();
	delete[] dists_1.ptr();
	delete[] dists_2.ptr();
	delete[] query.ptr();
	delete[] ground_truth.ptr();
	return;
}


/*
**@brief	LSH在GIST1M数据集下的测试函数
**
**@params	knn				-搜索的近邻个数
**			table_number	-用的哈希表个数
**			key_size		-关键字的个数。第一次哈希后得到的向量维数。
**			gap_w			-分割间隔
*/
void test_lsh_gist(
	size_t knn,
	unsigned int table_number,
	unsigned int key_size,
	float gap_w)
{
	
	Matrix<float> dataset;
	Matrix<float> dataset_1;
	Matrix<float> dataset_2;
	Matrix<float> query;
	Matrix<size_t> ground_truth;

	std::cout << std::endl;
	flann::fvecs_ivecs_read(query, "dataset/gist/gist_query.fvecs", 1, ANN_GIST1M_QUERY_NUM);
	flann::fvecs_ivecs_read(ground_truth, "dataset/gist/gist_groundtruth.ivecs", 1, ANN_GIST1M_GROUND_TRUTH_NUM);


	Matrix<size_t> indices(new size_t[query.rows*knn], query.rows, knn);
	Matrix<float> dists(new float[query.rows*knn], query.rows, knn);
	Matrix<size_t> indices_1(new size_t[query.rows*knn], query.rows, knn);
	Matrix<float> dists_1(new float[query.rows*knn], query.rows, knn);
	Matrix<size_t> indices_2(new size_t[query.rows*knn], query.rows, knn);
	Matrix<float> dists_2(new float[query.rows*knn], query.rows, knn);

	DWORD start_time, end_time, spend_time;
	DWORD buildtime, buildtime_1, buildtime_2;
	DWORD searchtime, searchtime_1, searchtime_2;
	int mem = 0, mem_1, mem_2;
	int dis_cnt = 0, dis_cnt_1 = 0, dis_cnt_2 = 0;
	float mean_precision = 0;


	std::cout << std::endl;
	std::cout << "LSH方式检索" << std::endl;
	std::cout << "哈希表数:" << table_number << "\t";
	std::cout << "关键字个数:" << key_size << "\t";
	std::cout << std::setprecision(4) << "分割间隔:" << gap_w << std::endl;


	/*部分搜索1*/
	flann::fvecs_ivecs_read(dataset, "dataset/gist/gist_base.fvecs", 1, ANN_GIST1M_BASE_NUM_1);

	start_time = GetTickCount();
	unsigned int table_size = dataset.rows;
	LshIndexParams lsh_param(table_number, key_size, gap_w, table_size);
	LshIndex<L2<float>> lsh_index(dataset, lsh_param, flann::L2<float>());
	lsh_index.buildIndex();
	end_time = GetTickCount();

	buildtime = end_time - start_time;
	start_time = GetTickCount();
	lsh_index.knnSearch(query, indices, dists, knn);
	end_time = GetTickCount();
	searchtime = end_time - start_time;
	mem = lsh_index.usedMemory();
	dis_cnt = lsh_index.get_m_distance_cnt();

	delete[] dataset.ptr();

	/*部分搜索2*/
	flann::fvecs_ivecs_read(dataset_1, "dataset/gist/gist_base.fvecs", ANN_GIST1M_BASE_NUM_1, ANN_GIST1M_BASE_NUM_2);
	start_time = GetTickCount();
	lsh_index.buildIndex(dataset_1);									//重建
	end_time = GetTickCount();
	buildtime_1 = end_time - start_time;

	start_time = GetTickCount();
	lsh_index.knnSearch(query, indices_1, dists_1, knn);
	end_time = GetTickCount();
	searchtime_1 = end_time - start_time;
	mem_1 = lsh_index.usedMemory();
	dis_cnt_1 = lsh_index.get_m_distance_cnt();

	delete[] dataset_1.ptr();

	/*部分搜索3*/
	flann::fvecs_ivecs_read(dataset_2, "dataset/gist/gist_base.fvecs", ANN_GIST1M_BASE_NUM_2, ANN_GIST1M_BASE_NUM_3);
	start_time = GetTickCount();
	lsh_index.buildIndex(dataset_2);									//第二次重建索引
	end_time = GetTickCount();
	buildtime_2 = end_time - start_time;

	start_time = GetTickCount();
	lsh_index.knnSearch(query, indices_2, dists_2, knn);				//重搜
	end_time = GetTickCount();
	searchtime_2 = end_time - start_time;
	mem_2 = lsh_index.usedMemory();
	dis_cnt_2 = lsh_index.get_m_distance_cnt();

	delete[] dataset_2.ptr();

	/*输出结果信息*/
	std::cout << "建立索引完成" << std::endl;
	std::cout << "建立索引所用时间: " << (buildtime + buildtime_1 + buildtime_2) << "ms" << std::endl;
	std::cout << "查询完成" << std::endl;
	std::cout << "查询所用时间: " << (searchtime + searchtime_1 + searchtime_2) << "ms" << std::endl;
	std::cout << "查询结果的平均精度：";
	gist_combine(indices, dists, indices_1, dists_1, indices_2, dists_2);
	mean_precision = flann::mean_precision(indices, ground_truth);
	std::cout << std::setprecision(4) << mean_precision * 100 << "%" << std::endl;
	std::cout << "内存开销:" << (mem + mem_1 + mem_2) << "Bytes" << std::endl;
	std::cout << "计算距离次数:" << (dis_cnt + dis_cnt_1 + dis_cnt_2) << "次" << std::endl;
	std::cout << std::endl;

	delete[] indices.ptr();
	delete[] indices_1.ptr();
	delete[] indices_2.ptr();
	delete[] dists.ptr();
	delete[] dists_1.ptr();
	delete[] dists_2.ptr();
	delete[] query.ptr();
	delete[] ground_truth.ptr();
	return;
	
}


//int main(int argc, char** argv)
int main( void )
{
    int nn = 3;
	//DWORD start_time, end_time, spend_time;
	//float mean_precision = 0;

	if (LSH_INDEX_ENABLE != 1)															//采用非LSH方法
	{

		/*测试组1-测试SIFT1M数据集下kdtree搜索的表现,变量 - 树的个数*/
		
		//test_kdtree_sift(nn, 2);
		//test_kdtree_sift(nn, 4);
		//test_kdtree_sift(nn, 8);
		//test_kdtree_sift(nn, 12);
		//test_kdtree_sift(nn, 16);
		
		
		/*测试组2-测试GIST1M数据集下kdtree搜索的表现，变量 - 树的个数*/
		
		/**
		test_kdtree_gist(nn, 2);
		test_kdtree_gist(nn, 4);
		test_kdtree_gist(nn, 8);
		test_kdtree_gist(nn, 12);
		test_kdtree_gist(nn, 16);
		**/
		
	}
	else
	{

		/*测试组1-测试SIFT1M数据集下LSH搜索的表现,变量 - 分割间隔*/
		test_lsh_sift(nn, 5, 4, 50);
		test_lsh_sift(nn, 5, 4, 100);
		test_lsh_sift(nn, 5, 4, 150);
		test_lsh_sift(nn, 5, 4, 200);
		//test_lsh_sift(nn, 5, 4, 250);
		//test_lsh_sift(nn, 5, 4, 300);
		//test_lsh_sift(nn, 5, 4, 350);
		
		/*测试组2-测试SIFT1M数据集下LSH搜索的表现,变量 - 哈希表数目*/
		/**
		test_lsh_sift(nn, 2, 4, 200);
		test_lsh_sift(nn, 4, 4, 200);
		test_lsh_sift(nn, 6, 4, 200);
		test_lsh_sift(nn, 8, 4, 200);
		test_lsh_sift(nn,10, 4, 200);
		**/
	
		/*测试组3-测试SIFT1M数据集下LSH搜索的表现,变量 - 关键字长度*/
		/**
		test_lsh_sift(nn, 4, 4, 210);
		test_lsh_sift(nn, 4, 8, 210);
		test_lsh_sift(nn, 4, 16, 210);
		test_lsh_sift(nn, 4, 32, 210);
		test_lsh_sift(nn, 4, 64, 210);
		**/

		/*测试组4-测试GIST1M数据集下LSH搜索的表现,变量 - 分割间隔*/

		//test_lsh_gist(nn, 3, 4, 1);
		//test_lsh_gist(nn, 3, 4, 1.5);
		//test_lsh_gist(nn, 3, 4, 2);
		//test_lsh_gist(nn, 3, 4, 3);
		//test_lsh_gist(nn, 3, 4, 4);



	}																			
	
	std::cout << "检索完毕" << std::endl;
	std::system("pause");
    return 0;
}

