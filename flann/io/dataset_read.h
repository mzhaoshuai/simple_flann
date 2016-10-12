/*
**@file  dataset_read.h		读取数据集的源文件
**
**@author zhaoshuai 
*/


#ifndef DATASET_READ_H_
#define DATASET_READ_H_

#include<fstream>
#include "../util/matrix.h"

namespace flann
{
	/*
	**@brief	将filename中从start到end的数据读取到到dataset中
	**			start从1开始
	**			Read a set of vectors stored in the fvec format (int + n * float)
	**			The function returns a set of output vector (one vector per row)暂定
	**
	**@tips		参考链接 http://corpus-texmex.irisa.fr/
	*/

	/*
	template<typename T>
	void fvecs_read(std::vector<T>& data,flann::Matrix<T>& dataset, const std::string& filename, int start, int end)
	{
		size_t a = 0, b = 0, bmax = 0;
		std::ifstream fid;
		fid.open(filename, std::ios::in|std::ios::binary);				//二进制，只读
		//fid.open(filename, std::ios::in);
		if (!fid)
		{
			std::cerr << "failed to open fvecs file" << std::endl;		//打开失败
			return;
		}
		//
		int d = 0;
		fid.seekg(0, std::ios::beg);
		char *point = reinterpret_cast<char*>(&d);
		fid.read(point, 4);												//读取向量维数
		//fid >> *point;
		int vecsizeof = 1 * 4 + d * 4;									//取得单个向量的大小，单位Bytes
		std::cout << "数据集向量维数为" << d << std::endl;
		//
		fid.seekg(0, std::ios::end);
		a = 1;
		bmax = fid.tellg() / vecsizeof;
		b = bmax;														//文件的总向量个数
		std::cout << "数据集向量个数为" << bmax << std::endl;
		//输入参数检查
		if (start < 1)
		{
			std::cerr << "failed to open fvecs file,parameter start wrong " << std::endl;
			return;
		}
		if ( end > b )
		{
			end = bmax;
		}
		if (end<1 || start>end)
		{
			std::cerr << "failed to open fvecs file,parameter end wrong " << std::endl;
			return;
		}
		//分配空间并传入数据
		int i, j;
		int n = end - start + 1;
		fid.seekg((start - 1)*vecsizeof, std::ios::beg);
		T temp;
		point = reinterpret_cast<char*>(&temp);
		//std::vector<T>one_vector(d);
		//char* p_vector = reinterpret_cast<char*>(&one_vector[0]);
		//
		//dataset = flann::Matrix<T>(new T[n*d], n, d);
		for (i = start - 1; i < end; i++)
		{
			for (j = 1; j < d + 1; j++)
			{
				fid.seekg(i*vecsizeof + 4 * j, std::ios::beg);
				//fid.read(dataset + (i - start + 1)*vecsizeof + 4 * (j - 1), 4);
				//point = dataset.ptr();
				fid.read(point, 4);
				data.push_back(temp);
				//dataset[(i - start + 1)*d + j - 1] = temp;
			}
			//fid.seekg(i*vecsizeof + 4 , std::ios::beg);
			//fid.read(p_vector, 4*d);
			
		}
		dataset =flann::Matrix<T>(&data[0], n, d);
		//关闭文件并退出
		fid.close();
		return;
	}
	*/


	/*
	**@brief	将filename中从start到end的数据读取到到dataset中
	**			start从1开始
	**			Read a set of vectors stored in the fvec format (int + n * float)
	**			The function returns a set of output vector (one vector per row)暂定
	**
	**@tips		参考链接 http://corpus-texmex.irisa.fr/
	*/
	template<typename T>
	void fvecs_ivecs_read(Matrix<T>& dataset, const std::string& filename, int start, int end)
	{
		size_t a = 0, b = 0, bmax = 0;
		std::ifstream fid;
		fid.open(filename, std::ios::in | std::ios::binary);			//二进制，只读
		if (!fid)
		{
			std::cerr << "failed to open fvecs file" << std::endl;		//打开失败
			return;
		}
		
		int d = 0;
		fid.seekg(0, std::ios::beg);
		char *point = reinterpret_cast<char*>(&d);
		fid.read(point, 4);												//读取向量维数
		int vecsizeof = 1 * 4 + d * 4;									//取得单个向量的大小，单位Bytes
		std::cout << "数据集向量维数为" << d << std::endl;
		

		fid.seekg(0, std::ios::end);
		a = 1;
		bmax = fid.tellg() / vecsizeof;
		b = bmax;														//文件的总向量个数
		std::cout << "数据集向量个数为" << bmax << std::endl;
		

		if (start < 1)													//输入参数检查
		{
			std::cerr << "failed to open fvecs file,parameter start wrong " << std::endl;
			return;
		}
		if (end > b)
		{
			end = bmax;
		}
		if (end<1 || start>end)
		{
			std::cerr << "failed to open fvecs file,parameter end wrong " << std::endl;
			return;
		}
		
		int i, j;														//分配空间并传入数据
		int n = end - start + 1;
		fid.seekg((start - 1)*vecsizeof, std::ios::beg);
		dataset = flann::Matrix<T>(new T[n*d], n, d);
		char* pointer = reinterpret_cast<char*>(dataset[0]);
		int stride = 4 * d;
		for (i = start - 1; i < end - 4; )								//一次复制5个向量
		{
			fid.seekg(i*vecsizeof + 4 , std::ios::beg);
			fid.read(pointer, stride);
			pointer += stride;
			fid.seekg((i + 1)*vecsizeof + 4, std::ios::beg);
			fid.read(pointer, stride);
			pointer += stride;
			fid.seekg((i + 2)*vecsizeof + 4, std::ios::beg);
			fid.read(pointer, stride);
			pointer += stride;
			fid.seekg((i + 3)*vecsizeof + 4, std::ios::beg);
			fid.read(pointer, stride);
			pointer += stride;
			fid.seekg((i + 4)*vecsizeof + 4, std::ios::beg);
			fid.read(pointer, stride);
			pointer += stride;
			i += 5;
		}
		while (i < end)
		{
			fid.seekg(i*vecsizeof + 4, std::ios::beg);
			fid.read(pointer, stride);
			pointer += stride;
			i++;
		}
		
		fid.close();													//关闭文件并退出
		point = NULL;
		pointer = NULL;
		return;
	}
}



#endif	//DATASET_READ_H_

