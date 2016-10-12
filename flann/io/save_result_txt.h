/*
**@file  save_result_txt.h
**
**@author zhaoshuai 
*/

#ifndef SAVE_RESULT_TXT_H_
#define SAVE_RESULT_TXT_H_


//#include "save_result_txt.h"
#include "../util/matrix.h"
#include <fstream>

namespace flann
{
	/*
	**@brief 将dataset存储到filename所指文件中，txt格式以便查看
	*/
	template<typename T>
	void save_result_txt(flann::Matrix<T>& dataset, const std::string& filename)
	{
		//size_t rows = 0, columns = 0;
		size_t i = 0, j = 0;

		//rows = dataset.rows;
		//columns = dataset.cols;

		std::ofstream result;
		result.open(filename, std::ios::trunc);							//打开待写入文件
		if (!result)
		{
			std::cerr << "failed to open result.txt" << std::endl;		//打开失败
			return;
		}
		T* pointer=dataset.ptr();										//some question,必须做一个转换
		for (i = 0; i < dataset.rows; i++)								//依次写入
		{
			for (j = 0; j < dataset.cols; j++)
			{
				result << *pointer << "\t";
				pointer += 1;											//这里+1,注意是T*
			}
			result << std::endl;
		}
		result.close();
		pointer=NULL;
		return;
	}
}



#endif //SAVE_RESULT_TXT_H_
