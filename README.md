/**********************************************************
*工程说明
**********************************************************/

simple_flann包括一个KDTrees检索算法和一个基于p稳定分布的LSH
算法。
KDTrees算法兼容于flann库原有的代码结构。
LSH算法则是从原来的flann库中独立出来的算法，不兼容flann库原有的代码结构。
只要包含必要的几个文件即可单独运行。

/**********************************************************
*文档结构
**********************************************************/

flann.vcxproj		VS2013建立的测试工程，用兼容的VS编译器
					打开即能运行。

test_main.cpp		这个源文件中包含了测试的主要函数，
					添加这个源文件到空白工程即可快速开始实验。
				
dataset/sift 		文件夹下面存放的是SIFT1M测试数据集。
dataset/gist 		文件夹下面存放的是GIST1M测试数据集。
				
flann/algorithms	包含LSH和KDTrees实现的具体源文件。	
flann/io		包含读取读取和存储数据的源文件。
flann/util		常用的一些数据结构等。

output				存储输出的文件

/**********************************************************
*快速开始
**********************************************************/

将test_main.cpp添加到新建立的空白工程中。不必改动其他文件的路径。

要求编译器支持中文输入输出，因为程序中有多处输出中文的地方。
要求将相应的测试数据集放到dataset文件夹中的两个文件夹下。否则要在程序中修改读取数据的路径。

flann/io文件夹下的dataset_read.h文件提供了一个函数用以读取ivecs和fvecs格式的文件。
fvecs_ivecs_read(Matrix<T>& dataset, const std::string& filename, int start, int end)
dataset			-	Matrix类型的类。
filename		-	文件路径
start			-	开始的位置。从1开始计算。
end			-	结束的位置。不能超过文档本身包含的数据个数。

test_main.cpp中提供了4个测试函数，用以测试不同数据集下不同的算法的性能。
摘录如下。

kd-tree在SIFT1M数据集下的测试函数
void test_kdtree_sift(size_t knn,int tree_num)
knn			-	搜索的近邻个数
tree_num		-	搜索时要使用的树的个数

LSH在SIFT1M数据集下的测试函数
void test_lsh_sift(size_t knn,unsigned int table_number,unsigned int key_size,float gap_w)
knn			-	搜索的近邻个数
table_number		-	用的哈希表个数
key_size		-	关键字的个数。第一次哈希后得到的向量维数。
gap_w			-	分割间隔

kd-tree在GIST1M数据集下的测试函数
void test_kdtree_gist(size_t knn,int tree_num)
knn			-	搜索的近邻个数
tree_num		-	搜索时要使用的树的个数

LSH在GIST1M数据集下的测试函数
void test_lsh_gist(size_t knn,unsigned int table_number,unsigned int key_size,float gap_w)
knn			-	搜索的近邻个数
table_number		-	用的哈希表个数
key_size		-	关键字的个数。第一次哈希后得到的向量维数。
gap_w			-	分割间隔

每个函数根据输入的参数，均独立的完成了读取数据，建立索引，查询数据，输出信息的功能。

在test_main.cpp的开头有一个宏定义
#define LSH_INDEX_ENABLE (1)	
为0选择kdtree方式检索，为1选择lsh方式检索。
可以修改代码忽略掉这个定义。

在test_main.cpp的int main( void )函数中，
有几组已经存在和测试过的测试组。以供参考。

更多的详细信息见注释。

/**********************************************************
*注意事项
**********************************************************/

虽然GIST1M数据集在测试时分为三次读取、查询，但有时仍可能分配内存失败，
请保证系统有足够的空间进行分配，否则程序可能会运行失败。

目前只支持fvecs和ivecs文件的读取，若是其他格式请修改读取数据的代码。

当测试的参数不当时可能程序的运行时间会非常长。

源文件采用的都是相对路径，改动文件路径请修改对应的代码。







