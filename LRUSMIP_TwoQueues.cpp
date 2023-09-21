// LRUSMIP_TwoQueues.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include "HeadOfHead.h"

int main(int argc, char** argv) {


#ifdef DEBUG
	//base_test();
	//test_NormOrder_Cover_Correct();
	test_run(std::cout);
#else

	if (argc < 6) {
		std::cout << "wrong num of args\n";
		return -1;
	}

	std::string path(argv[1]);
	std::string data_path = path + "_arma.bin";
	std::string query_path, groundtruth_path;
	int flag = atoi(argv[2]);
	if (flag) {//query inside
		query_path = path + "_inside_query.fvecs";
		groundtruth_path = path + "_inside_groundtruth.ivecs";
	} else {//query random
		query_path = path + "_random_query.fvecs";
		groundtruth_path = path + "_random_groundtruth.ivecs";
	}

	int k = atoi(argv[3]);
	double eps = atof(argv[4]);
	std::ofstream out(argv[5], std::ios::app);

	if (argc == 7) {
		int min_scale = atof(argv[6]);
		experiment_run(data_path, query_path, groundtruth_path, k, eps, out, min_scale);
	}
	if (argc == 6) {
		experiment_run(data_path, query_path, groundtruth_path, k, eps, out);
	}


	std::cout << "Hello World!\n";

#endif
}
// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
// 调试程序: F5 或调试 >“开始调试”菜单

// 入门使用技巧: 
//   1. 使用解决方案资源管理器窗口添加/管理文件
//   2. 使用团队资源管理器窗口连接到源代码管理
//   3. 使用输出窗口查看生成输出和其他消息
//   4. 使用错误列表窗口查看错误
//   5. 转到“项目”>“添加新项”以创建新的代码文件，或转到“项目”>“添加现有项”以将现有代码文件添加到项目
//   6. 将来，若要再次打开此项目，请转到“文件”>“打开”>“项目”并选择 .sln 文件
