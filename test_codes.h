#pragma once

void base_test() {
	std::map<double, std::vector<int>> x;
	std::vector<int> v1 = { 1,2,3 };
	std::vector<int> v2 = { 1,2,3,4 };
	x.emplace(3.0, v1);
	x.emplace(2.0, v2);
	std::cout << x.begin()->second.size() << std::endl;
}




void test_run(std::ostream& out) {
	int min_scale = -3; double eps = 0.8;
	arma::mat data; int k = 50;
	data.load("E:\\datasets\\sift\\sift_arma.bin",arma::arma_binary);
	out << "num\tdim\ttree_build_time\tsearch_time\tnaive_time\trecall\toverall_ratio\n";
	out << data.n_cols << " " << data.n_rows << " \n";

#ifdef _WIN32
	LARGE_INTEGER t1, t2, t3, t4, t5, t6, t7, tc;
	QueryPerformanceFrequency(&tc);
	QueryPerformanceCounter(&t1);
#else defined __linux__
	struct timeval t1, t2, t3, t4, t5, t6, t7, tc;
	gettimeofday(&t1, NULL);
#endif
	KMaxIP_LRUS xx(data, 1.3, min_scale);
#ifdef _WIN32
	QueryPerformanceCounter(&t2);//Tree build time;
#else defined __linux__
	gettimeofday(&t2, NULL);
#endif
	std::vector<double>  mip_naive, mip_sibling;
	std::vector<size_t> max_naive, max_sibling;

	double time_naive = 0, time_sibling = 0;
	double recall_sibling = 0;
	double overallratio_sibling = 0;
	arma::arma_rng::set_seed(time(0));
	arma::vec query;

	int numTests = 20;

	for (int i = 0; i < numTests; i++) {
		//query = arma::randu(data.n_rows);
		//query = arma::ones(data.n_rows);
		int ind = rand();
		//query = data.col(ind);
		query = data.col(i);
		double norm = arma::norm(query, 2);
		query = query / norm;
		//std::cout << arma::norm(query, 2) << "\n";

		//exact_maxInnerProduct_useInfo(xx, query, mip_info, max_info);
#ifdef _WIN32
		QueryPerformanceCounter(&t3);
#else defined __linux__
		gettimeofday(&t3, NULL);
#endif
		xx.naive(query, k, mip_naive, max_naive);
#ifdef _WIN32
		QueryPerformanceCounter(&t4);
#else defined __linux__
		gettimeofday(&t4, NULL);
#endif
		xx.search(query, k, eps, mip_sibling, max_sibling);
#ifdef _WIN32
		QueryPerformanceCounter(&t5);

#else defined __linux__
		gettimeofday(&t5, NULL);
#endif

		time_naive += (t4.QuadPart - t3.QuadPart) * 1.0 / tc.QuadPart;
		time_sibling += (t5.QuadPart - t4.QuadPart) * 1.0 / tc.QuadPart;
		
		recall_sibling += compute_recall(max_sibling, max_naive);
		overallratio_sibling += compute_overallratio(mip_sibling, mip_naive);
	}

	out << std::setw(12)<<(t2.QuadPart - t1.QuadPart) * 1.0 / tc.QuadPart << "\t";
	//out << (t3.QuadPart - t2.QuadPart) * 1.0 / tc.QuadPart << "\t";
	out << time_sibling / (numTests * 1.0) << "\t" << time_naive / (numTests * 1.0) << "\t";

	out << recall_sibling / (numTests * 1.0) << "\t";

	out << overallratio_sibling / (numTests * 1.0) << "\t";
}
