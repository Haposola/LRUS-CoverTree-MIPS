#pragma once

bool compare_results(std::vector<size_t>& res1, std::vector<size_t>& res2) {
	if (res1.size() != res2.size()) return false;
	bool res = true;
	for (size_t i = 0; i < res1.size(); i++) {
		if (res1[i] != res2[i])res = false;
	}
	return res;
}
bool varify_approximation(std::vector<double> approx, std::vector<double> exact, double eps) {
	assert(approx.size() == exact.size());
	std::sort(approx.begin(), approx.end());
	std::sort(exact.begin(), exact.end());
	if (approx[0] < exact[0] * eps)return false;
	else return true;
}
double compute_recall(std::vector<size_t>& approx, std::vector<size_t>& exact) {
	
	std::set<size_t> setoftruth;
	for (size_t i = 0; i < approx.size(); i++)
		setoftruth.insert(exact[i]);
	size_t numrecall = 0;
	for (size_t i = 0; i < approx.size(); i++) {
		if (setoftruth.contains(approx[i])) numrecall++;
	}
	return (double)numrecall / (double)(approx.size());
}
double compute_overallratio(std::vector<double>& approx, std::vector<double>& exact) {
	
	double sum_ratio = 0;
	for (size_t i = 0; i < approx.size(); i++)
		sum_ratio += approx[i] / exact[i];
	return sum_ratio / (double)approx.size();
}


void load_fvecs_query(std::string& query_path, arma::mat& query_mat) {
	std::ifstream ifile(query_path, std::ios::binary);
	int dim, num;
	ifile.read((char*)(&dim), sizeof(int));
	size_t size = std::filesystem::file_size(query_path);
	num = size / sizeof(float) / (dim + 1);
	std::cout << "query_dim " << dim << " query_num " << num << "\n";
	ifile.close();
	ifile.open(query_path, std::ios::binary);
	query_mat.set_size(dim, num);
	float ele; int iele;
	for (size_t i = 0; i < num; i++) {
		ifile.read((char*)(&iele), sizeof(iele));
		for (size_t j = 0; j < dim; j++) {
			ifile.read((char*)(&ele), sizeof(ele));
			query_mat(j, i) = ele;
		}
	}
}
size_t load_ives_groundtruth(std::string& groundtruth_path, std::vector<std::vector<size_t>>& groundtruth) {
	std::ifstream ifile(groundtruth_path, std::ios::binary);
	int groundtruth_k, num;
	ifile.read((char*)(&groundtruth_k), sizeof(int));
	ifile.close();
	size_t size = std::filesystem::file_size(groundtruth_path);
	num = size / sizeof(int) / (groundtruth_k + 1);
	std::cout << "groundtruth_k " << groundtruth_k << " num " << num << "\n";
	
	ifile.open(groundtruth_path, std::ios::binary);
	groundtruth.resize(num);
	int iele;
	for (size_t i = 0; i < num; i++) {
		ifile.read((char*)(&iele), sizeof(iele));
		for (size_t j = 0; j < groundtruth_k; j++) {
			ifile.read((char*)(&iele), sizeof(iele));
			groundtruth[i].push_back(iele);
		}
	}
	return groundtruth_k;

}

void experiment_run(std::string& data_path, std::string& query_path, std::string& groundtruth_path, int k, double eps, std::ostream& out = std::cout, int min_scale = INT_MIN) {
	out << query_path << " results:\n";
	out << "tree_build_time(s)\tsearch_time(ms)\trecall\toverall_ratio\teps\tmin_scale\tmemeory\n";
	arma::mat data, query_mat;
	data.load(data_path, arma::arma_binary);

	std::cout << "num " << data.n_cols << " dim " << data.n_rows << "\n";
	load_fvecs_query(query_path, query_mat);
	
	std::vector<std::vector<size_t>> groundtruth;
	size_t groundtruth_k = load_ives_groundtruth(groundtruth_path, groundtruth);

	std::vector<std::vector<double>> truth_ipvalues;
	truth_ipvalues.resize(query_mat.n_cols);
	for (int i = 0; i < query_mat.n_cols; i++) {
		for (int j = 0; j < groundtruth_k; j++) {
			arma::vec query = query_mat.col(i); 
			double norm = arma::norm(query);
			query = query / norm;
			truth_ipvalues[i].push_back(arma::dot(query, data.col(groundtruth[i][j])));
		}
	}

	auto tree_build_start = std::chrono::high_resolution_clock::now();
	KMaxIP_LRUS xx(data, 1.3, min_scale);
	auto tree_build_end = std::chrono::high_resolution_clock::now();

	std::vector<double>  mip_lrus;
	std::vector<size_t>  max_lrus;

	double time_lrus = 0;
	double recall_lrus = 0;
	double overallratio_lrus = 0;

	arma::vec query;

	for (int i = 0; i < query_mat.n_cols; i++) {
		query = query_mat.col(i);
		//query = data.col(i);
		double norm = arma::norm(query, 2);
		query = query / norm;

		auto search_start= std::chrono::high_resolution_clock::now();
		xx.search(query, k, eps, mip_lrus, max_lrus);
		auto search_end= std::chrono::high_resolution_clock::now();

		time_lrus += double(std::chrono::duration_cast<std::chrono::microseconds> (search_end - search_start).count())/double(1000);
		recall_lrus += compute_recall(max_lrus, groundtruth[i]);
		overallratio_lrus += compute_overallratio(mip_lrus, truth_ipvalues[i]);
	}

	out<< double(std::chrono::duration_cast<std::chrono::microseconds> (tree_build_end - tree_build_start).count())/double(1000000)<<"\t";
	out << time_lrus / (query_mat.n_cols * 1.0) << " \t";
	out << recall_lrus / (query_mat.n_cols * 1.0) << " \t";
	out << overallratio_lrus / (query_mat.n_cols * 1.0) << " \t";
	out << eps << "\t";
	out << min_scale << "\t";
	out << double(getCurrentRSS())/double(1000000) << "\n";
	//out<< "\n";

}

