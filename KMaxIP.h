#pragma once
#pragma once
//Priority: small first

class KMaxIP_LRUS {
public:

	KMaxIP_LRUS(arma::mat& data, double base = 2.0, int min_scale = INT_MIN) :
		base(base), min_scale(min_scale) {
		referenceTree = new LRUS_CoverTree(data, base, min_scale);
		squared_dist_close_desc = std::pow(base, min_scale) * std::pow(base, min_scale);
		cosine_close_desc = 1 - 0.5 * squared_dist_close_desc;
		sine_close_desc = std::sqrt(1 - cosine_close_desc * cosine_close_desc);
	}



	void search(arma::vec& query, int k, double eps, std::vector<double>& results_ipvalues, std::vector<size_t>& results_indexes) {
		assert(std::abs(arma::norm(query, 2) - 1) <= 1e-7);
		assert(eps <= 1);
		//resultsmap.clear();
		while (!results.empty())results.pop();
		results_ipvalues.clear(); results_indexes.clear();
		//const arma::mat& normalized = referenceTree->Normalized();
		//Info.type=1: search node; 0: closedesc
		std::priority_queue<MIPSearchInfo_Use_Sibling, std::vector<MIPSearchInfo_Use_Sibling>> searchQueue;
		std::priority_queue<MIPSearchInfo_Closedesc, std::vector<MIPSearchInfo_Closedesc>> closedescQueue;
		int tree_visited = 0;
		int close_visited = 0;

		double direcdist_squared_parent_query;
		double cosine_parent_query;
		double sine_parent_query;
		double ipvalue_parent;

		visit_node(referenceTree, query, k, eps, tree_visited, close_visited, closedescQueue, direcdist_squared_parent_query, cosine_parent_query, sine_parent_query, ipvalue_parent);

		bool parent_covered = (direcdist_squared_parent_query <= referenceTree->DirectionalFDDSquared());
		searchQueue.push(MIPSearchInfo_Use_Sibling(referenceTree->LongestDescendantNorm(),
												   direcdist_squared_parent_query,
												   cosine_parent_query, sine_parent_query,
												   parent_covered,
												   referenceTree));

		double direcdist_squared_child_query;
		double cosine_child_query;
		double ipvalue;
		double sine_child_query;
		double cosine_query_fd;
		bool parentCovered, siblingCovered, childCovered;
		double direcdist_gap_sibling_query, cosine_max_sibling;
		double close_desc_cosine_max;
		LRUS_CoverTree* child;
		LRUS_CoverTree* parent;

		while (searchQueue.size() > 0 || closedescQueue.size() > 0) {
			if ((results.size() >= k)
				&& (searchQueue.size() > 0 && eps * searchQueue.top().score <= candidate_ip_threshold())
				&& (closedescQueue.size() > 0 && eps * closedescQueue.top().score <= candidate_ip_threshold())) {
				break;
			}// Not possible to find an element with larger IP.
			
			if ( closedescQueue.empty() || (searchQueue.size() > 0 && searchQueue.top().score > closedescQueue.top().score)) {
				direcdist_squared_parent_query = searchQueue.top().distSquared;
				cosine_parent_query = searchQueue.top().cosine;
				sine_parent_query = searchQueue.top().sine;
				parent = searchQueue.top().node;
				parentCovered = searchQueue.top().covered;
				searchQueue.pop();

				siblingCovered = false; cosine_max_sibling = 1; direcdist_gap_sibling_query = 0;
				for (size_t i = 0; i < parent->NumChildren(); i++) {
					//LRUS_CoverTree*  child = parent->Children()[i];
					child = parent->Children()[i];
					if (results.size() >= k &&
						eps * ScoreNode_use_sibling(child,
													direcdist_squared_parent_query,
													cosine_parent_query,
													sine_parent_query,
													parentCovered,
													siblingCovered,
													direcdist_gap_sibling_query,
													cosine_max_sibling)
						<= candidate_ip_threshold()) continue;

					//visit_node(child, query, k, eps, tree_visited, close_visited,
						//	   direcdist_squared_child_query, cosine_child_query, sine_child_query, ipvalue);
					direcdist_squared_child_query = //arma::accu(arma::square(query - node->Normalized().unsafe_col(node->Point())));
						mlpack::LMetric<2, false>::Evaluate(query, child->Normalized().unsafe_col(child->Point()));
					cosine_child_query = 1 - 0.5 * direcdist_squared_child_query;
					sine_child_query = std::sqrt(1 - cosine_child_query * cosine_child_query);
					ipvalue = child->Norm() * cosine_child_query;

					add_results(k, child->Point(), ipvalue);
#ifdef OUT_VISITED
					tree_visited++;
#endif
					//visit close descendants
					if (!child->CloseDescendants().empty()) {
						double close_desc_cosine_max;
						if (direcdist_squared_child_query <= squared_dist_close_desc)
							close_desc_cosine_max = 1;
						else
							close_desc_cosine_max = cosine_child_query * cosine_close_desc + sine_child_query * sine_close_desc;
						//close_desc_cosine_max = 1;
						double score;
						if (close_desc_cosine_max > 0) {
							score = child->Norms().at(child->CloseDescendants()[0]) * close_desc_cosine_max;
						} else {
							score = child->Norms().at(child->CloseDescendants()[child->CloseDescendants().size() - 1]) * close_desc_cosine_max;
						}

						closedescQueue.push(MIPSearchInfo_Closedesc(score,
																	close_desc_cosine_max,
																	child));
					}
					// If covered, D(child,sibling)-child->fdd is used to set direcdist_gap_sibling_query
					// cosine<query,fd> is used to set cosine_max_sibling
					// If not covered, cosine<query,fd> is used to set score
					//cosine_query_fd = (cosine_child_query * child->FurthestDescendantCosine() + sine_child_query * child->FurthestDescendantSine());
					//childCovered = (!siblingCovered) && (direcdist_squared_child_query <= child->DirectionalFDDSquared());


					if (direcdist_squared_child_query <= child->PackingDistSquared()) {
						//compute cosine<query, fd>, which can be executed for at most once in a while loop.
						siblingCovered = true;
						direcdist_gap_sibling_query = child->PackingDistance() - std::sqrt(direcdist_squared_child_query);// 0;
						cosine_max_sibling = cosine_child_query * child->PackingCosine() + sine_child_query * child->PackingSine();// 1;
					}
					//add_results(k, child->Point(), ipvalue);

					if (!child->IsLeaf()) {
						double childScore = ScoreDescendants(child, ipvalue, direcdist_squared_child_query, cosine_child_query, sine_child_query);
						//double childScore = ScoreDescendants_use_sibling(child, ipvalue, childCovered, cosine_query_fd);
						if (results.size() < k || eps * childScore >= candidate_ip_threshold())
							searchQueue.push(MIPSearchInfo_Use_Sibling(childScore,
																	   direcdist_squared_child_query,
																	   cosine_child_query, sine_child_query,
																	   (direcdist_squared_child_query <= child->DirectionalFDDSquared()),
																	   child));
					}
				}


			} else if (closedescQueue.size() > 0) {
				child = closedescQueue.top().node;
				close_desc_cosine_max = closedescQueue.top().cosmax;
				closedescQueue.pop();
				if (close_desc_cosine_max > 0) {
					for (auto iter = child->CloseDescendants().begin(); iter != child->CloseDescendants().end(); iter++) {
						size_t cdi = *iter;
						if (results.size() < k || eps * child->Norms().at(cdi) * close_desc_cosine_max > candidate_ip_threshold()) {

							add_results(k, cdi, arma::dot(query, child->Dataset().col(cdi)));
#ifdef OUT_VISITED
							close_visited++;
#endif
						} else break;
					}
				} else {
					for (auto riter = child->CloseDescendants().rbegin(); riter != child->CloseDescendants().rend(); riter++) {
						size_t cdi = *riter;
						if (results.size() < k || child->Norms().at(cdi) * close_desc_cosine_max > candidate_ip_threshold()) {
							add_results(k, cdi, arma::dot(query, child->Dataset().col(cdi)));
#ifdef OUT_VISITED
							close_visited++;
#endif
						} else break;
					}
				}
			}
		}
		get_results(results_ipvalues, results_indexes);
#ifdef OUT_VISITED
		std::cout << tree_visited << ", " << close_visited << ", " << tree_visited + close_visited << " for sibling\n";
#endif

	}



	void naive(arma::vec& query, int k, std::vector<double>& results_ipvalues, std::vector<size_t>& results_indexes) {
		results_ipvalues.clear(); results_indexes.clear();
		for (size_t i = 0; i < referenceTree->Dataset().n_cols; i++) {
			add_results(k, i, arma::dot(query, referenceTree->Dataset().col(i)));
		}
		get_results(results_ipvalues, results_indexes);
	}
	inline double candidate_ip_threshold() {
		return results.top().second;
	}

private:
	LRUS_CoverTree* referenceTree;
	double base;
	int min_scale;
	typedef std::pair<size_t, double> Candidate;
	double squared_dist_close_desc;
	double cosine_close_desc;
	double sine_close_desc;

	void visit_node(LRUS_CoverTree* node, arma::vec& query, int k, double eps,
					int& tree_visited, int& close_visited,
					std::priority_queue<MIPSearchInfo_Closedesc, std::vector<MIPSearchInfo_Closedesc>>& closedescQueue,
					double& direcdist_squared_node_query,
					double& cosine_node_query,
					double& sine_node_query,
					double& ipvalue_node) {
		direcdist_squared_node_query = //arma::accu(arma::square(query - node->Normalized().unsafe_col(node->Point())));
			mlpack::LMetric<2, false>::Evaluate(query, node->Normalized().unsafe_col(node->Point()));
		cosine_node_query = 1 - 0.5 * direcdist_squared_node_query;
		sine_node_query = std::sqrt(1 - cosine_node_query * cosine_node_query);
		ipvalue_node = node->Norm() * cosine_node_query;

		add_results(k, node->Point(), ipvalue_node);
#ifdef OUT_VISITED
		tree_visited++;
#endif
		//visit close descendants
		if (!node->CloseDescendants().empty()) {

			double close_desc_cosine_max;
			if (direcdist_squared_node_query <= squared_dist_close_desc)
				close_desc_cosine_max = 1;
			else
				close_desc_cosine_max = cosine_node_query * cosine_close_desc + sine_node_query * sine_close_desc;
			double score;
			if (close_desc_cosine_max > 0) {
				score = node->Norms().at(node->CloseDescendants()[0]) * close_desc_cosine_max;
			} else {
				score = node->Norms().at(node->CloseDescendants()[node->CloseDescendants().size() - 1]) * close_desc_cosine_max;
			}

			closedescQueue.push(MIPSearchInfo_Closedesc(score,close_desc_cosine_max,node));


		}
		
	}

	struct paircomp {
		//priority queue is large first by default. Use > to ensure small first
		bool operator()(const Candidate& c1, const Candidate& c2) const {
			return c1.second > c2.second;
		}
	};
	std::priority_queue < std::pair<size_t, double>, std::vector<std::pair<size_t, double>>, paircomp> results;
	std::map<double, std::vector<size_t>> resultsmap;
	inline void add_results(int k, size_t index, double ipvalue) {
		if (results.size() < k) results.push(std::make_pair(index, ipvalue));
		else if (results.top().second < ipvalue) {
			//size_t topindex = results.top().first;
			//double topip= results.top().second;

			results.push(std::make_pair(index, ipvalue));
			results.pop();
			//if (topip == results.top().second) {
				//results.push(std::make_pair(topindex, topip));
			//}
		}

		//if (resultsmap.contains(ipvalue)) {
		//	resultsmap[ipvalue].push_back(index);
		//} else {
		//	std::vector<size_t> vec;
		//	vec.push_back(index);
		//	resultsmap.emplace(ipvalue, vec);
		//}
		//if (resultsmap.size() > k) {
		//	resultsmap.erase(resultsmap.begin());//map is in ascending order by default
		//}
	}
	inline void get_results(std::vector<double>& results_ipvalues, std::vector<size_t>& results_indexes) {
		//results is small first

		results_ipvalues.resize(results.size());
		results_indexes.resize(results.size());
		size_t index = results.size() - 1;
		while (results.size() > 0) {
			results_indexes[index] = (results.top().first);
			results_ipvalues[index] = (results.top().second);
			results.pop();
			index--;
		}
		/*for (auto riter = resultsmap.rbegin(); riter != resultsmap.rend();riter++) {
			for (auto index : riter->second) {
				results_indexes.push_back(index);
				results_ipvalues.push_back(riter->first);
			}
		}*/
	}
};