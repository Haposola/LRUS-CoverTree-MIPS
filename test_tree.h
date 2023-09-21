#pragma once

int test_NormOrder_Cover_Correct_dotest(LRUS_CoverTree& node) {
	size_t root = node.Point();
#ifdef TEST_COUT
	std::cout << "root: " << root << " at scale " << node.Scale() << " ";
	std::cout << "descendants: ";
#endif
	//close_descendants_correct
	for (size_t cdi = 0; cdi < node.CloseDescendants().size(); cdi++) {
		if (mlpack::EuclideanDistance::Evaluate(node.Normalized().col(node.Point()), node.Normalized().col(node.CloseDescendants()[cdi])) > std::pow(node.Base(), node.MinScale())) {
			std::cout << "close_desc not close\n";
		}
		if (node.Norms().at(node.CloseDescendants()[cdi]) < node.Norms().at(node.CloseDescendants()[cdi])) {
			std::cout << "close_desc not sorted\n";
		}
	}
	//desc_in_tree_match_desc_in_vec
	std::stack<LRUS_CoverTree*> stack;
	std::set<size_t> descendants_in_tree;
	stack.push(&node);
	while (!stack.empty()) {
		LRUS_CoverTree* current = stack.top();
		stack.pop();
		for (size_t i = 0; i < current->NumChildren(); i++) {
			stack.push(current->Children()[i]);
			descendants_in_tree.insert(current->Children()[i]->Point());
			for (size_t j = 0; j < current->Children()[i]->CloseDescendants().size(); j++) {
				descendants_in_tree.insert(current->Children()[i]->CloseDescendants()[j]);
			}
		}

	}
	if (descendants_in_tree.size() != node.NumDescendants()) {
		std::cout << "desc_in_tree no_match desc_in_vec, ";
		if (descendants_in_tree.size() > node.NumDescendants()) std::cout << "desc_in_tree more, dif " << descendants_in_tree.size() - node.NumDescendants() << "\n";
		else std::cout << "desc_in_vec more\n";
	}
	for (size_t i = 0; i < node.NumDescendants(); i++) {

		size_t desc = node.Descendant(i);
#ifdef TEST_COUT
		//std::cout << desc << " ";
#endif
		if (!descendants_in_tree.contains(desc))std::cout << "desc_in_vec not_in desc_in_tree\n";
		if (node.Norms().at(desc) > node.Norms().at(root)) {
			std::cout << "root: " << root << " norm: " << node.Norms().at(root) << " scale: " << node.Scale();
			std::cout << " point: " << node.Descendant(i) << " norm: " << node.Norms().at(node.Descendant(i)) << "\n";
			return -1;
		}
	}
#ifdef TEST_COUT
	std::cout << ". child list: ";
#endif
	for (size_t i = 0; i < node.NumChildren(); i++) {
		size_t i2 = node.Child(i).Point();
#ifdef TEST_COUT
		std::cout << i2 << " ";
#endif
		if (mlpack::EuclideanDistance::Evaluate(node.Normalized().col(root), node.Normalized().col(i2)) > pow(node.Base(), node.Scale())) {
			return -2;
		}
	}
#ifdef TEST_COUT
	std::cout << ".\n";
#endif

	for (size_t i = 0; i < node.NumChildren(); i++) {
		int res = test_NormOrder_Cover_Correct_dotest(node.Child(i));
		if (res != 1) return res;
	}
	return 1;
}



void test_NormOrder_Cover_Correct() {
	std::string ifname("E:\\datasets\\audio\\audio.data.txt");
	//std::string ifname("test_data.txt");
	arma::mat data;
	mlpack::data::Load(ifname, data);

	LRUS_CoverTree xx(data, 1.3, -6);
	std::cout << "tree built.\n";


	//std::cout << xx.Norms()<<"\n";
	//std::cout << xx.Normalized() << "\n";
	//std::cout << xx.NumChildren() << "\n";
	//std::cout << xx.Child(1).Norms() << "\n";
	int res = test_NormOrder_Cover_Correct_dotest(xx);
	if (res == -1) std::cout << "norm wrong" << "\n";
	if (res == -2)std::cout << "cover wrong" << "\n";
	if (res == 1)std::cout << "correct" << "\n";


}