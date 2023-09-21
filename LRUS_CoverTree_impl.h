#pragma once


LRUS_CoverTree::LRUS_CoverTree(
	const arma::mat& dataset,
	const double base,
	const int min_scale) :
	dataset(&dataset),
	base(base),
	scale(std::ceil( std::log(2) / std::log(base))),
	packing_scale(std::ceil(std::log(2) / std::log(base))),
	min_scale(min_scale),
	parent(NULL),
	directionalParentDistance(0),
	directionalParentDistanceSquared(0)
	//parentAngle(0),
	//parentDistance(0) 
{
	normalized = new arma::mat(dataset.n_rows, dataset.n_cols);
	norms = new arma::vec(dataset.n_cols);

	for (int i = 0; i < dataset.n_cols; i++) {
		double normi = norm(dataset.col(i), 2);
		normalized->col(i) = dataset.col(i) / normi;
		norms->at(i) = normi;
	}
	//RootPointPolicy: the point with max norm is root.

	point = norms->index_max();//Use armadillo function
	l2Norm = norms->at(point);

	for (size_t i = 0; i < dataset.n_cols; i++) {
		if (i == point) continue;
		double dist = mlpack::EuclideanDistance::Evaluate(normalized->col(point), normalized->col(i));
		add_descendant(i, dist);
	}

	setup_info();
	//Expand the root on initialization.
	Expand();
	std::stack<LRUS_CoverTree*> stack;
	stack.push(this);
	while (!stack.empty()) {
		LRUS_CoverTree* current = stack.top();
		//current->setup_info();
		std::sort(current->close_descendants.begin(), current->close_descendants.end(),
				  [&](size_t i1, size_t i2) {return norms->at(i1) > norms->at(i2); });
		stack.pop();
		for (size_t i = 0; i < current->NumChildren(); i++) {
			stack.push(current->Children()[i]);
		}
	}

}
LRUS_CoverTree::LRUS_CoverTree(
	const arma::mat* dataset,
	arma::mat* normalized,
	arma::vec* norms,
	size_t index,
	int scale,
	LRUS_CoverTree* parent,
	const double base,
	const int min_scale) :
	dataset(dataset),
	normalized(normalized),
	norms(norms),
	point(index),
	furthestDescendant(index),
	scale(scale),
	min_scale(min_scale),
	packing_scale(scale),
	parent(parent),
	base(base),
	directionalFurthestDescendantDistance(0),
	directionalFDDSquared(0),
	furthestDescendantDistance(0),
	longestDescendantNorm(0) {
	l2Norm = norms->at(point);
	if (parent != nullptr) {
		directionalParentDistanceSquared = mlpack::LMetric<2, false>::Evaluate(normalized->unsafe_col(point), normalized->unsafe_col(parent->Point()));
		directionalParentDistance = std::sqrt(directionalParentDistanceSquared);
		parentCosine = 1 - 0.5 * directionalParentDistanceSquared;
		parentSine = std::sqrt(1 - parentCosine * parentCosine);
		//parentAngle = std::acos(parentCosine);
		//parentDistance = mlpack::EuclideanDistance::Evaluate(dataset->unsafe_col(point), dataset->unsafe_col(parent->Point()));
	} else {
		directionalParentDistanceSquared = 0;
		directionalParentDistance = 0;
		parentCosine = 0;
		parentSine = 0;
		//parentAngle = 0;
		//parentDistance = 0;
	}
}
void LRUS_CoverTree::setup_info() {
	//packing_scale = scale;
	scale = ceil(log(directionalFurthestDescendantDistance) / log(base));
	numDescendants = descendants.size();
	directionalFDDSquared = directionalFurthestDescendantDistance * directionalFurthestDescendantDistance;
	furthestDescendantCosine = 1 - 0.5 * directionalFDDSquared;
	furthestDescendantSine = std::sqrt(1 - furthestDescendantCosine * furthestDescendantCosine);
	packing_dist = std::pow(base, packing_scale);
	packing_dist_squared = packing_dist * packing_dist;
	packing_cosine = 1 - 0.5 * packing_dist_squared;
	packing_sine = std::sqrt(1 - packing_cosine * packing_cosine);
	if (parent != nullptr) {
		direcdist_squared_parent_fd = mlpack::LMetric<2, false>::Evaluate(normalized->col(parent->Point()), normalized->col(furthestDescendant));
		direcdist_parent_fd = sqrt(direcdist_squared_parent_fd);
		if (direcdist_parent_fd > directionalFurthestDescendantDistance + directionalParentDistance) {
			std::cout << direcdist_parent_fd << " " << directionalFurthestDescendantDistance << " " << directionalParentDistance << " " << "\n";
		}
		cosine_parent_fd = 1 - 0.5 * direcdist_squared_parent_fd;
		sine_parent_fd = std::sqrt(1 - cosine_parent_fd * cosine_parent_fd);
	}
	//std::sort(close_descendants.begin(), close_descendants.end(),
		//	  [&](size_t i1, size_t i2) {return norms->at(i1) > norms->at(i2); });

}

bool LRUS_CoverTree::add_descendant(size_t point, double dist) {
	if (dist > directionalFurthestDescendantDistance) {
		directionalFurthestDescendantDistance = dist;
		furthestDescendant = point;
	}
	if (norms->at(point) > longestDescendantNorm) longestDescendantNorm = norms->at(point);
	if (norms->at(point) < shortestDescendantNorm) shortestDescendantNorm = norms->at(point);
	double fdd = mlpack::EuclideanDistance::Evaluate(dataset->col(this->point), dataset->col(point));
	if (fdd > furthestDescendantDistance) furthestDescendantDistance = fdd;
	if (dist > std::pow(base, min_scale)) {
		descendants.push_back(point);
		return false;
	} else {
		//if (point == 0) {
		//	std::cout << dist << " " << std::pow(base, min_scale) << "\n";
		//}
		close_descendants.push_back(point);
		return true;
	}
}

void LRUS_CoverTree::add_child(size_t index, size_t& pointLeft) {
	LRUS_CoverTree* child = new LRUS_CoverTree(dataset, normalized, norms, index, scale - 1, this, base, min_scale);
	children.push_back(child);
	for (size_t iter = 0; iter < pointLeft;) {
		double dist = mlpack::EuclideanDistance::Evaluate(normalized->col(index), normalized->col(descendants[iter]));
		if (dist <= std::pow(base, scale - 1)) {
			child->add_descendant(descendants[iter], dist);
			std::swap(descendants[iter], descendants[pointLeft - 1]);
			pointLeft--;
		} else iter++;
	}
	child->setup_info();
}

size_t LRUS_CoverTree::find_norm_max(const arma::mat* norms, const std::vector<size_t>& indices, size_t pointLeft) {
	if (pointLeft == 0)std::cout << "WRONG, FindMaxNorm in 0 size\n";
	if (pointLeft == 1) return 0;
	size_t resind = 0; double maxnorm = norms->at(indices[0]);
	for (size_t i = 1; i < pointLeft; i++) {
		if (norms->at(indices[i]) > maxnorm) {
			maxnorm = norms->at(indices[i]);
			resind = i;
		}
	}
	return resind;
}

void LRUS_CoverTree::Expand() {
	// Since we perform full expand for each expand, we do not need to delete an arbitrary element from sescendants. 
	// Thus, it need not to be a std::set, std::vector is enough
	size_t pointLeft = numDescendants;

	while (pointLeft > 0) {
		size_t index = find_norm_max(norms, descendants, pointLeft);//LongIsRoot
		size_t point = descendants[index];
		std::swap(descendants[index], descendants[pointLeft - 1]);
		pointLeft--;
		add_child(point, pointLeft);
	}
	//After the children are created, expand each of them.
	for (auto& child : children) child->Expand();
}


