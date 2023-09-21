#pragma once


class LRUS_CoverTree {
public:
	LRUS_CoverTree(const arma::mat& dataset,
				   const double base = 2.0,
				   const int min_scale = INT_MIN);
	LRUS_CoverTree(const arma::mat* dataset,
				   arma::mat* normalized,
				   arma::vec* norms,
				   size_t index,
				   int scale,
				   LRUS_CoverTree* parent,
				   const double base = 2.0,
				   const int min_scale = INT_MIN);
	//! Get a reference to the dataset.
	const arma::mat& Dataset() const {
		return *dataset;
	}
	const arma::mat& Normalized() const {
		return *normalized;
	}
	const arma::vec& Norms() const {
		return *norms;
	}
	//! Get the index of the point which this node represents. 
	size_t Point() const {
		return point;
	}
	//! For compatibility with other trees; the argument is ignored.
	size_t Point(const size_t) const {
		return point;
	}

	bool IsLeaf() const {
		return (children.size() == 0);
	}
	size_t NumPoints() const {
		return 1;
	}

	//! Get a particular child node.
	const LRUS_CoverTree& Child(const size_t index) const {
		return *children[index];
	}
	//! Modify a particular child node.
	LRUS_CoverTree& Child(const size_t index) {
		return *children[index];
	}

	LRUS_CoverTree*& ChildPtr(const size_t index) {
		return children[index];
	}

	//! Get the number of children.
	size_t NumChildren() const {
		return children.size();
	}

	//! Get the children.
	const std::vector<LRUS_CoverTree*>& Children() const {
		return children;
	}
	//! Modify the children manually (maybe not a great idea).
	std::vector<LRUS_CoverTree*>& Children() {
		return children;
	}
	const std::vector<size_t>& CloseDescendants() {
		return close_descendants;
	}
	//! Get the number of descendant points.
	size_t NumDescendants() const {
		return numDescendants;
	};

	//! Get the index of a particular descendant point.
	size_t Descendant(const size_t index) const {
		return descendants[index];
		//It is correct since we did not delete anything from @descendants, only swap used points to the rare parts.
	}

	//! Get the scale of this node.
	int Scale() const {
		return scale;
	}

	int PackingScale() const {
		return packing_scale;
	}
	double PackingDistSquared()const {
		return packing_dist_squared;
	}
	double PackingDistance()const {
		return packing_dist;
	}
	double PackingCosine() const {
		return packing_cosine;
	}
	double PackingSine()const {
		return packing_sine;
	}
	int MinScale() const {
		return min_scale;
	}
	//! Get the base.
	double Base() const {
		return base;
	}

	//! Get the parent node.
	LRUS_CoverTree* Parent() const {
		return parent;
	}
	//! Modify the parent node.
	LRUS_CoverTree*& Parent() {
		return parent;
	}
	double DirectionalParentDistance() const {
		return directionalParentDistance;
	}
	double DirectionalParentDistanceSquared() const {
		return directionalParentDistanceSquared;
	}
	double DirectionalFurthestDescendantDistance() const {
		return directionalFurthestDescendantDistance;
	}

	double DirectionalFDDSquared() const {
		return directionalFDDSquared;
	}
	double FurthestDescendantCosine() const {
		return furthestDescendantCosine;
	}
	double FurthestDescendantSine() const {
		return furthestDescendantSine;
	}
	double Norm() const {
		return l2Norm;
	}
	//double FurthestDescendantAngle() const {		return furthestDescendantAngle;	}
	double FurthestDescendantDistance() const {
		return furthestDescendantDistance;
	}
	double LongestDescendantNorm() const {
		return longestDescendantNorm;
	}
	double ShortestDescendantNorm() const {
		return shortestDescendantNorm;
	}
	double ParentFDCosine() const {
		return cosine_parent_fd;
	}
	double ParentFDSine() const {
		return sine_parent_fd;
	}
	double ParentFDDirectionalDistance() const {
		return direcdist_parent_fd;
	}
	double ParentFDDirectionalDistanceSquared() const {
		return direcdist_squared_parent_fd;
	}
	double ParentCosine() const {
		return parentCosine;
	}
	double ParentSine() const {
		return parentSine;
	}

	//double ParentDistance() const {		return parentDistance;}
	//double ParentAngle() const {		return parentAngle;	}
	void Expand();
	//The bool return value indicates whether the inserted point is a close_descendant
	bool add_descendant(size_t index, double dist);
	void setup_info();
private:

	std::vector<size_t> descendants;
	std::vector<size_t> close_descendants;

	//std::unordered_map<size_t, double> desc_distances;


	//! Reference to the matrix which this tree is built on.
	const arma::mat* dataset;
	arma::mat* normalized;
	arma::vec* norms;

	//! Index of the point in the matrix which this node represents.
	size_t point;
	//! The list of children; the first is the self-child.
	std::vector<LRUS_CoverTree*> children;
	//! Scale level of the node.
	int scale;
	//! the min_scale used to determine the close neighbors.
	int min_scale;
	int packing_scale;
	double packing_dist;
	double packing_dist_squared;
	double packing_cosine;
	double packing_sine;
	//! The base used to construct the tree.
	double base;

	//! The number of descendant points.
	size_t numDescendants;
	//! The parent node (NULL if this is the root of the tree).
	LRUS_CoverTree* parent;

	//! Information of descendants and parent.
	//! Distacn
	double directionalFurthestDescendantDistance;
	double directionalFDDSquared;

	double furthestDescendantSine;
	double furthestDescendantCosine;
	//double furthestDescendantAngle;


	double furthestDescendantDistance;
	//double parentDistance;
	double directionalParentDistance;
	double directionalParentDistanceSquared;
	double parentCosine;
	double parentSine;
	//double parentAngle;

	double cosine_parent_fd;
	double sine_parent_fd;
	double direcdist_parent_fd;
	double direcdist_squared_parent_fd;
	size_t furthestDescendant;
	double longestDescendantNorm;
	double shortestDescendantNorm;
	double l2Norm;

	//void compute_distances();
	//void delete_from_list(size_t index);


	void add_child(size_t index, size_t& pointLeft);
	size_t find_norm_max(const arma::mat* norms, const std::vector<size_t>& indices, size_t pointLeft);
	//insert a point into the constructed tree.
	// Used for insert points with shorter norm in construction.
	void insert(size_t index);
	//LRUS_CoverTree* create_childptr(size_t index);
};
#include "LRUS_CoverTree_impl.h"