#pragma once



inline double ScoreNode(LRUS_CoverTree* child, double direcdist_squared_parent_query, double cosine_parent_query, double sine_parent_query) {

	double score = child->Norm();
	if (direcdist_squared_parent_query > child->ParentFDDirectionalDistanceSquared())
		score *= (cosine_parent_query * child->ParentFDCosine() + sine_parent_query * child->ParentFDSine());
	return score;
}

inline double ScoreDescendants(LRUS_CoverTree* child, double dotvalue, double childDirecDistSquared, double cosine_child_query, double sine_child_query) {
	//return child->LongestDescendantNorm();
	double score1 = dotvalue + child->Norm() * child->DirectionalFurthestDescendantDistance();
	double score3 = dotvalue + child->FurthestDescendantDistance();
	
	
	double cosmax;
	if (childDirecDistSquared > child->DirectionalFDDSquared()) {
		cosmax= (cosine_child_query * child->FurthestDescendantCosine() + sine_child_query * child->FurthestDescendantSine());
	} else cosmax = 1;
	double score2;
	if (cosmax > 0) score2 = child->LongestDescendantNorm() * cosmax;
	else score2 = child->ShortestDescendantNorm() * cosmax;

	return std::min(score1, std::min(score2, score3));
}


inline double ScoreNode_use_sibling(
	LRUS_CoverTree* child,
	double direcdist_squared_parent_query,
	double cosine_parent_query,
	double sine_parent_query,
	bool parentCover,
	bool siblingCover,
	double siblingGap,
	double cosineMax) {
	//return child->Norm();
	if (!parentCover) {//This node can not cover  the query. Use parent-fd bound
		double cos = (cosine_parent_query * child->ParentFDCosine() + sine_parent_query * child->ParentFDSine());
		if (cos > 0) return child->Norm() * cos;
		else return child->ShortestDescendantNorm() * cos;
	} else {/*This node may cover the query.*/
		if (siblingCover) {
		//	// If sibling has covered, cosineMax would be less than 1.
		//	// Use D(child->fd, sibling->fd) as one score
		//	// Take maximum in case the estimated gap is less than 0

			double gap = std::max(siblingGap - child->DirectionalFurthestDescendantDistance(), 0.0);
			double cosine_sibling = std::min(cosineMax, 1 - 0.5 *gap * gap);
			if (cosine_sibling > 0)return child->Norm() * cosine_sibling;
			else return child->ShortestDescendantNorm() * cosine_sibling;

		} else {// Compare the parent-fd distance and parent-query distance, to see if this node may cover the query.
			//Which is the same with no-sibling score.
			double scoreParent = child->Norm();
			double cosmax;

			if (direcdist_squared_parent_query > child->ParentFDDirectionalDistanceSquared()) {
				cosmax = (cosine_parent_query * child->ParentFDCosine() + sine_parent_query * child->ParentFDSine());
			} else cosmax = 1;
			if (cosmax > 0)return child->Norm() * cosmax;
			else return child->ShortestDescendantNorm()*cosmax;
		}
	
	} 
}
