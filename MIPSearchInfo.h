#pragma once


class MIPSearchInfo {
public:
	MIPSearchInfo(double score,  double distSquared,double cosine, double sine,  LRUS_CoverTree* node):
		score(score), distSquared(distSquared), cosine(cosine),sine(sine), node(node)
	{
		//Nothing more to be initiallized
	}
	double score;
	double distSquared;
	double cosine;
	double sine;
	LRUS_CoverTree* node;
	//Priotity: large first
	friend bool operator<(const MIPSearchInfo& info1, const MIPSearchInfo& info2) {
		return info1.score< info2.score;
	}
};


class MIPSearchInfo_Use_Sibling {
public:
	MIPSearchInfo_Use_Sibling(double score, double distSquared, double cosine, double sine,bool covered, LRUS_CoverTree* node) :
		score(score), distSquared(distSquared), cosine(cosine), sine(sine),covered(covered), node(node) 
	{
		//Nothing more to be initiallized
	}
	double score;
	double distSquared;
	double cosine;
	double sine;
	bool covered;
	LRUS_CoverTree* node;
	//Priotity: large first
	friend bool operator<(const MIPSearchInfo_Use_Sibling& info1, const MIPSearchInfo_Use_Sibling& info2) {
		return info1.score < info2.score;
	}
};



class MIPSearchInfo_Composite {
public:
	MIPSearchInfo_Composite(double score, bool type, 
							double distSquared, double cosine, double sine, bool covered, 
							double cosmax,
							LRUS_CoverTree* node) :
							
		score(score), type(type),
		distSquared(distSquared), cosine(cosine), sine(sine), covered(covered), 
		cosmax(cosmax),
		node(node) {
		//Nothing more to be initiallized
	}
	double score;
	bool type;
	double distSquared;
	double cosine;
	double sine;
	bool covered;
	double cosmax;
	LRUS_CoverTree* node;
	//Priotity: large first
	friend bool operator<(const MIPSearchInfo_Composite& info1, const MIPSearchInfo_Composite& info2) {
		return info1.score < info2.score;
	}
};


class MIPSearchInfo_Closedesc {
public:
	MIPSearchInfo_Closedesc(double score, 
							double cosmax,
							LRUS_CoverTree* node) :

		score(score), 
		cosmax(cosmax),
		node(node) {
		//Nothing more to be initiallized
	}
	double score;
	double cosmax;
	LRUS_CoverTree* node;
	//Priotity: large first
	friend bool operator<(const MIPSearchInfo_Closedesc& info1, const MIPSearchInfo_Closedesc& info2) {
		return info1.score < info2.score;
	}
};