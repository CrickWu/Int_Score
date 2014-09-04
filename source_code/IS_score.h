#pragma once
#include <vector>
#include "TM_score.h"

//==== class: IS_score =====//
class IS_score : public TM_score
{
public:
	IS_score(int num=3000);
	~IS_score(void);

//---- variables ----//
public:
	vector < double > overlap_factor;  //length should be ali_orin

//---- functions ----//
public:
	void Calc_Overlap_Factor(vector < vector <int> > &in1, vector < vector <int> > &int2, 
		vector <pair<int,int> > &alignment ,vector < double > &out);
	double Calc_TM_Score_Single(XYZ *mol1,XYZ *mol2,int lali,double *rotmat_,double d0,double d8,int TM8orTM=0, double *ALLSCO=0);
};
