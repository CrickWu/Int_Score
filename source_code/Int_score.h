#pragma once
#include <vector>
#include "TM_score.h"

//==== class: Int_score =====//
class Int_score : public TM_score
{
public:
	Int_score(int num=3000);
	~Int_score(void);

//---- variables ----//
public:
	vector < double > overlap_factor;  //length should be ali_orin
	vector<vector<int> > contact_map; //the contact map which can be based on the results of a linked list

//---- functions ----//
public:
	void Calc_Overlap_Factor_and_Contact_Map(vector < vector <int> > &in1, vector < vector <int> > &int2, 
		vector <pair<int,int> > &alignment ,vector < double > &out);
	double Calc_TM_Score_Single(XYZ *mol1,XYZ *mol2,int lali,double *rotmat_,double d0,double d8,int TM8orTM=0, double *ALLSCO=0);
};
