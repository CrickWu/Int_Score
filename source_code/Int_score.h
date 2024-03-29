#pragma once
#include <vector>
#include "TM_score.h"
#include "IntConstants.h"

//==== class: Int_score =====//
class Int_score : public TM_score
{
public:
	Int_score(int num=3000);
	~Int_score(void);

//---- variables ----//
public:
	double dist_overlap_diff; // overlapped dist difference / min(total_contacts[x])
	double dist_overlap_ratio; // overlapped dist difference / min(total_contacts[x])
	double contact_overlap_ratio; // overlapped contacts / min(total_contacts[x])

	vector < double > overlap_factor;  //length should be ali_orin
	vector<vector<vector<int> > > contact_map; //the contact map which can be based on the results of a linked list: contact_map[0], contact_map[1] is the one for mol1 and mol2 respectively
	vector<vector<double> > blos; // the BLOSUM matrix related items for multiplying
	vector<double> distlap_factor; // the distance-related matrix which provides a more careful mutiplier
	vector<vector<vector<double> > > distance_matrix; // the distance matrix which denotes the distances between residues within a protein

//---- functions ----//
public:
	void Calc_Distlap_Factor(XYZ const *coords1, XYZ const *coords2//, const vector<char> &chain_labels1, const vector<char> &chain_labels2
		, const vector < vector <int> > &in1, const vector < vector <int> > &in2, const vector<pair<int, int> >&alignment, vector<double> &distlap_factor);
	void Calc_BLOSUM_Matrix(char const* ami1, char const* ami2, int moln1, int moln2, vector<vector<double> > &blos);
	void Calc_Contact_Map(const vector<vector<int> > &in, vector<vector<int> > &map);
	void Calc_Overlap_Factor(vector < vector <int> > &in1, vector < vector <int> > &int2, 
		vector <pair<int,int> > &alignment ,vector < double > &out);
	void Calc_Distance_Matrix(XYZ const* coords, int size, vector<vector<double> > &dist);
	double Calc_TM_Score_Single(XYZ *mol1,XYZ *mol2,int lali,double *rotmat_,double d0,double d8,int TM8orTM=0, double *ALLSCO=0);
};
