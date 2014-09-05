#include "Int_score.h"

//------ constructor -------//
Int_score::Int_score(int num)
:TM_score(num)
{
	// default value shows that the complex a dimer
	contact_map.resize(2);
}
Int_score::~Int_score(void)
{
}


//--------- Calc_Overlap_Factor -------//
//-> note: the size of in1 and in2 could NOT be identical !!
//   actually, they are moln1 and moln2, according to the alignment.
void Int_score::Calc_Overlap_Factor(vector < vector <int> > &in1, vector < vector <int> > &in2, 
		vector <pair<int,int> > &alignment ,vector < double > &out)
{
	int i,j,k;
	int ii,jj;
	int size=(int)alignment.size();
	int moln1=in1.size();
	int moln2=in2.size();
	vector <int> ali1 (moln1,-1);
	vector <int> ali2 (moln2,-2);

	//assign alignment
	for(i=0;i<size;i++)
	{
		ii=alignment[i].first;
		jj=alignment[i].second;
		if(ii>0 && jj>0 && ii<=moln1 && jj<=moln2)
		{
			ali1[ii-1]=i;
			ali2[jj-1]=i;
		}
	}

	//calculate overlap
	out.clear();
	for(i=0;i<size;i++)
	{
		ii=alignment[i].first;
		jj=alignment[i].second;
		if(ii>0 && jj>0)
		{
			//init
			int size1=(int)in1[ii-1].size();
			int size2=(int)in2[jj-1].size();
			//calc
			int col=0;
			for(j=0;j<size1;j++)
			{
				int pos1=in1[ii-1][j];
				for(k=0;k<size2;k++)
				{
					int pos2=in2[jj-1][k];
					if(ali1[pos1]==ali2[pos2])
					{
						col++;
						break;
					}
				}
			}
			//assign
			double fcol=0;
			if( size1>0 && size2>0 ) fcol=(1.0*col/size1 + 1.0*col/size2)/2.0;
			out.push_back(fcol);
		}
	}
}

// calculate the contact map
void Int_score::Calc_Contact_Map(const vector < vector <int> > &in, vector<vector<int> > &map) {
	int size = in.size();
	//map initialization
	map.resize(size);
	for (int i = 0; i < size; i++) {
		map[i].resize(size);
		map[i].assign(size, 0);
	}

	// calculate contact map from the linked list
	for (int i = 0; i < size; i++) {
		for (unsigned j = 0; j < in[i].size(); j++) {
			int index[2] = {i, in[i][j]};
			map[index[0]][index[1]] = 1;
			map[index[1]][index[0]] = 1;
		}
	}
}

void Int_score::Calc_BLOSUM_Matrix(char const * ami1, char const* ami2, int moln1, int moln2, vector<vector<double> > &blos) {
// initialize the blosum-related matrix
	blos.resize(moln1);
	for (int i = 0; i < moln1; ++i) {
		blos[i].resize(moln2);
		blos[i].assign(moln2, 0);
	}

	for (int i = 0; i < moln1; ++i) {
		for (int j = 0; j < moln2; ++j) {
			blos[i][j] = max(1.0, IntConstants::BLOSUM62_Calc(ami1[i], ami2[j]) / 5.0);
		}
	}
}

//--------- Calc_TM_Score_Single -------//
double Int_score::Calc_TM_Score_Single(XYZ *mol1,XYZ *mol2,int lali,double *rotmat_,double d0,double d8,
	int TM8orTM,double *ALLSCO)
{
	XYZ tempi, tempj;
	double dist2i, dist2j;
	double ori_d;
	ori_d=d0*d0;
	//calculate
	double score=0.0;

	for (int j = 0; j < lali; ++j)
	{
		TMs_Cache_Point(mol1,j,tempj,rotmat_);
		dist2j=mol2[j].distance_square(tempj);

		for(int i=0;i<lali;i++)
		{
			TMs_Cache_Point(mol1,i,tempi,rotmat_);
			dist2i=mol2[i].distance_square(tempi);
			// score+=overlap_factor[i]/(1.0+dist2i/ori_d) / (1.0+dist2j/ori_d);
			score+=blos[i][j]*contact_map[0][i][j]*contact_map[1][i][j]/(1.0+dist2i/ori_d) / (1.0+dist2j/ori_d);
		}
	}
	return score;
}
