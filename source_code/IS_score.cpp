#include "IS_score.h"

//------ constructor -------//
IS_score::IS_score(int num)
:TM_score(num)
{
}
IS_score::~IS_score(void)
{
}


//--------- Calc_Overlap_Factor -------//
//-> note: the size of in1 and in2 could NOT be identical !!
//   actually, they are moln1 and moln2, according to the alignment.
void IS_score::Calc_Overlap_Factor(vector < vector <int> > &in1, vector < vector <int> > &in2, 
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

//--------- Calc_TM_Score_Single -------//
double IS_score::Calc_TM_Score_Single(XYZ *mol1,XYZ *mol2,int lali,double *rotmat_,double d0,double d8,
	int TM8orTM,double *ALLSCO)
{
	int i;
	XYZ temp;
	double dist2;
	double ori_d;
	ori_d=d0*d0;
	//calculate
	double score=0.0;
	double factor;
	for(i=0;i<lali;i++)
	{
		TMs_Cache_Point(mol1,i,temp,rotmat_);
		dist2=mol2[i].distance_square(temp);
		factor=overlap_factor[i];
		score+=factor/(1.0+dist2/ori_d);
	}
	return score;
}
