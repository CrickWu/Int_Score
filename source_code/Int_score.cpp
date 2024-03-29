#include "Int_score.h"

//------ constructor -------//
Int_score::Int_score(int num)
:TM_score(num)
{
	// default value shows that we compare between two proteins
	contact_map.resize(2);
	distance_matrix.resize(2);
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

	int total_overlap = 0;
	// molx_contact = sum(len(inx[i]))
	int mol1_contact = 0, mol2_contact = 0;
	for (i = 0; i < moln1; i++)
		mol1_contact += in1[i].size();
	for (i = 0; i < moln2; i++)
		mol2_contact += in2[i].size();

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
			total_overlap += col;
			double fcol=0;
			if( size1>0 && size2>0 ) fcol=(1.0*col/size1 + 1.0*col/size2)/2.0;
			out.push_back(fcol);
		}
	}

	// calculate contact overlap ratio
	contact_overlap_ratio = (double) total_overlap / min(mol1_contact, mol2_contact);
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

// calculate the distance matrix
void Int_score::Calc_Distance_Matrix(XYZ const *coords, int size, vector<vector<double> > &dist) {
	// initialize dist
	dist.resize(size);
	for(int i = 0; i < size; ++i) {
		dist[i].resize(size);
		dist[i].assign(size, -1);
	}

	// calculate distance matrix from the coords
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			dist[i][j] = dist[j][i] = coords[i].distance(coords[j]);
		}
	}
}

// calculate the distance-related factor
// and the function should be called after Calc_Distance_Matrix (since it needs to utilize the information)
void Int_score::Calc_Distlap_Factor(XYZ const *coords1, XYZ const *coords2//, const vector<char> &chain_labels1, const vector<char> &chain_labels2
	, const vector < vector <int> > &in1, const vector < vector <int> > &in2, const vector<pair<int, int> >&alignment, vector<double> &distlap_factor) {
	int size = alignment.size();

	// calculate the factor, the function can be of two forms
	// 1. 1/1+(d_1-d_2)^2/c^2
	// 2. inner product

	// method 1. with only contact residues taken into consideration


	// note that coords has to be in agreement with ins (about the size)
	int i,j,k;
	int ii,jj;
	int moln1=in1.size();
	int moln2=in2.size();
	vector <int> ali1 (moln1,-1);
	vector <int> ali2 (moln2,-2);

	double total_overlap = 0.0;
	double mol1_diff = 0.0, mol2_diff = 0.0;
	// molx_contact = sum(len(inx[i]))
	int mol1_contact = 0, mol2_contact = 0;
	for (i = 0; i < moln1; i++)
		mol1_contact += in1[i].size();
	for (i = 0; i < moln2; i++)
		mol2_contact += in2[i].size();

	// first calculate the distant matrix
	vector<vector<double> > dist1, dist2;
	Calc_Distance_Matrix(coords1, moln1, dist1);
	Calc_Distance_Matrix(coords2, moln2, dist2);

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
	// initialize distlap_factor
	// cout << "size: " << size << endl;
	distlap_factor.clear();
	for(i=0;i<size;i++)
	{
		// cout << i << ": ";
		ii=alignment[i].first;
		jj=alignment[i].second;
		if(ii>0 && jj>0)
		{
			//init
			int size1=(int)in1[ii-1].size();
			int size2=(int)in2[jj-1].size();
			//calc
			double col=0;
			for(j=0;j<size1;j++)
			{
				int pos1=in1[ii-1][j];
				for(k=0;k<size2;k++)
				{
					int pos2=in2[jj-1][k];
					if(ali1[pos1]==ali2[pos2])
					{
						// distance difference
						double diff = dist1[ii-1][pos1] - dist2[jj-1][pos2];
						// or d0
						double denom = 4.0;
						// double denom = d0;
						// the score is changed
						col += 1.0/(1.0 + diff * diff / denom);
						break;
					}
				}
			}
			//assign
			double fcol=0;
			double denom1 = 0.0, denom2 = 0.0;
			for (j = 0; j < size1; j++) {
				// cout << dist1[ii-1][ in1[ii-1][j] ] << "--\n";
				denom1 += 1.0 / (1.0 + (dist1[ii-1][ in1[ii-1][j] ] - 3) * (dist1[ii-1][ in1[ii-1][j] ] - 3) / 2);
			}
			// cout << size2 << endl;
			for (j = 0; j < size2; j++) {
				denom2 += 1.0 / (1.0 + (dist2[jj-1][ in2[jj-1][j] ] - 3) * (dist2[jj-1][ in2[jj-1][j] ] - 3) / 2);
			}
			total_overlap += col;
			mol1_diff += denom1;
			mol2_diff += denom2;
			if( size1>0 && size2>0 ) fcol=(1.0*col/denom1 + 1.0*col/denom2)/2.0;
			// cout << fcol<<" "<< col << " " << denom1 << " " << denom2 << endl;
			// if( size1>0 && size2>0 ) fcol=1.0 * col / sqrt(1.0 * denom2 * denom1);
			distlap_factor.push_back(fcol);
		}
	}

	// calculate contact overlap ratio
	dist_overlap_ratio = total_overlap / min(mol1_contact, mol2_contact);
	dist_overlap_diff = total_overlap / min(mol1_diff, mol2_diff);
	// cout << "mol_diff:" << total_overlap << " " << mol1_diff << " " << mol2_diff << endl;
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
			// score+=blos[i][j]*contact_map[0][i][j]*contact_map[1][i][j]/(1.0+dist2i/ori_d) / (1.0+dist2j/ori_d);
			// score+=blos[i][j] / (1.0+(distance_matrix[0][i][j]-distance_matrix[1][i][j])*(distance_matrix[0][i][j]-distance_matrix[1][i][j])/4.0)
			score+=blos[i][j] * overlap_factor[i] * overlap_factor[j]
			// score+=blos[i][j] * distlap_factor[i] * distlap_factor[j]
			/(1.0+dist2i/ori_d) / (1.0+dist2j/ori_d); // norm factos is 4.0
		}
	}
	return score;
}
