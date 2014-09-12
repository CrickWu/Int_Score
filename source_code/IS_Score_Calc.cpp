#include <string>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <map>
#include "IS_score.h"
#include "Int_score.h"
#include "iAlign_Output.h"

using namespace std;

// tools for test
template<class T>
void printPairs(const vector<pair<T, T> > &chars) {
	for(unsigned i = 0; i < chars.size(); ++i) {
		cout << "(" <<chars[i].first << "," << chars[i].second << ")";
	}
	cout << endl;
}

template<class T>
void printVector(const vector<T> &chars) {
	for(unsigned i = 0; i < chars.size(); ++i) {
		cout << chars[i] << " ";
	}
	cout << endl;
}

template<class T>
void printVector(const vector<vector<T> > &chars) {
	for(unsigned i = 0; i < chars.size(); ++i) {
		cout << i << " - len " << chars[i].size() << ":";
		printVector<T>(chars[i]);
	}
	cout << endl;
}


void printContactMap(const vector<vector<int> > &map) {
	for (unsigned i = 0; i < map.size(); ++i)
	{
		cout << i << ": ";
		for (unsigned j = 0; j < map[i].size(); ++j)
		{
			cout << map[i][j] << " ";
		}
		cout << endl;
	}
}

//========== process PDB ============//
char WWW_Three2One_III(const char *input)
{
	int i;
	int len;
	int result;
	//encoding
	len=(int)strlen(input);
	if(len!=3)return 'X';
	result=0;
	for(i=0;i<len;i++)result+=(input[i]-'A')*(int)pow(1.0*26,1.0*i);
	//switch
	switch(result)
	{
		case 286:return 'A';
		case 4498:return 'R';
		case 9256:return 'N';
		case 10608:return 'D';
		case 12794:return 'C';
		case 9080:return 'Q';
		case 13812:return 'E';
		case 16516:return 'G';
		case 12383:return 'H';
		case 2998:return 'I';
		case 13635:return 'L';
		case 12803:return 'K';
		case 12960:return 'M';
		case 2901:return 'F';
		case 9921:return 'P';
		case 11614:return 'S';
		case 11693:return 'T';
		case 10601:return 'W';
		case 12135:return 'Y';
		case 7457:return 'V';
		default:return 'X';
	}
}

//--------- PDB_To_XYZ ----------//
int PDB_To_XYZ_NoChain(string &pdb,char * ami,XYZ *xyz)//, std::vector<int> resi_index)
{
	//--- list for mapping ---//
	map<string, int > ws_mapping;
	map<string, int>::iterator iter;
	ws_mapping.clear();
	ifstream fin;
	string buf,temp,name;
	//read
	fin.open(pdb.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"pdb_file %s not found!!\n",pdb.c_str());
		return -1;
	}
	int len;
	int count=0;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		len=(int)buf.length();
		if(len<3)continue;
		//check TER
		temp=buf.substr(0,3);
		if(temp=="END")break;
		//check ATOM
		if(len<4)continue;
		temp=buf.substr(0,4);
		if(temp!="ATOM"&&temp!="HETA")continue;
		//check CA
		temp=buf.substr(13,2);
		if(temp!="CA")continue;
		//record name
		name=buf.substr(21,6);
		iter = ws_mapping.find(name);
		if(iter != ws_mapping.end())continue;
		count++;
		ws_mapping.insert(map < string, int >::value_type(name, count));
		//record ami
		temp=buf.substr(17,3);
		ami[count-1]=WWW_Three2One_III(temp.c_str());
		//record xyz
		temp=buf.substr(30,8);
		xyz[count-1].X=atof(temp.c_str());
		temp=buf.substr(38,8);
		xyz[count-1].Y=atof(temp.c_str());
		temp=buf.substr(46,8);
		xyz[count-1].Z=atof(temp.c_str());

	}
	ami[count]='\0';
	return count;
}

//============== load iAlign interface file ============//
// one example as follows:
/*
Total residue-residue contacts found  :: 20
 ResIndex  NumCont Contact_Residues (starts 1 from the first interface residue)
RES     1   3   40   41   39
RES     2   4   41   39   40   24
RES     3   0
RES     4   0
RES     5   1   25
RES     6   1   41
RES     7   0
RES     8   3   40   38   39
RES     9   3   38   40   39
RES    10   3   38   39   40
RES    11   0
RES    12   1   40
....
RES    35   0
RES    36   0
RES    37   0
RES    38   3    8    9   10
RES    39   5    1    2    8    9   10
RES    40   6    1    2    8    9   10   12
RES    41   3    1    2    6
RES    42   1   13
*/

int Load_iAlign_Interface_File(string &int_file, vector < vector <int> > &output)
{
	ifstream fin;
	string buf,temp;
	//read
	fin.open(int_file.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"interface_file %s not found!!\n",int_file.c_str());
		return -1;
	}
	output.clear();
	//skip
	int found=0;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		istringstream www(buf);
		www>>temp;
		if(temp=="ResIndex")
		{
			found=1;
			break;
		}
	}
	if(found==0)
	{
		fprintf(stderr,"interface_file %s format bad!!\n",int_file.c_str());
		return -1;
	}
	//proc
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		istringstream www(buf);
		vector <int> tmp_rec;
		int num;
		www>> temp >>temp>>num;
		if(num==0)
		{
			output.push_back(tmp_rec);
			continue;
		}
		//load
		int pos;
		for(;;)
		{
			if( ! (www>>pos) )break;
			pos--;
			tmp_rec.push_back(pos);
		}
		output.push_back(tmp_rec);
	}
	//return
	return (int)output.size();
}
//========= modified main process ===========//
void Calc_IS_Score(iAlign_Output & ia)
{
	int maxnum=3000;
	int moln1=ia.coords[0].size();
	int moln2=ia.coords[1].size();
	XYZ *mol1=&(ia.coords[0][0]);
	XYZ *mol2=&(ia.coords[1][0]);
	char *ami1=&(ia.amis[0][0]);
	char *ami2=&(ia.amis[1][0]);
	vector<pair<int, int> > &alignment = ia.alignment;
	//------- suppose we're working on IS-score, then the size of each input should be EXACTLY the same !! -------//
	//-> if the alignment is not specified, generate alignment simply as one-to-one correspondence between two lines
	// and assume them to be of equal length
	int lali = 0;

	XYZ *effect_mol1 = new XYZ[maxnum];
	XYZ *effect_mol2 = new XYZ[maxnum];

	// values are ia.raw_labels[i].first/second
	vector<char> effect_chain_labels1, effect_chain_labels2;
	for(unsigned i = 0; i < ia.raw_labels.size(); ++i) {
		effect_chain_labels1.push_back(ia.raw_labels[i].first);
		effect_chain_labels2.push_back(ia.raw_labels[i].second);
	}

	// TODO may be modified to the second chain's length
	lali = 0;

	// this serves as extracting out the needed chains and 
	for (unsigned i = 0; i < alignment.size(); ++i) {
		if (alignment[i].first > 0 && alignment[i].second > 0) {
			effect_mol1[lali] = mol1[alignment[i].first - 1];
			effect_mol2[lali] = mol2[alignment[i].second - 1];

			lali++;
		}
	}

	vector<vector<int> > &int1 = ia.interfaces[0];
	vector<vector<int> > &int2 = ia.interfaces[1];

/*

	//-> calculate IS-score 
	IS_score is_score(max(moln1, moln2));
	// which one should be the valid normalization number ``L_Q'' (query length)
	is_score.Calc_TM_d0(moln2);
	vector <double> f_score;
	is_score.Calc_Overlap_Factor(int1,int2,alignment,f_score);

	is_score.overlap_factor=f_score;
	
	double isscore=is_score.Calc_TM_Score(effect_mol1,effect_mol2,lali,is_score.d0,is_score.d8,0,0)/moln2;

	//-> rescale IS-score
	// double f0=0.14-0.2*pow(1.0*min(moln1, moln2),-0.3);
	int scale = moln2;
	double f0=0.18-0.35*pow(1.0*scale,-0.3);

	double finscore=(isscore+f0)/(1.0+f0);


//--- test ---//
	printf("length_dep_score=%lf, raw_score=%lf, f0=%lf, d0=%lf, IS-score = %lf , lali = %d \n",isscore*moln2,isscore,f0,is_score.d0,finscore,lali);
//--- test ---//over
*/
	// here starts the scoring for my interface scoring function
	vector<vector<int > > contact_map;
	Int_score int_score(max(moln1, moln2));
	int_score.Calc_TM_d0(moln2);
	int_score.Calc_Overlap_Factor(int1,int2,alignment,int_score.overlap_factor);

	int_score.Calc_Contact_Map(int1, int_score.contact_map[0]);
	int_score.Calc_Contact_Map(int2, int_score.contact_map[1]);

	int_score.Calc_Distance_Matrix(effect_mol1, lali, int_score.distance_matrix[0]);
	int_score.Calc_Distance_Matrix(effect_mol2, lali, int_score.distance_matrix[1]);
	// printContactMap(int1);
	// printContactMap(contact_map); 
	int_score.Calc_BLOSUM_Matrix(ami1, ami2, moln1, moln2, int_score.blos);

	int_score.Calc_Distlap_Factor(mol1, mol2, int1, int2, alignment, int_score.distlap_factor);

	double intscore=int_score.Calc_TM_Score(effect_mol1,effect_mol2,lali,int_score.d0,int_score.d8,0,0)/moln2;
	cout << "intscore: " << intscore << endl;
	cout << "scale: " << moln2 << endl;
	cout << "contact_overlap_ratio:" << int_score.contact_overlap_ratio << endl;
	cout << "dist_overlap_ratio:" << int_score.dist_overlap_ratio << endl;
	cout << "dist_overlap_diff:" << int_score.dist_overlap_diff << endl;


	//--- delete ---//
	
	delete[] effect_mol1;
	delete[] effect_mol2;
}

//========= modified main process ===========//
void Calc_IS_Score(string &pdb1, string &interface1, string &pdb2, string &interface2, vector<pair<int, int> >&alignment)
{
	int i;
	int maxnum=3000;
	int moln1;
	int moln2;
	XYZ *mol1=new XYZ[maxnum];
	XYZ *mol2=new XYZ[maxnum];
	char *ami1=new char[maxnum+1];
	char *ami2=new char[maxnum+1];
	int retv;

	//load protein 1
	//-> load pdb
	retv=PDB_To_XYZ_NoChain(pdb1,ami1,mol1);
	if(retv<=0)exit(-1);
	//-> load interface
	vector < vector <int> > int1;
	moln1=Load_iAlign_Interface_File(interface1, int1);
	if(moln1<=0)exit(-1);
	//-> check length
	if(retv!=moln1)
	{
		fprintf(stderr,"protein_file_1 (%s,%s) -> pdb_size [%d] not equal to interface_size [%d] !!\n",
			pdb1.c_str(),interface1.c_str(),retv,moln1);
		exit(-1);
	}

	//load protein 2
	//-> load pdb
	retv=PDB_To_XYZ_NoChain(pdb2,ami2,mol2);
	if(retv<=0)exit(-1);
	//-> load interface
	vector < vector <int> > int2;
	moln2=Load_iAlign_Interface_File(interface2, int2);
	if(moln2<=0)exit(-1);
	//-> check length
	if(retv!=moln2)
	{
		fprintf(stderr,"protein_file_2 (%s,%s) -> pdb_size [%d] not equal to interface_size [%d] !!\n",
			pdb2.c_str(),interface2.c_str(),retv,moln2);
		exit(-1);
	}

	//------- suppose we're working on IS-score, then the size of each input should be EXACTLY the same !! -------//
	//-> if the alignment is not specified, generate alignment simply as one-to-one correspondence between two lines
	// and assume them to be of equal length
	int lali = 0;

	XYZ *effect_mol1 = new XYZ[maxnum];
	XYZ *effect_mol2 = new XYZ[maxnum];

	if(alignment.size() == 0) {
		lali = moln1;
		alignment.clear();
		for(i=0;i<lali;i++)
			alignment.push_back(pair<int,int>(i+1,i+1));

		for (unsigned i = 0; i < alignment.size(); ++i) {
			effect_mol1[i] = mol1[i];
			effect_mol2[i] = mol2[i];
		}
	}
	else {
		// TODO may be modified to the second chain's length
		lali = 0;

	// this serves as extracting out the needed chains and 
		for (unsigned i = 0; i < alignment.size(); ++i) {
			if (alignment[i].first > 0 && alignment[i].second > 0) {
				effect_mol1[lali] = mol1[alignment[i].first - 1];
				effect_mol2[lali] = mol2[alignment[i].second - 1];
				lali++;

				// cout << "(" << alignment[i].first - 1 << " " << alignment[i].second - 1 << ")";
			}
		}

	}

	//-> calculate IS-score 
	IS_score is_score(max(moln1, moln2));
	// which one should be the valid normalization number ``L_Q'' (query length)
	is_score.Calc_TM_d0(moln2);
	vector <double> f_score;
	is_score.Calc_Overlap_Factor(int1,int2,alignment,f_score);

/*
//--- test ---//
for(int k=0;k<(int)f_score.size();k++)
{
	printf("%lf %d %d \n",f_score[k],int2[k].size(),int1[k].size());
}
//--- test ---//over
*/
	is_score.overlap_factor=f_score;
	
	
	// for(unsigned i = 0; i < 5; ++i) {
	// 	cout << effect_mol2[i].X << " ";
	// }
	// cout << endl;
	
	double isscore=is_score.Calc_TM_Score(effect_mol1,effect_mol2,lali,is_score.d0,is_score.d8,0,0)/moln2;

/*
//--- test ---//
double rotmat[12];
for(int k=0;k<12;k++)rotmat[k]=is_score.finmat[k];
for(int k=0;k<(int)f_score.size();k++)
{
	XYZ xyz;
	is_score.rot_point(mol1[k],xyz,rotmat);
	double dist=xyz.distance(mol2[k]);
	printf("%lf \n",dist);
}
//--- test ---//over
*/

	//-> rescale IS-score
	// double f0=0.14-0.2*pow(1.0*min(moln1, moln2),-0.3);
	int scale = moln2;
	double f0=0.18-0.35*pow(1.0*scale,-0.3);

	double finscore=(isscore+f0)/(1.0+f0);


//--- test ---//
	printf("length_dep_score=%lf, raw_score=%lf, f0=%lf, d0=%lf, IS-score = %lf , lali = %d \n",isscore*moln2,isscore,f0,is_score.d0,finscore,lali);
//--- test ---//over

	// here starts the scoring for my interface scoring function
	vector<vector<int > > contact_map;
	Int_score int_score(max(moln1, moln2));
	int_score.Calc_TM_d0(moln2);
	int_score.Calc_Overlap_Factor(int1,int2,alignment,int_score.overlap_factor);

	int_score.Calc_Contact_Map(int1, int_score.contact_map[0]);
	int_score.Calc_Contact_Map(int2, int_score.contact_map[1]);
	// printContactMap(int1);
	// printContactMap(contact_map); 
	int_score.Calc_BLOSUM_Matrix(ami1, ami2, moln1, moln2, int_score.blos);
	double intscore=int_score.Calc_TM_Score(effect_mol1,effect_mol2,lali,int_score.d0,int_score.d8,0,0)/moln2;
	cout << "intscore: " << intscore << endl;


	//--- delete ---//
	
	delete[] effect_mol1;
	delete[] effect_mol2;

	delete [] mol1;
	delete [] mol2;
	delete [] ami1;
	delete [] ami2;
}
//default setting with no alignment
void Calc_IS_Score(string &pdb1, string &interface1, string &pdb2, string &interface2) {
	vector<pair<int, int> > empty_vector;
	empty_vector.clear();
	Calc_IS_Score(pdb1, interface1, pdb2, interface2, empty_vector);
}

// void Parse_iAlign_File(string filename, std::vector<pair<int, int> > &alignment;) {

// }

//---------- main ----------//
int main(int argc,char **argv)
{
	//---- Calc_IS_Score ----//
	{
		if(argc<6)
		{
			fprintf(stderr,"Calc_IS_Score <pdb1> <interface1> <pdb2> <interface2> <ialign_output_file>\n");
			exit(-1);
		}
		string pdb1=argv[1];
		string interface1=argv[2];
		string pdb2=argv[3];
		string interface2=argv[4];
		string ialign_out = argv[5];

		iAlign_Output ia;
/*
		vector<vector<string> > chain_numbers;
		vector<char> chain_types;
		vector<vector<char> > chain_amis;
		vector<vector<XYZ> > chain_coords;

		vector<pair<string, string> > raw_alignment;
		vector<pair<char, char> > raw_labels;

		ia.parsePDB(pdb1, chain_numbers, chain_types, chain_coords, chain_amis);
		ia.printChainResidues(chain_numbers, chain_types);
		// printVector<XYZ>(chain_coords[0]);
		printVector<char>(chain_amis);
		ia.parsePDB(pdb2, chain_numbers, chain_types, chain_coords, chain_amis);
		ia.printChainResidues(chain_numbers, chain_types);

		ia.parseiAlignRawAlignment(ialign_out, raw_alignment, raw_labels);
		printPairs(raw_alignment), printPairs(raw_labels);
*/
		ia.parseAllFiles(pdb1, interface1, pdb2, interface2, ialign_out);

		// note: that parseFile() function has to be written after parsePDB functions
		// since the alignment need information from pdb files (from Residue numbers to indexes in the contact files)
		
		// ia.parseFile(ialign_out);

		//process

		// this function is the extension of original Calc_IS_Score
		// Calc_IS_Score(pdb1, interface1, pdb2, interface2, ia.alignment);
		Calc_IS_Score(ia);

		// cout << ia.alignment.size() << endl;
		// ia.outputData();

		//exit
		exit(0);
	}
}
