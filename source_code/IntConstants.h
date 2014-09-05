#pragma once
#include <iostream>

using namespace std;

class IntConstants {
//variables
public:
	static const int AMI_NUM = 21;

	static const int BLOSUM62[AMI_NUM][AMI_NUM];
	static const int AA_Map[26];	
	static const int Blo_AA_Map[21];
//methods
public:
	//------ calculate -------//
	static int BLOSUM62_Calc(char a,char b);
};
