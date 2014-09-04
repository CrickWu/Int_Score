#pragma once
#include <map>
#include <vector>
#include <string>
using namespace std;

class iAlign_Output {
//functions
public:
	iAlign_Output();
	~iAlign_Output();
	// this function can only be called after the two PDB files have been parsed since it translates the original alignment into the region within the interface (which uses new sequential numbers for indexing)
	void parseFile(const string &filename);
	// complex_number is the number for which complex the results should be sent
	int parsePDB(const string &filename, const int complex_number);
	// complex_number is the number for which complex should the results be sent
	int parseInterface(const string &filename, const int complex_number);
	int size();
	// for test
	void outputData();
	void printAlignment(vector<pair<int, int> > alignment);
	bool is_number(const std::string& s);
// variables
public:
	// converted during the phase of reading parseFile
	vector<pair<int, int> > alignment;

private:
	// map from residue number to 
	vector<map<pair<char, int>, int> > number_to_index_map;
	// residues of each for comparision, which is read from PDB files
	vector<vector<int> > residue_number;
	vector<vector<char> > residue_chain;
};