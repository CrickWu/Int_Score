#pragma once
#include <map>
#include <vector>
#include <string>

#include "XYZ.h"

using namespace std;

class iAlign_Output {
//functions
public:
	iAlign_Output();
	~iAlign_Output();
	// parse the iAlign output file for the label orders of two proteins' chains, label[0], and label[1]s
	static void parseiAlignRawAlignment(const string &filename, vector<pair<string, string> > &raw_alignment,
		vector<pair<char, char> > &raw_labels);
	// parse the results as a multimer
	static int parsePDB(const string &filename, vector<vector<string> > &chain_numbers, vector<char> &chain_types
		, vector<vector<XYZ> >&chain_coords, vector<vector<char> >&chain_amis);
	// realign the interface list according to the iAlign's labels (note the labels are all compressed)
	// residue_number are for the convinience of counting the number of residues in a label (assoc with ori_labels)
	static void realignInterface(const vector<vector<int > >& raw_interface, const vector<char> &ori_labels, 
		const vector<char> &ialign_labels, const vector<vector<string> > &residue_number, vector<vector<int> >& interface);
	// fill the gaps in the (raw) alignment
	static void fillGapInAlignment(vector<pair<int, int> >&alignment, int residue_number_size1, int residue_size_number2);
	// parse the interface file as a list, return the size of the interface
	static int parseInterface(const string &interface_file, vector<vector<int> > &intList);
	// this function can only be called after the two PDB files have been parsed since it translates the original alignment into the region within the interface (which uses new sequential numbers for indexing)
	void parseFile(const string &filename);
	// complex_number is the number for which complex the results should be sent
	int parsePDB(const string &filename, const int complex_number);
	// parse the entire files
	void parseAllFiles(const string &pdb1, const string &interface1, 
		const string &pdb2, const string &interface2, const string &ialign_out);
	// Fill in the map based on sep_residue_number, used after gaining the raw_alignment and raw_label
	// contact map
	void generateFullAlignment();
	/*const vector<pair<string, string> > &raw_alignment, 
		const vector<pair<string, string> > &raw_labels, // the paras above are for alignment
		const vector<vector<XYZ> > &sep_residue_coord, 
		const vector<vector<vector<string> > > &sep_residue_number, 
		const vector<vector<char> > &sep_reisdue_label, 
		const vector<vector<char> > &sep_reisdue_ami, // the paras above are for raw separate chains' inputs
		vector<vector<XYZ> > &protein_coords,
		vector<vector<char> > &protein_amis,
		);
		*/
	int size();
	// for test
	void outputData();
	static void printChainResidues(const vector<vector<string> >& chain_numbers, const vector<char> &chain_types);
	static void printAlignment(const vector<pair<int, int> >& alignment);
	// tools
	static bool is_number(const std::string& s);
	static char WWW_Three2One_III(const string &input);
// variables
public:
	// converted during the phase of reading parseFile
	vector<pair<int, int> > alignment;
	// aminos returned by joining all proteins
	vector<vector<char> > amis;
	// residue coordinates returned by joining all proteins
	vector<vector<XYZ> > coords;
	// the realigned interfaces based on the results of iAlign result file
	vector< vector<vector<int> > > interfaces;

private:
	// map from residue number to 
	vector<map<pair<char, string>, int> > number_to_index_map;
	// the final residues joined by multiple chains and the results from iAlign result file
	vector<vector<string> > residue_number;
	vector<vector<char> > residue_label;
	// the parsed linked results of interfaces
	vector< vector<vector<int> > > raw_interfaces;
	// the results of directly parsing the alignment in iAlign result files
	vector<pair<string, string> > raw_alignment;
	vector<pair<char, char> > raw_labels;

	//separate reisidue_chains 
	// for vector<vector<vector<T> > > **[0][0][*] **[0][1][*] ... are the chains for dimer[0]
	// for vector<vector<T> > **[0][0] **[0][1] ... are the chains for dimer[0]
	vector< vector<vector<char> > > sep_amis;
	vector< vector<vector<XYZ> > > sep_coords; 
	vector< vector<vector<string> > > sep_residue_number;
	vector< vector<char> > sep_residue_label; // e.g. sep_reisdue_label[0][0] sep_reisdue_label[0][1] are the markers of which chains is it
};