#include "iAlign_Output.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>

using namespace std;

int iAlign_Output::size() {

	return max(residue_number[0].size(), residue_number[1].size());
}

bool iAlign_Output::is_number(const std::string& s) {
    std::string::const_iterator it = s.begin();
    while (it != s.end() && std::isdigit(*it)) ++it;
    return !s.empty() && it == s.end();
}

void iAlign_Output::parseFile(const string &filename) {
	ifstream fin;
	fin.open(filename.c_str(), ios::in);

	if (!fin.is_open()) {
		cout << "Error opening file." << endl;
		return;
	}

	string line;
	while (! fin.eof()) {
		getline(fin, line);
		// starting point of the alignment
		// note that the result is unsigned in which we can not use >= 0 comparisons
		if (line.find("Index Ch1") != string::npos) {
			// cout << line << " " << line.find("Index Ch1") << endl;
			break;
		}
	}

	stringstream str;
	string temp;
	// int index = 0;
	string index_string;


	while(! fin.eof()) {
		getline(fin, line);
		// first end point of the alignment if the length is way too few to parse
		if (line.length() < 2)
			break;

		str.str(string());
		str << line;
		char c0, c1;
		int res0, res1;
		str >> index_string;
		// end point of the alignment
		// if the first entry is not a number
		if (!is_number(index_string))
			break;
		// index = atoi(index_string.c_str());
		str >> c0 >> res0 >> temp >> c1 >> res1;
		// cout << index << " " << c0 << " " << " " << res0 << " " << " " << temp << " " << c1 << " " << res1 << endl;

		// use the pdb files' output to determine the mapping
		int index0 = number_to_index_map[0][pair<char, int>(c0, res0)];
		int index1 = number_to_index_map[1][pair<char, int>(c1, res1)];

		alignment.push_back(make_pair(index0 + 1, index1 + 1));
	}
	fin.close();

//	cout << "alignment without inserting gaps: ";
//	printAlignment(alignment);

	//-- fill the alignment with gaps, the second chain would be gapped ealier than the first chain --
	vector<pair<int, int> > tmp_alignment = alignment;

	alignment.insert(alignment.begin(), make_pair(0, 0));
	// ``+1'' since the alignment begins from 1 ends with ``size()''
	alignment.push_back(make_pair(residue_number[0].size() + 1, residue_number[1].size() + 1));
	// generate vectors to insert at each intervel
	vector<vector<pair<int, int> > > insert_vector(alignment.size() - 1);
	for (unsigned i = 0; i < alignment.size() - 1; ++i) {
		int start[2] = {alignment[i].first, alignment[i].second};
		int end[2] = {alignment[i + 1].first, alignment[i + 1].second};

		for (int j = start[0] + 1; j <= end[0] - 1; ++j) {
			insert_vector[i].push_back(make_pair( j, -start[1] ));
		}

		for (int j = start[1] + 1; j <= end[1] - 1; ++j) {
			insert_vector[i].push_back(make_pair( -(end[0] - 1), j ));
		}
	}

	// combine the intervel vectors and original vectors together
	alignment.clear();
	for (unsigned i = 0; i < tmp_alignment.size(); i++) {
		alignment.insert(alignment.end(), insert_vector[i].begin(), insert_vector[i].end());
		alignment.push_back(tmp_alignment[i]);
	}
	alignment.insert(alignment.end(), insert_vector[tmp_alignment.size()].begin(), insert_vector[tmp_alignment.size()].end());
	//-- filling completes --
}

int iAlign_Output::parsePDB(const string &filename, const int complex_number) {
		//--- list for mapping ---//
	map<string, int > ws_mapping;
	map<string, int>::iterator iter;
	ws_mapping.clear();
	ifstream fin;
	string buf,temp,name;
	//read
	fin.open(filename.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"pdb_file %s not found!!\n", filename.c_str());
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
		ws_mapping.insert(map < string, int >::value_type(name, count));
		
		char tmp_chain = buf[21];
		string tmp_str = buf.substr(22, 5);
		int tmp_number = atoi(tmp_str.c_str());

		residue_number[complex_number].push_back(tmp_number);
		residue_chain[complex_number].push_back(tmp_chain);

		// cout << complex_number << " " << tmp_chain << " " << tmp_number << endl;

		number_to_index_map[complex_number][pair<char, int>(tmp_chain, tmp_number)] = count;
		count++;
	}
	fin.close();

	return count;
}

void iAlign_Output::printAlignment(vector<pair<int, int> > alignment) {
	int size = alignment.size();

	cout << "Size: " << size << "\t";
	for (int i = 0; i < size; ++i) {
		cout << "(" << alignment[i].first << "," << alignment[i].second << ")" << " ";
	}
	cout << endl;

}

void iAlign_Output::outputData() {
	for (int j = 0; j < 2; j++) {
		int size = residue_number[j].size();

		cout << "Size: " << size << "\t";
		for (int i = 0; i < size; ++i) {
			cout << residue_number[j][i] << " ";
		}	
		cout << endl;
	}

/*
	int size = alignment.size();

	cout << "Size: " << size << "\t";
	for (int i = 0; i < size; ++i) {
		cout << "(" << alignment[i].first << "," << alignment[i].second << ")" << " ";
	}
	cout << endl;
*/
}

iAlign_Output::iAlign_Output() {
	residue_chain.clear();
	residue_chain.resize(2);

	residue_number.clear();
	residue_number.resize(2);

	number_to_index_map.clear();
	number_to_index_map.resize(2);
}

iAlign_Output::~iAlign_Output() {

}