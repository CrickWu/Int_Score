#include "iAlign_Output.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>

using namespace std;

// tools for test
// remember to delete otherwise causing double definitions?
template<class T>
void printPairs(const vector<pair<T, T> > &chars) {
	for(unsigned i = 0; i < chars.size(); ++i) {
		cout << "(" <<chars[i].first << "," << chars[i].second << ")";
	}
	cout << endl;
}

template<class T>
void printVector(const vector<T> &chars) {
	cout << "len " << chars.size() << ":";
	for(unsigned i = 0; i < chars.size(); ++i) {
		cout << " "<< chars[i];
	}
	cout << endl;
}

template<class T>
void printVector(const vector<vector<T> > &chars) {
	for(unsigned i = 0; i < chars.size(); ++i) {
		cout << i << ":";
		printVector<T>(chars[i]);
	}
	cout << endl;
}

//------- tools -------
//========== process PDB ============//
char iAlign_Output::WWW_Three2One_III(const string &input)
{
	int i;
	int len;
	int result;
	//encoding
	len=input.length();
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

int iAlign_Output::size() {

	return max(residue_number[0].size(), residue_number[1].size());
}

bool iAlign_Output::is_number(const std::string& s) {
    std::string::const_iterator it = s.begin();
    while (it != s.end() && std::isdigit(*it)) ++it;
    return !s.empty() && it == s.end();
}

void iAlign_Output::parseiAlignRawAlignment(const string &filename, vector<pair<string, string> > &raw_alignment, vector<pair<char, char> > &raw_labels) {
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
			break;
		}
	}

	stringstream str;
	string temp;
	string index_string;

	raw_alignment.clear();

	while(! fin.eof()) {
		getline(fin, line);
		// first end point of the alignment if the length is way too few to parse
		if (line.length() < 2)
			break;

		str.str(string());
		str << line;
		char c0, c1;
		string res0, res1;
		str >> index_string;
		// end point of the alignment
		// if the first entry is not a number
		if (!is_number(index_string))
			break;
		str >> c0 >> res0 >> temp >> c1 >> res1;

		raw_alignment.push_back(make_pair(res0, res1));
		raw_labels.push_back(make_pair(c0, c1));
	}
	fin.close();
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
		string res0, res1;
		str >> index_string;
		// end point of the alignment
		// if the first entry is not a number
		if (!is_number(index_string))
			break;
		// index = atoi(index_string.c_str());
		str >> c0 >> res0 >> temp >> c1 >> res1;
		// cout << index << " " << c0 << " " << " " << res0 << " " << " " << temp << " " << c1 << " " << res1 << endl;

		// use the pdb files' output to determine the mapping
		int index0 = number_to_index_map[0][pair<char, string>(c0, res0)];
		int index1 = number_to_index_map[1][pair<char, string>(c1, res1)];

		alignment.push_back(make_pair(index0 + 1, index1 + 1));
	}
	fin.close();

	cout << "alignment without inserting gaps: ";
	printAlignment(alignment);
	outputData();


}

// note that this function can read in multiple chains (>= 2)
int iAlign_Output::parsePDB(const string &filename, vector<vector<string> > &chain_number, vector<char> &chain_type
	, vector<vector<XYZ> >&chain_coord, vector<vector<char> >&chain_ami) {
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

	unsigned which_chain = 0;
	chain_type.clear();
	chain_number.clear();
	chain_coord.clear();
	chain_ami.clear();

	vector<string> empty_vec_str;
	vector<XYZ> empty_vec_XYZ;
	vector<char> empty_vec_char;

	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		len=(int)buf.length();
		if(len<3)continue;
		//check TER
		temp=buf.substr(0,3);
		if(temp=="END")break;
		if(temp=="TER") which_chain++;
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

		// record chain labels
		char tmp_chain = buf[21];
		// when encountered "TER" append the new chain
		if (which_chain == chain_type.size()) {
			chain_type.push_back(tmp_chain);
			chain_number.push_back(empty_vec_str); // append a new empty chain
			chain_ami.push_back(empty_vec_char);
			chain_coord.push_back(empty_vec_XYZ);
		}

		// record residue number
		string tmp_str = buf.substr(22, 6);
		// trim the white spaces TODO but not considering all whitespaces
		tmp_str = tmp_str.substr(tmp_str.find_first_not_of(" "));
		tmp_str = tmp_str.substr(0, tmp_str.find_last_not_of(" ")+1);

		chain_number[which_chain].push_back(tmp_str);

		// record ami
		temp=buf.substr(17,3);
		chain_ami[which_chain].push_back( WWW_Three2One_III(temp) );
		// record xyz
		XYZ xyz;
		temp=buf.substr(30,8);
		xyz.X=atof(temp.c_str());
		temp=buf.substr(38,8);
		xyz.Y=atof(temp.c_str());
		temp=buf.substr(46,8);
		xyz.Z=atof(temp.c_str());
		chain_coord[which_chain].push_back(xyz);

		count++;
	}
	fin.close();

	// total number of all chains
	return count;
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
		string tmp_str = buf.substr(22, 6);
		// trim the white spaces TODO but not considering all whitespaces
		tmp_str = tmp_str.substr(tmp_str.find_first_not_of(" "));
		tmp_str = tmp_str.substr(0, tmp_str.find_last_not_of(" ")+1);

		residue_number[complex_number].push_back(tmp_str);
		residue_label[complex_number].push_back(tmp_chain);

		number_to_index_map[complex_number][pair<char, string>(tmp_chain, tmp_str)] = count;
		count++;
	}
	fin.close();

	return count;
}
void iAlign_Output::fillGapInAlignment(vector<pair<int, int> >&alignment, int residue_number_size1, int residue_number_size2) {
	//-- fill the alignment with gaps, the second chain would be gapped ealier than the first chain --
	vector<pair<int, int> > tmp_alignment = alignment;

	alignment.insert(alignment.begin(), make_pair(0, 0));
	// ``+1'' since the alignment begins from 1 ends with ``size()''
	alignment.push_back(make_pair(residue_number_size1 + 1, residue_number_size2 + 1));
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
void iAlign_Output::generateFullAlignment() {
	// labels which show the sequence of appearing chains
	vector<vector<char> > compressed_labels;
	compressed_labels.resize(2);
	pair<char, char> pred = make_pair('0', '0');
	for(unsigned i = 0; i < raw_labels.size(); ++i) {
		if (raw_labels[i] != pred) {
			pred = raw_labels[i];
			compressed_labels[0].push_back(raw_labels[i].first);
			compressed_labels[1].push_back(raw_labels[i].second);
		}
	}

	// join chains to a whole string based on the order of compressed_labels
	residue_number.clear();
	residue_label.clear();
	amis.clear();
	coords.clear();
	for(unsigned i = 0; i < 2; ++i) {
		for(unsigned j = 0; j < compressed_labels[i].size(); ++j) {
			// find the proper, join when sep_residue_label[i(0 or 1)][?] == compressed_labels[0 or 1][j]
			for(unsigned k = 0; k < sep_residue_label[i].size(); ++k) {
				// join
				if (sep_residue_label[i][k] == compressed_labels[i][j]) {
					residue_number[i].insert(residue_number[i].end(), sep_residue_number[i][k].begin(), sep_residue_number[i][k].end());
					residue_label[i].insert(residue_label[i].end(), sep_residue_number[i][k].size(), sep_residue_label[i][k]);
					amis[i].insert(amis[i].end(), sep_amis[i][k].begin(), sep_amis[i][k].end());
					coords[i].insert(coords[i].end(), sep_coords[i][k].begin(), sep_coords[i][k].end());

					break;
				}
			}
		}
	}

	// set up the (number, label)->index map
	for(unsigned i = 0; i < 2; ++i) {
		for(unsigned j = 0; j < residue_number[i].size(); ++j) {
			number_to_index_map[i][make_pair(residue_label[i][j], residue_number[i][j])] = j;
		}
	}

	// compute the alignment based on the generated map
	for(unsigned i = 0; i < raw_alignment.size(); ++i) {
		int index0 = number_to_index_map[0][make_pair(raw_labels[i].first, raw_alignment[i].first)];
		int index1 = number_to_index_map[1][make_pair(raw_labels[i].second, raw_alignment[i].second)];

		alignment.push_back(make_pair(index0 + 1, index1 + 1));
	}

	// cout << "alignment without inserting gaps: ";
	// printAlignment(alignment);

	fillGapInAlignment(alignment, residue_number[0].size(), residue_number[1].size());

	// realign the interface list
	interfaces.resize(2);
	realignInterface(raw_interfaces[0], sep_residue_label[0], compressed_labels[0], sep_residue_number[0], interfaces[0]);
	realignInterface(raw_interfaces[1], sep_residue_label[1], compressed_labels[1], sep_residue_number[1], interfaces[1]);
}

void iAlign_Output::realignInterface(const vector<vector<int > >& raw_interface, const vector<char> &ori_labels, 
	const vector<char> &ialign_labels, const vector<vector<string> > &residue_number, vector<vector<int> >& interface) {
	int total_size = raw_interface.size();
	// vectors recording the size of each chain
	vector<int> residue_number_size(residue_number.size(), 0);
	for(unsigned i = 0; i < residue_number.size(); ++i) {
		residue_number_size[i] = residue_number[i].size();
	}
	// vector each denoting the aggregated size
	vector<int> aggregated_residue_number_size(residue_number.size(), 0);
	for(unsigned i = 1; i < residue_number.size(); ++i) {
		aggregated_residue_number_size[i] = aggregated_residue_number_size[i-1] + residue_number_size[i-1];
	}
	// construct the mapping
	vector<int> index_map(total_size, 0);

	int current_filling_number = 0;
	for(unsigned i = 0; i < residue_number.size(); ++i) {
		// find the original order of ialign_labels
		// index of j
		int ori_index = 0;
		for(unsigned j = 0; j < residue_number.size(); ++j) {
			if (ialign_labels[i] == ori_labels[j]) {
				ori_index = j;
				break;
			}
		}
		// number of elements before the chain in the original labels
		int elements_before = aggregated_residue_number_size[ori_index];
		for(unsigned j = 0; j < residue_number[ori_index].size(); ++j) {
			index_map[elements_before + j] = current_filling_number;
			current_filling_number++;
		}
	}

	// generate invert mapping
	vector<int> invert_map(total_size, 0);
	for(int i = 0; i < total_size; ++i) {
		invert_map[ index_map[i] ] = i;
	}

	printVector(index_map);
	printVector(invert_map);
	
	// generate the new interface
	interface.resize(total_size);
	for(int i = 0; i < total_size; ++i) {
		// interface[i]  invert_map[i]'s entries, convert with index_map
		int ori_index = invert_map[i];
		interface[i].clear();
		for(unsigned j = 0; j < raw_interface[ori_index].size(); ++j) {
			interface[i].push_back(index_map[ raw_interface[ori_index][j] ]);
		}
	}

	printVector(interface);
}


int iAlign_Output::parseInterface(const string &interface_file, vector<vector<int> > &intList) {
	ifstream fin;
	string buf,temp;
	//read
	fin.open(interface_file.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"interface_file %s not found!!\n",interface_file.c_str());
		return -1;
	}
	intList.clear();
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
		fprintf(stderr,"interface_file %s format bad!!\n",interface_file.c_str());
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
			intList.push_back(tmp_rec);
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
		intList.push_back(tmp_rec);
	}
	//return
	return (int)intList.size();	
}

void iAlign_Output::parseAllFiles(const string &pdb1, const string &interface1,
	const string &pdb2, const string &interface2, const string &ialign_out) {

	parsePDB(pdb1, sep_residue_number[0], sep_residue_label[0], sep_coords[0], sep_amis[0]);
	// printChainResidues(sep_residue_number[0], sep_residue_label[0]);
	// printVector(sep_amis[0][0]);printVector(sep_amis[0][1]);

	parsePDB(pdb2, sep_residue_number[1], sep_residue_label[1], sep_coords[1], sep_amis[1]);
	// printChainResidues(sep_residue_number[1], sep_residue_label[1]);

	parseiAlignRawAlignment(ialign_out, raw_alignment, raw_labels);
	// printPairs(raw_alignment), printPairs(raw_labels);

	parseInterface(interface1, raw_interfaces[0]);
	parseInterface(interface2, raw_interfaces[1]);

	generateFullAlignment();
}

void iAlign_Output::printChainResidues(const vector<vector<string> > &chain_numbers, const vector<char> &chain_types) {
	int size = chain_numbers.size();
	for (int i = 0; i < size; i++) {
		cout << chain_types[i] << " - len " << chain_numbers[i].size() << ":";
		for(unsigned j = 0; j < chain_numbers[i].size(); ++j) {
			cout << " " << chain_numbers[i][j];
		}
		cout << endl;
	}
}

void iAlign_Output::printAlignment(const vector<pair<int, int> > &alignment) {
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
	residue_label.resize(2);

	residue_number.resize(2);

	number_to_index_map.resize(2);

	raw_interfaces.resize(2);

	sep_amis.resize(2);
	sep_coords.resize(2);
	sep_residue_number.resize(2);
	sep_residue_label.resize(2);

	//output
	amis.resize(2);
	coords.resize(2);
}

iAlign_Output::~iAlign_Output() {

}
