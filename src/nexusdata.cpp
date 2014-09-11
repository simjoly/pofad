//
//  nexusdata.cpp
//  
//
//  Created by Simon Joly on 13-11-27.
//  Copyright (c) 2013 __MyCompanyName__. All rights reserved.
//

#include "nexusdata.h"

using namespace std;

nexusdata::nexusdata() {
    TaxaRead = 0;
    diagonal = 1;
    triangle = 0;
    Labels = 1;
    interleave = 0;
    datatype = "DNA";
	IsDistanceMatrix = false;
}
nexusdata::~nexusdata(){
    delete [] dist_matrix;
    delete [] char_matrix;
}
void nexusdata::InitializeMember(){
    NTaxa = 0;
    NChars = 0;
    Taxa.clear();
    TaxaRead = 0;
    diagonal = 1;
    Labels = 1;
    interleave = 0;
    missingchar.clear();
    gapchar.clear();
    triangle = 0;
    delete [] dist_matrix;
}
void nexusdata::GetNTaxa(int nb_taxa) {
    this->NTaxa = nb_taxa;
}
int nexusdata::ReturnNTax() {
    return NTaxa;
}
void nexusdata::GetNChars(int nb_chars) {
    this->NChars = nb_chars;
}
int nexusdata::ReturnNChars() {
    return NChars;
}
void nexusdata::ReadTaxa(string a_taxa) {
    Taxa.push_back(a_taxa);
}
string nexusdata::ReturnTaxa(int number) {
    return Taxa[number];
}
int nexusdata::NumberTaxaLabels() {
    return Taxa.size();
}
void nexusdata::isTaxa() {
    this->TaxaRead = 1;
}
void nexusdata::IsDiagonal(int diag) {
    this->diagonal = diag;
}
int nexusdata::ReturnDiagonal() {
    return diagonal;
}
void nexusdata::isLabels(int number) {
    this->Labels = number;
}
int nexusdata::ReturnIsLabels() {
    return Labels;
}
void nexusdata::IsInterleave(int a_number) {
    this->interleave = a_number;
}
int nexusdata::ReturnInterleave() {
    return interleave;
}
void nexusdata::MissingChar(string a_character) {
    this->missingchar = a_character;
}
string nexusdata::ReturnMissingChar() {
    return missingchar;
}
void nexusdata::GetGapChar(string a_character) {
    this->gapchar = a_character;
}
string nexusdata::ReturnGapChar() {
    return gapchar;
}
void nexusdata::InitMatrix() {
    int x,y;
    dist_matrix = new double* [NTaxa];   //Allocate place for matrix
    for (y=0; y < NTaxa; y++)
        dist_matrix[y] = new double[NTaxa];
    if (!dist_matrix)
    {
        cerr << "Can't allocate space for creating the distance matrix of haplotype" << endl;
        exit(1);
    }
    for (x=0; x < NTaxa; x++) //Initialize the matrix with 0s
    {
        for (y=0; y < NTaxa; y++)
        {
            dist_matrix[x][y] = 0;
        }
    }
}
void nexusdata::IsTriangle(int a_number) {
    this->triangle = a_number;
}
int nexusdata::ReturnTriangle(void){
    return triangle;
}
void nexusdata::EnterMatrix(double distance, int i, int j) {
    dist_matrix[i][j] = distance;
}
double nexusdata::ReturnDist(int i, int j) {
    
#ifdef DEBUG
    cout << "nex_data:" << i << "," << j << ":";
    cout << dist_matrix[i][j] << endl;
#endif
    
    return dist_matrix[i][j];
}

/**** New Functions ****/

bool nexusdata::ReturnIsDistanceMatrix(void) {
	return IsDistanceMatrix;
}
void nexusdata::SetIsDistanceMatrix(bool a_bool) {
	IsDistanceMatrix = a_bool;
}
int nexusdata::GetPositionofTaxa(string a_taxa) {
    unsigned int i;
    int position;
    vector<string> a_vector;   //TODO: essayer a_vector = Taxa
    for (i=0; i < Taxa.size() ; i++)
    {
        a_vector.push_back(Taxa[i]);
    }
    vector<string>::iterator an_iterator;
    an_iterator = find(a_vector.begin(), a_vector.end(), a_taxa);
    if (an_iterator == a_vector.end())
    {
        cout << a_taxa << endl;
        cerr << "Taxa: \"" << a_taxa << "\" was not found in the taxa block";
        exit(1);
    }
    else
    {
        a_vector.erase(an_iterator, a_vector.end());
        position = a_vector.size();
    } 
    return position;
}

int nexusdata::ReturnisTaxa() {
    return TaxaRead;
}

void nexusdata::GetDatatype(string the_datatype) {
    this->datatype = the_datatype;
}
string nexusdata::ReturnDatatype() {
    return datatype;
}
void nexusdata::InitDataMatrix() {
    char_matrix = new string[NTaxa];
    if (!char_matrix) {
        cerr << "Can't allocate space for creating the character matrix" << endl;
        exit(1);
    }
    int i;
    /* Initialise the matrix */
    for (i=0; i<NTaxa; i++) {
        char_matrix[i] = "";
    }
}
void nexusdata::AddCharacter(int i, string characters) {
    char_matrix[i] += characters;
}
int nexusdata::ReturnCharacterForTaxa(int i) {
    return char_matrix[i].size();
}
char nexusdata::ReturnChar(int i, int j) {
    return char_matrix[i][j];
}
string nexusdata::ReturnSequence(int i) {
    return char_matrix[i];
}
void nexusdata::ComputeDistanceMatrix(void) {
	int i,j,k;
	for (i=0;i<NTaxa;i++)
    {
        cout << "  -> Computing allelic distances: " << ((i)*100/NTaxa) << " % \r";
        cout.flush();
		for (j=0;j<i;j++)
        {
		    double finaldistance;
		    double totaldistance = 0;
		    double numbernucleotides = 0;
			for (k=0; k < NChars; k++)
            {
				string char_i="", char_j="";
				char_i += char_matrix[i][k];
				char_j += char_matrix[j][k];
				if (info_data.ReturnIsIgnoreMissingData()) {
					if ( (char_i == "?") || (char_j == "?") ) continue;
					if ( (info_data.ReturnGapHandling() == 0) && ( (char_i == "-") || (char_j == "-") ) ) continue;
                }
				else if (!info_data.ReturnIsIgnoreMissingData()) {
					if ((char_i == "?") || (char_j == "?")) {
						numbernucleotides += 1;
						continue;
                    }
					else if ((info_data.ReturnGapHandling() == 0) && ( (char_i == "-") || (char_j == "-") )) {
						numbernucleotides += 1;
						continue;
                    }
                }
				if (char_i != char_j) totaldistance += 1;
				numbernucleotides += 1;
            }
			finaldistance = (totaldistance / numbernucleotides);
			dist_matrix[i][j] = finaldistance;
			dist_matrix[j][i] = finaldistance;
        }
		dist_matrix[i][j] = 0;
    }
	cout << endl;
}
