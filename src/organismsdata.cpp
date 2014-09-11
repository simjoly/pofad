//
//  organismsdata.cpp
//  
//
//  Created by Simon Joly on 13-11-27.
//  Copyright (c) 2013 __MyCompanyName__. All rights reserved.
//

#include <string>
#include <iostream>
#include <vector>
#include "organismsdata.h"

using namespace std;

organismsdata::organismsdata() {
    IsOrganismsRead = false;
}
organismsdata::~organismsdata(){
    delete [] dist_matrix;
    delete [] char_matrix;
	delete [] FRQMatrix;
	delete [] RelFRQMatrix;
}
void organismsdata::InitializeMember(){
    NOrg = 0;
    NChar = 0;
    Organisms.clear();
    IsOrganismsRead = false;
    datatype.clear();
    missingchar.clear();
    gapchar.clear();
}
void organismsdata::DeleteCharMatrix() {
    delete [] char_matrix;
}
void organismsdata::DeleteFRQMatrix() {
    delete [] FRQMatrix;
    delete [] RelFRQMatrix;
}
void organismsdata::DeleteDistanceMatrix(){
    delete [] dist_matrix;
}
void organismsdata::GetNOrg(int nb_taxa) {
    this->NOrg = nb_taxa;
}
int organismsdata::ReturnNOrg() {
    return NOrg;
}
void organismsdata::GetNChar(int nb_chars) {
    this->NChar = nb_chars;
}
int organismsdata::ReturnNChars() {
    return NChar;
}
void organismsdata::ReadOrganism(string an_organism) {
    Organisms.push_back(an_organism);
}
string organismsdata::ReturnOrganism(int number) {
    return Organisms[number];
}
int organismsdata::NumberOrgLabels() {
    return Organisms.size();
}
void organismsdata::OrganismsRead() {
    this->IsOrganismsRead = true;
}
bool organismsdata::IsOrganisms() {
    return IsOrganismsRead;
}
void organismsdata::MissingChar(string a_character) {
    this->missingchar = a_character;
}
string organismsdata::ReturnMissingChar() {
    return missingchar;
}
void organismsdata::GetGapChar(string a_character) {
    this->gapchar = a_character;
}
string organismsdata::ReturnGapChar() {
    return gapchar;
}
void organismsdata::InitDistMatrix() {
    int x,y;
    dist_matrix = new double* [NOrg];   //Allocate place for matrix
    for (y=0; y < NOrg; y++)
        dist_matrix[y] = new double[NOrg];
    if (!dist_matrix)
    {
        cerr << "Can't allocate space for creating the distance matrix of organisms" << endl;
    }
    for (x=0; x < NOrg; x++) //Initialize the matrix with 0s
    {
        for (y=0; y < NOrg; y++)
        {
            dist_matrix[x][y] = 0.;
        }
    }
}
void organismsdata::EnterDist(double distance, int i, int j) {
    dist_matrix[i][j] = distance;
}
double organismsdata::ReturnDist(int i, int j) {
    
#ifdef DEBUG
    cout << "nex_data:" << i << "," << j << ";";
    cout << dist_matrix[i][j];
#endif
    
    return dist_matrix[i][j];
}
void organismsdata::GetDatatype(string the_datatype) {
    this->datatype = the_datatype;
}
string organismsdata::ReturnDatatype() {
    return datatype;
}
void organismsdata::InitDataMatrix() {
    char_matrix = new string[NOrg];
    if (!char_matrix) {
        cerr << "Can't allocate space for creating the character matrix" << endl;
    }
    int i;
    /* Initialise the matrix */
    for (i=0; i<NOrg; i++) {
        char_matrix[i] = "";
    }
}
void organismsdata::AddCharacter(int i, string characters) {
    char_matrix[i] += characters;
}
int organismsdata::ReturnCharacterForOrg(int i) {
    return char_matrix[i].size();
}
int organismsdata::ReturnLengthofOrg(int i) {
    return Organisms[i].size();
}
char organismsdata::ReturnChar(int i, int j) {
    return char_matrix[i][j];
}
string organismsdata::ReturnSequence(int i) {
    return char_matrix[i];
}

void organismsdata::CalculateDistances(string method) {
    int i=0;
    int j;
    
    if (method == "genpofad") {
        //#pragma omp parallel for schedule(dynamic) collapse(2)
	    for (i=0;i<NOrg;i++) {
            //#ifndef _OPENMP
            cout << "  -> Calculating distances: " << ((i)*100/NOrg) << " %\r";
            cout.flush();
            //#endif
	        for (j=0;j<=i;j++) {
	            if (i==j) {
	                dist_matrix[i][j] = 0;
                    }
	            if ( (char_matrix[i] == "") || (char_matrix[j] == "") ) {//Data is missing for this organisms for this dataset
	                dist_matrix[i][j] = -999;
                    }
	            else {
	                dist_matrix[i][j] = organismsdata::getdistance(i,j,"genpofad");
	                dist_matrix[j][i] = dist_matrix[i][j];
                    }
                }
            }
        }
    
	if (method == "mrca")
    {
	    for (i=0;i<NOrg;i++)
        {
            cout << "  -> Calculating distances: " << ((i)*100/NOrg) << " %\r";
            cout.flush();
	        for (j=0;j<=i;j++)
            {
	            if (i==j)
                {
	                dist_matrix[i][j] = 0;
                }
	            if ( (char_matrix[i] == "") || (char_matrix[j] == "") ) //Data is missing for this organisms for this dataset
                {
	                dist_matrix[i][j] = -999;
                }
	            else
                {
	                dist_matrix[i][j] = organismsdata::getdistance(i,j,"mrca");
	                dist_matrix[j][i] = dist_matrix[i][j];
                }
            }
        }
    }
    
	if (method == "matchstates")
    {
	    for (i=0;i<NOrg;i++)
        {
            cout << "  -> Calculating distances: " << ((i)*100/NOrg) << " %\r";
            cout.flush();
	        for (j=0;j<=i;j++)
            {
	            if (i==j)
                {
	                dist_matrix[i][j] = 0;
                }
	            if ( (char_matrix[i] == "") || (char_matrix[j] == "") ) //Data is missing for this organisms for this dataset
                {
	                dist_matrix[i][j] = -999;
                }
	            else
                {
	                dist_matrix[i][j] = organismsdata::getdistance(i,j,"matchstates");
	                dist_matrix[j][i] = dist_matrix[i][j];
                }
            }
        }
    }
	if (method == "PP")
    {
	    for (i=0;i<NOrg;i++)
        {
            cout << "  -> Calculating distances: " << ((i)*100/NOrg) << " %\r";
            cout.flush();
	        for (j=0;j<=i;j++)
            {
	            if (i==j)
                {
	                dist_matrix[i][j] = 0;
                }
	            if ( (char_matrix[i] == "") || (char_matrix[j] == "") ) //Data is missing for this organisms for this dataset
                {
	                dist_matrix[i][j] = -999;
                }
	            else
                {
	                dist_matrix[i][j] = organismsdata::getdistance(i,j,"PP");
	                dist_matrix[j][i] = dist_matrix[i][j];
                }
            }
        }
    }
	if (method == "FRQ")
    {
		IsMissingOrg.clear();
		// Build the FREQ matrix
		BuildFRQMatrix();
		int w,x,y,z;
		int alleles;
		double distance;
		
		if (info_data.ReturnIsIgnoreMissingData()) {
			//Calculate relative frequencies in the FRQ matrix
			x=0;
			for (x=0;x<NOrg;x++) {
                cout << "  -> Building Relative FRQ matrix: " << ((x)*100/NOrg) << " %\r";
                cout.flush();
				if (IsMissingOrg[x]) continue;
				for (w=0;w<NChar;w++) {
					alleles=0;
					for (y=0;y<4;y++) {
						alleles += FRQMatrix[x][w][y];
                    }
					if (alleles == 0) { //If the organism has missing states at this nucleotide position
						for (y=0;y<4;y++) {
							RelFRQMatrix[x][w][y] = 0.;
                        }
						continue;
                    }
					for (y=0;y<4;y++) {
						RelFRQMatrix[x][w][y] = FRQMatrix[x][w][y]/alleles;
                    }
                }
            }
			cout << endl;
			//Calculate organism distances
			x=0;
			for (x=0;x<NOrg;x++) {
                cout << "  -> Calculating distances: " << ((x)*100/NOrg) << " %\r";
                cout.flush();
				for (y=0;y<NOrg;y++) {
					distance=0.;
					if (IsMissingOrg[x] || IsMissingOrg[y]) {
						dist_matrix[i][j] = -999;
						dist_matrix[j][i] = -999;
						continue;
                    }
					for (w=0;w<NChar;w++) {
						for (z=0;z<4;z++) {
							distance += pow( (RelFRQMatrix[x][w][z] - RelFRQMatrix[y][w][z]) , 2);
                        }
                    }
					dist_matrix[x][y] = sqrt(distance);
					dist_matrix[y][x] = dist_matrix[x][y];
                }			
				dist_matrix[x][x] = 0.;
            }
        } //if (info_data.ReturnIsIgnoreMissingData())
		else {
			//Calculate relative frequencies in the FRQ matrix
			x=0;
			for (x=0;x<NOrg;x++) {
                cout << "  -> Building Relative FRQ matrix: " << ((x)*100/NOrg) << " %\r";
                cout.flush();
				if (IsMissingOrg[x]) continue;
				for (w=0;w<NChar;w++) {
					alleles=0;
					for (y=0;y<5;y++) {
						alleles += FRQMatrix[x][w][y];
                    }
					if (alleles == 0) { //If the organism has missing states at this nucleotide position
						for (y=0;y<5;y++) {
							RelFRQMatrix[x][w][y] = 0.;
                        }
						continue;
                    }
					for (y=0;y<5;y++) {
						RelFRQMatrix[x][w][y] = FRQMatrix[x][w][y]/alleles;
                    }
                }
            }
			cout << endl;
			//Calculate organisms distances
			x=0;
			for (x=0;x<NOrg;x++) {
                cout << "  -> Calculating distances: " << ((x)*100/NOrg) << " %\r";
                cout.flush();
				for (y=0;y<x;y++) {
					distance=0.;
					if (IsMissingOrg[x] || IsMissingOrg[y]) {
						dist_matrix[x][y] = -999;
						dist_matrix[y][x] = -999;
						continue;
                    }
					for (w=0;w<NChar;w++) {
						for (z=0;z<5;z++) {
							distance += pow( (RelFRQMatrix[x][w][z] - RelFRQMatrix[y][w][z]) , 2);
                        }
                    }
					dist_matrix[x][y] = sqrt(distance);
					dist_matrix[y][x] = dist_matrix[x][y];
                }
				dist_matrix[x][x] = 0.;
            }
        } // else
    }
	if (method == "nei") {
		IsMissingOrg.clear();
		// Build the FREQ matrix
		BuildFRQMatrix();
		int w,x,y,z;
		int alleles;
		double distance,distance_char;
		
		if (info_data.ReturnIsIgnoreMissingData()) {
			//Calculate relative frequencies in the FRQ matrix
			x=0;
			for (x=0;x<NOrg;x++) {
                cout << "  -> Building Relative FRQ matrix: " << ((x)*100/NOrg) << " %\r";
                cout.flush();
				if (IsMissingOrg[x]) continue;
				for (w=0;w<NChar;w++) {
					alleles=0;
					for (y=0;y<4;y++) {
						alleles += FRQMatrix[x][w][y];
                    }
					if (alleles == 0) { //If the organism has missing states at this nucleotide position
						for (y=0;y<4;y++) {
							RelFRQMatrix[x][w][y] = 0.;
                        }
						continue;
                    }
					for (y=0;y<4;y++) {
						RelFRQMatrix[x][w][y] = FRQMatrix[x][w][y]/alleles;
                    }
                }
            }
			cout << endl;
			//Calculate organism distances
			x=0;
			for (x=0;x<NOrg;x++) {
                cout << "  -> Calculating distances: " << ((x)*100/NOrg) << " %\r";
                cout.flush();
				for (y=0;y<NOrg;y++) {
					distance=0.;
					if (IsMissingOrg[x] || IsMissingOrg[y]) {
						dist_matrix[i][j] = -999;
						dist_matrix[j][i] = -999;
						continue;
                    }
					for (w=0;w<NChar;w++) {
                        distance_char=0.;
						for (z=0;z<4;z++) {
							distance_char +=  sqrt(RelFRQMatrix[x][w][z] * RelFRQMatrix[y][w][z]);
                        }
                        distance += (1-distance_char);
                    }
					dist_matrix[x][y] = distance/NChar;
					dist_matrix[y][x] = dist_matrix[x][y];
                }			
				dist_matrix[x][x] = 0.;
            }
        } //if (info_data.ReturnIsIgnoreMissingData())
		else {
			//Calculate relative frequencies in the FRQ matrix
			x=0;
			for (x=0;x<NOrg;x++) {
                cout << "  -> Building Relative FRQ matrix: " << ((x)*100/NOrg) << " %\r";
                cout.flush();
				if (IsMissingOrg[x]) continue;
				for (w=0;w<NChar;w++) {
					alleles=0;
					for (y=0;y<5;y++) {
						alleles += FRQMatrix[x][w][y];
                    }
					if (alleles == 0) { //If the organism has missing states at this nucleotide position
						for (y=0;y<5;y++) {
							RelFRQMatrix[x][w][y] = 0.;
                        }
						continue;
                    }
					for (y=0;y<5;y++) {
						RelFRQMatrix[x][w][y] = FRQMatrix[x][w][y]/alleles;
                    }
                }
            }
			cout << endl;
			//Calculate organisms distances
			x=0;
			for (x=0;x<NOrg;x++) {
                cout << "  -> Calculating distances: " << ((x)*100/NOrg) << " %\r";
                cout.flush();
				for (y=0;y<x;y++) {
					distance=0.;
					if (IsMissingOrg[x] || IsMissingOrg[y]) {
						dist_matrix[x][y] = -999;
						dist_matrix[y][x] = -999;
						continue;
                    }
					for (w=0;w<NChar;w++) {
                        distance_char=0.;
						for (z=0;z<5;z++) {
							distance_char += sqrt(RelFRQMatrix[x][w][z] * RelFRQMatrix[y][w][z]);
                        }
                        distance += (1-distance_char);
                    }
					dist_matrix[x][y] = distance/NChar;
					dist_matrix[y][x] = dist_matrix[x][y];
                }
				dist_matrix[x][x] = 0.;
            }
        } // else
    }
	cout << endl;
}

double organismsdata::getdistance(int i, int j, string method) {
    double finaldistance;
    double totaldistance = 0;
    double numbercomparisons = 0;
    int k;
    
    for (k=0; k < NChar; k++)
    {
        string char_i, char_j;
        char_i += char_matrix[i][k];
        char_j += char_matrix[j][k];
        if (info_data.ReturnIsIgnoreMissingData()) {
        	if ( (char_i == "?") || (char_j == "?") ) continue;
	        if ( (info_data.ReturnGapHandling() == 0) && ( (char_i == "-") || (char_j == "-") ) ) continue;
        }
        if (method == "genpofad") totaldistance += organismsdata::getpofaddist(char_i,char_j);
        if (method == "mrca") totaldistance += organismsdata::getmrcadist(char_i,char_j);
        if (method == "matchstates") totaldistance += organismsdata::getmatchstatesdist(char_i,char_j);
        if (method == "PP") totaldistance += organismsdata::getPPdist(char_i,char_j);
        numbercomparisons += 1;
    }
    
    finaldistance = (totaldistance / numbercomparisons);
    return finaldistance;
}

double organismsdata::getpofaddist(string firstchar, string secondchar) {
    
	int c;
	int PositionOfFirstChar=-999;
	int PositionOfSecondChar=-999;
    
    // If both characters are identicals
    if (firstchar == secondchar) return 0;
    
	if ((firstchar == "?") || (secondchar == "?")) return 0;
    
    // Handling of gaps if treated as missing data (only option presently)
	if ( (info_data.ReturnGapHandling() == 0) && ((firstchar == "-" ) || (secondchar == "-" )) ) return 0;
    
	for (c=0; c<=31;c++)
    {
		if (Characters[c] == firstchar)
        {
			PositionOfFirstChar = c;
			break;
        }
    }
    
	for (c=0; c<=31;c++)
    {
		if (Characters[c] == secondchar)
        {
			PositionOfSecondChar = c;
			break;
        }
    }
    
	if (PositionOfFirstChar == -999)
    {
        cout << endl << endl << "Unknown character: " << firstchar << endl;
		cerr << "FATAL!: Look at the manual for allowed characters" << endl;
        exit(1);
    }
	if (PositionOfSecondChar == -999)
    {
        cout << endl << endl << "Unknown character: " << secondchar << endl;
		cerr << "FATAL!: Look at the manual for allowed characters" << endl;
        exit(1);
    }
	
	return DistGenpofad[PositionOfFirstChar][PositionOfSecondChar];
	
}


double organismsdata::getmrcadist(string firstchar, string secondchar) {
    
	int c;
	int PositionOfFirstChar=-999;
	int PositionOfSecondChar=-999;
    
    // If both characters are identicals
    if (firstchar == secondchar) return 0;
    
	if ((firstchar == "?") || (secondchar == "?")) return 0;
    
	for (c=0; c<=31;c++)
    {
		if (Characters[c] == firstchar)
        {
			PositionOfFirstChar = c;
			break;
        }
    }
    
	for (c=0; c<=31;c++)
    {
		if (Characters[c] == secondchar)
        {
			PositionOfSecondChar = c;
			break;
        }
    }
    
	if (PositionOfFirstChar == -999)
    {
        cout << endl << endl << "Unknown character: " << firstchar << endl;
		cerr << "FATAL!: Look at the manual for allowed characters" << endl;
        exit(1);
    }
	if (PositionOfSecondChar == -999)
    {
        cout << endl << endl << "Unknown character: " << secondchar << endl;
		cerr << "FATAL!: Look at the manual for allowed characters" << endl;
        exit(1);
    }
	
	return DistMRCA[PositionOfFirstChar][PositionOfSecondChar];
    
}


double organismsdata::getmatchstatesdist(string firstchar, string secondchar) {
    
	int c;
	int PositionOfFirstChar=-999;
	int PositionOfSecondChar=-999;
    
    // If both characters are identicals
    if (firstchar == secondchar) return 0;
    
	if ((firstchar == "?") || (secondchar == "?")) return 0;
    
	for (c=0; c<=31;c++)
    {
		if (Characters[c] == firstchar)
        {
			PositionOfFirstChar = c;
			break;
        }
    }
    
	for (c=0; c<=31;c++)
    {
		if (Characters[c] == secondchar)
        {
			PositionOfSecondChar = c;
			break;
        }
    }
    
	if (PositionOfFirstChar == -999)
    {
        cout << endl << endl << "Unknown character: " << firstchar << endl;
		cerr << "FATAL!: Look at the manual for allowed characters" << endl;
        exit(1);
    }
	if (PositionOfSecondChar == -999)
    {
        cout << endl << endl << "Unknown character: " << secondchar << endl;
		cerr << "FATAL!: Look at the manual for allowed characters" << endl;
        exit(1);
    }
	
	return DistMatchstates[PositionOfFirstChar][PositionOfSecondChar];
    
}


double organismsdata::getPPdist(string firstchar, string secondchar) {
    
	int c;
	int PositionOfFirstChar=-999;
	int PositionOfSecondChar=-999;
    
    // If both characters are identicals
    if (firstchar == secondchar) return 0;
    
	if ((firstchar == "?") || (secondchar == "?")) return 0;
    
	for (c=0; c<=31;c++)
    {
		if (Characters[c] == firstchar)
        {
			PositionOfFirstChar = c;
			break;
        }
    }
    
	for (c=0; c<=31;c++)
    {
		if (Characters[c] == secondchar)
        {
			PositionOfSecondChar = c;
			break;
        }
    }
    
	if (PositionOfFirstChar == -999)
    {
        cout << endl << endl << "Unknown character: " << firstchar << endl;
		cerr << "FATAL!: Look at the manual for allowed characters" << endl;
        exit(1);
    }
	if (PositionOfSecondChar == -999)
    {
        cout << endl << endl << "Unknown character: " << secondchar << endl;
		cerr << "FATAL!: Look at the manual for allowed characters" << endl;
        exit(1);
    }
    
    // If character have gaps...
    if ( (PositionOfFirstChar > 16) || (PositionOfSecondChar > 16)) {
		cerr << "Method PP does not allow the inclusion of gap characters in polymorphic bases. Exiting..." << endl;
        exit(1);       
    }
	
	return DistPP[PositionOfFirstChar][PositionOfSecondChar];
    
}

//-------------------------------------------
// InitFRQMatrix
//-------------------------------------------

void organismsdata::InitFRQMatrix(void) {
	int x,y,z;
	FRQMatrix = new int**[NOrg];
    for (y=0; y < NOrg; y++) {
        FRQMatrix[y] = new int*[NChar];
		for (z=0; z < NChar; z++) {
	        FRQMatrix[y][z] = new int[5];
        }
    }
	if (!FRQMatrix) {
		cerr << "Can't allocate space for creating the FRQ matrix" << endl;
        exit(1);
    }
	/* Initialise the matrix */
    for (x=0; x < NOrg; x++) //Initialize the matrix with 0s
    {
        for (y=0; y < NChar; y++)
        {
			for (z=0; z < 5; z++) {
			    FRQMatrix[x][y][z] = 0;
            }
        }
    }
    
	//Initialise the RelFRQMatrix
	RelFRQMatrix = new double**[NOrg];
    for (y=0; y < NOrg; y++) {
        RelFRQMatrix[y] = new double*[NChar];
		for (z=0; z < NChar; z++) {
	        RelFRQMatrix[y][z] = new double[5];
        }
    }
	if (!RelFRQMatrix) {
		cerr << "Can't allocate space for creating the relative FRQ matrix" << endl;
        exit(1);
    }
	/* Initialise the matrix */
    for (x=0; x < NOrg; x++) //Initialize the matrix with 0s
    {
        for (y=0; y < NChar; y++)
        {
			for (z=0; z < 5; z++) {
			    RelFRQMatrix[x][y][z] = 0.;
            }
        }
    }
}


//-----------------------------------------------------------------
//  BuildFRQMatrix
//-----------------------------------------------------------------

void organismsdata::BuildFRQMatrix(void)
{
	int i=0,j;
	for (i=0; i<NOrg; i++) {
        cout << "  -> Building FRQ Matrix: " << ((i)*100/NOrg) << " % \r";
        cout.flush();
		//If this organisms is not represented in the matrix, ignore it
        if (AllelesInOrganisms[i].size() == 0) {
            IsMissingOrg.push_back(true);
            continue;
        }
        else {
            IsMissingOrg.push_back(false);
            for (j=0; j<nexus_data.ReturnNChars(); j++) {
                get_FRQ_char(i,j);
            }            
        }
    }
	cout << endl;
}


//-----------------------------------------------------------------
// Function : get_FRQ_char
//-----------------------------------------------------------------
 
int organismsdata::get_FRQ_char(int organism_nb, int position)
{
    int i;
    for (i=0;i<AllelesInOrganisms[organism_nb].size();i++) {
        FRQMatrix_AddCharacter(organism_nb,
            position, 
            nexus_data.ReturnChar(get_position_nexus(AllelesInOrganisms[organism_nb][i]),position) );
    }
    return(0);
    // Check if we need to use To_Uppercase for the allele name
    //astring = To_Uppercase(Organisms[organism_nb]);
}

//--------------------------------------------------
//  FRQMatrix_AddCharacter
//--------------------------------------------------

void organismsdata::FRQMatrix_AddCharacter(int organism_nb, int position, char a_character)
{
	if (a_character == 'A') FRQMatrix[organism_nb][position][0]++;
	if (a_character == 'C') FRQMatrix[organism_nb][position][1]++;
	if (a_character == 'T') FRQMatrix[organism_nb][position][2]++;
	if (a_character == 'G') FRQMatrix[organism_nb][position][3]++;
	if (a_character == '-') FRQMatrix[organism_nb][position][4]++;
    //Add the possibility to add polymorphisms
}


//--------------------------------------------------
//  InitializeAllelesInOrganisms
//--------------------------------------------------

void organismsdata::InitializeAllelesInOrganisms(void)
{
    int i;
    for (i=0;i<AllelesInOrganisms.size();i++){
        AllelesInOrganisms[i].clear();
    }
    AllelesInOrganisms.clear();
    vector<string> p;
    for (i=0;i<NOrg;i++) {
        AllelesInOrganisms.push_back(p);
    }
}


//--------------------------------------------------
//  AddAllele
//--------------------------------------------------

int organismsdata::AddAllele(string an_organism,string an_allele)
{
    int i;
    bool OrganismFound=false;
    for (i=0;i<NOrg;i++) {
        if (an_organism != Organisms[i]) continue;
        else {
            AllelesInOrganisms[i].push_back(an_allele);
            OrganismFound=true;
            //cout << an_allele << " in " << an_organism << endl;
            break;
        }
    }
    if (OrganismFound) return(0);
    else return(1);
}


//--------------------------------------------------
//  IsOrgPresent
//--------------------------------------------------

bool organismsdata::IsOrgPresent(int an_org) {
    if (AllelesInOrganisms[an_org].size() == 0) return(true);
    else return(false);
}


//--------------------------------------------------
//  ReturnNbAllelesforOrganism
//--------------------------------------------------

int organismsdata::ReturnNbAllelesforOrganism(int an_org) {
    return (AllelesInOrganisms[an_org].size());
}


//--------------------------------------------------
//  ReturnAlleleFromOrganism
//--------------------------------------------------

string organismsdata::ReturnAlleleFromOrganism(int an_org, int an_allele) {
    return (AllelesInOrganisms[an_org][an_allele]);
}

