/**********************************************\
 |                                              |
 |                    POFAD                     |
 |                                              |
 |   Phylogeny of Organisms From Allelic Data   |
 |                                              |
 |                                              |
 |                version 1.06                  |
 |                                              |
 |                                              |
 |----------------------------------------------|
 |                                              |
 |              pre-release version             |
 |                                              |
 |                                              |
 |   (c) Copyright 2005-2014 by Simon Joly      |
 |                                              |
 \**********************************************/

/*
 Pofad is a free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, version 2.
 
 Pofad is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/*
 
 Instructions
 
 Description:
 This program computes the genetic distance between organisms from the distance between their
 alleles. The exact description of the method can be found in Joly and Bruneau, in prep. A
 maximum of two alleles per individual is acepted in the current version.
 
 Usage:
 pofad [-i] [file_containning_the_organims]
 [-b] [Run in batch file mode]
 [-d] [name of file containing datasets names]
 [-a] [distance method: 0=genpofad, 1=matchstates (default),
 2=mrca, 3=classic pofad, 4=FRQ, 5=PBC,
 6=Output consensus sequences, 8=MIN, 9=Nei, 10=2ISP]
 [-c] [name of file for consensus sequence]
 [-m] [multiple hits correction: 0=none, 1=JC model]
 [-w] [distance standardization: 0=false, 1=true]
 [-g] [Gives gap handling: 0=missing, 1=fifth state]
 [-?] [missing nucleotides: 0=ignore, 1=infer]
 [-z] [missing distances: 0=leave them, 1=infer]
 [-o] [name_for_the_output_file]
 [-v] [sets the program in verbose mode]
 
 Warnings:
 - The programs does not check if the number of haplotypes in the haplotype file correspond to the size of the matrix. It is assume that the matrix will have the size corresponding to the number of haplotypes.
 - The maximal length of the organisms or haplotype is of 99 characters.
 
 Notes:
 - The maximum length of each distance output by the program is set to "8" by default. If you
 want to modify this number, change the number attributed to 'PRECISION' in the definitions at
 the beginning of the program.
 
 */


#ifndef _pofad_h
#define _pofad_h

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <time.h>                      // define time()
#include "GenFunctions.h"
#include "information.h"
#include "organismsdata.h"
#include "nexusdata.h"
#include "ReadNexusFile.h"
#include "information.h"

#define	MAX_BUFFER_SIZE	20000000
#define VERSION 1.07
#define PRECISION 8
//#define DEBUG                   /* To turn on debug, remove the "//" bofore "#define" */

using namespace std;

const static string Consensus[31][31]= { \
    /*		A		C		T		G		R		Y		K		M		S		W		B		D		H		V		N		-	a	c	t	g	r	y	k	m	s	w	b	d	h	v	n */
    /*A*/	{"A",	"M",	"W",	"R",	"R",	"H",	"D",	"M",	"V",	"W",	"N",	"D",	"H",	"V",	"N",	"a",	"a",	"m",	"w",	"g",	"r",	"h",	"d",	"m",	"v",	"w",	"n",	"d",	"h",	"v",	"n"},
    /*C*/	{"M",	"C",	"Y",	"S",	"V",	"Y",	"B",	"M",	"S",	"H",	"B",	"N",	"H",	"V",	"N",	"c",	"m",	"c",	"y",	"s",	"v",	"y",	"b",	"m",	"s",	"h",	"b",	"n",	"h",	"v",	"n"},
    /*T*/	{"W",	"Y",	"T",	"K",	"D",	"Y",	"K",	"H",	"B",	"W",	"B",	"D",	"H",	"N",	"N",	"t",	"w",	"y",	"t",	"k",	"d",	"y",	"k",	"h",	"b",	"w",	"b",	"d",	"h",	"n",	"n"},
    /*G*/	{"R",	"S",	"K",	"G",	"R",	"B",	"K",	"V",	"S",	"D",	"B",	"D",	"N",	"V",	"N",	"g",	"r",	"s",	"k",	"g",	"r",	"b",	"k",	"v",	"s",	"d",	"b",	"d",	"n",	"v",	"n"},
    /*R*/	{"R",	"V",	"D",	"R",	"R",	"N",	"D",	"V",	"V",	"D",	"N",	"D",	"N",	"V",	"N",	"r",	"r",	"v",	"d",	"r",	"r",	"n",	"d",	"v",	"v",	"d",	"n",	"d",	"n",	"v",	"n"},
    /*Y*/	{"H",	"Y",	"Y",	"B",	"N",	"Y",	"B",	"H",	"B",	"H",	"B",	"N",	"H",	"N",	"N",	"y",	"h",	"y",	"y",	"b",	"n",	"y",	"b",	"h",	"b",	"h",	"b",	"n",	"h",	"n",	"n"},
    /*K*/	{"D",	"B",	"K",	"K",	"D",	"B",	"K",	"N",	"B",	"D",	"B",	"D",	"N",	"N",	"N",	"k",	"d",	"b",	"k",	"k",	"d",	"b",	"k",	"n",	"b",	"d",	"b",	"d",	"n",	"n",	"n"},
    /*M*/	{"M",	"M",	"H",	"V",	"V",	"H",	"N",	"M",	"V",	"H",	"N",	"N",	"H",	"V",	"N",	"m",	"m",	"m",	"h",	"v",	"v",	"h",	"n",	"m",	"v",	"h",	"n",	"n",	"h",	"v",	"n"},
    /*S*/	{"V",	"S",	"B",	"S",	"V",	"B",	"B",	"V",	"S",	"N",	"B",	"N",	"N",	"T",	"N",	"s",	"v",	"s",	"b",	"s",	"v",	"b",	"b",	"v",	"s",	"n",	"b",	"n",	"n",	"t",	"n"},
    /*W*/	{"W",	"H",	"W",	"D",	"D",	"H",	"D",	"H",	"N",	"W",	"N",	"D",	"H",	"N",	"N",	"w",	"w",	"h",	"w",	"d",	"d",	"h",	"d",	"h",	"n",	"w",	"n",	"d",	"h",	"n",	"n"},
    /*B*/	{"N",	"B",	"B",	"B",	"N",	"B",	"B",	"N",	"B",	"N",	"B",	"N",	"N",	"N",	"N",	"b",	"n",	"b",	"b",	"b",	"n",	"b",	"b",	"n",	"b",	"n",	"b",	"n",	"n",	"n",	"n"},
    /*D*/	{"D",	"N",	"D",	"D",	"D",	"N",	"D",	"N",	"N",	"D",	"N",	"D",	"N",	"N",	"N",	"d",	"d",	"n",	"d",	"d",	"d",	"n",	"d",	"n",	"n",	"d",	"n",	"d",	"n",	"n",	"n"},
    /*H*/	{"H",	"H",	"H",	"N",	"N",	"H",	"N",	"H",	"N",	"H",	"N",	"N",	"H",	"N",	"N",	"h",	"h",	"h",	"h",	"n",	"n",	"h",	"n",	"h",	"n",	"h",	"n",	"n",	"h",	"n",	"n"},
    /*V*/	{"V",	"V",	"N",	"V",	"V",	"N",	"N",	"V",	"T",	"N",	"N",	"N",	"N",	"V",	"N",	"v",	"v",	"v",	"n",	"v",	"v",	"n",	"n",	"v",	"t",	"n",	"n",	"n",	"n",	"v",	"n"},
    /*N*/	{"N",	"N",	"N",	"N",	"N",	"N",	"N",	"N",	"N",	"N",	"N",	"N",	"N",	"N",	"N",	"n",	"n",	"n",	"n",	"n",	"n",	"n",	"n",	"n",	"n",	"n",	"n",	"n",	"n",	"n",	"n"},
    /*-*/		{"a",	"c",	"t",	"g",	"r",	"y",	"k",	"m",	"s",	"w",	"b",	"d",	"h",	"v",	"n",	"-",	"a",	"c",	"t",	"g",	"r",	"y",	"k",	"m",	"s",	"w",	"b",	"d",	"h",	"v",	"n"},
    /*a*/	{"a",	"m",	"w",	"g",	"r",	"h",	"d",	"m",	"v",	"w",	"n",	"d",	"h",	"v",	"n",	"a",	"a",	"m",	"w",	"g",	"r",	"h",	"d",	"m",	"v",	"w",	"n",	"d",	"h",	"v",	"n"},
    /*c*/	{"m",	"c",	"y",	"s",	"v",	"y",	"b",	"m",	"s",	"h",	"b",	"n",	"h",	"v",	"n",	"c",	"m",	"c",	"y",	"s",	"v",	"y",	"b",	"m",	"s",	"h",	"b",	"n",	"h",	"v",	"n"},
    /*t*/		{"w",	"y",	"t",	"k",	"d",	"y",	"k",	"h",	"b",	"w",	"b",	"d",	"h",	"n",	"n",	"t",	"w",	"y",	"t",	"k",	"d",	"y",	"k",	"h",	"b",	"w",	"b",	"d",	"h",	"n",	"n"},
    /*g*/	{"r",	"s",	"k",	"g",	"r",	"b",	"k",	"v",	"s",	"d",	"b",	"d",	"n",	"v",	"n",	"g",	"r",	"s",	"k",	"g",	"r",	"b",	"k",	"v",	"s",	"d",	"b",	"d",	"n",	"v",	"n"},
    /*r*/		{"r",	"v",	"d",	"r",	"r",	"n",	"d",	"v",	"v",	"d",	"n",	"d",	"n",	"v",	"n",	"r",	"r",	"v",	"d",	"r",	"r",	"n",	"d",	"v",	"v",	"d",	"n",	"d",	"n",	"v",	"n"},
    /*y*/	{"h",	"y",	"y",	"b",	"n",	"y",	"b",	"h",	"b",	"h",	"b",	"n",	"h",	"n",	"n",	"y",	"h",	"y",	"y",	"b",	"n",	"y",	"b",	"h",	"b",	"h",	"b",	"n",	"h",	"n",	"n"},
    /*k*/	{"d",	"b",	"k",	"k",	"d",	"b",	"k",	"n",	"b",	"d",	"b",	"d",	"n",	"n",	"n",	"k",	"d",	"b",	"k",	"k",	"d",	"b",	"k",	"n",	"b",	"d",	"b",	"d",	"n",	"n",	"n"},
    /*m*/	{"m",	"m",	"h",	"v",	"v",	"h",	"n",	"m",	"v",	"h",	"n",	"n",	"h",	"v",	"n",	"m",	"m",	"m",	"h",	"v",	"v",	"h",	"n",	"m",	"v",	"h",	"n",	"n",	"h",	"v",	"n"},
    /*s*/		{"v",	"s",	"b",	"s",	"v",	"b",	"b",	"v",	"s",	"n",	"b",	"n",	"n",	"t",	"n",	"s",	"v",	"s",	"b",	"s",	"v",	"b",	"b",	"v",	"s",	"n",	"b",	"n",	"n",	"t",	"n"},
    /*w*/	{"w",	"h",	"w",	"d",	"d",	"h",	"d",	"h",	"n",	"w",	"n",	"d",	"h",	"n",	"n",	"w",	"w",	"h",	"w",	"d",	"d",	"h",	"d",	"h",	"n",	"w",	"n",	"d",	"h",	"n",	"n"},
    /*b*/	{"n",	"b",	"b",	"b",	"n",	"b",	"b",	"n",	"b",	"n",	"b",	"n",	"n",	"n",	"n",	"b",	"n",	"b",	"b",	"b",	"n",	"b",	"b",	"n",	"b",	"n",	"b",	"n",	"n",	"n",	"n"},
    /*d*/	{"d",	"n",	"d",	"d",	"d",	"n",	"d",	"n",	"n",	"d",	"n",	"d",	"n",	"n",	"n",	"d",	"d",	"n",	"d",	"d",	"d",	"n",	"d",	"n",	"n",	"d",	"n",	"d",	"n",	"n",	"n"},
    /*h*/	{"h",	"h",	"h",	"n",	"n",	"h",	"n",	"h",	"n",	"h",	"n",	"n",	"h",	"n",	"n",	"h",	"h",	"h",	"h",	"n",	"n",	"h",	"n",	"h",	"n",	"h",	"n",	"n",	"h",	"n",	"n"},
    /*v*/	{"v",	"v",	"n",	"v",	"v",	"n",	"n",	"v",	"t",	"n",	"n",	"n",	"n",	"v",	"n",	"v",	"v",	"v",	"n",	"v",	"v",	"n",	"n",	"v",	"t",	"n",	"n",	"n",	"n",	"v",	"n"},
    /*n*/	{"n",	"n",	"n",	"n",	"n",	"n",	"n",	"n",	"n",	"n",	"n",	"n",	"n",	"n",	"n",	"n",	"n",	"n",	"n",	"n",	"n",	"n",	"n",	"n",	"n",	"n",	"n",	"n",	"n",	"n",	"n"}
};


/***** Functions declarations *****/

void read_organisms(string name);
static void doOrgDimensions(void);
static void doOrgLabels(void);
int get_position_nexus(string the_string);
double get_nexus_distance(int m, int n);
double get_PBC_distance(int m, int n);
double get_MIN_distance(int m, int n);
double get_pofadX_distance(int m, int n);
void read_nexus_matrix(double ***dist_matrix, double ***stand_dist_matrix, int set);
void AttributeSequencesToOrganisms(int set);
void read_nexus_file(string name);
char get_org_char(int organism, int position);
string char_consensus(char char1, char char2);
void output_file(double **final_matrix, double **final_stand_matrix);
void doDatasetDimensions(void);
void doDatasetLabels(void);
void consensus(string FileName);
void output_file(void);

// Variables and functions used for tree reconstrcution from incomplete matrices  

double	**Dist, **Da, *Long, Dmax=0., eps=0.001;
char	**Et, Nom[20], car;
int		N, Ns, NbInc=0, *Tree, *Ch, *Pl, *Clx, *Cly;
int miss=-99;

/* Function prototypes */

static void EstimateMissingDist(double **DD);
static void alea (int T, int *aa);
static void Ultra1(int *sum, int *sum2, int T, int **B2, double **m, double *max, int *aa);
static void Ultra2(int *sum, int *sum2, int T, int **B2, double **m, double *max);
static void Additif(int *sum, int *sum2, int T, int **B2, double **m, double *max);		
static int floor1(double x);


/***** Private Functions *****/

static void read_datasets(string file_name);
static string nextToken(void);
static string nextToken2(void);
static int parse_assignment2(string target);
static void doUnrecognizedBlock(void);



/***** Global Variables declarations *****/

int output_cons_seq=0;           //A flag that indicates is the program should output a consensus sequence matrix: 0=No, 1=Yes
string buffer;                   //string that will serve to store various entry data
organismsdata org_data;
information info_data;
nexusdata nexus_data;
int i,j,k,z;
string name_of_output_file;      //Name of the output file entered as argument
int NumberOfAlleles=0;
ofstream logfile("logfile.txt", ios::trunc);



/***** Local variables *****/

static char *bufPtr;
static string aTokenPtr;
static string  LocalToken;
int nbultra;



#endif
