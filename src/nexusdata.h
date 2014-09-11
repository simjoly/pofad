/***************************\
*                           *
*      nexusdata class      *
*                           *
\***************************/

/*

    (c) Copyright 2005-2007 by Simon Joly

    This file is part of Pofad.

    Pofad is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    Pofad is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Foobar; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

*/

#ifndef _NEXUSDATA_CLASS_DEFINED_

#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include "GenFunctions.h"
#include "information.h"

extern information info_data;

using namespace std;

class nexusdata {
    private:
        int NTaxa;            // number of taxa in the dataset
        int NChars;           // The number of characters
        vector<string> Taxa;  // A vector containning the taxa's name
        int TaxaRead;         // A flag to show if taxa labels are read: 0=not read; 1=read.
        int diagonal;         // 0 = no diagonal, 1 = diagonal
		int Labels;           // Flag to indicate is labels are in the distance matrix: 0=NoLabels, 1=Labels
        int interleave;       // Flag to indicate if matrix is interleaved: 0=no; 1=yes
        string datatype;      // Indicate the dataype
        string missingchar;   // Character indicating the symbol for missing characters
        string gapchar;       // Character for gaps
        int triangle;         // 0=Lower, 1=Upper, 2=both
        double **dist_matrix; // The distance matrix
		bool IsDistanceMatrix;

/*New variables*/
        string *char_matrix;  // The character matrix

    public:
        nexusdata();
		~nexusdata();
        void InitializeMember();
        void GetNTaxa(int nb_taxa);
        int ReturnNTax();
        void GetNChars(int nb_chars);
        int ReturnNChars();
        void ReadTaxa(string a_taxa);
		string ReturnTaxa(int number);
        int NumberTaxaLabels();
        void isTaxa();
        void IsDiagonal(int diag);
        int ReturnDiagonal();
        void isLabels(int number);
        int ReturnIsLabels();
        void IsInterleave(int a_number);
        int ReturnInterleave();
        void MissingChar(string a_character);
        string ReturnMissingChar();
        void GetGapChar(string a_character);
        string ReturnGapChar();
        void InitMatrix();
        void IsTriangle(int a_number);
        int ReturnTriangle(void);
        void EnterMatrix(double distance, int i, int j);
        double ReturnDist(int i, int j);

/*New functions */
		bool ReturnIsDistanceMatrix(void); //Returns whether a distance matrix was provided in the nexus file
		void SetIsDistanceMatrix(bool);
        int GetPositionofTaxa(string a_taxa);
        int ReturnisTaxa();
        void GetDatatype(string the_datatype);
        string ReturnDatatype();
        void InitDataMatrix();
        void AddCharacter(int i, string characters);
        int ReturnCharacterForTaxa(int i);
        char ReturnChar(int i, int j);                    //Returns character at position j+1 for Taxa[i]
        string ReturnSequence(int i);                     //Returns the sequence of Taxa[i]
		void ComputeDistanceMatrix(void);
};


#define _NEXUSDATA_CLASS_DEFINED_
#endif /* _NEXUSDATA_CLASS_DEFINED_ */
