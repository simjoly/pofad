/***************************\
*                           *
*     information class     *
*                           *
\***************************/

/*

    (c) Copyright 2006-2008 by Simon Joly

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


#ifndef _INFORMATION_CLASS_DEFINED_

#include <string>
#include <iostream>
#include <vector>
#include "GenFunctions.h"

using namespace std;

class information {
    private:
        bool Verbose;
        bool OrganismsFile;
        bool DatasetFile;
        bool OutputFile;
		bool BatchFile;
		bool SequenceFile;
		bool JukesCantor;
		int DistanceMethod;  // 0=genpofad (default), 1=matchstates, 2=mrca, 3=pofad, 4=Output consensus sequences
        string file_of_organisms;
        string file_of_datasets;
        string file_for_sequences;
        string output_file;
        int GapHandling;           //0=missing, 1=5th state
        int Number_of_datasets;
        vector<string> datasets;
        bool StandDistanceFileOutput; //true = standardize output file
        int **AllelicInfo;			  //This contains how many alleles were found for each ind for all datasets
        bool IgnoreMissingData;       //true=ignore missing data, false=optimize to minimize distance
        bool EstimateMissingDist;     //true=estimate missing data, false=leave them in the final matrix (-999)

    public:
        information();
		~information();
        void SetVerbose();
        bool IsVerbose();
        void SetIsBatchFile();
        bool ReturnIsBatchFile();
        void IsOrganismsFile();
        bool ReturnIsOrganismsFile(); 
        void IsDatasetFile();
        bool ReturnIsDatasetFile();
        void IsOutputFile();
        bool ReturnIsOutputFile();
        void IsSequenceFile();  	 		/* Use for consensus function */
        bool ReturnIsSequenceFile();  		/* Use for consensus function */
        void SetIsJukesCantor(bool abool);
        bool ReturnIsJukesCantor();
        void IsDistanceFileOutput();
        bool ReturnIsDistanceFileOutput();
        void SetDistanceMethod(int method);
        int ReturnDistanceMethod();
        void SetOrganismsFile(string a_file);
        string ReturnOrganismsFile();
        void SetDatasetsFile(string a_file);
        string ReturnDatasetsFile();
        void SetOutputFile(string a_file);
        string ReturnOutputFile();
        void SetSequenceFile(string a_file);
        string ReturnSequenceFile();
        void SetGapHandling(int a_number);
        int ReturnGapHandling();
        void GetDatasetName(string a_string);
        string ReturnDatasetName(int a_number);
        void GetDatasetsNumber(int a_number);
        int ReturnDatasetsNumber();
        int ReturnDatasetsSize();
        void SetStandDistanceFileOutput();
        void SetRawDistanceFileOutput();
        bool ReturnIsStandDistanceFileOutput();
        void InitAllelicInfo(int NOrg, int NDatasets);
		void EnterNbAllelesforOrganism(int Organism, int dataset, int NbAlleles);
		int ReturnNbAllelesforOrganism(int Organism, int dataset);
        void SetIsIgnoreMissingData(bool boolvalue);
        bool ReturnIsIgnoreMissingData();
        void SetIsEstimateMissingDist(bool boolvalue);
        bool ReturnIsEstimateMissingDist();

};

#define _INFORMATION_CLASS_DEFINED_
#endif /* _INFORMATION_CLASS_DEFINED_ */
