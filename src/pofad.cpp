//
//  pofad.cpp
//
//  This file is part of POFAD
//  Created by Simon Joly on 13-11-27.

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include "pofad.h"

using namespace std;

/***** Main program *****/

int main(int argc, char **argv)
{
    char *cur;
    char ch;
    char character,selection;
	int w;
	
    for(w=1;w<argc;w++)
        {
        cur=argv[w];
        ch=cur[0];
      
        if (ch == '-')
            {
            ch=toupper(cur[1]);
            switch(ch)
                {
                case 'V':
	            	info_data.SetVerbose();
                    continue;
                case 'B':
	            	info_data.SetIsBatchFile();
                    continue;
                case 'A':
                    w++;
		            if ( w<argc && (argv[w][0] != '-'))    
		                {
		                info_data.SetDistanceMethod(atoi(argv[w]));
			        	}
		            else w--;
			        continue;
                case 'M':
                    w++;
                    if ( w<argc && (argv[w][0] != '-'))
                    {
                        if (atoi(argv[w]) == 0) info_data.SetIsJukesCantor(false);
                        else if (atoi(argv[w]) == 1) info_data.SetIsJukesCantor(true);
                    }
                    else w--;
                    continue;
                case 'W':
                    w++;
                    if ( w<argc && (argv[w][0] != '-'))
                    {
                        if (atoi(argv[w]) == 0) info_data.SetRawDistanceFileOutput();
                        else if (atoi(argv[w]) == 1) info_data.SetStandDistanceFileOutput();
                    }
                    else w--;
                    continue;
                case '?':
                    w++;
                    if ( w<argc && (argv[w][0] != '-'))
                    {
                        if (atoi(argv[w]) == 0) info_data.SetIsIgnoreMissingData(true);
                        else if (atoi(argv[w]) == 1) info_data.SetIsIgnoreMissingData(false); //infer as to minimize distances
                    }
                    else w--;
                    continue;
                case 'Z':
                    w++;
                    if ( w<argc && (argv[w][0] != '-'))
                    {
                        if (atoi(argv[w]) == 0) info_data.SetIsEstimateMissingDist(false);
                        else if (atoi(argv[w]) == 1) info_data.SetIsEstimateMissingDist(true);
                    }
                    else w--;
                    continue;
                case 'I':
                    w++;
		            if ( w<argc && (argv[w][0] != '-'))    
		                {
		                info_data.SetOrganismsFile(argv[w]);
	                    read_organisms(info_data.ReturnOrganismsFile());
			        	}
		            else {w--;}	
		            continue;
                case 'D':
                    w++;
		            if(w<argc && argv[w][0] != '-')    
		                {
                        info_data.SetDatasetsFile(argv[w]);
                        read_datasets(info_data.ReturnDatasetsFile());
                        info_data.IsDatasetFile();
			        	}
		            else {w--;}	
		            continue;
                case 'O':
                    w++;
		            if(w<argc && argv[w][0] != '-')    
		                {
                        info_data.SetOutputFile(argv[w]);
                        info_data.IsOutputFile();
				        }
		            else {w--;}	
		            continue;
                case 'C':
					info_data.SetDistanceMethod(6);
                    w++;
		            if(w<argc && argv[w][0] != '-')    
		                {
                        info_data.SetSequenceFile(argv[w]);
	                    info_data.IsSequenceFile();
				        }
		            else w--;
	            	continue;
                case 'G':
                    w++;
		            if(w<argc && argv[w][0] != '-')    
		                {
                        if (atoi(argv[w]) > 1)
                            {
                            cout << endl << "*** invalid value for gap handling: " << argv[w] << " ***" << endl;
                            cout << "Gaps will be treated as missing" << endl;
                            info_data.SetGapHandling(0);
                            pauseprogram();
                            continue;
                            }
                        info_data.SetGapHandling(atoi(argv[w]));
				        }
		            else w--;
			        continue;
                default:
                    continue;
                }
            }
        }

    string file_of_datasets, name_of_output_file, file_of_organisms, sequence_file;

	if (!info_data.ReturnIsBatchFile())  /* If running batch file, don't do this */
		{
	    for (;;)
	        {
	
			clearscreen();
	
	        cout <<                                                                            endl;
	        cout << " ************************  POFAD version "<<VERSION<<" ***********************" << endl;
	        cout <<                                                                            endl;
	        cout << " Select one of the following options:"                                 << endl;
	        cout <<                                                                            endl;
	        cout << " A     Distance algorithm                  ";
		    if (info_data.ReturnDistanceMethod() == 0) cout <<"genpofad"                    << endl;
		    else if (info_data.ReturnDistanceMethod() == 1) cout <<"matchstates"            << endl;
		    else if (info_data.ReturnDistanceMethod() == 2) cout <<"mrca"                   << endl;
		    else if (info_data.ReturnDistanceMethod() == 3) cout <<"classic pofad"          << endl;
		    else if (info_data.ReturnDistanceMethod() == 4) cout <<"FRQ distance"           << endl;
		    else if (info_data.ReturnDistanceMethod() == 5) cout <<"Phylogenetic Bray-Curtis (PBC)"<< endl;
		    else if (info_data.ReturnDistanceMethod() == 6) cout <<"Consensus sequences"    << endl;
		    else if (info_data.ReturnDistanceMethod() == 7) cout <<"pofadX"                 << endl;
		    else if (info_data.ReturnDistanceMethod() == 8) cout <<"MIN"                    << endl;
		    else if (info_data.ReturnDistanceMethod() == 9) cout <<"Nei's genetic distance" << endl;
		    else if (info_data.ReturnDistanceMethod() == 10) cout <<"PP Distance" << endl;
	        cout << " I     File containing the organisms       ";
	        if (info_data.ReturnIsOrganismsFile()) cout << info_data.ReturnOrganismsFile() << endl;
	        else cout << "(no file entered)" << endl;

			if (info_data.ReturnDistanceMethod() == 6)
				{
	        	cout << " S     Sequence file                       ";
		        if (info_data.ReturnIsSequenceFile()) cout << info_data.ReturnSequenceFile() << endl;
	    	    else cout << "(no file entered)" << endl;
				}
			else
				{
		        cout << " D     File containing the datasets        ";
		        if (info_data.ReturnIsDatasetFile()) cout << info_data.ReturnDatasetsFile() << endl;
	    	    else cout << "(no file entered)" << endl;
	    	    }
	        cout << " O     Name of the output file             ";
	        if (info_data.ReturnIsOutputFile()) {
	            cout << info_data.ReturnOutputFile() << endl;
	            }
	        else {
	            cout << "(no file entered)" << endl;
	            }
	        cout << " W     Distance standardization            ";
	        if (info_data.ReturnIsOutputFile())
	            {
	            if (info_data.ReturnIsStandDistanceFileOutput())
	                {
	                cout << "Standardized distances" << endl;
	                }
	            else cout << "No" << endl;
	            }
	        else {
	            cout << "standardized distances" << endl;
	            }
			if ( (info_data.ReturnDistanceMethod() != 4)  && (info_data.ReturnDistanceMethod() != 6) ) { //No Jukes and Cantor correction for FRQ distance and consensus
				cout << " M     Multiple hits correction            ";
				if (info_data.ReturnIsJukesCantor()) {
					cout << "Jukes and Cantor" << endl;
					}
				else if (!info_data.ReturnIsJukesCantor()) {
					cout << "No" << endl;
					}
				}
			if (info_data.ReturnDistanceMethod() != 3) {
		        cout << " G     Gap handling                        ";
		        if (info_data.ReturnGapHandling() == 0) {
		            cout << "missing" << endl;
		            }
		        else if (info_data.ReturnGapHandling() == 1){
		            cout << "5th state" << endl;
		            }
		        cout << " ?     Missing nucleotides                 ";
		        if (info_data.ReturnIsIgnoreMissingData()) {
		            cout << "Ignore" << endl;
		            }
		        else if (!info_data.ReturnIsIgnoreMissingData()){
		            cout << "Infer" << endl;
		            }
				}
	        cout << " Z     Missing distances                   ";
	        if (info_data.ReturnIsEstimateMissingDist()) {
	            cout << "Estimate" << endl;
	            }
	        else if (!info_data.ReturnIsEstimateMissingDist()){
	            cout << "Leave them (-999)" << endl;
	            }
	
	
	        cout <<                                                                            endl;
	        cout << " Other options:"                                                       << endl;
	        cout <<                                                                            endl;
	        cout << " Q     Quit the program"                                               << endl;
	        cout <<                                                                            endl;
	        cout << " ********************************************************************" << endl;
	        cout <<                                                                            endl;
	        cout << " Type \"Y\" to accept these settings or the letter for one to change: ";
	
	        cin >> character;

            clearscreen();
	
	        switch (toupper(character))
	            {
	            case 'A':
	                for (;;)
	                    {
	                    cout <<                                                        endl;
	                    cout <<                                                        endl;
	                    cout << " Choose one distance method                      " << endl;
	                    cout <<                                                        endl;
	                    cout << " (a)   genpofad                                  " << endl;
	                    cout << " (b)   matchstates                               " << endl;
	                    cout << " (c)   mrca                                      " << endl;
	                    cout << " (d)   classic pofad                             " << endl;
	                    cout << " (e)   FRQ distance                              " << endl;
	                    cout << " (f)   Phylogenetic Bray-Curtis                  " << endl;
	                    cout << " (g)   MIN"                                        << endl;
                        cout << " (h)   Nei's genetic distance"                     << endl;
                        cout << " (i)   PP distance                               " << endl;
	                    cout << " (j)   Create consensus sequences                " << endl;
//	                    cout << " (k)   pofadX                                    " << endl;
	                    cout <<                                                        endl;
	                    cout << " Enter your selection: ";
	
	                    cin >> selection;
	
	                    switch (toupper(selection))
	                        {
	                        case 'A':
	                            info_data.SetDistanceMethod(0); //genpofad
	                            break;
	                        case 'B':
	                            info_data.SetDistanceMethod(1); //matchstates
	                            break;
	                        case 'C':
	                            info_data.SetDistanceMethod(2); //mrca
	                            break;
	                        case 'D':
	                            info_data.SetDistanceMethod(3); //pofad
	                            break;
	                        case 'E':
	                            info_data.SetDistanceMethod(4); //FRQ
	                            break;
	                        case 'F':
	                            info_data.SetDistanceMethod(5); //PBC
	                            break;
	                        case 'K':
	                            info_data.SetDistanceMethod(7); //pofadX
	                            break;
	                        case 'G':
	                            info_data.SetDistanceMethod(8); //MIN
	                            break;
                            case 'H':
                                info_data.SetDistanceMethod(9); //Nei
                                break;
                            case 'I':
                                info_data.SetDistanceMethod(10); //PP
                                break;
	                        case 'J':
	                            info_data.SetDistanceMethod(6); //consensus
	                            break;
	                        default:
	                            cout << endl << " *** Not a valid selection ***" << endl;
	                            cout << " Enter a valid choice: a, b, c, d, e, f, g, h, i, j, or k  " << endl;
	                            continue;
	                        }
	                    break;
	                    }
	                continue;
	            case 'W':
	                for (;;)
	                    {
	                    cout <<                                                        endl;
	                    cout <<                                                        endl;
	                    cout << " How should distance be scaled in the final matrix?" << endl;
	                    cout <<                                                        endl;
	                    cout << " (a)   Standardized per dataset                  " << endl;
	                    cout << " (b)   Use raw distances                         " << endl;
	                    cout <<                                                        endl;
	                    cout << " Enter your selection: ";
	
	                    cin >> selection;
	
	                    switch (toupper(selection))
	                        {
	                        case 'A':
	                            info_data.SetStandDistanceFileOutput();
	                            break;
	                        case 'B':
	                            info_data.SetRawDistanceFileOutput();
	                            break;
	                         default:
	                            cout << endl << " *** Not a valid selection ***" << endl;
	                            cout << " Enter a valid choice: a or b" << endl;
	                            continue;
	                        }
	                    break;
	                    }
	                continue;
	            case 'G':
	                for (;;)
	                    {
	                    cout <<                                                        endl;
	                    cout <<                                                        endl;
	                    cout << " Choose one method for gap handling              " << endl;
	                    cout <<                                                        endl;
	                    cout << " (a)   gap = missing                             " << endl;
	                    cout << " (b)   gap = 5th state                           " << endl;
	                    cout <<                                                        endl;
	                    cout << " Enter your selection: ";
	
	                    cin >> selection;
	
	                    switch (toupper(selection))
	                        {
	                        case 'A':
	                            info_data.SetGapHandling(0);
	                            break;
	                        case 'B':
	                            info_data.SetGapHandling(1);
	                            break;
	                        default:
	                            cout << endl << " *** Not a valid selection ***" << endl;
	                            cout << " Enter a valid choice: a or b" << endl;
	                            continue;
	                        }
	                    break;
	                    }
	                continue;
	            case 'M':
	                for (;;)
	                    {
	                    cout <<                                                        endl;
	                    cout <<                                                        endl;
	                    cout << " Choose a correction method for multiple hits    " << endl;
	                    cout <<                                                        endl;
	                    cout << " (a)   Jukes and Cantor                          " << endl;
	                    cout << " (b)   No correction method                      " << endl;
	                    cout <<                                                        endl;
	                    cout << " Enter your selection: ";
	
	                    cin >> selection;
	
	                    switch (toupper(selection))
	                        {
	                        case 'A':
	                            info_data.SetIsJukesCantor(true);
	                            break;
	                        case 'B':
	                            info_data.SetIsJukesCantor(false);
	                            break;
	                         default:
	                            cout << endl << " *** Not a valid selection ***" << endl;
	                            cout << " Enter a valid choice: a or b" << endl;
	                            continue;
	                        }
	                    break;
	                    }
	                continue;
	            case 'I':
	                cout << endl << " Enter the name of the file containing the organisms block: ";
	                cin >> file_of_organisms;
	                info_data.SetOrganismsFile(file_of_organisms);
	                read_organisms(info_data.ReturnOrganismsFile());
	                pauseprogram();
	                continue;
	            case 'S':
	                cout << endl << " Enter the name of the sequence file: ";
	                cin >> sequence_file;
	                info_data.SetSequenceFile(sequence_file);
					info_data.IsSequenceFile();
	                pauseprogram();
	                continue;
	            case 'D':
	                cout << endl << " Enter the name of the file containing the datasets block: ";
	                cin >> file_of_datasets;
	                info_data.SetDatasetsFile(file_of_datasets);
	                read_datasets(info_data.ReturnDatasetsFile());
	                pauseprogram();
	                continue;
	            case 'O':
	                cout << endl << " Enter the name of the output file: ";
	                cin >> name_of_output_file;
	                info_data.SetOutputFile(name_of_output_file);
	                info_data.IsOutputFile();
	                continue;
	            case '?':
	                for (;;)
	                    {
	                    cout <<                                                        endl;
	                    cout <<                                                        endl;
	                    cout << " Choose one method handling missing data         " << endl;
	                    cout <<                                                        endl;
	                    cout << " (a)   Ignore                                    " << endl;
	                    cout << " (b)   Infer                                     " << endl;
	                    cout <<                                                        endl;
	                    cout << " Enter your selection: ";
	
	                    cin >> selection;
	
	                    switch (toupper(selection))
	                        {
	                        case 'A':
	                            info_data.SetIsIgnoreMissingData(true);
	                            break;
	                        case 'B':
	                            info_data.SetIsIgnoreMissingData(false);
	                            break;
	                        default:
	                            cout << endl << " *** Not a valid selection ***" << endl;
	                            cout << " Enter a valid choice: a or b" << endl;
	                            continue;
	                        }
	                    break;
	                    }
	                continue;
	            case 'Z':
	                for (;;)
	                    {
	                    cout <<                                                        endl;
	                    cout <<                                                        endl;
	                    cout << " Do you want to estimate the missing distances in the final matrix?" << endl;
	                    cout <<                                                        endl;
	                    cout << " (a)   Yes (using additive criterion)         " << endl;
	                    cout << " (b)   No                                     " << endl;
	                    cout <<                                                        endl;
	                    cout << " Enter your selection: ";
	
	                    cin >> selection;
	
	                    switch (toupper(selection))
	                        {
	                        case 'A':
	                            info_data.SetIsEstimateMissingDist(true);
	                            break;
	                        case 'B':
	                            info_data.SetIsEstimateMissingDist(false);
	                            break;
	                        default:
	                            cout << endl << " *** Not a valid selection ***" << endl;
	                            cout << " Enter a valid choice: a or b" << endl;
	                            continue;
	                        }
	                    break;
	                    }
	                continue;
	            case 'Q':
	                exit(1);
	            case 'Y': /* Make sure everything is OK before continuing */
	            	if (!info_data.ReturnIsOrganismsFile()) { /* Make sure the is an organisms file */
	            		cout << endl << " You need to enter a file of organisms!" << endl;
		                pauseprogram();
		                continue;
		                }
		            if ((info_data.ReturnDistanceMethod() == 6) && (!info_data.ReturnIsSequenceFile()))
		            	{
		            	cout <<  endl << " You need to enter a sequence file of haplotypes!" << endl;
		                pauseprogram();
		                continue;
		            	}
	                break;
	            default:
	                cout << endl << " *** Not a valid selection ***" << endl << endl;
	                pauseprogram();
	                continue;
	            }
	        cout << endl;
	        break;
	        }
		}

    while ( !info_data.ReturnIsOrganismsFile() )  //If no organism file was entered as arguments in batchfile mode
        {
        cout << endl << " Enter the name of the file containning the organisms: ";
        cin >> file_of_organisms; 
        info_data.SetOrganismsFile(file_of_organisms);
        read_organisms(info_data.ReturnOrganismsFile());
        }

    // If creating a consensus file
    if (info_data.ReturnDistanceMethod() == 6) /* Create consensus of a file */
    	{
       	consensus(info_data.ReturnSequenceFile());
    	output_file();
    	exit(1);
    	}
	
    if (!info_data.ReturnIsDatasetFile())
        {
        int sets;
        cout << endl << " How many datasets do you have ? ";
        cin >> sets;
        info_data.GetDatasetsNumber(sets);
        for (z=0; z < info_data.ReturnDatasetsNumber(); z++)
            {
            cout << " Enter the filename (no white space) of the dataset #" << (z+1) << ": ";
            cin >> buffer;
            info_data.GetDatasetName(buffer);
            }
        }

	/* Create the distance matrices that will store distance information for all datasets */

    double ***dist_matrix = new double** [info_data.ReturnDatasetsNumber()];          //Allocate place for matrix. The 1st dimension is used to refer 
    for (i=0; i < info_data.ReturnDatasetsNumber(); i++) {                             //to the datasets, the 2nd is for the matrix' rows whereas the
        dist_matrix[i] = new double* [org_data.ReturnNOrg()];                           //third is for the column of the matrix.
        for (j=0; j < org_data.ReturnNOrg(); j++) {
            dist_matrix[i][j] = new double[org_data.ReturnNOrg()];
            }
        }

    if (!dist_matrix) {
        cout << "Can't allocate space for creating the distance matrix of organisms" << endl;
        exit(1);
        }

    double ***stand_dist_matrix = new double** [info_data.ReturnDatasetsNumber()];        //Allocate place for matrix. The 1st dimension is used to refer 
    for (i=0; i < info_data.ReturnDatasetsNumber(); i++) {                                 //to the datasets, the 2nd is for the matrix' rows whereas the
        stand_dist_matrix[i] = new double* [org_data.ReturnNOrg()];                         //third is for the column of the matrix.
        for (j=0; j < org_data.ReturnNOrg(); j++) {
            stand_dist_matrix[i][j] = new double[org_data.ReturnNOrg()];
            }
        }

    if (!stand_dist_matrix) {
        cout << "Can't allocate space for creating the standardized distance matrix of organisms" << endl;
        exit(1);
        }


	/* Create an array that will contain information about how many alleles were found
	   for each organisms in each dataset */
	
	info_data.InitAllelicInfo(org_data.ReturnNOrg(),info_data.ReturnDatasetsNumber());

	/* Writing stuff to log File */
    if (logfile.bad()) {
        cout << "Error in opening output file" << endl;
        }
	time_t tim=time(NULL);
	tm *now=localtime(&tim);
    logfile << "POFAD, version " << VERSION << endl << endl;
	logfile << " log file written on the " << (now->tm_mday) << "/" << (now->tm_mon+1) << "/" << (now->tm_year+1900);
	logfile << " at " << (now->tm_hour) << "h" << (now->tm_min) << "min" << (now->tm_sec) << endl << endl;
    logfile << " Settings:" << endl << endl;
    if (info_data.ReturnIsBatchFile()) logfile << " -> batch file mode" << endl;
    logfile << "  -> Organisms file: " << '\t' << info_data.ReturnOrganismsFile() << endl;
	logfile << "       - Number of organisms read: " << org_data.ReturnNOrg() << endl;
    if (info_data.ReturnIsDatasetFile()) {
    	logfile << "  -> Dataset file: " << '\t' << info_data.ReturnDatasetsFile() << endl;
		logfile << "       - Number of datasets: " << info_data.ReturnDatasetsNumber() << endl;
		}
    logfile << "  -> Distance algorithm used: " << '\t';
    if (info_data.ReturnDistanceMethod() == 0) logfile <<"genpofad"<< endl;
    else if (info_data.ReturnDistanceMethod() == 1) logfile <<"matchstates"<< endl;
    else if (info_data.ReturnDistanceMethod() == 2) logfile <<"mrca"<< endl;
    else if (info_data.ReturnDistanceMethod() == 3) logfile <<"classic pofad"<< endl;
    else if (info_data.ReturnDistanceMethod() == 4) logfile <<"FRQ"<< endl;
    else if (info_data.ReturnDistanceMethod() == 5) logfile <<"Phylogenetic Bray-Curtis"<< endl;
    else if (info_data.ReturnDistanceMethod() == 7) logfile <<"POFADX"<< endl;
    else if (info_data.ReturnDistanceMethod() == 8) logfile <<"MIN"<< endl;
    else if (info_data.ReturnDistanceMethod() == 9) logfile <<"Nei"<< endl;
    else if (info_data.ReturnDistanceMethod() == 10) logfile <<"PP distance"<< endl;
    if (info_data.ReturnIsOutputFile()) logfile << "  -> Output file: " << '\t' << info_data.ReturnOutputFile() << endl;
    logfile << "  -> Distance standardization: " << '\t';
    if (info_data.ReturnIsStandDistanceFileOutput()) logfile <<"Standardized distances"<< endl;
    else if (!info_data.ReturnIsStandDistanceFileOutput()) logfile <<"Non-standardized distances"<< endl;
    logfile << "  -> Multiple hits correction: " << '\t';
    if (info_data.ReturnIsJukesCantor()) logfile <<"Jukes and Cantor"<< endl;
    else if (!info_data.ReturnIsJukesCantor()) logfile <<"No"<< endl;
	if (info_data.ReturnDistanceMethod() != 3) {
	    logfile << "  -> Gap handling: " << '\t';
	    if (info_data.ReturnGapHandling() == 0) logfile <<"missing"<< endl;
	    else if (info_data.ReturnGapHandling() == 1) logfile <<"5th state"<< endl;
	    logfile << "  -> Missing nucleotides: " << '\t';
	    if (info_data.ReturnIsIgnoreMissingData()) logfile <<"ignore"<< endl;
	    else if (!info_data.ReturnIsIgnoreMissingData()) logfile <<"Infer"<< endl;
		}
    logfile << "  -> Missing distances: " << '\t';
    if (info_data.ReturnIsEstimateMissingDist()) logfile <<"Estimate using the additive procedure"<< endl;
    else logfile <<"Leave the missing data (-999)"<< endl;
	logfile << endl << endl;

    /*** Now we have to enter each dataset and make a distance matrix of organisms for each one of them ***/
    for (z=0; z < info_data.ReturnDatasetsNumber(); z++) {
        read_nexus_matrix(dist_matrix, stand_dist_matrix, z);
        }

    /*** Now let's compute the average matrices of all the input matrices ***/
    int missing_cells;         //Number of datasets that lack a given between individual distance
    int nb_sets_contributing;  //Number of datasets that have a given between individual distance

    double **final_matrix = new double* [org_data.ReturnNOrg()];   //Allocate place for matrix
    for (i=0; i < org_data.ReturnNOrg(); i++)
        final_matrix[i] = new double[org_data.ReturnNOrg()];

    if (!final_matrix) {
        cout << "Can't allocate space for creating the distance matrix of organisms" << endl;
        exit(1);
        }

    #ifdef DEBUG
        cout << "space allocated for final matrix !" << endl;
    #endif

    for (i=0; i < org_data.ReturnNOrg(); i++)   //Initialyse the matrix with zeros
        {
        for (j=0; j < org_data.ReturnNOrg(); j++)
            {
            final_matrix[i][j] = 0;
            }
        }

    for (i=0; i < org_data.ReturnNOrg(); i++)    //Adds values of the different datasets
        {
        for (j=0; j <= i; j++)
            {
            missing_cells = 0;
            for (k=0; k < info_data.ReturnDatasetsNumber(); k++)
                {
                if (dist_matrix[k][i][j] == -999)
                    {
                    missing_cells += 1;
                    continue;
                    }
                final_matrix[i][j] += dist_matrix[k][i][j];
                }
            if (missing_cells == info_data.ReturnDatasetsNumber())
                {
                final_matrix[i][j] = -999;
                final_matrix[j][i] = -999;
                }
            else
                {
                nb_sets_contributing = ( info_data.ReturnDatasetsNumber() - missing_cells );
                final_matrix[i][j] /= nb_sets_contributing;
                final_matrix[j][i] = final_matrix[i][j];
                }
            }
        }

    /*** Calculate a final matrix where each individual matrix has the same weight ***/

    double **final_stand_matrix = new double* [org_data.ReturnNOrg()];   //Allocate place for matrix
    for (i=0; i < org_data.ReturnNOrg(); i++)
        final_stand_matrix[i] = new double[org_data.ReturnNOrg()];

    if (!final_stand_matrix) {
        cout << "Can't allocate space for creating the standardized distance matrix of organims" << endl;
        exit(1);
        }

    #ifdef DEBUG
        cout << "space allocated for final matrix !" << endl;
    #endif

    for (i=0; i < org_data.ReturnNOrg(); i++)   //Initialyse the matrix with zeros
        {
        for (j=0; j < org_data.ReturnNOrg(); j++)
            {
            final_stand_matrix[i][j] = 0;
            }
        }

    for (i=0; i < org_data.ReturnNOrg(); i++)    //Adds values of the different datasets
        {
        for (j=0; j <= i; j++)
            {
            missing_cells = 0;
            for (k=0; k < info_data.ReturnDatasetsNumber(); k++)
                {
                if (stand_dist_matrix[k][i][j] == -999)
                    {
                    missing_cells += 1;
                    continue;
                    }
                final_stand_matrix[i][j] += stand_dist_matrix[k][i][j];
                }
            if (missing_cells == info_data.ReturnDatasetsNumber())
                {
                final_stand_matrix[i][j] = -999;
                final_stand_matrix[j][i] = -999;
                }
            else
                {
                nb_sets_contributing = ( info_data.ReturnDatasetsNumber() - missing_cells );
                final_stand_matrix[i][j] /= nb_sets_contributing;
                final_stand_matrix[j][i] = final_stand_matrix[i][j];
                }
            }
        }

    #ifdef DEBUG
        cout << "Average matrix constructed !" << endl;
    #endif

	/* Calculate the number of missing cells in the final matrix */
    double final_missing_cells = 0;         //Number of datasets that lack a given between individual distance
    double total_number_of_cells = 0;       //Number of datasets that have a given between individual distance
    for (i=0; i < org_data.ReturnNOrg(); i++)   //Adds values of the different datasets
		{
        for (j=0; j < i; j++)
            {
            total_number_of_cells += 1;
            if (final_matrix[i][j] == -999) final_missing_cells += 1;
            }
        }

	/* Estimate missing cells in the final matrices*/
	if( (final_missing_cells > 0) && (info_data.ReturnIsEstimateMissingDist()) )
		{
		EstimateMissingDist(final_matrix);
		EstimateMissingDist(final_stand_matrix);
		}

    /*** Calculate several statistics of the final matrix ***/

    double final_percent_missing;
    double max_distance = 0;            //To find the maximum distance of the matrix
    double max_stand_distance = 0;            //To find the maximum distance of the matrix
    double mean_stand_distance;
    double mean_distance;
    double total_distance = 0;
    double total_stand_distance = 0;

    for (i=0; i < org_data.ReturnNOrg(); i++)    //Adds values of the different datasets
        {
        for (j=0; j < i; j++)
            {
            if (final_stand_matrix[i][j] > max_distance) max_stand_distance = final_stand_matrix[i][j];
            if (final_stand_matrix[i][j] != -999) total_stand_distance += final_stand_matrix[i][j];
            }
        }

    for (i=0; i < org_data.ReturnNOrg(); i++)    //Adds values of the different datasets
        {
        for (j=0; j < i; j++)
            {
            if (final_matrix[i][j] > max_distance) max_distance = final_matrix[i][j];
            if (final_matrix[i][j] != -999) total_distance += final_matrix[i][j];
            }
        }

    final_percent_missing = (final_missing_cells / total_number_of_cells) * 100;
    mean_distance = total_distance / (total_number_of_cells - final_missing_cells);
    mean_stand_distance = total_stand_distance / (total_number_of_cells - final_missing_cells);

    cout << " Final matrix" << endl << endl;
    cout << "  -> Number of organisms = ";
    cout << org_data.ReturnNOrg() << endl;
    cout << "  -> Percentage of missing cells in the final matrix: ";
    cout << setprecision(2) << final_percent_missing << "%" << endl;
	if (info_data.ReturnIsEstimateMissingDist() && (final_missing_cells>0)){
		cout << "  -> Number of cells estimated using the ultrametric property: ";
		cout << setprecision(2) << nbultra << " (" << (nbultra/final_missing_cells) << "%)" << endl;
		}
	if (info_data.ReturnIsStandDistanceFileOutput()) {
	    cout << "  -> Average distance (standardized matrix, ignoring missing cells): ";
	    cout << setprecision(5) << mean_stand_distance << endl;
	    cout << "  -> Maximum pair-wise distance (standardized matrix): ";
	    cout << setprecision(5) << max_stand_distance << endl << endl;
		}
	else {
		cout << "  -> Average distance (unstandardized matrix, ignoring missing cells): ";
    	cout << setprecision(5) << mean_distance << endl;
	    cout << "  -> Maximum pair-wise distance (unstandardized matrix): ";
	    cout << setprecision(5) << max_distance << endl << endl;
    	}

	/* Let's print this also in the log file */
    logfile << " Final matrix" << endl << endl;
    logfile << "  -> Number of organisms = ";
    logfile << org_data.ReturnNOrg() << endl;
    logfile << "  -> Percentage of missing cells in the final matrix: ";
    logfile << setprecision(2) << final_percent_missing << "%" << endl;
	if (info_data.ReturnIsEstimateMissingDist() && (final_missing_cells>0)){
		logfile << "  -> Number of cells estimated using the ultrametric property: ";
		logfile << setprecision(2) << nbultra << " (" << (nbultra/final_missing_cells) << "%)" << endl;
		}
	if (info_data.ReturnIsStandDistanceFileOutput()) {
	    logfile << "  -> Average distance (standardized matrix, ignoring missing cells): ";
	    logfile << setprecision(5) << mean_stand_distance << endl;
	    logfile << "  -> Maximum pair-wise distance (standardized matrix): ";
	    logfile << setprecision(5) << max_stand_distance << endl << endl;
		}
	else {
		logfile << "  -> Average distance (unstandardized matrix, ignoring missing cells): ";
    	logfile << setprecision(5) << mean_distance << endl;
	    logfile << "  -> Maximum pair-wise distance (unstandardized matrix): ";
	    logfile << setprecision(5) << max_distance << endl << endl;
    	}

    /*** Finally, let's print the final average matrix to a nexus file ***/

    output_file(final_matrix, final_stand_matrix);


	/* let's print a matrix containning the number of alleles for each individuals for 
	   each dataset in the logfile */
	
    logfile << endl << "Matrix of the number of alleles found for each organisms for every dataset" << endl << endl;
	for (i=0;i<info_data.ReturnDatasetsNumber();i++) {
		logfile << '\t' << info_data.ReturnDatasetName(i);
		}
	logfile << endl;
	for (i=0; i<org_data.ReturnNOrg();i++) {
		logfile << org_data.ReturnOrganism(i);
		for (j=0; j<info_data.ReturnDatasetsNumber();j++) {
			logfile << '\t' << info_data.ReturnNbAllelesforOrganism(i, j); //TODO: implement this also for the classic pofad
			}
		logfile << endl;
		}
	logfile.close();

    delete [] dist_matrix;
    delete [] stand_dist_matrix;
    delete [] final_matrix;
    delete [] final_stand_matrix;

    return 0;
}




/***** Functions *****/


/************************************************************************************

Function : read_organisms

Description: The function reads the organisms from a file and places them in a vector of organisms

Arguments:   - the name of a file containning the organisms 

************************************************************************************/


void read_organisms(string name)    //When using the batch mode
{
    string stemp;
    char *BigBuffer;

    ifstream inorganisms(name.c_str());

    while (!inorganisms)   //Check if the input file can be found
        {
        cerr << endl << "Error with input organisms file: " << name << endl;
        cerr << "Make sure the file exist and that it contains a valid organisms block" << endl;
        exit(1);
        }

    BigBuffer = new char[MAX_BUFFER_SIZE];
    if (!BigBuffer)
    {
        cerr << "Could not allocate file buffer";
        exit(1);		
    }

    for (i=0; !(inorganisms.eof()); i++)
        {
        if (i >= MAX_BUFFER_SIZE-1)  //have to save room for terminating null
            {
            cerr << "Nexus file exceeds 500k maximum";
            exit(1);
            }
        inorganisms.get(BigBuffer[i]);
        }        

    inorganisms.close();
    BigBuffer[i]='\0';
    bufPtr=BigBuffer;

    if ( bufPtr != NULL )
        {
        while ((aTokenPtr=nextToken()) != "") //This has to be checked
            {
            if (aTokenPtr == "BEGIN")
                {
                stemp = nextToken();   /* get the block name and store in 'stemp'*/
                if (stemp == "")
                    {
                    cerr << "Error reading block name";
                    exit(1);
                    }
                if ((aTokenPtr=nextToken()) == ";") /* pop the terminating semicolon */
                    {
                    if (stemp == "ORGANISMS")
                        {
                        do            /* need to put in error checking in case no DIMENSIONS statement */
		            {
                            aTokenPtr=nextToken();
                            if (aTokenPtr == "DIMENSIONS")  
                                doOrgDimensions();
                            if (aTokenPtr == "ORGLABELS")
                                doOrgLabels();
                            }
                        while ( (aTokenPtr != "END") && (aTokenPtr != "ENDBLOCK") );
                        aTokenPtr=nextToken();
                            if (aTokenPtr != ";") {
                                cerr << "Block not terminated with semicolon";
                                exit(1);
                            }
                        }
                    else  /* token is not a recognized block */
                        {  
                        doUnrecognizedBlock();
                        }
                    }
                stemp.clear();				
                }
            }
        }

    if (info_data.IsVerbose())
        {
        cout << endl << "File containning the organisms: " << name << endl;
        cout << "Number of organisms read: " << org_data.ReturnNOrg() << endl;
        cout << "Content of file: " << endl;
        for (i=0; i < org_data.ReturnNOrg(); i++)
            {
            cout << org_data.ReturnOrganism(i) << '\t';
            }
        cout << endl << endl;
        }

    info_data.IsOrganismsFile();

    delete [] BigBuffer;
}



/************************************************************************************

Function : read_datasets

The function reads the organisms from a file and places them in a vector of organisms

This function is for the batch mode

************************************************************************************/

void read_datasets(string file_name)
{
    string stemp;
    char *BigBuffer;

    ifstream indatasetfile(file_name.c_str());

    while (!indatasetfile)   //Check if the input file can be found
        {
        cerr << endl << "Error in opening the file containning the datasets: " << file_name << "." << endl;
        cerr << "Make sure the name of the file is correct." << endl;
        exit(1);
        }

    BigBuffer = new char[MAX_BUFFER_SIZE];
    if (!BigBuffer)
    {
        cerr << "Could not allocate file buffer";
        exit(1);		
    }

    for (i=0; !(indatasetfile.eof()); i++)
        {
        if (i >= MAX_BUFFER_SIZE-1)  //have to save room for terminating null
            {
            cerr << "Nexus file exceeds 500k maximum";
            exit(1);
            }
        indatasetfile.get(BigBuffer[i]);
        }        

    indatasetfile.close();
    BigBuffer[i]='\0';
    bufPtr=BigBuffer;

    if ( bufPtr != NULL )
        {
        while ((aTokenPtr=nextToken()) != "")
            {
            if (aTokenPtr == "BEGIN")
                {
                stemp = nextToken();   /* get the block name and store in 'stemp'*/
                if (stemp == "")
                    {
                    cerr << "Error reading block name";
                    exit(1);
                    }
                if ((aTokenPtr=nextToken()) == ";") /* pop the terminating semicolon */
                    {
                    if (stemp == "DATASETS")
                        {
                        do
		            		{
                            aTokenPtr=nextToken();
                            if (aTokenPtr == "DIMENSIONS")  
                                doDatasetDimensions();
                            if (aTokenPtr == "DATASETSLABELS")
                                doDatasetLabels();
                            }
                            while ( (aTokenPtr != "END") && (aTokenPtr != "ENDBLOCK") );
                        aTokenPtr=nextToken();
                        if (aTokenPtr != ";") {
                            cerr << "Block not terminated with semicolon";
                            exit(1);
                            }
                        }
                    else  /* token is not a recognized block */
                        {  
                        doUnrecognizedBlock();
                        }
                    }
                stemp.clear();				
                }
            }
        }

    info_data.IsDatasetFile();

    if (info_data.IsVerbose())
        {
        cout << endl << "File containning the name of the datasets: " << file_name << endl;
        cout << "Number of dataset read: " << info_data.ReturnDatasetsNumber() << endl;
        cout << "Content of file: " << endl;
        for (i=0; i < info_data.ReturnDatasetsNumber(); i++)
            {
            cout << info_data.ReturnDatasetName(i) << '\t';
            }
        cout << endl << endl;
        }

    delete [] BigBuffer;

}



/***********************************

Function : doDatasetDimensions

***********************************/

void doDatasetDimensions(void)
{
    while ( (aTokenPtr=nextToken()) != ";")
        {
        if (parse_assignment2("NDATASETS"))
            {
            info_data.GetDatasetsNumber( atoi(LocalToken.c_str()) );
            }
        }
    return;
}


/***********************************

Function : doDatasetLabels

***********************************/


void doDatasetLabels(void)
{
    while ( (aTokenPtr=nextToken2()) != ";") {
        info_data.GetDatasetName(aTokenPtr);
    }
    if ( info_data.ReturnDatasetsSize() < info_data.ReturnDatasetsNumber()) {
        cerr << "Too few Dataset labels";
        exit(1);
    }
    else if ( info_data.ReturnDatasetsSize() > info_data.ReturnDatasetsNumber()) {
        cerr << "Too many Dataset labels";
        exit(1);
    }
    return;
}


/***********************************

Function : consensus

***********************************/


void consensus(string FileName)
{

    nexus_data.InitializeMember();   //initialize the object nexus_data so new informations can be put in
    char *BigBuffer;
    int c;
    long i=0;
    ifstream infile;

    infile.open(FileName.c_str());

    while (!infile)   //Check if the input file can be found
        {
        cerr << "Error in opening nexus file: " << FileName << endl;
        cerr << "Make sure the name is correct and the file exists" << endl << endl;
        exit(1);
        }

    /* Scan the content of the file and transfer it to a buffer */

    BigBuffer = new char[MAX_BUFFER_SIZE];
    if (!BigBuffer)
        {
        cerr << "Could not allocate file buffer";
        exit(1);		
        }

    c = infile.get();

    while (c != EOF)	// have to define c as int so that EOF can be detected
        {
        if (i >= MAX_BUFFER_SIZE-1)  //have to save room for terminating null
            {
            cerr << "Nexus file exceeds 500k maximum";
            exit(1);
            }
        BigBuffer[i]=c;
        ++i;
        c = infile.get();	
        }

    infile.close();

    BigBuffer[i]='\0';

    readNexusFile(BigBuffer);  //Extract the information from the nexus file contained in the buffer

    delete [] BigBuffer;

    org_data.GetNChar(nexus_data.ReturnNChars());
    org_data.MissingChar(nexus_data.ReturnMissingChar());
    org_data.GetGapChar(nexus_data.ReturnGapChar());
    org_data.GetDatatype(nexus_data.ReturnDatatype());

    /*** Compute consensus sequences for each organisms ***/
    
    if ((To_Uppercase(nexus_data.ReturnDatatype()) != "DNA") || (To_Uppercase(nexus_data.ReturnDatatype()) != "NUCLEOTIDE"))
        {
        cerr << "Only DNA characters are allowed";
        exit(1);
        }

    org_data.DeleteCharMatrix();
    org_data.InitDataMatrix();

    string temp;

    AttributeSequencesToOrganisms(1);

    for (i=0; i<org_data.ReturnNOrg(); i++) {
        cout << "  -> Computing consensus: " << ((i+1)*100/org_data.ReturnNOrg()) << " % \r";
        cout.flush();
        //Check if organisms is in the matrix
        if (org_data.ReturnNbAllelesforOrganism(i) == 0) continue;
        for (j=0; j<nexus_data.ReturnNChars(); j++) {
            temp = get_org_char(i,j);
            org_data.AddCharacter(i,temp);
        }
    }
    cout << endl;

}


/*********************************************

Function : read_nexus_matrix

The function read a nexus file that contains a distance matrix and translate the distance matrix in the big matrix containning 
all the datasets.

variables:
    -***dist_matrix: The distance matrix where the distance matrix will be stored
    -***stand_dist_matrix: The standardized distance matrix where the standardized distance matrix will be stored
    -set: The number of the dataset used

*********************************************/

void read_nexus_matrix(double ***dist_matrix, double ***stand_dist_matrix, int set)
{

    string input_file;
    nexus_data.InitializeMember();   //initialize the object nexus_data so new informations can be put in
    char *BigBuffer;
    int c;
    long i=0;
    ifstream infile;

    cout << " Openning file " << info_data.ReturnDatasetName(set) << endl << endl;
    input_file = info_data.ReturnDatasetName(set);
    infile.open(input_file.c_str());

	/* Write in log file */
    logfile << " Openning file " << info_data.ReturnDatasetName(set) << endl << endl;

    while (!infile)   //Check if the input file can be found
        {
        cerr << " Error in opening nexus file: " << input_file << endl;
        cerr << " Make sure the name is correct and the file exists" << endl << endl;
        exit(1);
        }

    /* Scan the content of the file and transfer it to a buffer */

    BigBuffer = new char[MAX_BUFFER_SIZE];
    if (!BigBuffer)
    {
        cerr<<"Could not allocate file buffer";
        exit(1);		
    }

    c = infile.get();

    while (c != EOF)	// have to define c as int so that EOF can be detected
    {
        if (i >= MAX_BUFFER_SIZE-1)  //have to save room for terminating null
        {
            cerr <<"Nexus file exceeds 500k maximum";
            exit(1);
        }
        BigBuffer[i]=c;
        ++i;
        c = infile.get();	
    }

    infile.close();

    BigBuffer[i]='\0';

    readNexusFile(BigBuffer);  //Extract the information from the nexus file contained in the buffer

    delete [] BigBuffer;

    org_data.DeleteCharMatrix();
    AttributeSequencesToOrganisms(set);
    org_data.DeleteDistanceMatrix();
    
    //Obtain consensus sequences for relevant methods
    if ((info_data.ReturnDistanceMethod() == 0) || (info_data.ReturnDistanceMethod() == 1) || (info_data.ReturnDistanceMethod() == 2) || (info_data.ReturnDistanceMethod() == 10)) {
	    org_data.GetNChar(nexus_data.ReturnNChars());
	    org_data.MissingChar(nexus_data.ReturnMissingChar());
	    org_data.GetGapChar(nexus_data.ReturnGapChar());
	    org_data.GetDatatype(nexus_data.ReturnDatatype());
        
	    /*** Compute consensus sequences for each organism ***/

	    if ( (To_Uppercase(nexus_data.ReturnDatatype()) != "DNA") && (To_Uppercase(nexus_data.ReturnDatatype()) != "NUCLEOTIDE")) {
            cerr << "Only DNA characters are allowed";
            exit(1);
            }
        
	    org_data.InitDataMatrix();
            	
	    string temp;
	    i=0;
	    for (i=0; i<org_data.ReturnNOrg(); i++) {
	        info_data.EnterNbAllelesforOrganism(i, set, org_data.ReturnNbAllelesforOrganism(i));
            cout << "  -> Computing consensus: " << ((i+1)*100/org_data.ReturnNOrg()) << " % \r";
            cout.flush();
            //cout << org_data.ReturnOrganism(i) << " has " << org_data.ReturnNbAllelesforOrganism(i) << " alleles" << endl;
            /*cout << org_data.ReturnOrganism(i) << " has alleles ";
            for (j=0;j<org_data.ReturnNbAllelesforOrganism(i);j++) {
                cout << org_data.ReturnAlleleFromOrganism(i,j) << " ";
            }
            cout << endl;
            */
            //Check if organisms is in the matrix
            if (org_data.ReturnNbAllelesforOrganism(i) == 0) continue;
	        for (j=0; j<nexus_data.ReturnNChars(); j++) {
	            temp = get_org_char(i,j);
	            org_data.AddCharacter(i,temp);
	            }
	        }
        cout << endl;
	
	    /*** Calculate the distance between the consensus sequences ***/
	
	    org_data.InitDistMatrix();
	    
	    if (info_data.ReturnDistanceMethod() == 0) org_data.CalculateDistances("genpofad");
	    if (info_data.ReturnDistanceMethod() == 1) org_data.CalculateDistances("matchstates");
	    if (info_data.ReturnDistanceMethod() == 2) org_data.CalculateDistances("mrca");
	    if (info_data.ReturnDistanceMethod() == 10) org_data.CalculateDistances("PP");
		} /* end if (polymorphisms SNP based methods) */

	if (info_data.ReturnDistanceMethod() == 4 || info_data.ReturnDistanceMethod() == 9) {
	    org_data.GetNChar(nexus_data.ReturnNChars());
	    org_data.MissingChar(nexus_data.ReturnMissingChar());
	    org_data.GetGapChar(nexus_data.ReturnGapChar());
	    org_data.GetDatatype(nexus_data.ReturnDatatype());
	    org_data.DeleteFRQMatrix();
		org_data.InitFRQMatrix();
		org_data.InitDistMatrix();
		if (info_data.ReturnDistanceMethod() == 4) org_data.CalculateDistances("FRQ");  //FRQ distance of Goker and Grimm
        if (info_data.ReturnDistanceMethod() == 9) org_data.CalculateDistances("nei"); //Nei's genetic distance
		}
	
	/*** Now we can construct the matrix of pairwise organisms distances ***/
    
    double raw_distance;
    
    if (info_data.ReturnDistanceMethod() == 3) //Original POFAD distance
    	{
        if (!nexus_data.ReturnIsDistanceMatrix()) {
            nexus_data.InitMatrix();
            nexus_data.ComputeDistanceMatrix();
			}
	    for (i=0; i < org_data.ReturnNOrg(); i++)
		    {
            cout << "  -> Calculating distances: " << ((i+1)*100/org_data.ReturnNOrg()) << " % \r";
            cout.flush();
	        for (j=0; j <=i ; j++)
		        {
				raw_distance = get_nexus_distance(i, j);
				if (info_data.ReturnIsJukesCantor()) { 		/* Jukes and Cantor distance */
					if (raw_distance == -999) { /* No data */
		            	dist_matrix[set][i][j] = raw_distance;
		            	dist_matrix[set][j][i] = raw_distance;						
						}
					else {
						dist_matrix[set][i][j] = (-(3./4)*log(1-((4./3)*raw_distance)));
						dist_matrix[set][j][i] = dist_matrix[set][i][j];
						}
					}
				else if (!info_data.ReturnIsJukesCantor()){ /* Hamming distance */	
	            	dist_matrix[set][i][j] = raw_distance;
	            	dist_matrix[set][j][i] = raw_distance;
	            	}
		        }
		    }
		cout << endl;
		}
    else if (info_data.ReturnDistanceMethod() == 5) //Phylogenetic Bray-Curtis distance
    	{
		if (!nexus_data.ReturnIsDistanceMatrix()) {
			nexus_data.InitMatrix();
			nexus_data.ComputeDistanceMatrix();
			}
	    for (i=0; i < org_data.ReturnNOrg(); i++) {
            cout << "  -> Calculating organisms distances: " << ((i+1)*100/org_data.ReturnNOrg()) << " % \r";
            cout.flush();
	        for (j=0; j <= i ; j++) {
                if (i==j) {
                    dist_matrix[set][i][j] = 0;
                    continue;
                    }
				raw_distance = get_PBC_distance(i, j);
				if (info_data.ReturnIsJukesCantor()) { 		/* Jukes and Cantor distance */
					if (raw_distance == -999) { /* No data */
		            	dist_matrix[set][i][j] = raw_distance;
		            	dist_matrix[set][j][i] = raw_distance;						
						}
					else {
						dist_matrix[set][i][j] = (-(3./4)*log(1-((4./3)*raw_distance)));
						dist_matrix[set][j][i] = dist_matrix[set][i][j];
						}
					}
				else if (!info_data.ReturnIsJukesCantor()){ /* Hamming distance */	
	            	dist_matrix[set][i][j] = raw_distance;
	            	dist_matrix[set][j][i] = raw_distance;
	            	}
		        }
		    }
		cout << endl;
		}
    else if (info_data.ReturnDistanceMethod() == 8) //MIN distance of Goker and Grimm
    	{
		if (!nexus_data.ReturnIsDistanceMatrix()) {
			nexus_data.InitMatrix();
			nexus_data.ComputeDistanceMatrix();
			}
    	i=0;
	    for (i=0; i < org_data.ReturnNOrg(); i++)
		    {
            cout << "  -> Calculating organisms distances: " << ((i+1)*100/org_data.ReturnNOrg()) << " % \r";
            cout.flush();
	        for (j=0; j <= i ; j++)
		        {
				raw_distance = get_MIN_distance(i, j);
				if (info_data.ReturnIsJukesCantor()) { 		/* Jukes and Cantor distance */
					if (raw_distance == -999) { /* No data */
		            	dist_matrix[set][i][j] = raw_distance;
		            	dist_matrix[set][j][i] = raw_distance;						
						}
					else {
						dist_matrix[set][i][j] = (-(3./4)*log(1-((4./3)*raw_distance)));
						dist_matrix[set][j][i] = dist_matrix[set][i][j];
						}
					}
				else if (!info_data.ReturnIsJukesCantor()){ /* Hamming distance */	
	            	dist_matrix[set][i][j] = raw_distance;
	            	dist_matrix[set][j][i] = raw_distance;
	            	}
		        }
		    }
		cout << endl;
		}
    else if (info_data.ReturnDistanceMethod() == 7) //pofadX
    	{
		if (!nexus_data.ReturnIsDistanceMatrix()) {
			nexus_data.InitMatrix();
			nexus_data.ComputeDistanceMatrix();
			}
	    for (i=0; i < org_data.ReturnNOrg(); i++)
		    {
            cout << "  -> Calculating organisms distances: " << ((i+1)*100/org_data.ReturnNOrg()) << " % \r";
            cout.flush();
	        for (j=0; j <= i ; j++)
		        {
				raw_distance = get_pofadX_distance(i, j);
				if (info_data.ReturnIsJukesCantor()) { 		// Jukes and Cantor distance
					if (raw_distance == -999) { // No data
		            	dist_matrix[set][i][j] = raw_distance;
		            	dist_matrix[set][j][i] = raw_distance;						
						}
					else {
						dist_matrix[set][i][j] = (-(3./4)*log(1-((4./3)*raw_distance)));
						dist_matrix[set][j][i] = dist_matrix[set][i][j];
						}
					}
				else if (!info_data.ReturnIsJukesCantor()){ // Hamming distance
	            	dist_matrix[set][i][j] = raw_distance;
	            	dist_matrix[set][j][i] = raw_distance;
	            	}
		        }
		    }
		cout << endl;
		}
    else if ( info_data.ReturnDistanceMethod() != 4 ) {   /* Either genpofad, matchstates, mrca, or PP */
	    for (i=0; i < org_data.ReturnNOrg(); i++) {
	        for (j=0; j <= i; j++) {
				raw_distance = org_data.ReturnDist(i, j);
				if (info_data.ReturnIsJukesCantor()) { 		/* Jukes and Cantor distance */
					if (raw_distance == -999) { /* No data */
		            	dist_matrix[set][i][j] = raw_distance;
		            	dist_matrix[set][j][i] = raw_distance;						
						}
					else {
						dist_matrix[set][i][j] = (-(3./4)*log(1-((4./3)*raw_distance)));
						dist_matrix[set][j][i] = dist_matrix[set][i][j];
						}
					}
				else if (!info_data.ReturnIsJukesCantor()){ /* Hamming distance */	
	            	dist_matrix[set][i][j] = raw_distance;
	            	dist_matrix[set][j][i] = raw_distance;
	            	}
	            }
		    }
		}

    /*** Calculate the average distance of the matrix of organisms ***/

    int matrix_size = 0;
    double max_distance = 0;         //To find the maximum distance of the matrix
    double mean_distance;
    double total_distance = 0;
    double total_cells = 0;
    double missing_distances = 0;
    double percent_missing;

    for (i=0; i < org_data.ReturnNOrg(); i++)
        {
        for (j=0; j < i; j++)
            {
            total_cells += 1;
            if ( dist_matrix[set][i][j] == -999 )
               {
               missing_distances += 1;
               continue;
               }
            total_distance += dist_matrix[set][i][j];
            matrix_size++;
            if (dist_matrix[set][i][j] > max_distance) {
                max_distance = dist_matrix[set][i][j];
                }
            }
        }

    mean_distance = (total_distance / matrix_size ) * 100;  //Will give the mean distance in %
    percent_missing = (missing_distances / total_cells) * 100;

    cout << "  -> Average distance of organisms for dataset " << info_data.ReturnDatasetName(set) << " is: ";
    cout << setprecision(2) << mean_distance << "%" << endl;
    cout << "  -> Maximum pair-wise distance between organisms is: ";
    cout << setprecision(5) << max_distance << endl;
    cout << "  -> Percentage of organisms distances missing in the dataset: ";
    cout << setprecision(2) << percent_missing << "%" << endl << endl;


	/* Write in log file */
	logfile << "  -> Number of taxa in matrix = " << nexus_data.ReturnNTax() << endl;
	logfile << "  -> Number of characters in the matrix = " << nexus_data.ReturnNChars() << endl;
    logfile << "  -> Average distance of organisms for dataset " << info_data.ReturnDatasetName(set) << " is: ";
    logfile << setprecision(2) << mean_distance << "%" << endl;
    logfile << "  -> Maximum pair-wise distance between organisms is: ";
    logfile << setprecision(5) << max_distance << endl;
    logfile << "  -> Percentage of organisms distances missing in the dataset: ";
    logfile << setprecision(2) << percent_missing << "%" << endl << endl << endl;

	
    /*** Construct the standardized matrix of organisms ***/

    #ifdef DEBUG
    cout << endl << "standardized distances: " << endl;
    #endif

    for (i=0; i < org_data.ReturnNOrg(); i++)
        {
        for (j=0; j < org_data.ReturnNOrg(); j++)
            {
            if ( dist_matrix[set][i][j] == -999 )
               {
               stand_dist_matrix[set][i][j] = dist_matrix[set][i][j];
               continue;
               }
            stand_dist_matrix[set][i][j] = ( dist_matrix[set][i][j] / max_distance );

            #ifdef DEBUG
            cout << stand_dist_matrix[set][i][j] << '\t';
            #endif

            }
        }


    #ifdef DEBUG
        cout << " distance matrices of dataset " << info_data.ReturnDatasetName(set) << " constructed !" << endl;
    #endif
}



/**************************************************************************************************************************
 
Function : AttributeSequencesToOrganisms
 
Read the taxa names and attributes them to an organism from the organism file. Organisms need to start the name of the
allele and is separated from the rest of the allele name by a '_'. ex: Rnitida761_A, Rnitida761_B, etc.
 
variables :
 - dataset number being treated
 
**************************************************************************************************************************/
 
void AttributeSequencesToOrganisms(int set)
{
    int i,x;
    string temp,an_organism;
    org_data.InitializeAllelesInOrganisms();
    for (i=0; i < nexus_data.ReturnNTax(); i++) {
        an_organism="";
        temp = nexus_data.ReturnTaxa(i);
        for (x=0;x<temp.length();x++) {
            if (temp[x] == '_') break;
            else an_organism+=temp[x];
        }
        if (org_data.AddAllele(an_organism,nexus_data.ReturnTaxa(i)) == 1) {
            cout << " Organism " << an_organism << " in dataset " << info_data.ReturnDatasetName(set) << " is not present in organisms list." << endl;
            exit(1);
        }
    }
}


/**************************************************************************************************************************

Function : get_position_nexus

The function get_position returns the position "i" in a haplotype vector (vector[i]) that has been read from the nexus file

variables : 
	- the name of the string to search for

**************************************************************************************************************************/

int get_position_nexus(string the_string)
{
    int i;
    vector<string> a_vector;
    for (i=0; i < nexus_data.ReturnNTax(); i++)
    {
        a_vector.push_back(nexus_data.ReturnTaxa(i));
    }

    #ifdef DEBUG
    cout << "Vec" << '\t';
    #endif

    int position;
    vector<string>::iterator an_iterator;
    an_iterator = find(a_vector.begin(), a_vector.end(), the_string);
    if (an_iterator == a_vector.end())
    {
        return -1;
    }
    else
    {
        a_vector.erase(an_iterator, a_vector.end());
        position = a_vector.size();
        #ifdef DEBUG
        cout << "Pos: " << position << '\t';
        #endif
        return (position);
    }
}



/************************************************************************************************************

Function: get_nexus_distance

************************************************************************************************************/

double get_nexus_distance(int org1, int org2)
{

    int i,j;
    double distance;  //We define the double distance that will be return by the function

    //If we compare the individual with itself
    if (org1 == org2) return (0);
    
    //If one fo the organism is absent from the matrix
    if (org_data.ReturnNbAllelesforOrganism(org1)==0 || org_data.ReturnNbAllelesforOrganism(org2)==0) return(-999);

    //If one organism has more than 2 alleles
    if (org_data.ReturnNbAllelesforOrganism(org1)>2) {
        cerr << " Organism " << org_data.ReturnOrganism(org1) << " has more than two alleles." << endl;
        cerr << " This is not permitted with the standard pofad method." << endl << endl;
        exit(1);
    }
    if (org_data.ReturnNbAllelesforOrganism(org2)>2) {
        cerr << " Organism " << org_data.ReturnOrganism(org2) << " has more than two alleles." << endl;
        cerr << " This is not permitted with the standard pofad method." << endl << endl;
        exit(1);
    }
    
    //If both organisms have a single haplotype
    if (org_data.ReturnNbAllelesforOrganism(org1)==1 || org_data.ReturnNbAllelesforOrganism(org2)==1) {
        distance = nexus_data.ReturnDist(get_position_nexus(org_data.ReturnAlleleFromOrganism(org1,0))
                                         , get_position_nexus(org_data.ReturnAlleleFromOrganism(org2,0)));
        return distance;
    }
    
    //If both organisms have two haplotypes
    if (org_data.ReturnNbAllelesforOrganism(org1)==2 || org_data.ReturnNbAllelesforOrganism(org2)==2) {
        double distance1;
        double distance2;
        distance1 = (nexus_data.ReturnDist(get_position_nexus(org_data.ReturnAlleleFromOrganism(org1,0))
                                        , get_position_nexus(org_data.ReturnAlleleFromOrganism(org2,0))))
                     + (nexus_data.ReturnDist(get_position_nexus(org_data.ReturnAlleleFromOrganism(org1,1))
                                        , get_position_nexus(org_data.ReturnAlleleFromOrganism(org2,1))));
        distance2 = (nexus_data.ReturnDist(get_position_nexus(org_data.ReturnAlleleFromOrganism(org1,0))
                                        , get_position_nexus(org_data.ReturnAlleleFromOrganism(org2,1))))
                    + (nexus_data.ReturnDist(get_position_nexus(org_data.ReturnAlleleFromOrganism(org1,1))
                                        , get_position_nexus(org_data.ReturnAlleleFromOrganism(org2,0))));
        if (distance1 < distance2) distance = (distance1/2);
        else distance = distance2/2;
        return distance;
    }
    
    // One organism has one allele and the other 2
    else {
        distance = .0;
        for (i=0;i<org_data.ReturnNbAllelesforOrganism(org1);i++) {
            for (j=0;j<org_data.ReturnNbAllelesforOrganism(org2);j++) {
                distance += nexus_data.ReturnDist(get_position_nexus(org_data.ReturnAlleleFromOrganism(org1,i)),
                                                  get_position_nexus(org_data.ReturnAlleleFromOrganism(org2,j)));
            }
        }
        distance = distance / 2;
        return distance;
    }
}

//------------------------------------------
// get_PBC_distance
//------------------------------------------

double get_PBC_distance(int org1, int org2)
{
	int x,j;
    //If we compare the individual with itself
    if (org1 == org2) return (0);
    //If one fo the organism is absent from the matrix
    if (org_data.ReturnNbAllelesforOrganism(org1)==0 || org_data.ReturnNbAllelesforOrganism(org2)==0) return(-999);

    double distance=0.;
    double tempdist1,tempdist2;

    for (i=0;i<org_data.ReturnNbAllelesforOrganism(org1);i++) {
        tempdist2=9999999.0;
        for (j=0;j<org_data.ReturnNbAllelesforOrganism(org2);j++) {
            tempdist1 = nexus_data.ReturnDist(get_position_nexus(org_data.ReturnAlleleFromOrganism(org1,i)), 
                                              get_position_nexus(org_data.ReturnAlleleFromOrganism(org2,j)));
            if (tempdist1 < tempdist2) tempdist2=tempdist1;
        }
        distance+=tempdist2;
    }
    for (j=0;j<org_data.ReturnNbAllelesforOrganism(org2);j++) {
        tempdist2=9999999.0;
        for (i=0;i<org_data.ReturnNbAllelesforOrganism(org1);i++) {
            tempdist1 = nexus_data.ReturnDist(get_position_nexus(org_data.ReturnAlleleFromOrganism(org1,i)), 
                                              get_position_nexus(org_data.ReturnAlleleFromOrganism(org2,j)));
            if (tempdist1 < tempdist2) tempdist2=tempdist1;
        }
        distance+=tempdist2;
    }    
    distance = distance/(org_data.ReturnNbAllelesforOrganism(org1)+org_data.ReturnNbAllelesforOrganism(org2));
    return distance;
}


//------------------------------------------
// get_MIN_distance
//------------------------------------------

double get_MIN_distance(int org1, int org2)
{
	int x,j;
    //If we compare the individual with itself
    if (org1 == org2) return (0);
    //If one fo the organism is absent from the matrix
    if (org_data.ReturnNbAllelesforOrganism(org1)==0 || org_data.ReturnNbAllelesforOrganism(org2)==0) return(-999);
    
    double distance=0.;
    double tempdist1,tempdist2;

    tempdist2=9999999.0;
    for (i=0;i<org_data.ReturnNbAllelesforOrganism(org1);i++) {
        for (j=0;j<org_data.ReturnNbAllelesforOrganism(org2);j++) {
            tempdist1 = nexus_data.ReturnDist(get_position_nexus(org_data.ReturnAlleleFromOrganism(org1,i)), 
                                              get_position_nexus(org_data.ReturnAlleleFromOrganism(org2,j)));
            if (tempdist1 < tempdist2) tempdist2=tempdist1;
        }
    }
    distance=tempdist2;
    return(distance);
}


//------------------------------------------
// get_pofadX_distance
//------------------------------------------

//TODO: modify this function

double get_pofadX_distance(int m, int n)
{
	int x,y,k;
    //If we compare the individual with itself
    if (m == n) return (0);

	//vector that contain the names of the alleles in each individual
	vector<string> IndividualOne,IndividualTwo;	
	
    string string_1, string_1a, string_1b, string_2, string_2a, string_2b;
    double distance;  //We define the double distance that will be return by the function
	
    //Following are the different possibilities for the haplotype names to be found in the haplotype vector
    string_1 = To_Uppercase(org_data.ReturnOrganism(m));
    string_1a = string_1 + "_A";
    string_1b = string_1 + "_B";
    string_2 = To_Uppercase(org_data.ReturnOrganism(n));
    string_2a = string_2 + "_A";
    string_2b = string_2 + "_B";

    //If one of the organisms is absent from the dataset, return the flag "-999".
    if ( ( (get_position_nexus(string_1) == -1) && (get_position_nexus(string_1a) == -1) ) || ( (get_position_nexus(string_2) == -1) && (get_position_nexus(string_2a) == -1) ) )
    {
    return (-999);
    }
    //If both organisms compared have only a single haplotype
    if ( (get_position_nexus(string_1a) == -1) && (get_position_nexus(string_2a) == -1) )
    {
        distance = nexus_data.ReturnDist(get_position_nexus(string_1), get_position_nexus(string_2));
        return distance;
    }

	IndividualOne.push_back(string_1a);
	IndividualTwo.push_back(string_2a);

	//Calculate the number of alleles in each individual and get the names of the alleles
	int AllelesIn1=0;
	int AllelesIn2=0;
	string temp = string_1a;
	while (get_position_nexus(temp) != -1) {
		AllelesIn1++;
		IndividualOne.push_back(temp);
		if (AllelesIn1 == 1) temp=(string_1 + "_B");
		else if (AllelesIn1 == 2) {temp=(string_1 + "_C");}
		else if (AllelesIn1 == 3) {temp=(string_1 + "_D");}
		else if (AllelesIn1 == 4) {temp=(string_1 + "_E");}
		else if (AllelesIn1 == 5) {temp=(string_1 + "_F");}
		else if (AllelesIn1 == 6) {temp=(string_1 + "_G");}
		else if (AllelesIn1 == 7) {temp=(string_1 + "_H");}
		else if (AllelesIn1 == 8) {temp=(string_1 + "_I");}
		else if (AllelesIn1 == 9) {temp=(string_1 + "_J");}
		else if (AllelesIn1 == 10) {temp=(string_1 + "_K");}
		else if (AllelesIn1 == 11) {temp=(string_1 + "_L");}
		else if (AllelesIn1 == 12) {temp=(string_1 + "_M");}
		else if (AllelesIn1 == 13) {temp=(string_1 + "_N");}
		else if (AllelesIn1 == 14) {temp=(string_1 + "_O");}
		else if (AllelesIn1 == 15) {temp=(string_1 + "_P");}
		else if (AllelesIn1 == 16) {temp=(string_1 + "_Q");}
		else if (AllelesIn1 == 17) {temp=(string_1 + "_R");}
		else if (AllelesIn1 == 18) {temp=(string_1 + "_S");}
		else if (AllelesIn1 == 19) {temp=(string_1 + "_T");}
		else if (AllelesIn1 == 20) {temp=(string_1 + "_U");}
		else if (AllelesIn1 == 21) {temp=(string_1 + "_V");}
		else if (AllelesIn1 == 22) {temp=(string_1 + "_W");}
		else if (AllelesIn1 == 23) {temp=(string_1 + "_X");}
		else if (AllelesIn1 == 24) {temp=(string_1 + "_Y");}
		else if (AllelesIn1 == 25) {temp=(string_1 + "_Z");}
		else {cerr << " Too many alleles (<26) in the dataset" << string_1 << "." << endl;exit(1);}
		}
	temp = string_2a;
	while (get_position_nexus(temp) != -1) {
		AllelesIn2++;
		IndividualTwo.push_back(temp);
		if (AllelesIn2 == 1) temp=(string_2 + "_B");
		else if (AllelesIn2 == 2) {temp=(string_2 + "_C");}
		else if (AllelesIn2 == 3) {temp=(string_2 + "_D");}
		else if (AllelesIn2 == 4) {temp=(string_2 + "_E");}
		else if (AllelesIn2 == 5) {temp=(string_2 + "_F");}
		else if (AllelesIn2 == 6) {temp=(string_2 + "_G");}
		else if (AllelesIn2 == 7) {temp=(string_2 + "_H");}
		else if (AllelesIn2 == 8) {temp=(string_2 + "_I");}
		else if (AllelesIn2 == 9) {temp=(string_2 + "_J");}
		else if (AllelesIn2 == 10) {temp=(string_2 + "_K");}
		else if (AllelesIn2 == 11) {temp=(string_2 + "_L");}
		else if (AllelesIn2 == 12) {temp=(string_2 + "_M");}
		else if (AllelesIn2 == 13) {temp=(string_2 + "_N");}
		else if (AllelesIn2 == 14) {temp=(string_2 + "_O");}
		else if (AllelesIn2 == 15) {temp=(string_2 + "_P");}
		else if (AllelesIn2 == 16) {temp=(string_2 + "_Q");}
		else if (AllelesIn2 == 17) {temp=(string_2 + "_R");}
		else if (AllelesIn2 == 18) {temp=(string_2 + "_S");}
		else if (AllelesIn2 == 19) {temp=(string_2 + "_T");}
		else if (AllelesIn2 == 20) {temp=(string_2 + "_U");}
		else if (AllelesIn2 == 21) {temp=(string_2 + "_V");}
		else if (AllelesIn2 == 22) {temp=(string_2 + "_W");}
		else if (AllelesIn2 == 23) {temp=(string_2 + "_X");}
		else if (AllelesIn2 == 24) {temp=(string_2 + "_Y");}
		else if (AllelesIn2 == 25) {temp=(string_2 + "_Z");}
		else {cerr << " Too many alleles (<26) in individual " << string_2 << "." << endl;exit(1);}
		}
	if (AllelesIn1==0) {IndividualOne.push_back(string_1);AllelesIn1++;}
	if (AllelesIn2==0) {IndividualTwo.push_back(string_2);AllelesIn2++;}

	//Bool vector that indicate whether the alelles have been picked
	vector<bool> IsPickedAlleleInOne (AllelesIn1,false);
	vector<bool> IsPickedAlleleInTwo (AllelesIn2,false);

	//Calculate minimum distance between the alleles of the two species
	double tempdist1,tempdist2;
	int TempAllele1,TempAllele2;
	distance=0;
	//TODO: The following is not correct!
	if (AllelesIn1 >= AllelesIn2) {
		for(k=0;k < AllelesIn2;k++) {
			tempdist2=9999999.0;
			for (x=0;x<AllelesIn1;x++) {
				if (IsPickedAlleleInOne[x]) continue;
				for (y=0;y<AllelesIn2;y++) {
					if (IsPickedAlleleInTwo[y]) continue;
					tempdist1 = nexus_data.ReturnDist(get_position_nexus(IndividualOne[x]), get_position_nexus(IndividualTwo[y]));
					if (tempdist1 < tempdist2) {tempdist2=tempdist1; TempAllele1=x; TempAllele2=y;}
					}
				}
			IsPickedAlleleInOne[TempAllele1]=true;
			IsPickedAlleleInTwo[TempAllele2]=true;
			distance += tempdist2;
			}
		for(k=0;k<(AllelesIn1-AllelesIn2);k++) {
			tempdist2=9999999.0;
			for (x=0;x<AllelesIn1;x++) {
				if (IsPickedAlleleInOne[x]) continue;
				for (y=0;y<AllelesIn2;y++) {
					tempdist1 = nexus_data.ReturnDist(get_position_nexus(IndividualOne[x]), get_position_nexus(IndividualTwo[y]));
					if (tempdist1 < tempdist2) {tempdist2=tempdist1; TempAllele1=x;}
					}
				}
			IsPickedAlleleInOne[TempAllele1]=true;
			distance += tempdist2;
			}
		distance = (distance / AllelesIn1);
		return(distance);
		} // if (AllelesIn1 >= AllelesIn2)
	else {
		for(k=0;k < AllelesIn1;k++) {
			tempdist2=9999999.0;
			for (x=0;x<AllelesIn2;x++) {
				if (IsPickedAlleleInTwo[x]) continue;
				for (y=0;y<AllelesIn1;y++) {
					if (IsPickedAlleleInOne[y]) continue;
					tempdist1 = nexus_data.ReturnDist(get_position_nexus(IndividualTwo[x]), get_position_nexus(IndividualOne[y]));
					if (tempdist1 < tempdist2) {tempdist2=tempdist1; TempAllele2=x; TempAllele1=y;}
					}
				}
			IsPickedAlleleInOne[TempAllele1]=true;
			IsPickedAlleleInTwo[TempAllele2]=true;
			distance += tempdist2;
			}
		for(k=0;k<(AllelesIn2-AllelesIn1);k++) {
			tempdist2=9999999.0;
			for (x=0;x<AllelesIn2;x++) {
				if (IsPickedAlleleInTwo[x]) continue;
				for (y=0;y<AllelesIn1;y++) {
					tempdist1 = nexus_data.ReturnDist(get_position_nexus(IndividualTwo[x]), get_position_nexus(IndividualOne[y]));
					if (tempdist1 < tempdist2) {tempdist2=tempdist1; TempAllele2=x;}
					}
				}
			IsPickedAlleleInTwo[TempAllele2]=true;
			distance += tempdist2;
			}
		distance = (distance / AllelesIn2);
		return(distance);
		}
}


/*********************************************

Function : get_org_char

*********************************************/


char get_org_char(int organism_nb, int position) {
    int i;
    string an_allele;
    string character = "";
    string temp_character = "";
    
    for (i=0;i<org_data.ReturnNbAllelesforOrganism(organism_nb);i++) {
        an_allele = To_Uppercase(org_data.ReturnAlleleFromOrganism(organism_nb,i));
        if (i==0) {
            character = nexus_data.ReturnChar( get_position_nexus(an_allele), position);
        }
        else {
            temp_character = char_consensus(character[0],nexus_data.ReturnChar( get_position_nexus(an_allele), position));
            character = temp_character;            
        }
    }
    return character[0];
}



/*********************************************

Function : char_consensus

Returns the consensus base of two nucleic acids

Variables: - two nucleic acids


*********************************************/


string char_consensus(char char1, char char2)
{

    string firstchar;
    string secondchar;

    firstchar = char1;
    secondchar = char2;

	int c;
	int PositionOfFirstChar=-999;
	int PositionOfSecondChar=-999;

    // If both characters are identicals
    if (firstchar == secondchar) return firstchar;

	if (firstchar == "?") return secondchar;
	if (secondchar == "?") return firstchar;

    // Handling of gaps if treated as missing data
	if (info_data.ReturnGapHandling() == 0) {
	    if (firstchar == "-" ) return secondchar;
    	if (secondchar == "-" ) return firstchar;
		}

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

    // Handling 0-1 characters
    if (firstchar == "0"  || secondchar == "1") {
    	cerr << "Binary characters are not allowed! Try converting them to nucleotide data...";
        exit(1);
    	}
	if (PositionOfFirstChar == -999) {
        cout << endl << endl << "Unknown character: " << firstchar << endl;
        cerr << "Look at the manual for allowed characters";
        exit(1);
        }
	if (PositionOfSecondChar == -999) {
        cout << endl << endl << "Unknown character: " << secondchar << endl;
        cerr << "Look at the manual for allowed characters";
        exit(1);
        }
	
	return Consensus[PositionOfFirstChar][PositionOfSecondChar];
}



/*****************************************************************************

Function: output_file

Description: This function outputs a distance matrix in a file in nexus format

Variables:
    - the name of the 2-dimensions matrix to be output

******************************************************************************/

void output_file(double **final_matrix, double **final_stand_matrix)
{
    string outputfile;

    if (info_data.ReturnIsOutputFile())
        {
        outputfile = info_data.ReturnOutputFile();
        }
    else
        {
        cout << endl << "Enter the name (no space) of the output file: ";
        cin >> outputfile;
        }

    ofstream outfile(outputfile.c_str());

    if (outfile.bad())
        {
        cout << "Error in opening output file" << endl;
        cout << "Verify that there are no white space in the name" << endl;
        }

	logfile << " Output file written to: " <<  outputfile << endl << endl << endl;

    if (info_data.ReturnIsStandDistanceFileOutput()) // output is standardized
        {
        outfile << "#nexus" << endl << endl;
        outfile << "[!********************************************************************]" << endl;
        outfile << "[!*                                                                  *]" << endl;
        outfile << "[!* File generated by POFAD "<<VERSION<<"                                     *]" << endl;
        outfile << "[!*                                                                  *]" << endl;

    if (info_data.ReturnDistanceMethod() == 0) outfile << "[!* Distance matrix generated using the genpofad algorithm           *]" << endl;
    if (info_data.ReturnDistanceMethod() == 1) outfile << "[!* Distance matrix generated using the matchstates distance         *]" << endl;
    if (info_data.ReturnDistanceMethod() == 2) outfile << "[!* Distance matrix generated using the mrca distance                *]" << endl;
    if (info_data.ReturnDistanceMethod() == 3) outfile << "[!* Distance matrix generated using the original POFAD algorithm     *]" << endl;
    if (info_data.ReturnDistanceMethod() == 4) outfile << "[!* Distance matrix generated using the FRQ distance                 *]" << endl;
    if (info_data.ReturnDistanceMethod() == 5) outfile << "[!* Distance matrix generated using the PBC distance                 *]" << endl;
    if (info_data.ReturnDistanceMethod() == 8) outfile << "[!* Distance matrix generated using the MIN distance                 *]" << endl;
    if (info_data.ReturnDistanceMethod() == 9) outfile << "[!* Distance matrix generated using the Nei distance                 *]" << endl;
    if (info_data.ReturnDistanceMethod() == 10) outfile << "[!* Distance matrix generated using the PP distance               *]" << endl;

        outfile << "[!*                                                                  *]" << endl;
        outfile << "[!* refs.: -> S. Joly and A. Bruneau (2006) Syst. Biol. 55:623-636   *]" << endl;
        outfile << "[!*        -> S. Joly D. Bryant, and P.J. Lockhart. Unpublished data *]" << endl;
        outfile << "[!*                                                                  *]" << endl;
        outfile << "[!* Final matrix obtained from standardized distance matrices of     *]" << endl;
        outfile << "[!* organisms                                                        *]" << endl;
        outfile << "[!*                                                                  *]" << endl;
        outfile << "[!********************************************************************]" << endl;
        outfile << endl << endl;
        outfile << "Begin taxa;" << endl;
        outfile << "    Dimensions ntax=" << org_data.ReturnNOrg() << ";" << endl;
        outfile << "    Taxlabels" << endl;

        for (i=0; i < org_data.ReturnNOrg(); i++)
            {
            outfile << org_data.ReturnOrganism(i) << endl;
            }

        outfile << ";" << endl;
        outfile << "end;" << endl << endl;
        outfile << "Begin distances;" << endl;
        outfile << "    Format triangle=both labels diagonal;" << endl;
        outfile << "matrix" << endl << endl;

	    //Claculate the length of the longest taxa to make nice output

	    int max_length=0;
	    int current_length;

 	   	for (i=0; i < org_data.ReturnNOrg(); i++)
  	  		{
   	     	current_length=org_data.ReturnLengthofOrg(i);
  	      	if (current_length > max_length)
            	{
            	max_length = current_length;
            	}
    		}

    	int size_of_current_taxa;

    	for (i=0; i < org_data.ReturnNOrg(); i++)
        	{
            outfile << org_data.ReturnOrganism(i); 
        	size_of_current_taxa = org_data.ReturnLengthofOrg(i);
        	for (j=0; j < ( (max_length) - size_of_current_taxa); j++)
            	{
            	outfile << " ";
            	} 
			outfile << '\t';
            for (j=0; j < org_data.ReturnNOrg(); j++)
                {
                outfile << setprecision(PRECISION) << final_stand_matrix[i][j] << '\t';
                }
            outfile << endl;
        	}

        outfile << ";" << endl;
        outfile << "end;" << endl;
        }

    else // The output matrix is not standardized
        {
        outfile << "#nexus" << endl << endl;
        outfile << "[!********************************************************************]" << endl;
        outfile << "[!*                                                                  *]" << endl;
        outfile << "[!* File generated by POFAD "<<VERSION<<"                                     *]" << endl;
        outfile << "[!*                                                                  *]" << endl;

        if (info_data.ReturnDistanceMethod() == 0) outfile << "[!* Distance matrix generated using the genpofad algorithm           *]" << endl;
        if (info_data.ReturnDistanceMethod() == 1) outfile << "[!* Distance matrix generated using the matchstates distance         *]" << endl;
        if (info_data.ReturnDistanceMethod() == 2) outfile << "[!* Distance matrix generated using the mrca distance                *]" << endl;
        if (info_data.ReturnDistanceMethod() == 3) outfile << "[!* Distance matrix generated using the original POFAD algorithm     *]" << endl;
        if (info_data.ReturnDistanceMethod() == 4) outfile << "[!* Distance matrix generated using the FRQ distance                 *]" << endl;
        if (info_data.ReturnDistanceMethod() == 5) outfile << "[!* Distance matrix generated using the PBC distance                 *]" << endl;
        if (info_data.ReturnDistanceMethod() == 8) outfile << "[!* Distance matrix generated using the MIN distance                 *]" << endl;
        if (info_data.ReturnDistanceMethod() == 9) outfile << "[!* Distance matrix generated using the Nei distance                 *]" << endl;
        if (info_data.ReturnDistanceMethod() == 10) outfile << "[!* Distance matrix generated using the PP distance               *]" << endl;

        outfile << "[!*                                                                  *]" << endl;
        outfile << "[!* refs.: -> S. Joly and A. Bruneau (2006) Syst. Biol. 55:623-636   *]" << endl;
        outfile << "[!*        -> S. Joly D. Bryant, and P.J. Lockhart. Unpublished data *]" << endl;
        outfile << "[!*                                                                  *]" << endl;
        outfile << "[!* Final matrix obtained from non-standardized distance matrices    *]" << endl;
        outfile << "[!* of organisms                                                     *]" << endl;
        outfile << "[!*                                                                  *]" << endl;
        outfile << "[!********************************************************************]" << endl;
        outfile << endl << endl;
        outfile << "Begin taxa;" << endl;
        outfile << "    Dimensions ntax=" << org_data.ReturnNOrg() << ";" << endl;
        outfile << "    Taxlabels" << endl;

        for (i=0; i < org_data.ReturnNOrg(); i++)
            {
            outfile << org_data.ReturnOrganism(i) << endl;
            }

        outfile << ";" << endl;
        outfile << "end;" << endl << endl;
        outfile << "Begin distances;" << endl;
        outfile << "    Format triangle=both labels diagonal;" << endl;
        outfile << "matrix" << endl << endl;

	    //Claculate the length of the longest taxa to make nice output

	    int max_length=0;
	    int current_length;

 	   	for (i=0; i < org_data.ReturnNOrg(); i++)
  	  	{
   	     	current_length=org_data.ReturnLengthofOrg(i);
  	      	if (current_length > max_length)
            	{
            	max_length = current_length;
            	}
    	}

    	int size_of_current_taxa;

    	for (i=0; i < org_data.ReturnNOrg(); i++)
        	{
            outfile << org_data.ReturnOrganism(i); 
        	size_of_current_taxa = org_data.ReturnLengthofOrg(i);
        	for (j=0; j < ( (max_length) - size_of_current_taxa); j++)
            	{
            	outfile << " ";
            	} 
			outfile << '\t';
            for (j=0; j < org_data.ReturnNOrg(); j++)
                {
                outfile << setprecision(PRECISION) << final_matrix[i][j] << '\t';
                }
            outfile << endl;
        	}

        outfile << ";" << endl;
        outfile << "end;" << endl;
        }

    outfile.close();
}



/*****************************************************************************

Function: output_file

Description: This function outputs a distance matrix in a file in nexus format


******************************************************************************/

void output_file(void)
{

    string outputfile;

    if (info_data.ReturnIsOutputFile())
        {
        outputfile = info_data.ReturnOutputFile();
        }
    else
        {
        cout << endl << "Enter the name (no space) of the output file: ";
        cin >> outputfile;
        }

    ofstream outfile(outputfile.c_str());

    if (outfile.bad())
        {
        cout << "Error in opening output file" << endl;
        cout << "Verify that there are no white space in the name" << endl;
        }

    outfile << "#nexus" << endl << endl;
    outfile << "[!********************************************************************]" << endl;
    outfile << "[!*                                                                  *]" << endl;
    outfile << "[!* File generated by POFAD "<<VERSION<<"                                     *]" << endl;
    outfile << "[!*                                                                  *]" << endl;
    outfile << "[!* refs.: -> S. Joly and A. Bruneau (2006) Syst. Biol. 55:623-636   *]" << endl;
    outfile << "[!*        -> S. Joly D. Bryant, and P.J. Lockhart. Unpublished data *]" << endl;
    outfile << "[!*                                                                  *]" << endl;
    outfile << "[!********************************************************************]" << endl;
    outfile << endl << endl;
    outfile << "Begin taxa;" << endl;
    outfile << "    Dimensions ntax=" << org_data.ReturnNOrg() << ";" << endl;
    outfile << "    Taxlabels" << endl;

    for (i=0; i < org_data.ReturnNOrg(); i++)
        {
        outfile << org_data.ReturnOrganism(i) << endl;
        }

    outfile << ";" << endl;
    outfile << "end;" << endl << endl;
    outfile << "Begin data;" << endl;
    outfile << "    Dimensions ntax=" << org_data.ReturnNOrg() << " nchar=" << org_data.ReturnNChars() << " ;" << endl;
    outfile << "    Format datatype=" << org_data.ReturnDatatype() << " gap=" << org_data.ReturnGapChar() << " ;" << endl;
    outfile << "matrix" << endl << endl;

    //Claculate the length of the longest taxa to make nice output

    int max_length=0;
    int current_length;

   	for (i=0; i < org_data.ReturnNOrg(); i++)
  		{
     	current_length=org_data.ReturnLengthofOrg(i);
      	if (current_length > max_length)
           	{
           	max_length = current_length;
           	}
   		}

   	int size_of_current_taxa;

   	for (i=0; i < org_data.ReturnNOrg(); i++)
       	{
        outfile << org_data.ReturnOrganism(i); 
       	size_of_current_taxa = org_data.ReturnLengthofOrg(i);
       	for (j=0; j < ( (max_length) - size_of_current_taxa); j++)
           	{
           	outfile << " ";
           	} 
		outfile << '\t';
        outfile << org_data.ReturnSequence(i) << endl;
	   	}

    outfile << ";" << endl;
    outfile << "end;" << endl;

    outfile.close();
}



/***********************************

Function : EstimateMissingDist

void GL_Main(double **DD, int nbObj, int methode, double **distanceArbre,
				 double *RESULTATS, long int *ARETE, double *LONGUEUR, double **W)

***********************************/

// Main function calling the tree reconstruction method from incomplete matrices

void EstimateMissingDist(double **DD)
{

	bool add;
//	bool complete;	   	  
	int i, j, T, nbmiss, sum, sum2=0, **B2, *aa;
	int  nbcycles;
    double  *max, **m;     	
//	int tt, k; 
//    int NbInc=0;
	double **DI;
//	int error, Iternum=50;
	int methode=3;               //addtitive estimation
	int n=org_data.ReturnNOrg();

	if (methode!=4)
		{
		for (i=0;i<n;i++)
		    {
			for (j=i+1;j<n;j++)
				{
				if ((DD[i][j]<0.0)&&(DD[i][j]!=-999)) 
				{cerr << "Dissimilarity matrix contains negative values not equal to -999 (value marking missing data). Computation is not possible for this method."; exit(1); }	  			
				}
	    	}
		}


	DI=(double **) malloc((n)*sizeof(double*)); 
    for (i=0;i<n;i++)
		{
    	DI[i]=(double*)malloc((n)*sizeof(double));
    	if (DI[i]==NULL) { cerr << "Data matrix is too large"; exit(1);}
		}

/*	{ 		
		Dist=(double **) malloc((n+1)*sizeof(double*)); 
		for (i=0;i<=n;i++)
		{
			Dist[i]=(double*)malloc((n+1)*sizeof(double));
			if (Dist[i]==NULL)
			{
				printf("Data matrix is too large"); exit(1); 
			}
		}	
	}
*/
/*		
	for (i=1;i<=n;i++)
    {
		for (j=i+1;j<=n;j++)
		{
//			k=1+floor1((n-0.5*i)*(i-1)+j-i-1);    			
			DI[i][j]=DD[i][j];
			DI[j][i]=DI[i][j];
			if (Dist[i-1][j-1]<0)	NbInc++;
			else if (Dist[i-1][j-1]>Dmax) Dmax=Dist[i-1][j-1];
		}
		DI[i][i]=0.0;
    }
*/
    // Test for avoiding entire rows and columns with missing entries
	int FLAG1;
	for (i=0;i<n;i++)
    {
		FLAG1=1;
		for (j=0;j<n;j++)
		{			   
			if ((i!=j)&&(DI[i][j]!=-999)) { FLAG1=0; break; }
		}
		if (FLAG1==1) { cerr <<"Too many missing values. You have to have at least one value by row and by colomn.";
		                exit(1); } 
    }

    //Selection of a tree reconstruction method
	if ((methode==2)||(methode==3)) //Ultrametric and Additive procedures
	{

	T=org_data.ReturnNOrg();	
	max = (double*)malloc((T) * sizeof(double));
	m = (double**)malloc((T) * sizeof(double *));
	B2 = (int**)malloc((T) * sizeof(int *));		
    if ((max== NULL) ||(m == NULL) || (B2 == NULL))
	{ printf("Data matrix is too large"); exit(1);} 	
	
	for (i=0;i<T;i++)
	{ 			
		m[i] = (double*)malloc((T) * sizeof(double));
		B2[i] = (int*)malloc((T) * sizeof(int));
		if ((m[i] == NULL) || (B2[i] == NULL))
		{ printf("Data matrix is too large"); exit(1);} 
	}				

    for (i=0;i<n;i++)
    {
		for (j=i+1;j<n;j++)
		{
//			k=1+floor1((n-0.5*i)*(i-1)+j-i-1);    
			m[i][j]=m[j][i]=DD[i][j];
         }			
		m[i][i]=0.0;
    }
	
    nbmiss=0;
//    complete=true;	
//	tt=(T * (T - 1))/2;
	for (i=0; i<T; i++)
	{
		for (j=0; j<T; j++)
		{
			if (m[i][j]!=-999) B2[i][j]=1;
			else { B2[i][j]=0; nbmiss++;}
		}
	}
//    if (nbmiss>0) complete=false;    
	  	 
	aa=(int*)malloc((T) * sizeof(int));
	  	  
	if (methode==2) add=false; 
	if (methode==3) add=true;  
	 
	if (nbmiss > 0)
	{
		nbultra=0;
		nbcycles=1;
	    if (add) Additif(&sum, &sum2, T, B2, m, max);
	    else Ultra2(&sum, &sum2, T, B2, m, max);
	}
	
//	cout << "OK!" << endl << endl;
//	cout << "Number of cycles: " << nbcycles << '\t' << "sum2: " << sum2 << endl;
		 
	while (sum2>0)
	{
		nbcycles++;
		if (!add) Ultra2(&sum, &sum2, T, B2, m, max);
		if (add)
		{
			if (sum2==sum)
			{
				 nbultra=nbultra+1;
			  	 alea(T, aa);
			     Ultra1(&sum, &sum2, T, B2, m, max, aa);
			}
			Additif(&sum, &sum2, T, B2, m, max);
		}
	    cout << "Number of cycles: " << nbcycles << '\t' << "sum2: " << sum2 << endl;
	 }
	 
/*	   
	    //method MW to complete Additive and Ultrametric
	    for (i=1; i<=n; i++)
		{
  		 for (j=i+1; j<=n; j++)
		 { 
		  W[i][j]=W[j][i]=1.0;
		 }   
		 W[i][i]=1.0;
		} 
		W[0][0]=1.0;

		for (i=1; i<=n; i++)
		{
  		  for (j=1; j<=n; j++)
			if (m[i][j] <= 0.0) m[i][j] = 0.0005;

		}

		error=parcour211(m,W,Dist,&Iternum,ARETE,LONGUEUR); 
		if (error==-1) {exit(1);}
		
		for (i=1;i<=n;i++)
		{
			for (j=i+1;j<=n;j++)
			{
				k=1+floor1((n-0.5*i)*(i-1)+j-i-1); 
				distanceArbre[i][j]=distanceArbre[j][i]=Dist[i][j];
			}
			distanceArbre[i][i]=0.0;	     
		}		
		compute_criteres11(DI, Dist, RESULTATS, LONGUEUR, n);
*/	 }

/*     
	 else if ((methode==4)||(methode==5)) // MW-modified method or MW*
	 {

		int FLAG=0;
		for (i=1; i<=n; i++)
		{
  			for (j=i+1; j<=n; j++)
			{  		  		
  				if ((DI[i][j]==-99.00)&&(DI[j][i]==-99.00)) {W[i][j]=W[j][i]=0.0;FLAG=1;}
  				else if ((DI[i][j]==-99.00)&&(DI[j][i]!=-99.00)) {DI[i][j]=DI[j][i]; W[i][j]=W[j][i]=1.0;}
  				else if ((DI[i][j]!=-99.00)&&(DI[j][i]==-99.00)) {DI[j][i]=DI[i][j]; W[i][j]=W[j][i]=1.0;}
  				else {W[i][j]=W[j][i]=1.0;}
			}   
			DI[i][i]=0.0;
			W[i][i]=1.0;
		}
		DI[0][0]=0; W[0][0]=0;

		if ((FLAG==1)&&(methode==5)) FLAG=tryout(DI, W);     
  
		if (n==2)
		{
			Dist[1][2]=DI[1][2];
			Dist[2][1]=DI[1][2];
			Dist[1][1]=0.0;
			Dist[2][2]=0.0;  
		}      
		else 
		{    
			Iternum=50; 
			DI[0][0]=1;

		    if (FLAG==0) error=parcour211(DI,W,Dist,&Iternum,ARETE,LONGUEUR); 						   
			else error=parcour2MisVal(DI,W,Dist,&Iternum,ARETE,LONGUEUR);
           
			if (error==-1) { exit(1);}
		}   
        
		for (i=1;i<=n;i++)
		{
			for (j=i+1;j<=n;j++)
			{
				k=1+floor1((n-0.5*i)*(i-1)+j-i-1); 
				distanceArbre[i][j]=distanceArbre[j][i]=Dist[i][j];
			}	     
		}		
		compute_criteres11(DD, Dist, RESULTATS, LONGUEUR, n);
	}
*/

	for (i=0;i<n;i++) {
		for (j=0;j<n;j++) {
			DD[i][j]=m[i][j];
			}
		}

  // memory liberation
	for (i=0;i<n;i++)
	{ 
		free(DI[i]);
//		free(Dist[i]);
		if ((methode==2)||(methode==3)) free(m[i]);
        if ((methode==2)||(methode==3)) free(B2[i]);
	}
	free(DI);
//	free(Dist);
	if ((methode==2)||(methode==3)) free(m);
    if ((methode==2)||(methode==3)) free(B2);
	if ((methode==2)||(methode==3)) free(max);

	return;				
}


// returns the closest integer value of x
int floor1(double x)
{
  int i;
  
  if (ceil(x)-floor(x)==2.0) i=(int)x; 
  else if (fabs(x-floor(x)) > fabs(x-ceil(x))) i=(int)ceil(x);
  else i=(int)floor(x);
  return i;
} 



void alea (int T, int *aa)
{  
	  int i, j, k;
	  double x, *P;
	  double XY=0.98674532;
	  
	  P = (double*)malloc((T+1) * sizeof(double));
//	  assert(P!=NULL);
			  			
      for (i=0;i<T;i++)
      {
		double D2P31M=2147483647.0, D2P31=2147483648.0;                                                      
        double a;                                                                   
        a=16807 * XY;                           
        XY=a - floor1(a/D2P31M) * D2P31M;
	    P[i]=XY/D2P31; 
	  }
	    
      x=1;
      for (i=0;i<T;i++)
      {
        for (j=0;j<T;j++)
	    {
            if (P[j] < x) 
               { x = P[j]; aa[i] = j; }		        		     
        }       
      for (k=0;k<T;k++)
      {
	    if (P[k] == x) P[k] = 1;
      }
      x=1;
     } 
}// end alea

 void Ultra1(int *sum, int *sum2, int T, int **B2, double **m, double *max, int *aa)
 {			
	    int i, j, k, a;
	    double min, val1, val2;
	    bool done;
				
		done = false;
		*sum = 0;
	    *sum2 = 0;
		min = 1000;
		for (i=0;i<T;i++)
        {
          for (j=0;j<T;j++)
	      {
			 if (B2[aa[i]][aa[j]]==0) *sum=*sum+1;
		   } 
		}   
		for (i=0;i<T;i++)
        {
           for (j=0;j<T;j++)
	       {
			  if (B2[aa[i]][aa[j]]==1) goto label2;
						for (k=0;k<T;k++)
                        {
						  if ((k==i)||(k==j)) goto label1;									
								if ((B2[aa[i]][aa[k]]==0)||(B2[aa[j]][aa[k]]==0)) goto label1;																	
								val1=(m[aa[i]][aa[k]] + m[aa[k]][aa[i]])/2.0;
								val2=(m[aa[j]][aa[k]] + m[aa[k]][aa[j]])/2.0;
              					if (val1>val2) max[aa[k]]=val1;          			   
                                else max[aa[k]]=val2;
                                   								
								if (max[aa[k]] < min)
								{
							  min=max[aa[k]];
								  m[aa[i]][aa[j]]=min;
								  m[aa[j]][aa[i]]=min;     
								  done=true;
								}
label1:                      a=1;
						}
						min=1000;
						if (done) goto label3;
label2:                  a=1;
			}			
		}
        
label3:     a=1;
       for (i=0;i<T;i++)
       {
          for (j=0;j<T;j++)
	      {
			 if (m[aa[i]][aa[j]]!=miss) B2[aa[i]][aa[j]]=1;			
		  }
	   }  			
				     
        for (i=0;i<T;i++)
        {
          for (j=0;j<T;j++)
	      {
			if (B2[aa[i]][aa[j]]==0) *sum2=*sum2+1;
		  }
		} 		
	} // end Ultra1


void Ultra2(int *sum, int *sum2, int T, int **B2, double **m, double *max)
{			
	    int i, j, k;
	    double min, val1, val2;

		*sum = 0;
		*sum2 = 0;
		min = 1000;
		for (i=0;i<T;i++)
        {
          for (j=0;j<T;j++)
	      {
			 if (B2[i][j]==0) *sum=*sum+1;
		   } 
		}   
		for (i=0;i<T;i++)
        {
           for (j=0;j<T;j++)
	       {
			  if (B2[i][j]==1) continue;
						for (k=0;k<T;k++)
                        {
						  if ((k==i)||(k==j)) continue;									
								if ((B2[i][k]==0)||(B2[j][k]==0)) continue;																	
								val1=(m[i][k] + m[k][i])/2.0;
								val2=(m[j][k] + m[k][j])/2.0;
              					if (val1>val2) max[k]=val1;          			   
                                else max[k]=val2;
                                   								
								if (max[k] < min)
								{
							      min=max[k];
								  m[i][j]=min;
						
								}
						}
						min=1000;
			}
		}
        
       for (i=0;i<T;i++)
       {
          for (j=0;j<T;j++)
	      {
			 if (m[i][j]!=miss) B2[i][j]=1;			
		  }
	   }  			
				     
        for (i=0;i<T;i++)
        {
          for (j=0;j<T;j++)
	      {
			if (B2[i][j]==0) *sum2=*sum2+1;
		  }
		} 		
} // end Ultra2
	
void Additif(int *sum, int *sum2, int T, int **B2, double **m, double *max)
{	
	int i, j, k, l;
	double min, val1, val2, val3, newval;
	*sum=0;
	*sum2=0;
	min=1000;

	for (i=0;i<T;i++)
    	{
    	for (j=0;j<T;j++)
    		{
			if (B2[i][j]==0) *sum=*sum+1;
			} 
		}   
	
	for (i=0;i<T;i++)
    	{
    	for (j=0;j<T;j++)
    		{
			if (B2[i][j]==1) continue;
			for (k=0;k<T;k++)
				{									
				if ((k==i)||(k==j)) continue;
				if ((B2[i][k]==0)||(B2[j][k]==0)) continue;
				for (l=0;l<T;l++)
			        {
			        if ((l==i)||(l==j)||(l==k)) continue;
			        if ((B2[i][l]==0)||(B2[j][l]==0)||(B2[k][l]==0)) continue;
					val1=(m[k][l]+m[l][k])/2.0;
					val2=(m[i][l]+m[l][i])/2.0 + (m[j][k]+m[k][j])/2.0;
					val3=(m[i][k]+m[k][i])/2.0 + (m[j][l]+m[l][j])/2.0;
					if (val2>val3) max[l]=val2;
			        else max[l]=val3;
					newval=max[l]-val1;											
		            if (newval<min)
				    	{
					    min=newval;
					    m[i][j]=min;
				    	}
					}
				}
			min=1000;
			}
		}				
   
    for (i=0;i<T;i++)
    	{
    	for (j=0;j<T;j++)
    		{
			if (m[i][j]!=-999) B2[i][j]=1;			
			}
    	}  							     
    for (i=0;i<T;i++)
    	{
      	for (j=0;j<T;j++)
      		{
			if (B2[i][j]==0) *sum2=*sum2+1;
	  		}
		} 										
} // end Additif


/***********************************

Function : doOrgDimensions

***********************************/

void doOrgDimensions(void)
{
    while ( (aTokenPtr=nextToken()) != ";")
        {
        if (parse_assignment2("NORG"))
            {
            org_data.GetNOrg( atoi(LocalToken.c_str()) );
            }
        }
    return;
}


/***********************************

Function : doOrgLabels

***********************************/


void doOrgLabels(void)
{
    int size;
    while ( (aTokenPtr=nextToken()) != ";")
        {
        org_data.ReadOrganism(aTokenPtr);
        }
    size = org_data.NumberOrgLabels();       //Number of organisms that were read
    if ( size < org_data.ReturnNOrg()) {
        cerr << "Too few Organism labels";
        exit(1);
    }
    else if ( size > org_data.ReturnNOrg()) {
        cerr << "Too many Organism labels";
        exit(1);
    }
    org_data.OrganismsRead();                //Flag to indicate the organisms were read
    return;
}

/*******************************

Function : doUnrecognizedBlock

*******************************/

void doUnrecognizedBlock(void)
{
	do 				
		{
		aTokenPtr=nextToken();
		}  while ( (aTokenPtr != "END")  &&
						(aTokenPtr != "ENDBLOCK" ) );
	aTokenPtr=nextToken();  /* Pop the terminating semicolon */
	if (aTokenPtr != ";") {
        cerr << "Block not terminated with semicolon";
        exit(1);
    }
	return;
}





/****************************************************************/
int parse_assignment2(string target)

/* on entry 'aTokenPtr' contains the putative first word of a three token
assignment statement of the form 'word1=word2'.  This function checks to see
if word1 is the same as 'target' and if so, it returns the address of a string
containing 'word2' or NULL if an error occurs.  aTokenPtr is
set to the last token in the assignment statement
If no match, aTokenPtr is left unchanged!! */

{
	if (aTokenPtr == target)
			{
			aTokenPtr=nextToken();
			if (aTokenPtr != "=")
				{
				cerr << "Bad assignment statement= " << aTokenPtr << endl;
				exit(1);
				}
			aTokenPtr=nextToken();
			LocalToken = aTokenPtr;
			return 1;
			}
	else return 0;
}


/***************************************************

Function: nextToken

***************************************************/

/*
	Function originally written by Michael J. Sanderson
	original name: nextToken2.c
	included in the program r8s
*/


/*	
	Gets the next token from input stream 'fpointer', and copies it onto the global
	buffer pointed to by 'aTokenPtr'.  If there is NO next token, we copy a null
	string onto that buffer.  That's a signal for the main caller routine...

	If the global variable gNewLine=1 then the newline characters, '\n' and '\r'
	ARE returned as individual tokens,  when encountered.  The normal state is
	gNewLine=0,  which treats these as white space delimiters too.  The only time
	NEXUS file needs to think about newlines is when reading interleaved matrices!

*/

string nextToken(void)
	{
	string TempPtr;		// a pointer to manipulate aTokenPtr
	char c;
	
	if  ((c=*bufPtr++) == '\0') NULL_RETURN;
	
	while (( isNEXUSwhiteSpace(c) ) || (c=='['))  
			    /* this whole loop is in case multiple comments separated by whitespace */
		{
		while ( isNEXUSwhiteSpace(c) )  /* skip white space and land on next 'c'*/
			{
			c=*bufPtr++;
			if (c=='\0')  NULL_RETURN;/* check for embedded EOF */
			}
			    
		if (c=='[')		/* skip the comment and land on next 'c' after comment */
			{
			while (c !=']')
				{
				c=*bufPtr++;
				if (c=='\0') NULL_RETURN;  /* check for embedded EOF */
				}
			c=*bufPtr++;	/* get next char after ']' */
			if (c=='\0') NULL_RETURN;
			}
		}	    
    
	if (c=='\'')		/* deal with single-quoted tokens */
		{
		TempPtr.push_back(toupper(c));
		while (  (c=*bufPtr++) != '\'')
			{
			if (c=='\0') NULL_RETURN;  /* check for embedded EOF */
			TempPtr.push_back(toupper(c));	/* this is a valid character in the word, add to token */
			}
		TempPtr.push_back(toupper(c));	/* add the terminating quote too */
		return(TempPtr);
		}	/* return everything between single quotes, including the quotes, as a token*/		

	TempPtr.push_back(toupper(c));		/* char is either punctuation or part of word, so add it to token */
	    
	if (!isNEXUSpunct(c))	/* next char is part of word, so add all of word until white,punct,eof,
								    or Token size exceeded */
	    {
	    for (;;)
		    {
		    c=*bufPtr++;
		    if (  isNEXUSpunct(c) || isNEXUSwhiteSpace(c) ||  (c == '\0') )
			    {
			    --bufPtr; /* word is terminated by some c that is not part of word;
							     push c back into stream and deal with it on
							    next call to this function; meantime, break out, 
							    and return this token*/
			    break;
			    }
		    TempPtr.push_back(toupper(c));	/* this is a valid character in the word, add to token */
		    }
	    }
	return(TempPtr);
	}

/***************************************************

Function: nextToken2

***************************************************/

/*
	Function originally written by Michael J. Sanderson
	original name: nextToken2.c
	included in the program r8s
*/


/*	
	Sames as nextToken but it does not change the case of character (important when reading files
        in UNIX).

*/

string nextToken2(void)
	{
	string TempPtr;		// a pointer to manipulate aTokenPtr
	char c;
	
	if  ((c=*bufPtr++) == '\0') NULL_RETURN
	
	while (( isNEXUSwhiteSpace(c) ) || (c=='['))  
			    /* this whole loop is in case multiple comments separated by whitespace */
		{
		while ( isNEXUSwhiteSpace(c) )  /* skip white space and land on next 'c'*/
			{
			c=*bufPtr++;
			if (c=='\0')  NULL_RETURN;/* check for embedded EOF */
			}
			    
		if (c=='[')		/* skip the comment and land on next 'c' after comment */
			{
			while (c !=']')
				{
				c=*bufPtr++;
				if (c=='\0') NULL_RETURN;  /* check for embedded EOF */
				}
			c=*bufPtr++;	/* get next char after ']' */
			if (c=='\0') NULL_RETURN;
			}
		}	    
    
	if (c=='\'')		/* deal with single-quoted tokens */
		{
		TempPtr.push_back(c);
		while (  (c=*bufPtr++) != '\'')
			{
			if (c=='\0') NULL_RETURN;  /* check for embedded EOF */
			TempPtr.push_back(c);	/* this is a valid character in the word, add to token */
			}
		TempPtr.push_back(c);	/* add the terminating quote too */
		return(TempPtr);
		}	/* return everything between single quotes, including the quotes, as a token*/		

	TempPtr.push_back(c);		/* char is either punctuation or part of word, so add it to token */
	    
	if (!isNEXUSpunct(c))	/* next char is part of word, so add all of word until white,punct,eof,
								    or Token size exceeded */
	    {
	    for (;;)
		    {
		    c=*bufPtr++;
		    if (  isNEXUSpunct(c) || isNEXUSwhiteSpace(c) ||  (c == '\0') )
			    {
			    --bufPtr; /* word is terminated by some c that is not part of word;
							     push c back into stream and deal with it on
							    next call to this function; meantime, break out, 
							    and return this token*/
			    break;
			    }
		    TempPtr.push_back(c);	/* this is a valid character in the word, add to token */
		    }
	    }
	return(TempPtr);
	}
