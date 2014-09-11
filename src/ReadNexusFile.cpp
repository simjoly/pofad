#include "ReadNexusFile.h"


//TODO: Need to check whether a taxa name appears more than once in the nexus file!

/* 
	Read a NEXUS file buffer and set up a global data structure containing everything. 
	See nexus.h for that data structure.
	Returns NULL on error. 
*/

using namespace std;


/*****  private functions ******/

static int parse_assignment2(string target);
static void doTaxaBlock(void);
static void doDataBlock(void);
static void doCharsBlock(void);
static void doDistBlock(void);
//static void doSpeciesBlock(void);
static void doDataFormat(void);
//static void doSpecies(void);
static void doMatrix(void);
static void doTaxLabels(void);
static void doDistFormat(void);
static void doDistMatrix(void);
static void doUnrecognizedBlock(void);
static void doDimensions(void);
static string nextToken(void);
static string nextToken2(void);
static string nextToken3(void);

/**** Global variables ****/

string  LocalToken;
char 	*bufPtr;
string aTokenPtr;



/***** Main Function *****/

void readNexusFile(char *theNexusFileBuffer)
{
	string stemp;

	bufPtr=theNexusFileBuffer;	/* Initialize this global to beginning of the 
					buffer and will sweep through it until end of buffer */
	
	if ( bufPtr != NULL )
		{
		aTokenPtr=nextToken();
		if(aTokenPtr != "#NEXUS")
			{
			cerr << "Error: Not a NEXUS file";
            exit(1);
			}
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
					if      (stemp == "TAXA") 	doTaxaBlock();
//					else if (stemp == "SPECIES") 	doSpeciesBlock();
					else if (stemp == "CHARACTERS")	doCharsBlock();
					else if (stemp == "DATA")	doDataBlock();	
					else if (stemp == "DISTANCES")	doDistBlock();
                                        else  /* token is not a recognized block */
                                            {  
                                            doUnrecognizedBlock();
                                            }
					}
				stemp.clear();				
				}
			}
		cout << "  -> nexus file processed" << endl;
		return;	/* normal return */
		}
	else
		{
		cerr << "Error opening nexus file";
        exit(1);
		}
	
}


/****************************************************************/
/****************   BLOCK PROCESSING FUNCTIONS  *****************/
/****************************************************************/


/****************************************************************/
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

void doTaxaBlock(void)
{
	do 				/* need to put in error checking in case no DIMENSIONS statement */
		{
		aTokenPtr=nextToken();
		if (aTokenPtr == "DIMENSIONS")  
			doDimensions();
		if (aTokenPtr == "TAXLABELS")
			doTaxLabels();
		}  while ( (aTokenPtr != "END") && (aTokenPtr != "ENDBLOCK") );
	aTokenPtr=nextToken();
	if (aTokenPtr != ";") {
        cerr << "Block not terminated with semicolon";
        exit(1);
    }
	return;
}

/**************************************************************/

void doDataBlock(void)
{
	do
		{
		aTokenPtr=nextToken();
		if (aTokenPtr == "DIMENSIONS")
			doDimensions();
		if (aTokenPtr == "FORMAT")
			doDataFormat();
		if (aTokenPtr == "MATRIX")
			doMatrix();

		}  while ( (aTokenPtr != "END") && (aTokenPtr != "ENDBLOCK") );
	aTokenPtr=nextToken();
	if (aTokenPtr != ";") {
        cerr << "Block not terminated with semicolon";
        exit(1);
    }
	return;
}

/**************************************************************/

void doCharsBlock(void)
{
	do
		{
		aTokenPtr=nextToken();
		if (aTokenPtr == "DIMENSIONS")
			doDimensions();
		if (aTokenPtr == "FORMAT")
			doDataFormat();
		if (aTokenPtr == "MATRIX")
			doMatrix();

		}  while ( (aTokenPtr != "END") && (aTokenPtr != "ENDBLOCK") );
	aTokenPtr=nextToken();
	if (aTokenPtr != ";") {
        cerr << "Block not terminated with semicolon";
        exit(1);
    }
	return;
}


/****************************************************************/

void doDistBlock(void)
{
	do
		{
		aTokenPtr=nextToken();
		if (aTokenPtr == "DIMENSIONS")  
			doDimensions();
		if (aTokenPtr == "FORMAT")  
			doDistFormat();
		if (aTokenPtr == "TAXLABELS")
			doTaxLabels();
		if (aTokenPtr == "MATRIX")
			doDistMatrix();
		}  while ( (aTokenPtr != "END") && (aTokenPtr != "ENDBLOCK") );
	aTokenPtr=nextToken();
	if (aTokenPtr != ";") {
        cerr << "Block not terminated with semicolon";
        exit(1);
    }
    return;
}


/*************************************************************/

/*
void doSpeciesBlock(void)
{
    do 
        {
        aTokenPtr=nextToken();
        if (aTokenPtr == "DIMENSIONS") doDimensions();
        if (aTokenPtr == "SPECIES") doSpecies();
        }  while ( (aTokenPtr != "END") && (aTokenPtr != "ENDBLOCK") );

	aTokenPtr=nextToken();
	if (aTokenPtr != ";")
		doGenericAlert("Block not terminated with semicolon");
	return;

}
*/
/****************************************************************/
/**************** COMMAND PROCESSING FUNCTIONS ******************/
/****************************************************************/


/*--------------------------------------------------------------*/
/*
void doSpecies(void)
{
    
    nexus_data.initializeSpeciesPosition();
    nexus_data.initNbLineages();
    int SpePosition = 0;
    int NbSpecies = 0;
    int NbLineagesInSpecies;

    while ( (aTokenPtr=nextToken()) != ";")
        {
        nexus_data.ReadSpecies(aTokenPtr);
        nexus_data.GetSpeciesPosition(NbSpecies, SpePosition);
        if ((aTokenPtr=nextToken()) != "=")
            {
            doGenericAlert("The species name is not followed by = ");
            }
        NbLineagesInSpecies = 0;
        while ( (aTokenPtr=nextToken()) != "," )
            {
            nexus_data.ReadLineage(aTokenPtr);
            SpePosition++;
            NbLineagesInSpecies++;
            }
        nexus_data.NbLineagesInSpecies(NbSpecies, NbLineagesInSpecies);

        NbSpecies++;
        }
    return;
}
*/


/*----------------------------------------------------------------*/
void doDimensions(void)
{
	while ( (aTokenPtr=nextToken()) != ";")
		{
		if (parse_assignment2("NCHAR"))
			{
			nexus_data.GetNChars(atoi(LocalToken.c_str()));
			cout << "  -> Number of characters in matrix = " << nexus_data.ReturnNChars() << endl;
			}
		if (parse_assignment2("NTAX"))
			{
			nexus_data.GetNTaxa(atoi(LocalToken.c_str()));
			cout << "  -> Number of taxa in matrix = " << nexus_data.ReturnNTax() << endl;
			}

/*		if (parse_assignment2("NSPE"))
			{
			nexus_data.GetNSpe(atoi(LocalToken.c_str()));
//			cout << "  -> Number of species in matrix = " << nexus_data.ReturnNSpe() << endl;
			}
*/
		}
	return;
}


/*----------------------------------------------------------------*/

void doTaxLabels(void)
{
	while ( (aTokenPtr=nextToken()) != ";")
		{
		nexus_data.ReadTaxa(aTokenPtr);
		}
	if (nexus_data.NumberTaxaLabels() < nexus_data.ReturnNTax()) {
        cerr << "Too few taxon labels";
        exit(1);
    }
	else if (nexus_data.NumberTaxaLabels() > nexus_data.ReturnNTax()) {
        cerr << "Too many taxa labels";
        exit(1);
    }
	else
		nexus_data.isTaxa();	/* Set flag showing that labels read */
	return;
}


/*----------------------------------------------------------------*/

void doDataFormat(void)
{
	while ( (aTokenPtr=nextToken()) != ";") {
		if(aTokenPtr == "INTERLEAVE")
			nexus_data.IsInterleave(1);

		if (aTokenPtr == "DATATYPE") {
			aTokenPtr=nextToken();
			if (aTokenPtr != "=") {
				cerr << "Bad assignment statement= " << aTokenPtr << endl;
                exit(1);
            }
			aTokenPtr=nextToken();
			nexus_data.GetDatatype(aTokenPtr);
        }
		if (aTokenPtr == "GAP") {
			aTokenPtr=nextToken();
			if (aTokenPtr != "=") {
				cerr << "Bad assignment statement= " << aTokenPtr << endl;
				exit(1);
            }
			aTokenPtr=nextToken();
			nexus_data.GetGapChar(aTokenPtr);
        }
		if (aTokenPtr == "MISSING") {
			aTokenPtr=nextToken();
			if (aTokenPtr != "=") {
				cerr << "Bad assignment statement= " << aTokenPtr << endl;
				exit(1);
            }
			aTokenPtr=nextToken();
			nexus_data.MissingChar(aTokenPtr);
        }
    }
    /* Report */
/*    cout << "  -> Datatype = " << nexus_data.ReturnDatatype() << endl;
    cout << "  -> Gap character = " << nexus_data.ReturnGapChar() << endl;
    cout << "  -> Missing character = " << nexus_data.ReturnMissingChar() << endl;
*/
    return;
}



/*----------------------------------------------------------------*/
void doMatrix(void)

/* reads a NEXUS data matrix and stores in global dmatrix as a set of strings corresponding 
to the rows of the matrix; interleaved and whitespace is allowed */

{

    int i;
    string name;  //To store the taxa name
    nexus_data.InitDataMatrix();

    if(nexus_data.ReturnisTaxa() == 0) //i.e., no taxa block
        {

/* This works OK */
        if (nexus_data.ReturnInterleave() == 0) // The dataset is not interleaved
            {
            i=0;
            while ( (aTokenPtr=nextToken()) != ";")            
                {
                nexus_data.ReadTaxa(aTokenPtr);
                while ( nexus_data.ReturnCharacterForTaxa(i) < nexus_data.ReturnNChars() )
                    {
                    aTokenPtr=nextToken();
                    nexus_data.AddCharacter(i, aTokenPtr);
                    }
                i++;
                }
            if ( i != nexus_data.ReturnNTax() )
                {
                cout << "WARNING: number of taxa in matrix is not the same as in the taxa block!" << endl << endl;
                pauseprogram();
                }        
            }


/* This part works well */
        else if (nexus_data.ReturnInterleave() == 1)  // The dataset is interleaved
            {
            for (i=0; i < nexus_data.ReturnNTax(); i++)
                {
                aTokenPtr=nextToken2();
                while (aTokenPtr == "eol")
                    {
                    aTokenPtr=nextToken2();
                    }
                if ( aTokenPtr == ";")
                    {
                    cerr << "Expecting characters, found ;";
                    exit(1);
                    }
                nexus_data.ReadTaxa(aTokenPtr);
                while ( (aTokenPtr=nextToken2()) != "eol" )
                    {
                    nexus_data.AddCharacter(i, aTokenPtr);
                    }
                }
            while ( nexus_data.ReturnCharacterForTaxa( nexus_data.ReturnNTax() - 1 ) < nexus_data.ReturnNChars() )
                {
                for (i=0; i < nexus_data.ReturnNTax(); i++)
                    {
                    aTokenPtr=nextToken2();  //Pop off the taxon name
                    while (aTokenPtr == "eol")
                        {
                        aTokenPtr=nextToken2();
                        }
                    if ( aTokenPtr == ";") 
                        {
                        cerr << "Expecting characters, found ;";
                        exit(1);
                        }
                    while ( (aTokenPtr=nextToken2()) != "eol" )
                        {
                        nexus_data.AddCharacter(i, aTokenPtr);
                        }
                    }
                }
            }
        }

/* This works OK */
    else if (nexus_data.ReturnisTaxa() == 1) // taxa block present
        {
        if (nexus_data.ReturnInterleave() == 0) //Dataset is not interleaved
            {
            i = 0;
            while ( (name=nextToken()) != ";")
                {
                while ( nexus_data.ReturnCharacterForTaxa(nexus_data.GetPositionofTaxa(name)) < nexus_data.ReturnNChars() )
                    {
                    aTokenPtr=nextToken();
                    nexus_data.AddCharacter(i, aTokenPtr);
                    }
                i++;
                }
            if ( i != nexus_data.ReturnNTax() )
                {
                cout << "WARNING: number of taxa in matrix is not the same as in the taxa block!" << endl << endl;
                pauseprogram();
                }        
            }

/* This works OK */
/* TODO: Needs to check in case there are too many characters */
        else if (nexus_data.ReturnInterleave() == 1) //Dataset is interleaved
            {
            while ( nexus_data.ReturnCharacterForTaxa( (nexus_data.ReturnNTax() - 1) ) < nexus_data.ReturnNChars() )
                {
                if ( (name=nextToken2()) == "eol") continue; //This pop off the taxon name
                if (name == ";") 
                    {
                    break;
                    }
                while ( (aTokenPtr=nextToken2()) != "eol" )
                    {
                    nexus_data.AddCharacter(nexus_data.GetPositionofTaxa(name), aTokenPtr);
                    }
                }
            if ( nexus_data.ReturnCharacterForTaxa( nexus_data.ReturnNTax() - 1 ) != nexus_data.ReturnNChars() )
                {
                cerr << "Problem in reading matrix: expecting characters, found ;";
                exit(1);
                }
            }
        }
    return;

}


/**********************************************/

void doDistFormat(void)
{
	while ( (aTokenPtr=nextToken()) != ";")
		{
		if(aTokenPtr == "TRIANGLE")
			{
			if((aTokenPtr=nextToken()) == "=")
				{
				aTokenPtr=nextToken();
            			if (aTokenPtr == "LOWER")
					{
					nexus_data.IsTriangle(0);
					}
            			else if (aTokenPtr == "UPPER")
					{
					nexus_data.IsTriangle(1);
					}
            			else if (aTokenPtr == "BOTH")
					{
					nexus_data.IsTriangle(2);
					}
				else
					{
					cerr << "Description of triangle not recognized in distance block; use lower, upper, or both";
                    exit(1);
					}
				}
			else
				{
				cerr << "No = sign after Triangle in the nexus file" << endl;
				exit(1);
				}
			continue;
			}
		if (aTokenPtr == "DIAGONAL")
			{
			nexus_data.IsDiagonal(1);
			continue;
			}
		if (aTokenPtr == "NODIAGONAL")
			{
			nexus_data.IsDiagonal(0);
			continue;
			}
		if (aTokenPtr == "LABELS")
			{
			nexus_data.isLabels(1);
			continue;
			}
		if (aTokenPtr == "NOLABELS")
			{
			nexus_data.isLabels(0);
			continue;
			}
		if (aTokenPtr == "INTERLEAVE")
			{
			cerr << "FATAL: Presently, Pofad does not accept interleaved distance matrices" << endl;
			exit(1);
			}
		}
	return;
}	



/**************************************************************/

void doDistMatrix(void)
{
    int i,j;
    double distance;

    nexus_data.InitMatrix();

    if (nexus_data.ReturnTriangle() == 0) //If the matrix is a lower triangle
        {

        /* If there is a diagonal and if there are labels in the data matrix */

	if ( (nexus_data.ReturnDiagonal() == 1) && (nexus_data.ReturnIsLabels() == 1) )
		{
		for (i=0; i < nexus_data.ReturnNTax(); i++)
			{
			aTokenPtr=nextToken();  // Get rid of the label
			for (j=0; j <= i; j++)
				{
				aTokenPtr=nextToken();
				distance = atof(aTokenPtr.c_str()); 	 //Transform the token in a float dummy
				nexus_data.EnterMatrix(distance, i, j);  //place the float in the dist_matrix
				nexus_data.EnterMatrix(distance, j, i);  //To make a square symetric matrix
				}
			}
		}

	/* If there is no diagonal and if there are labels in the data matrix */

	if ( (nexus_data.ReturnDiagonal() == 0) && (nexus_data.ReturnIsLabels() == 1) )
		{
		for (i=0; i < nexus_data.ReturnNTax(); i++)
			{
			aTokenPtr=nextToken();  // Get rid of the label
			for (j=0; j < i; j++)
				{
				aTokenPtr=nextToken();
				distance = atof(aTokenPtr.c_str()); 	 //Transform the token in a float dummy
				nexus_data.EnterMatrix(distance, i, j);  //place the float in the dist_matrix
				nexus_data.EnterMatrix(distance, j, i);  //To make a square symetric matrix
				}
			nexus_data.EnterMatrix(0, i, i);  //place the "0" diagonal in the dist_matrix
			}		
		}

	/* If there is a diagonal and if there are no labels in the data matrix */

	if ( (nexus_data.ReturnDiagonal() == 1) && (nexus_data.ReturnIsLabels() == 0) )
		{
		for (i=0; i < nexus_data.ReturnNTax(); i++)
			{
			for (j=0; j <= i; j++)
				{
				aTokenPtr=nextToken();
				distance = atof(aTokenPtr.c_str()); 	 //Transform the token in a float dummy
				nexus_data.EnterMatrix(distance, i, j);  //place the float in the dist_matrix
				nexus_data.EnterMatrix(distance, j, i);  //To make a square symetric matrix
				}
			}
		}

	/* If there is no diagonal and if there are no labels in the data matrix */

	if ( (nexus_data.ReturnDiagonal() == 0) && (nexus_data.ReturnIsLabels() == 0) )
		{
		for (i=0; i < nexus_data.ReturnNTax(); i++)
			{
			for (j=0; j < i; j++)
				{
				aTokenPtr=nextToken();
				distance = atof(aTokenPtr.c_str()); 	 //Transform the token in a float dummy
				nexus_data.EnterMatrix(distance, i, j);  //place the float in the dist_matrix
				nexus_data.EnterMatrix(distance, j, i);  //To make a square symetric matrix
				}
			nexus_data.EnterMatrix(0, i, i);  //place the "0" diagonal in the dist_matrix
			}		
		}
        } // if (nexus_data.ReturnTriangle() == 0)

    if (nexus_data.ReturnTriangle() == 1)  //If the matrix is a upper diagonal
        {

        /* If there is a diagonal and if there are labels in the data matrix */

	if ( (nexus_data.ReturnDiagonal() == 1) && (nexus_data.ReturnIsLabels() == 1) )
		{
		for (i=0; i < nexus_data.ReturnNTax(); i++)
			{
			aTokenPtr=nextToken();  // Get rid of the label
			for (j=i; j < nexus_data.ReturnNTax() ; j++)
				{
				aTokenPtr=nextToken();
				distance = atof(aTokenPtr.c_str()); 	 //Transform the token in a float dummy
				nexus_data.EnterMatrix(distance, i, j);  //place the float in the dist_matrix
				nexus_data.EnterMatrix(distance, j, i);  //To make a square symetric matrix
				}
			}
		}

	/* If there is no diagonal and if there are labels in the data matrix */

	if ( (nexus_data.ReturnDiagonal() == 0) && (nexus_data.ReturnIsLabels() == 1) )
		{
		for (i=0; i < nexus_data.ReturnNTax(); i++)
			{
			aTokenPtr=nextToken();  // Get rid of the label
			for (j=(i+1); j < nexus_data.ReturnNTax() ; j++)
				{
				aTokenPtr=nextToken();
				distance = atof(aTokenPtr.c_str()); 	 //Transform the token in a float dummy
				nexus_data.EnterMatrix(distance, i, j);  //place the float in the dist_matrix
				nexus_data.EnterMatrix(distance, j, i);  //To make a square symetric matrix
				}
			nexus_data.EnterMatrix(0, i, i);  //place the "0" diagonal in the dist_matrix
			}		
		}

	/* If there is a diagonal and if there are no labels in the data matrix */

	if ( (nexus_data.ReturnDiagonal() == 1) && (nexus_data.ReturnIsLabels() == 0) )
		{
		for (i=0; i < nexus_data.ReturnNTax(); i++)
			{
			for (j=i; j < nexus_data.ReturnNTax(); j++)
				{
				aTokenPtr=nextToken();
				distance = atof(aTokenPtr.c_str()); 	 //Transform the token in a float dummy
				nexus_data.EnterMatrix(distance, i, j);  //place the float in the dist_matrix
				nexus_data.EnterMatrix(distance, j, i);  //To make a square symetric matrix
				}
			}
		}

	/* If there is no diagonal and if there are no labels in the data matrix */

	if ( (nexus_data.ReturnDiagonal() == 0) && (nexus_data.ReturnIsLabels() == 0) )
		{
		for (i=0; i < nexus_data.ReturnNTax(); i++)
			{
			for (j=(i+1); j < nexus_data.ReturnNTax() ; j++)
				{
				aTokenPtr=nextToken();
				distance = atof(aTokenPtr.c_str()); 	 //Transform the token in a float dummy
				nexus_data.EnterMatrix(distance, i, j);  //place the float in the dist_matrix
				nexus_data.EnterMatrix(distance, j, i);  //To make a square symetric matrix
				}
			nexus_data.EnterMatrix(0, i, i);  //place the "0" diagonal in the dist_matrix
			}		
		}
        }  // if (nexus_data.ReturnTriangle() == 1)

    if (nexus_data.ReturnTriangle() == 2)  // if the distance matrix is a square matrix
        {

        if (nexus_data.ReturnDiagonal() == 0)
            {
            cerr << "Exiting: A square matrix must have a diagonal..." << endl;
            exit(1);
            }

        /* If there are labels in the data matrix */

	if ( nexus_data.ReturnIsLabels() == 1 )
		{
		for (i=0; i < nexus_data.ReturnNTax(); i++)
			{
			aTokenPtr=nextToken();  // Get rid of the label
			for (j=0; j < nexus_data.ReturnNTax(); j++)
				{
				aTokenPtr=nextToken();
				distance = atof(aTokenPtr.c_str()); 	 //Transform the token in a float dummy
				nexus_data.EnterMatrix(distance, i, j);  //place the float in the dist_matrix
				}
			}
		}

	/* If there are no labels in the data matrix */

	if ( (nexus_data.ReturnDiagonal() == 1) && (nexus_data.ReturnIsLabels() == 0) )
		{
		for (i=0; i < nexus_data.ReturnNTax(); i++)
			{
			for (j=0; j < nexus_data.ReturnNTax(); j++)
				{
				aTokenPtr=nextToken();
				distance = atof(aTokenPtr.c_str()); 	 //Transform the token in a float dummy
				nexus_data.EnterMatrix(distance, i, j);  //place the float in the dist_matrix
				}
			}
		}

        }   // if (nexus_data.ReturnTriangle() == 2)

    aTokenPtr=nextToken();
    if (aTokenPtr != ";")
	{
	cerr << "  -> Matrix does not end by a semicolumn" << endl;
	cerr << "     or the size of matrix does not fit the number of taxa" << endl;
	exit(1);
	}
	nexus_data.SetIsDistanceMatrix(true);
    return;
}



/****************************************************************/
/****************  MISCELLANEOUS FUNCTIONS **********************/
/****************************************************************/


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
	string TempPtr;	// a pointer to manipulate aTokenPtr
	char c;
	
	c = *bufPtr++;

	if  (c == '\0') NULL_RETURN
	
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
			for(;;)
				{
				c=*bufPtr++;
				if (c=='\0') NULL_RETURN;   /* check for embedded EOF */
				if (c == '[')   /* In case there are brackets within braquets! */
					{
					while (c !=']')
						{
						c=*bufPtr++;
						if (c=='\0') NULL_RETURN;  /* check for embedded EOF */
						}
					continue;
					}
				if (c == ']')
					{
					break;
					}
				}
			c=*bufPtr++;	/* get next char after ']' */
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


Function: nextToken2()

Description: is used when dealing with interleaved 
             matrices.


***************************************************/

string nextToken2(void)
	{
	string TempPtr;	// a pointer to manipulate aTokenPtr
	char c;
	
	c = *bufPtr++;

	if  (c == '\0') NULL_RETURN
	
	while (( isNL_NEXUSwhiteSpace(c) ) || (c=='['))  
			    /* this whole loop is in case multiple comments separated by whitespace */
		{
		while ( isNL_NEXUSwhiteSpace(c) )  /* skip white space and land on next 'c'*/
			{
			c=*bufPtr++;
			if (c=='\0')  NULL_RETURN;/* check for embedded EOF */
			}
			    
		if (c=='[')		/* skip the comment and land on next 'c' after comment */
			{
			for(;;)
				{
				c=*bufPtr++;
				if (c=='\0') NULL_RETURN;   /* check for embedded EOF */
				if (c == '[')   /* In case there are brackets within braquets! */
					{
					while (c !=']')
						{
						c=*bufPtr++;
						if (c=='\0') NULL_RETURN;  /* check for embedded EOF */
						}
					continue;
					}
				if (c == ']')
					{
					break;
					}
				}
			c=*bufPtr++;	/* get next char after ']' */
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

        if ( ((c)=='\n') || ((c)=='\r') )
            {
            return "eol";
            }

	TempPtr.push_back(c);		/* char is either punctuation or part of word, so add it to token */
	    
	if (!isNL_NEXUSpunct(c))	/* next char is part of word, so add all of word until white,punct,eof,
								    or Token size exceeded */
	    {
	    for (;;)
		    {
		    c=*bufPtr++;
		    if (  isNL_NEXUSpunct(c) || isNL_NEXUSwhiteSpace(c) ||  (c == '\0') )
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


/***************************************************

Function: nextToken3

Description: Like nextToken3, but do not capitalize 
             characters (for reading matrices).

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

string nextToken3(void)
	{
	string TempPtr;	// a pointer to manipulate aTokenPtr
	char c;
	
	c = *bufPtr++;

	if  (c == '\0') NULL_RETURN
	
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
			for(;;)
				{
				c=*bufPtr++;
				if (c=='\0') NULL_RETURN;   /* check for embedded EOF */
				if (c == '[')   /* In case there are brackets within braquets! */
					{
					while (c !=']')
						{
						c=*bufPtr++;
						if (c=='\0') NULL_RETURN;  /* check for embedded EOF */
						}
					continue;
					}
				if (c == ']')
					{
					break;
					}
				}
			c=*bufPtr++;	/* get next char after ']' */
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

