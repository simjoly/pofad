#include "GenFunctions.h"

using namespace std;

/****  Miscellaneous utility commands ****/


/**********************
*   Function: Fatal   *
**********************/

void fatal(char *s)
{
	long i;
	cout << endl << "!!!!!!!!!!!!!!!!!!! FATAL ERROR !!!!!!!!!!!!!!!!!!!!" << endl << endl;
	doGenericAlert(s);
	pauseprogram();
	for (i=1;i<100000;i++);
	exit(1);
	return;
}


/*********************************

Function: doGenericAlert

*********************************/

void doGenericAlert(char* errorMsg)
{
cout << "********************* WARNING **********************" << endl << endl;
cout << errorMsg << endl << endl;
cout << "****************************************************" << endl;
return;
}



/***********************************************************

Function : strtoupper

Definition: Convert a c-style string to a capitalized string

Variables:

***********************************************************/


void strtoupper(char *s)
{
	char *temps;
	temps=s; 
	while(  *temps ) 
		{
		*temps=toupper(*temps); 
		++temps;
		}
	/* converts string to upper case */
	return;
}



/******************************************************

Function: To_Uppercase()

Takes a string as input and outputs an uppercase string

******************************************************/

string To_Uppercase(string input_string)
{
    int i;
    char character;
    string output_string, temp_string;
    temp_string = input_string;
    int size_of_string = temp_string.size();

    for (i=0; i < size_of_string; i++)
    {
        character = temp_string[i];
        output_string += toupper(character);
    }
    return output_string;
}

/******************************************************

Function: clearscreen()

Takes a string as input and outputs an uppercase string

******************************************************/

void clearscreen()
{

#ifdef __WIN32__
	system("cls");
#else
	cout << string( 40, '\n' );
#endif
}



/******************************************************

Function: pauseprogram()

Takes a string as input and outputs an uppercase string

******************************************************/

void pauseprogram(void)
{
	cout << endl << " Press enter to continue...";
	cin.get();
}

