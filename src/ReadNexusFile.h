/*

    (c) Copyright 2002 by Michael J. Sanderson.

    Permission is granted to copy and use this program provided 
    no fee is charged for it and provided that this copyright 
    notice is not removed.


    This function was translated and modified from the function
    ReadNexusFile2 written by Michael J. Sanderson and distributed
    in the program r8s. All problems related to this function and
    program should be attributed to the present version, and not
    to the original one.

    (c) Copyright 2005-2008 by Simon Joly

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


#include <vector>
#include <string>
#include <string.h>
#include <ctype.h>
#include <iostream>
#include <stdio.h>

#include "nexusdata.h"
#include "GenFunctions.h"

//The following was in nextToken.h ...
#define punct "()[]{}/\\,;:=*\'\"`+"			/* these are NEXUS definitions */
#define NL_punct "()[]{}/\\,;:=*\'\"`+\r\n"             /*                             */

#define isNL_NEXUSwhiteSpace(c)  ( strchr(" \t\v\f", (c)) || (((c) <= 6) && ((c) >= 0)))
#define isNL_NEXUSpunct(c) ( strchr(NL_punct,(c)) )
#define isNEXUSpunct(c) ( strchr(punct,(c)) )
#define isNEXUSwhiteSpace(c)	( isspace((c)) )
	/* current NEXUS format does not exclude ASCII 0-6 */
#define NULL_RETURN {TempPtr=""; return(TempPtr);}
#define CHECK_OVERFLOW  if (cix>=MAX_TOKEN_SIZE-1) doGenericAlert("Token Size Exceeded in nextToken")

using namespace std;

extern nexusdata nexus_data; /* This is THE class for the NEXUS data */

/**** Main function ****/

void readNexusFile(char *);
