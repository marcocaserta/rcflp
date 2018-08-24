/***************************************************************************
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

/*! \file options.cpp 
  \brief Read options from command line.

    Command line options (see parseCommandLine):
    -i inputfile

    -t type:
    1 -> orLibrary
    2 -> Avella
    3 -> SS-orLibrary
    4 -> Robust orLibrary

    -v version (1 -> Single-source; 2 -> Multi-source)

*/

#include <iostream>
#include <cstdlib>
/**********************************************************/
#define   _TIMELIMITdef  3600   //!< default wall-clock time limit
#define   _VERSIONdef    1     //!< single source by default
/**********************************************************/

using namespace std;

extern char* _FILENAME; 	//!< name of the instance file
extern int timeLimit;		//!< wall-clock time limit
extern int fType;           //!< instance type (1-4)
extern int version;         //!< 1-SS; 2-MS
extern string instanceType;
extern string versionType;


/// Parse command line options
int parseOptions(int argc, char* argv[])
{
   bool setFile = false;
   bool setType = false;
   timeLimit    = _TIMELIMITdef;
   version      = _VERSIONdef;
   cout <<endl << "R-CLSP v1.0 " << endl;
   if (argc == 1)
   {
      cout << "No options specified. Try -h " << endl;
      return -1;
   }  
  
 
   int i = 0;
   while (++i < argc)
   {
      const char *option = argv[i];
      if (*option != '-')
	 return i;
      else if (*option == '\0')
	 return i;
      else if (*option == '-')
      {
	 switch (*++option)
	 {
	    case '\0':
	       return i + 1;
	    case 'i':
	       _FILENAME = argv[i+1];
	       setFile = true;
	       i++;
	       break;
	    case 't':
	       fType = atol(argv[i+1]);
           setType = true;
	       i++;
	       break;
	    case 'l':
	       timeLimit = atol(argv[i+1]);
	       i++;
	       break;
	    case 'v':
	       version = atol(argv[i+1]);
	       i++;
	       break;
	    case 'h':
	       cout << "OPTIONS :: " << endl;
	       cout << "-i : problem instance file" << endl;
	       cout << "-l : time limit (real)" << endl;
	       cout << "-v : problem version (1-SS; 2-MS; 3-SOCP; 4- Poly)" << endl;
	       cout << "-t : instance type (1-OR Library; 2-Avella)" << endl;
	       cout << endl;
	       return -1;
	 }
      }
   }
 
   if (setFile && setType)
   {
        if (fType == 1)
            instanceType = "OR Library";
        else if (fType == 2)
            instanceType = "Avella";
        else
            instanceType = "***";

        if (version == 1)
            versionType = "Single-source-Nominal";
        else if (version == 2)
            versionType = "Multi-source-Nominal";
        else if (version == 3)
            versionType = "Multi-source-Ellipsoidal";
        else if (version == 4)
            versionType = "Polyhedral Uncertainty (Wd <= h)";

        return 0;
   }
   else
   {
      cout <<"Options -i and -t are mandatory. Try ./rcflp -h" << endl;
      return -1;
   }
}
