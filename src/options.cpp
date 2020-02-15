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

    Command line options are:

    - **-h** : help (visualize the list of options)

    - **-i** : inputfile

    - **-t** : instance type:
            -# orLibrary
            -# Avella

    - **-v** : version:
            -# Single-source
            -# Multi-source
            -# Ellipsoidal
            -# Polyhedral

    - **-u** : uncertainty set
            -# Box uncertainty set
            -# Budget uncertanty set

    - **-r** : read from disk
            -# 0 No: A new Budget set \f$B_l\f$ is generated and stored
            -# 1 Yes: The Budget set is read from disk

    - **-e** : epsilon, i.e., the size of the uncertainty set used during the 
                optimization phase. Typically, \f$\varepsilon \in \left[0.1, 1\right]\f$.

    <b>Note:</b> The meaning of the flag `-r` is as follows. When a budget 
            uncertainty set is used, we randomly select which customers, i.e., 
            columns, will be included in the budget constraint. Next, this set, \f$B_l\f$,
            is used to find a robust optimal solution for a given \f$\varepsilon\f$.
            If we want to compare solutions obtained for different values of \f$\varepsilon\f$, 
            we need to ensure that the same budget constraint(s) is(are) used in the
            optimziation phase. Otherwise, the comparison would be meaningless.
*/

#include <iostream>
#include <cstdlib>
/**********************************************************/
#define   _TIMELIMITdef  3600   //!< default wall-clock time limit
#define   _VERSIONdef    1      //!< single source by default
#define   _FROMDISKdef    1      //!< single source by default
/**********************************************************/

using namespace std;

extern char* _FILENAME; 	//!< name of the instance file
extern int timeLimit;		//!< wall-clock time limit
extern int fType;           //!< instance type (1-2)
extern int version;         //!< 1-SS; 2-MS; 3-Ellipsoidal; 4-Polyhedral
extern int support;         //!< 1-Box; 2-Budget
extern int readFromDisk;    //!< 0-No; (Generate a new Budget set B_l); 1-Yes
extern string instanceType;
extern string versionType;
extern string supportType;
extern double _epsilon;


/// Parse command line options.
/** Use -h to visualize the list of options.
 */
int parseOptions(int argc, char* argv[])
{
   bool setFile = false;
   bool setType = false;
   timeLimit    = _TIMELIMITdef;
   version      = _VERSIONdef;
   readFromDisk= _FROMDISKdef;
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
        case 'u':
	       support = atol(argv[i+1]);
	       i++;
	       break;
        case 'r':
	       readFromDisk = atol(argv[i+1]);
	       i++;
	       break;
        case 'e':
	       _epsilon= atof(argv[i+1]);
	       i++;
	       break;
	    case 'h':
	       cout << "OPTIONS :: " << endl;
	       cout << "-i : problem instance file" << endl;
	       cout << "-l : time limit (real)" << endl;
	       cout << "-v : problem version (1-SS; 2-MS; 3-SOCP; 4- Poly; 5 - Bin Packing; 6- MMCF)" << endl;
	       cout << "-t : instance type (1-OR Library; 2-Avella). Relevant only if CFLP" << endl;
	       cout << "-u : support type (1-Box; 2-Budget)" << endl;
	       cout << "-r : read Budget support set from disk (0-No; 1-Yes)" << endl;
	       cout << "-e : epsilon (for box and budget support) " << endl;
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
        else if (fType == 0)
            instanceType = "Original OR Lib";
        else 
            instanceType = "***";

        // if (version == 1)
            // versionType = "Single-source-Nominal";
        // else if (version == 2)
            // versionType = "Multi-source-Nominal";
        // else if (version == 3)
            // versionType = "Multi-source-Ellipsoidal";
        // else if (version == 4)
            // versionType = "Polyhedral Uncertainty (Wd <= h)";
//
        // if (support == 1)
            // supportType = "Box Uncertainty Set";
        // else if (support == 2)
            // supportType = "Budget Uncertainty Set";

        return 0;
   }
   else
   {
      cout <<"Options -i and -t are mandatory. Try ./rcflp -h" << endl;
      return -1;
   }
}
