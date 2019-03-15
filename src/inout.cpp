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

/*! \file inout.cpp 
  \brief Manage input/output.

 * We manage here the following operations:
 * * Read nominal instance from disk (both OR Library and Avella). See the
 *   introduction part of rcflp.cpp to see how the costs \f$c_{ij}\f$ are
 *   managed in the two instance types.
 * * Read parameters for the different support sets. Currently, we read
 *   parameters for the following sets:
 *   * Ellipsoidal support set. See read_parameters_ellipsoidal()
 *   * Box support set. See read_parameters_box()
 *   * Budget support set. See read_parameters_budget()
 *

*/

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <fstream>

using namespace std;


struct INSTANCE { /// See same data structure define in rcflp.cpp
    int nF;
    int nC;
    double  *f;
    double  *s;
    double  *d;
    double **c;
    double   totS;
    double   totD;

    int nR;        // number of constraints polyhedron uncertainty set
    double  *h;    // rhs of polyhedron definining support
    int     *W;    // matrix W in column major format
    int *index;    // index of column major format for w
    int *start;    // starting position for elements of column j
};

extern string instanceType;
extern string versionType;
extern string supportType;
extern int support;


/// Read benchmark instances
/**
 * Currently, two types of instances can be imported:
 * type 1: OR Library
 * type 2: Avella (Test Bed 1, Test Bed A. Test Bed B)
 */
int readProblemData(char * _FILENAME, int fType, INSTANCE & inp)
{
    inp.totS = 0.0;
    inp.totD = 0.0;
    ifstream fReader(_FILENAME, ios::in);
    if (!fReader)
    {
        cout << "cannot open file " << _FILENAME << endl;
        exit(1);
    }
    // read OR Library instances
    if (fType == 1)
    {
        
        fReader >> inp.nF >> inp.nC;
        inp.s = new double[inp.nF];
        inp.f = new double[inp.nF];
        inp.d = new double[inp.nC];
        inp.c = new double*[inp.nF];
        for (int i = 0; i < inp.nF; i++)
            inp.c[i] = new double[inp.nC];

        for (int i = 0; i < inp.nF; i++)
        {
            fReader >> inp.s[i] >> inp.f[i];
            inp.totS += inp.s[i];
        }

        for (int j = 0; j < inp.nC; j++)
        {
            fReader >> inp.d[j];
            inp.totD += inp.d[j];
        }

        for (int i = 0; i < inp.nF; i++)
            for (int j = 0; j < inp.nC; j++)
            {
                fReader >> inp.c[i][j];
                inp.c[i][j] /= inp.d[j];
            }

    }
    // read Avella instances
    else if (fType == 2)
    {
        fReader >> inp.nC >> inp.nF;
        inp.f = new double[inp.nF];
        inp.s = new double[inp.nF];
        inp.d = new double[inp.nC];
        inp.c = new double*[inp.nF];
        for (int i = 0; i < inp.nF; i++)
            inp.c[i] = new double[inp.nC];
        for (int j = 0; j < inp.nC; j++)
        {
            fReader >> inp.d[j];
            inp.totD += inp.d[j];
        }
        for (int i = 0; i < inp.nF; i++)
        {
            fReader >> inp.s[i];
            inp.totS += inp.s[i];
        }
        for (int i = 0; i < inp.nF; i++)
            fReader >> inp.f[i];

        for (int i = 0; i < inp.nF; i++)
            for (int j = 0; j < inp.nC; j++)
                fReader >> inp.c[i][j];
    }
    else
        cout << "Problem type not defined (-t option). Use '-h' for help. " << endl;

    fReader.close();

    return 1;
}

/// Print instance info and algorithmic parameters.
void printOptions(char * _FILENAME, INSTANCE inp, int timeLimit)
{
   double epsilon = 0.0;
   double delta   = 0.0;
   double gamma   = 0.0;
   int    L       = 0;
   cout << "-------------------------------------" << endl;
   cout << "- OPTIONS : " << endl;
   cout << "-------------------------------------" << endl;
   cout << "  DATA FILE      = " << _FILENAME        << endl;
   cout << "  Instance type  = " << instanceType << endl;
   cout << "  Version        = " << versionType << endl;
   cout << "  Support        = " << supportType << endl;
   cout << "  Time Limit     = " << timeLimit << endl;
   cout << "  Nr. Facilities = " << inp.nF << endl;
   cout << "  Nr. Customers  = " << inp.nC << endl;
   cout << "-------------------------------------" <<  endl << endl;   
}


/// Read parameters to define the ellipsoidal uncertainty set.
/** We need to define two parameters:
 *  * `Omega` : It determines the size of the ellipsoidal, i.e., the size of
 *              the uncertainty set.
 *  * `epsilon`: It defines the variance of each variable \f$x_{ij}\f$. We
 *               assume that the covariance matrix is diagonal and with
 *               constant values \f$\epsilon\f$ for each element.
 */
void read_parameters_ellipsoidal(double & epsilon, double & Omega)
{
    ifstream fReader("parameters/paramsEllipsoidal.txt", ios::in);
    if (!fReader)
    {
        cout << "Cannot open file 'paramsEllipsoidal.txt'." << endl;
        exit(111);
    }
    cout << "[** Reading parameters from file 'paramsEllipsoidal.txt'.]" << endl;
    string line;
    fReader >> epsilon; 
    getline(fReader, line); // read comment on same line
    fReader >> Omega;

    cout << "[** Uncertainty Set Parameters :: epsilon = " << epsilon 
         << "; Omega = " <<  Omega << "]\n" << endl;

    fReader.close();
}
/// Read parameters to define the Box support.
/**
 *  The file 'paramsBox.txt' has the following format:
 *      - first row : `epsilon`, i.e., the parameter used to define the with of
 *                    the box. Given a nominal demand value \f$d_j\f$, we define 
 *                    the interval around \f$d_j\f$ as:
 *                    \f[ (1-\epsilon)*d_j <= d_j <= (1+\epsilon)d_j, \quad 0
 *                    \leq \epsilon \leq 1 \f]
 */
void read_parameters_box(double & epsilon)
{
    ifstream fReader("parameters/paramsBox.txt", ios::in);
    if (!fReader)
    {
        cout << "Cannot open file 'paramsBox.txt'." << endl;
        exit(111);
    }

    cout << "[** Reading parameters from file 'paramsBox.txt'.]" << endl;
    fReader >> epsilon;
    cout << "[** Uncertainty Set Parameters :: epsilon = " << epsilon << "]" << endl;

    fReader.close();
}

/// Read parameters to define the budget support.
/**
 * * File : `paramsBudget.txt`
 *      - first row : `epsilon`, as above
 *      - second row: `delta`, i.e., how the \f$b_l\f$ value (the r.h.s. value of each
 *                    budget constraint) is defined. We define \f$b_l\f$ as a
 *                    percentage of the total demand of customers included in
 *                    the budget constraint, i.e.:
 *                    \f[ b_l = \delta* (\sum_{j \in B_l} d_j), \quad 0 \leq
 *                    \delta \leq 1 \f]
 *      - third row : `gamma`, i.e., percentage of columns included in each
 *                    budget constraint:
                      \f[|B_l| = \lfloor \gamma n \rfloor, \quad 0 \leq \gamma \leq 1\f]
 *      - forth row : \f$L \geq 1\f$, i.e., number of budget constraints.
*/
void read_parameters_budget(double & epsilon, double & delta, double & gamma, 
                            int & L)
{
    ifstream fReader("parameters/paramsBudget.txt", ios::in);
    if (!fReader)
    {
        cout << "Cannot open file 'paramsBudget.txt'." << endl;
        exit(111);
    }
    cout << "[** Reading parameters from file 'paramsBudget.txt'.]" << endl;
    string line;
    fReader >> epsilon; 
    getline(fReader, line); // read comment on same line
    /* line.erase( find( line.begin(), line.end(), '#' ), line.end() ); */

    fReader >> delta;
    getline(fReader, line); // read comment on same line
    fReader >> gamma;
    getline(fReader, line); // read comment on same line
    fReader >> L;
    cout << "[** Uncertainty Set Parameters :: epsilon = " << epsilon 
         << "; delta = " << delta << "; gamma = " << gamma << "; L = " << L 
         << "]" << endl;

    fReader.close();
}
