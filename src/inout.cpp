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
#include <vector>
#include "options.h"
#include "binPacking.h"
#include "mmcf.h"

using namespace std;

struct INSTANCE {
    int nF;        //!< Number of facilities
    int nC;        //!< Number of customers
    double  *f;    //!< Fixed costs
    double  *s;    //!< Capacity
    double  *d;    //!< Demand
    double **c;    //!< Allocation costs
    double   totS; //!< Total supply
    double   totD; //!< Total demand

    int     nR;    //!< Number of constraints polyhedron uncertainty set
    double  *h;    //!< Rhs of polyhedron definining support
    int     *W;    //!< Matrix W in column major format
    int *index;    //!< Index of column major format for w
    int *start;    //!< Starting position for elements of column j

    std::vector <int> listB; //!< List of columns in budget constraints
    int * customer_type; //!< 1, 2, 3 (for Roberto B. instances, otherwise not used)
    int  *nType;
};
extern double _Omega;
extern double _epsilon;
extern double _delta;
extern double _gamma;
extern double _beta0;
extern double _beta1;
extern double _beta2;
extern int    L;      
double old_epsilon;  // _epsilon is now read from command line

// extern string instanceType;
// extern string versionType;
// extern string supportType;
extern int support;
extern int version;




int readMMCF(char* _FILENAME, InstanceMMCF& inpMMCF)
{
    int origin, destin, source, sink;
    double cost, cap, fixed, temp1, temp2, dem;

    ifstream fReader(_FILENAME, ios::in);
    if (!fReader) {
        cout << "cannot open file " << _FILENAME << endl;
        exit(1);
    }
    string ss;
    getline(fReader, ss);

    fReader >> inpMMCF.nNodes >> inpMMCF.nArcs >> inpMMCF.nK;

    std::vector<int> aux;
    for (int j = 0; j < inpMMCF.nNodes; j++) {
        inpMMCF.outbound.push_back(aux);
        inpMMCF.inbound.push_back(aux);
    }
    
    for (int j = 0; j < inpMMCF.nArcs; j++) {
        fReader >> origin >> destin >> cost >> cap >> fixed >> temp1 >> temp2;
        origin -= 1;
        destin -= 1;
        inpMMCF.origin.push_back(origin);
        inpMMCF.destin.push_back(destin);
        inpMMCF.c.push_back(cost);
        inpMMCF.cap.push_back(cap);
        inpMMCF.fixed.push_back(fixed);

        inpMMCF.outbound[origin].push_back(j);
        inpMMCF.inbound[destin].push_back(j);

    }

    for (int j = 0; j < inpMMCF.nK; j++) {
        fReader >> source >> sink >> dem;
        source -= 1;
        sink   -= 1;
        inpMMCF.source.push_back(source);
        inpMMCF.sink.push_back(sink);
        inpMMCF.d.push_back(dem);
    }

    /* for (int j = 0; j < inpMMCF.nArcs; j++) 
     *     cout << "Arc " << j << " from " << inpMMCF.origin[j] << " to " << inpMMCF.destin[j] << endl; */
    cout << "MMCF Instance " << endl;
    cout << "---------------------" << endl;
    cout << "Tot Nodes \t :  " << inpMMCF.nNodes << endl;
    cout << "Tot Arcs  \t :  " << inpMMCF.nArcs << endl;
    cout << "Tot Commodities\t :  " << inpMMCF.nK << endl;
    cout << "---------------------" << endl;
        
    return 1;

}

int readBinPacking(char * _FILENAME, InstanceBin& inpBin)
{
    ifstream fReader(_FILENAME, ios::in);
    if (!fReader)
    {
        cout << "cannot open file " << _FILENAME << endl;
        exit(1);
    }
    string ss;
    getline(fReader, ss);
    getline(fReader, ss);

    fReader >> inpBin.totS >> inpBin.nI >> inpBin.opt;
    double temp;
    for (int j = 0; j < inpBin.nI; j++)
    {
        fReader >> temp;
        inpBin.d.push_back(temp);
    }


    cout << "BIN Packing Instance " << endl;
    cout << "---------------------" << endl;
    cout << "Tot Capacity \t :  " << inpBin.totS << endl;
    cout << "Tot Items  \t :  " << inpBin.nI << endl;
    cout << "Optimum    \t :  " << inpBin.opt << endl;
    cout << "---------------------" << endl;

    return 1;

}


/// Read benchmark instances
/**
 * Currently, two types of instances can be imported:
 * type 1: OR Library
 * type 2: Avella (Test Bed 1, Test Bed A. Test Bed B)
 */
int readCFLP(char * _FILENAME, int fType, INSTANCE & inp)
{
    char * iType = "";
    double temp = 0.0;
    inp.totS = 0.0;
    inp.totD = 0.0;
    ifstream fReader(_FILENAME, ios::in);
    if (!fReader)
    {
        cout << "cannot open file " << _FILENAME << endl;
        exit(1);
    }
    // read OR Library instances (official format : email Roberto B., 20.08.19)
    // June 2020: Added part used to read customer type
    if (fType == 0)
    {
        iType = "Baldacci"; 
        fReader >> inp.nF >> inp.nC;
        inp.s = new double[inp.nF];
        inp.f = new double[inp.nF];
        inp.d = new double[inp.nC];
        inp.c = new double*[inp.nF];
        inp.customer_type = new int[inp.nC];
        for (int i = 0; i < inp.nF; i++)
            inp.c[i] = new double[inp.nC];

        for (int i = 0; i < inp.nF; i++)
        {
            fReader >> inp.s[i] >> inp.f[i];
            // inp.s[i] *= 2.5; // <------ REMOVE THIS !!!! Roberto Instances!!!
            inp.totS += inp.s[i];
            // cout << "s[" << i << "] = " <<  inp.s[i] << endl;
        }

        double maxD = 0.0;
        for (int j = 0; j < inp.nC; j++)
        {
            fReader >> inp.d[j];
            inp.totD += inp.d[j];
            if (inp.d[j] > maxD)
                maxD = inp.d[j];
            for (int i = 0; i < inp.nF; i++)
            {
                fReader >> inp.c[i][j];
                inp.c[i][j] /= inp.d[j];
            }
        }

        cout << "** Tot S vs Tot D " << inp.totS << ", " << inp.totD << endl;
        cout << "Max D " << maxD << endl;

        // skip coordinates (for each facility and each customer)
        for (int i = 0; i < inp.nF; i++)
            fReader >> temp >> temp;

        for (int j = 0; j < inp.nC; j++)
            fReader >> temp >> temp;

        // read customer type (1,2, or 3) and rescale to 0, 1, 2
        inp.nType = new int[3];
        for (int j = 0; j < inp.nC; j++) {
            fReader >> inp.customer_type[j];
            inp.customer_type[j] -= 1;
            inp.nType[inp.customer_type[j]]++;
        }
        cout << "Total number of customers in each budget constraint :: ";
        for (int l = 0; l < 3; l++)
            cout << " " << inp.nType[l];
        cout << endl;
        

        /* for (int j = 0; j < inp.nC; j++)
         * {
         *     cout << "Customer " << j << " d= " << inp.d[j] << " Costs = ";
         *     for (int i = 0; i < inp.nF; i++)
         *         cout << " " << inp.c[i][j]*inp.d[j];
         *     cout << endl;
         * } */
    }
    else if (fType == 1)
    {
        
        iType = "OR-Library";
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


        // 03.08.19 ::
        // This part is needed to save the instances in a slightly different
        // format. The demands and transportation costs are written using
        // end-of-line.
        // While in c++ the format is not a problem, if we want to read these
        // instances in python (for example, to run the scenario simulation
        // using the rKnap.py code), we need to have an "end-of-line" after
        // all the demands are written, and an end-of-line after each supplier
        // (see list for cycle below.)
        //
        // ofstream fWriter("temp.txt", ios::out);
        // fWriter << inp.nF << " " << inp.nC << "\n";
        // for (int i = 0; i < inp.nF; i++)
            // fWriter << inp.s[i] << " " << inp.f[i] << "\n";
//
        // for (int j = 0; j < inp.nC; j++)
            // fWriter << inp.d[j] << " ";
        // fWriter << "\n";
        // for (int i = 0; i < inp.nF; i++)
        // {
            // for (int j = 0; j < inp.nC; j++)
                // fWriter << inp.c[i][j] <<" ";
            // fWriter << "\n";
        // }
        // fWriter.close();
        // exit(123);



    }
    // read Avella instances
    else if (fType == 2)
    {
        iType = "Avella";
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

    cout << "CFLP Instance " << endl;
    cout << "---------------------" << endl;
    cout << "Tot Facilities \t :  " << inp.nF << endl;
    cout << "Tot Customers  \t :  " << inp.nC << endl;
    cout << "Instance Type  \t :  " << iType << endl;
    if (fType == 0) {
        cout << "Beta0 = \t : " << _beta0 << endl;
        cout << "Beta1 = \t : " << _beta1 << endl;
        cout << "Beta2 = \t : " << _beta2 << endl;
    }
        
    cout << "---------------------" << endl;

    return 1;
}

/// Print instance info and algorithmic parameters.
void printOptions(char * _FILENAME, int timeLimit)
{
   string versionType = versionLabel[version-1];
   string supportType  = supportLabel[support-1];
   cout << "-------------------------------------" << endl;
   cout << "- OPTIONS : " << endl;
   cout << "-------------------------------------" << endl;
   cout << "  DATA FILE      = " << _FILENAME        << endl;
   cout << "  Version        = " << versionType << endl;
   cout << "  Support        = " << supportType << endl;
   cout << "  Time Limit     = " << timeLimit << endl;
   cout << "-------------------------------------" <<  endl << endl;   
}


/// Read parameters to define the ellipsoidal uncertainty set.
/** We need to define two parameters:
 *  * `_Omega` : It determines the size of the ellipsoidal, i.e., the size of
 *              the uncertainty set.
 *  * `_epsilon`: It defines the variance of each variable \f$x_{ij}\f$. We
 *               assume that the covariance matrix is diagonal and with
 *               constant values \f$\epsilon\f$ for each element.
 */
// void read_parameters_ellipsoidal(double & _epsilon, double & _Omega)
void read_parameters_ellipsoidal()
{
    ifstream fReader("parameters/paramsEllipsoidal.txt", ios::in);
    if (!fReader)
    {
        cout << "Cannot open file 'paramsEllipsoidal.txt'." << endl;
        exit(111);
    }
    cout << "[** Reading parameters from file 'paramsEllipsoidal.txt'.]" << endl;
    string line;
    fReader >> old_epsilon; 
    getline(fReader, line); // read comment on same line
    fReader >> _Omega;

    cout << "[** Uncertainty Set Parameters :: _epsilon = " << _epsilon 
         << "; _Omega = " <<  _Omega << "]\n" << endl;

    fReader.close();
}
/// Read parameters to define the Box support.
/**
 *  The file 'paramsBox.txt' has the following format:
 *      - first row : `_epsilon`, i.e., the parameter used to define the with of
 *                    the box. Given a nominal demand value \f$d_j\f$, we define 
 *                    the interval around \f$d_j\f$ as:
 *                    \f[ (1-\_epsilon)*d_j <= d_j <= (1+\_epsilon)d_j, \quad 0
 *                    \leq \epsilon \leq 1 \f]
 *
 *  MARCO 11.02.20: This is no longer used, since the only parameter of a box
 *  set, i.e., epsilon, is now read from command line with flag -e.
 */
// void read_parameters_box(double & _epsilon)
void read_parameters_box()
{
    ifstream fReader("parameters/paramsBox.txt", ios::in);
    if (!fReader)
    {
        cout << "Cannot open file 'paramsBox.txt'." << endl;
        exit(111);
    }

    cout << "[** Reading parameters from file 'paramsBox.txt'.]" << endl;
    fReader >> old_epsilon;
    cout << "[** Uncertainty Set Parameters :: _epsilon = " << _epsilon << "]" << endl;

    fReader.close();
}

/// Read parameters to define the budget support.
/**
 * * File : `paramsBudget.txt`
 *      - first row : `_epsilon`, as above
 *      - second row: `_delta`, i.e., how the \f$b_l\f$ value (the r.h.s. value of each
 *                    budget constraint) is defined. We define \f$b_l\f$ as a
 *                    percentage of the total demand of customers included in
 *                    the budget constraint, i.e.:
 *                    \f[ b_l = \delta* (\sum_{j \in B_l} d_j), \quad 0 \leq
 *                    \delta \leq 1 \f]
 *      - third row : `_gamma`, i.e., percentage of columns included in each
 *                    budget constraint:
                      \f[|B_l| = \lfloor \_gamma n \rfloor, \quad 0 \leq \gamma \leq 1\f]
 *      - forth row : \f$L \geq 1\f$, i.e., number of budget constraints.
*/
// void read_parameters_budget(double & _epsilon, double & _delta, double & gamma,
void read_parameters_budget()
{
    ifstream fReader("parameters/paramsBudget.txt", ios::in);
    if (!fReader)
    {
        cout << "Cannot open file 'paramsBudget.txt'." << endl;
        exit(111);
    }
    cout << "[** Reading parameters from file 'paramsBudget.txt'.]" << endl;
    string line;
    // MARCO 11.02.20: The next two lines are commented out, since epsilon is
    // now read from command line with flag -e
    // fReader >> old_epsilon;
    // getline(fReader, line); // read comment on same line
    /* line.erase( find( line.begin(), line.end(), '#' ), line.end() ); */

    fReader >> _delta;
    getline(fReader, line); // read comment on same line
    fReader >> _gamma;
    getline(fReader, line); // read comment on same line
    fReader >> L;
    cout << "[** Uncertainty Set Parameters :: _epsilon = " << _epsilon 
         << "; _delta = " << _delta << "; _gamma = " << _gamma << "; L = " << L 
         << "]" << endl;

    fReader.close();
}
