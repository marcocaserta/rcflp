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

/*! \mainpage Robust Capacitated Facility Location Problem

  Description here.

  \authors
  \version v. 1.0.0
  \date Begins: 25.05.18
  \date Ends:

  The project is compiled using a makefile and run via command line:
  ~~~
  make
  ./bin/rcflp -h
  ~~~

  The flag`-h` allows to visualize the options available via command line.

  The basic organization of the files is as follows:
  - options.cpp: Managing the command line. 
  - inout.cpp: Managing the input/output. Here we both read the instance and 
               define the support sets.
  - rcflp.cpp: Main implementation of the CFLP models.

  \file rcflp.cpp
  \brief General Implementation of the compact formulations for the (R)-CFLP.

  ### Different Versions of the CFLP

  We define here the following variants of the CFLP. For a full description of 
  the different variants, along with other parameters, see options.cpp.

  Current versions of the CFLP implemented here are (controlled via the command 
  line flag **-v**):
  - Single Source Nominal: see define_SS_CFLP()
  - Multi Source Nominal: see define_MS_CFLP()
  - Ellipsoidal Support Set (multi-source only?): see define_SOCP_CFLP()
  - Polyhedra Support Set (both single and multi-source): see define_POLY_CFLP()

  ### Instance Types

  In addition, we need to define the type of instance used. We are currently
  working with the following instances (controlled via the command line 
  flag **-t**):
  - 1 : OR Library instances
  - 2 : Avella instances (Type 1, Type A, Type B)

  __Note__: These instances define the costs in different ways and, therefore
  the way in which the \f$x_{ij}\f$ variables are defined changes. More precisely:
  - For the OR Library instances, the \f$c_{ij}\f$ values are the cost of delivering
  the **entire demand** to the customer. Thus, we define \f$x_{ij}\f$ as fractional 
  and multiply \f$c_{ij}\times x_{ij}\f$ to get the real cost.
  - For the Avella instances, \f$c_{ij}\f$ is the cost of transporting **one unit**. 
  Thus, \f$x_{ij}\f$ is no longer the fraction of demand, but the actual number 
  of units transported.

  Obviously, this different treatment of the \f$x_{ij}\f$ variables has an effect on 
  the definition of the demand and capacity constraints as well.
  Therefore, to make the treatment of the two types of instances uniform, the
  cost \f$c_{ij}\f$ of the OR Library is divided by the total demand:
  \f[
    c_{ij} = \frac{c_{ij}}{d_j}
  \f]

  To summarize, after the transformation, for both instance types (OR Library and Avella), 
  we have:
  - \f$x_{ij} \in [0,1]\f$ : percentage of demand
  - \f$c_{ij}\f$ : cost of delivering one unit of demand
  - \f$d_j\f$  : total demand of customer \f$j\f$

  This allows to define the demand constraint as:
  \f[
    \sum_{i=1}^m x_{ij} = 1
  \f]
  and the capacity constraint as:
  \f[
    \sum_{j=1}^m d_j x_{ij} \leq s_iy_i
  \f]
  Finally, the transportation cost is computed as:
  \f[
    \sum_{i=1}^m \sum_{j=1}^n d_jx_{ij}c_{ij}
  \f]

  ### MIP Solver

  We use __cplex branch-and-cut__ to solve the different programs. 

  In addition, we can also use __cplex Benders__ as implemented in this version 
  of the solver. To activate Benders, see define_benders().
  
  See the makefile to determine how 
  to link the library to the code.


  03.08.19: Added code to solve the Instances provided by Roberto B. They are in
  the same format of the OR Library. However, to solve them in robust format, 
  i.e., with epsilon > 0, we need to increase the capacity of each facility.

  07.08.19: Modified code to read epsilon from command line. This is useful to
  make the tests automatic. Note that, in the case of the budget constraints, 
  all the other parameters are still read from disk.

*/

#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

#include <limits> 

#include <chrono>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <cstring>
#include <string>
#include <cmath>
#include <climits>
#include <cfloat>
#include <cassert>
#include <random>
#include <algorithm>
#include <cstdlib>
#include <sstream>

/* #include "timer.h" */
#include "options.h"

using namespace std;

double INFTY = std::numeric_limits<double>::infinity();
const long _MAXRANDOM  = 2147483647;
const double ZERO      = 0.0e0;
const double EPSI      = 0.00001;
const int seed = 27; //!< Initialization of Random Number Generator
mt19937_64 gen(seed); //!< 64-bit Mersenne Twister by Matsumoto and Nishimura, 2000
/* mt19937_64 gen(seed()); */

/****************** VARIABLES DECLARATION ***************************/
char * _FILENAME;		//!< Instance name file
int fType;              //!< instance type (1-4)
int version;            //!< 1-SS; 2-MS; 3-SOCP
int support;            //!< 1-Box; 2-Budget
int readFromDisk;       //!< 0-No; (Generate a new Budget set B_l); 1-Yes
double _epsilon;        //!< epsilon used in robust optimization (box and budget)
string instanceType;
string versionType;
string supportType;

typedef std::vector <int> MyVect;
double _Omega;
double _delta;
double _gamma;
int    L;      

/// Structure used to define the instance data
// NOTE: Change the same structure in the file inout.cpp !!!
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

};
INSTANCE inp; //!< Instance data

/// Optimal solution and obj function value
struct SOLUTION {
    int nOpen;
    int  *ySol;
    double **xSol;
    double zStar;
    IloAlgorithm::Status zStatus;
    IloNum startTime;
    IloNum cpuTime;
};
SOLUTION opt; //!< Solution data structure


/**** CPLEX DEFINITION ****/
typedef IloArray <IloNumVarArray> TwoD;
IloEnv env;
IloModel model(env, "cflp");
/* IloCplex cplex(model); */
TwoD x_ilo;
IloNumVarArray y_ilo;
IloNumVarArray q_ilo;
IloNumVar w_ilo;
TwoD psi_ilo;
IloNumVarArray u_ilo;
IloNumVarArray delta_ilo;
int solLimit     = 9999;
int displayLimit = 4;
int timeLimit;

/****************** FUNCTIONS DECLARATION ***************************/
int readProblemData(char * _FILENAME, int fType, INSTANCE & inp);
void printOptions(char * _FILENAME, INSTANCE & inp, int timeLimit);
void define_MS_CFLP(INSTANCE inp, int fType, IloModel & model, IloCplex & cplex);
void define_SS_CFLP(INSTANCE inp, int fType, IloModel & model, IloCplex & cplex);
void define_SOCP_CFLP(INSTANCE inp, int fType, IloModel & model, IloCplex & cplex);
void define_POLY_CFLP(INSTANCE inp, int fType, IloModel & model, IloCplex & cplex, int support);
int solveCplexProblem(IloModel model, IloCplex cplex, INSTANCE inp, int solLimit, int timeLimit, int displayLimit);
void getCplexSol(INSTANCE inp, IloCplex cplex, SOLUTION & opt);
void printSolution(char * _FILENAME, INSTANCE inp, SOLUTION opt, bool toDisk, int fullOutput);
/* void read_parameters_box(double & _delta); */
void read_parameters_box();
void read_parameters_ellipsoidal();
void read_parameters_budget();
void define_box_support(INSTANCE & inp);
void define_budget_support(INSTANCE & inp, bool fromDisk);
void save_instance_2_disk(double _epsilon, double _delta, double _gamma, int L, 
                          int nBl, int ** Bl, double * budget);
void read_instance_from_disk(double & _epsilon, double & _delta, double & _gamma, 
                             int & L, int & nBl, int ** Bl, double * budget);
void define_benders(IloModel & model, IloCplex & cplex, INSTANCE inp);
void define_new_budget_support(INSTANCE & inp, bool fromDisk);
void define_new2_budget_support(INSTANCE & inp, bool fromDisk);
void define_new_POLY_CFLP(INSTANCE inp, int fType, IloModel & model, IloCplex & cplex, int support);
void define_new2_POLY_CFLP(INSTANCE inp, int fType, IloModel & model, IloCplex & cplex, int support);
/****************** FUNCTIONS DECLARATION ***************************/

/************************ main program ******************************/
/// main program
/************************ main program ******************************/
int main(int argc, char *argv[])
{
    // _epsilon = 0.0;
    _delta   = 0.0;
    _gamma   = 0.0;
    L        = 0;
    std::srand ( unsigned ( std::time(0) ) );

    int err = parseOptions(argc, argv);
    if (err != 0) exit(1);

    readProblemData(_FILENAME, fType, inp);
    printOptions(_FILENAME, inp, timeLimit);

    auto start = chrono::system_clock::now();

    IloCplex cplex(model);
    switch(version)
    {
        case 1 :  // single source nominal
            define_SS_CFLP(inp, fType, model, cplex);
            break;
        case 2 : // multi source nominal
            define_MS_CFLP(inp, fType, model, cplex);
            break;
        case 3 : // multi source ellipsoidal
            define_SOCP_CFLP(inp, fType, model, cplex);
            break;
        case 4 : // robust polyhedral uncertainty set (both SS and MS)
            // define_POLY_CFLP(inp, fType, model, cplex, support);
            // define_new_POLY_CFLP(inp, fType, model, cplex, support);
            define_new2_POLY_CFLP(inp, fType, model, cplex, support);
            break;
        default :
            cout << "ERROR : Version type not defined.\n" << endl;
            exit(123);
    }

    // define_benders(model, cplex, inp);

    solveCplexProblem(model, cplex, inp, solLimit, timeLimit, displayLimit);

    opt.cpuTime = chrono::duration_cast<chrono::seconds>(chrono::system_clock::now()-start).count();

    getCplexSol(inp, cplex, opt);
    printSolution(_FILENAME, inp, opt, true, 1);


    env.end();
    return 0;
}
/************************ main program ******************************/
/// END main program
/************************ main program ******************************/

/****************** FUNCTIONS DEFINITION ***************************/
/// Get and store cplex solution in data structure opt
void getCplexSol(INSTANCE inp, IloCplex cplex, SOLUTION & opt)
{
    opt.nOpen = 0;
    opt.ySol = new int[inp.nF];
    opt.xSol = new double*[inp.nF];
    for (int i = 0; i < inp.nF; i++)
        opt.xSol[i] = new double[inp.nC];

    opt.zStar = cplex.getObjValue();
    opt.zStatus = cplex.getStatus();

    for (int i = 0; i < inp.nF; i++)
        if (cplex.getValue(y_ilo[i]) >= 1.0-EPSI)
        {
            opt.ySol[i] = 1;
            opt.nOpen++;
        }
        else
            opt.ySol[i] = 0;

    for (int i = 0; i < inp.nF; i++)
        for (int j = 0; j < inp.nC; j++)
            opt.xSol[i][j] = cplex.getValue(x_ilo[i][j]);

    switch(version)
    {
        case 1 :  // single source nominal
            versionType = "ss";
            break;
        case 2 : // multi source nominal
            versionType = "ms";
            break;
        case 3 : // multi source ellipsoidal
            versionType = "ellipsoidal";
            break;
        case 4 : // robust polyhedral uncertainty set (both SS and MS)
            if (support==1)
                versionType = "box";
            else
                versionType = "budget";
            break;
        default :
            cout << "ERROR in WRITING SOLUTION: Version type not defined.\n" << endl;
            exit(123);
    }

    string  s1      = string(_FILENAME);
    s1              = s1.substr(s1.find_last_of("\\/"), 100);
    ostringstream obj;
    obj << _epsilon << "-" <<  _delta << "-" << _gamma << "-" << L;
    string filename = "solutions" + s1 + "-" + versionType + "-" + obj.str();

    ofstream fWriter(filename, ios::out);
    fWriter << inp.nF << " " << inp.nC << endl;
    fWriter << setprecision(15) << opt.zStar << endl;
    fWriter << opt.zStatus << endl;
    fWriter << opt.nOpen << endl;
    for (int i = 0; i < inp.nF; i++)
        if (opt.ySol[i] == 1)
            fWriter << " " << i;
    fWriter << endl;
    for (int i = 0; i < inp.nF; i++)
        for (int j = 0; j < inp.nC; j++)
            if (opt.xSol[i][j] > 0.0)
                fWriter << i << " " << j << " " << opt.xSol[i][j] << endl;


    /* for (int i = 0; i < inp.nF; i++)
     *     for (int t = 0; t < inp.nR; t++)
     *         if (cplex.getValue(psi_ilo[i][t] > 0.0)
     *             fWriter << "psi " << i << " " << t << " " << cplex.getValue(psi_ilo[i][t]) << endl; */

    fWriter.close();
    cout << "Solution written to disk. ('" << filename <<"')" << endl;

    // print here capacity consumptions
    for (int i = 0; i < inp.s[i]; i++)
    {
        if (opt.ySol[i] == 1)
        {
            float tot = 0.0;
            for (int j = 0; j < inp.nC; j++)
                if (opt.xSol[i][j] > 0.0)
                    tot += opt.xSol[i][j]*inp.d[j];
            cout << "Facility " << i << " with capacity = " << inp.s[i] << " and consumption = " << tot << endl;
        }
        
    }


    
}


/// Print solution to screen (and disk, if required)
void printSolution(char * _FILENAME, INSTANCE inp, SOLUTION opt, bool toDisk, 
int fullOutput)
{
    cout << endl << "** SOLUTION **" << endl;
    cout << " ..z*     \t= " << setprecision(15) << opt.zStar << endl;
    cout << " ..time   \t= " << opt.cpuTime << endl;
    cout << " ..status \t= " << opt.zStatus << endl;

    ofstream fWriter("solution.txt", ios::out);
    fWriter << _FILENAME << "\t" << instanceType << "\t"
            << versionType << "\t" << supportType << "\t" 
            << setprecision(15) << opt.zStar << "\t" 
            << opt.zStatus << "\t" << opt.cpuTime <<"\t"
            << _Omega << "\t" << _epsilon << "\t" << _delta << "\t" 
            << _gamma << "\t" << L << endl;
    fWriter.close();

    if (fullOutput >= 1)
    {
        cout << "Open Facilities: ";
        for (int i = 0; i < inp.nF; i++)
            if (opt.ySol[i] == 1)
                cout << setw(4) << i;
        cout << endl;
        if (fullOutput >= 2)
        {
            cout << "Allocation variables :: " << endl;
            for (int i = 0; i < inp.nF; i++)
                for (int j = 0; j < inp.nC; j++)
                    if (opt.xSol[i][j] >= EPSI)
                        cout << "x(" << i << "," << j << ") = " << setprecision(3) 
                        << opt.xSol[i][j] << endl;
        }
    }
}

/// Define the Multi-source Capacitated Facility Location Model [Nominal]
void define_MS_CFLP(INSTANCE inp, int fType, IloModel & model, IloCplex & cplex)
{
    IloEnv env = model.getEnv();

    char varName[100];

    // location variables 
    y_ilo = IloNumVarArray(env, inp.nF, 0, 1, ILOINT);

    // allocation variables
    x_ilo = TwoD(env, inp.nF);
    for (int i = 0; i < inp.nF; i++)
        x_ilo[i] = IloNumVarArray(env, inp.nC, 0.0, 1.0, ILOFLOAT);

    // set vars name
    for (int i = 0; i < inp.nF; i++)
    {
        /* sprintf(varName, "y.%d", (int)i);
         * y_ilo[i].setName(varName); */
        for (int j = 0; j < inp.nC; j++)
        {
            sprintf(varName, "x.%d.%d", (int)i, (int) j);
            x_ilo[i][j].setName(varName);

        }
    }

    // customers demand
    for (int j = 0; j < inp.nC; j++)
    {
        IloExpr sum(env);
        for (int i = 0; i < inp.nF; i++)
            sum += x_ilo[i][j];
        model.add(sum == 1.0);
    }

    // facility capacity 
    for (int i = 0; i < inp.nF; i++)
    {
        IloExpr sum(env);
        for (int j = 0; j < inp.nC; j++)
            sum += x_ilo[i][j]*inp.d[j];
        sum -= y_ilo[i]*inp.s[i];

        model.add(sum <= 0.0);
    }

    // thightening the model (does not seem to be beneficial)
    /* for (int i = 0; i < inp.nF; i++)
     *     for (int j = 0; j < inp.nC; j++)
     *         model.add(x_ilo[i][j] - y_ilo[i] <= 0.0); */

    // objective function
    IloExpr totCost(env);
    for (int i = 0; i < inp.nF; i++)
    {
        totCost += y_ilo[i]*inp.f[i];
        for (int j = 0; j < inp.nC; j++)
            totCost += inp.d[j]*x_ilo[i][j]*inp.c[i][j];
    }

    model.add(IloMinimize(env,totCost));
}

/// Define the Single-source Capacitated Facility Location Model [Nominal] 
void define_SS_CFLP(INSTANCE inp, int fType, IloModel & model, IloCplex & cplex)
{

    char varName[100];
    IloEnv env = model.getEnv();

    // location variables 
    y_ilo = IloNumVarArray(env, inp.nF, 0, 1, ILOINT);

    // allocation variables
    x_ilo = TwoD(env, inp.nF);
    // single source case
    for (int i = 0; i < inp.nF; i++)
        x_ilo[i] = IloNumVarArray(env, inp.nC, 0.0, 1.0, ILOINT);

    // set var names
    for (int i = 0; i < inp.nF; i++)
    {
        sprintf(varName, "y.%d", (int)i);
        y_ilo[i].setName(varName);
        for (int j = 0; j < inp.nC; j++)
        {
            sprintf(varName, "x.%d.%d", (int)i, (int) j);
            x_ilo[i][j].setName(varName);
        }
    }

    // customers demand
    for (int j = 0; j < inp.nC; j++)
    {
        IloExpr sum(env);
        for (int i = 0; i < inp.nF; i++)
            sum += x_ilo[i][j];
        model.add(sum == 1.0);
    }

    // facility capacity 
    for (int i = 0; i < inp.nF; i++)
    {
        IloExpr sum(env);
        for (int j = 0; j < inp.nC; j++)
            sum += x_ilo[i][j]*inp.d[j];
        sum -= y_ilo[i]*inp.s[i];

        model.add(sum <= 0.0);
    }

    // does it tighten the formulation?
    /* for (int i = 0; i < inp.nF; i++)
     *     for (int j = 0; j < inp.nC; j++)
     *         model.add(x_ilo[i][j] - y_ilo[i] <= 0.0); */

    // objective function
    IloExpr totCost(env);
    for (int i = 0; i < inp.nF; i++)
    {
        totCost += y_ilo[i]*inp.f[i];
        for (int j = 0; j < inp.nC; j++)
            totCost += x_ilo[i][j]*inp.c[i][j]*inp.d[j];
    }

    model.add(IloMinimize(env,totCost));
}

/// Define the Multi-source Capacitated Facility Location Model [Ellipsoidal]
void define_SOCP_CFLP(INSTANCE inp, int fType, IloModel & model, IloCplex & cplex)
{
    read_parameters_ellipsoidal();

    IloEnv env = model.getEnv();

    // location variables 
    y_ilo = IloNumVarArray(env, inp.nF, 0, 1, ILOINT);

    // allocation variables
    x_ilo = TwoD(env, inp.nF);
    for (int i = 0; i < inp.nF; i++)
        x_ilo[i] = IloNumVarArray(env, inp.nC, 0.0, 1.0, ILOFLOAT);

    // Q_i vars
    q_ilo = IloNumVarArray(env, inp.nC, 0.0, IloInfinity, ILOFLOAT);

    // W var
    w_ilo = IloNumVar(env, 0.0, IloInfinity, ILOFLOAT);

    // customers demand
    for (int j = 0; j < inp.nC; j++)
    {
        IloExpr sum(env);
        for (int i = 0; i < inp.nF; i++)
            sum += x_ilo[i][j];
        model.add(sum == 1.0);
    }

    // second order cone W 
    IloExpr sum(env);
    // sum = -w_ilo;
    sum = -w_ilo*w_ilo;
    for (int i = 0; i < inp.nF; i++)
        for (int j = 0; j < inp.nC; j++)
            sum += x_ilo[i][j]*x_ilo[i][j]*inp.c[i][j]*_epsilon*inp.c[i][j]*_epsilon;
    model.add(sum <= 0.0);

    // Q conic constraints
    for (int i = 0; i < inp.nF; i++)
    {
        IloExpr sum(env);
        // sum = -q_ilo[i];
        sum = -q_ilo[i]*q_ilo[i];
        for (int j = 0; j < inp.nC; j++)
            sum += x_ilo[i][j]*x_ilo[i][j]*_epsilon*_epsilon;
        model.add(sum <= 0.0);
    }

    // capacity constraints QUADRATIC
    for (int i = 0; i < inp.nF; i++)
    {
        IloExpr sum(env);
        for (int j = 0; j < inp.nC; j++)
            sum += x_ilo[i][j]*inp.d[j];

        sum += _Omega*q_ilo[i];
        sum -= y_ilo[i]*inp.s[i];

        model.add(sum <= 0.0);
    }

    // objective function
    IloExpr totCost(env);
    for (int i = 0; i < inp.nF; i++)
    {
        totCost += y_ilo[i]*inp.f[i];
        for (int j = 0; j < inp.nC; j++)
            totCost += x_ilo[i][j]*inp.c[i][j]*inp.d[j];
    }
    totCost += _Omega*w_ilo;

    model.add(IloMinimize(env,totCost));
}

/// Define robust model based on polyhedral uncertainty set
/**
 *
 * We implement here both the SINGLE SOURCE and the MULTI SOURCE versions. Just
 * change the way in which the allocation variables \c x_ilo are defined, to 
 * switch between the SS and the MS versions. The notation here follows the one
 * presented in the paper.
 *
 * Each support type is build from the nominal instance. That is, we read the
 * instance from disk, we read a few parameters (depending on the support type)
 * and we then build the support set around the nominal values provided by the
 * instance.
 *
 * The polyhedral uncertainty set is defined in the form:
 * \f[ W\mathbf{d} \leq \mathbf{h} \f]
 * Therefore, before we can solve the corresponding model, we need to define
 * the structure of \f$W\f$ and \f$\mathbf{h}\f$:
 *
 * * \f$W\f$ is of size \f$r \times n\f$, where \f$r\f$ is the number of budget
 *   constraints we want to define;
 * * \f$\mathbf{h}\f$ is of size \f$r \times 1\f$.
 *
 * The support set is defined in a different function, depending on the type of 
 * support we want to define (box, budget, etc.) The functions
 * `define_type_support()` is used to create \f$W\f$ and \f$\mathbf{h}\f$.
 * Each support type requires a set of parameters. These parameters are read
 * from disk and are stored in a file saved in the folder "parameters". The
 * file name is of type `paramsSupportType.txt`. Each file contains comments
 * w.r.t. how to interpret the parameters. However, the interpretation is as
 * follows:
 *
 * * File : `paramsBox.txt`
 *      - first row : `_epsilon`, i.e., the parameter used to define the with of
 *                    the box. Given a nominal demand value \f$d_j\f$, we define 
 *                    the interval around \f$d_j\f$ as:
 *                    \f[ (1-\epsilon)*d_j <= d_j <= (1+\epsilon)d_j, \quad 0
 *                    \leq \epsilon \leq 1 \f]
 *
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
                      \f[|B_l| = \lfloor \_gamma n \rfloor, \quad 0 \leq \_gamma \leq 1\f]
 *      - forth row : \f$L \geq 1\f$, i.e., number of budget constraints.
 *
 * To define which columns are included in each budget constraint, we first
 * compute the cardinality of each set \f$B_l\f$ (using `_gamma`); next, we randomly
 * select \f$\_gamma \times n\f$ columns from the set \f$N = \left\{0, ..., n-1\right\}\f$.
 *
 * **Note**: Sets \f$B_l\f$, with \f$l = 1, ..., L\f$, are NOT disjoint, i.e., the same
 * customer j can appear in more than one budget constraint.
 *
 * To test it, use the toy problem, whose instance values are defined in the
 * file 'data/toy.txt'. These values are the one reported in the paper, under the
 * example section. Changing the value of the parametes above, we can transform
 * a robust problem into a nomimal one (e.g., setting \f$\epsilon = 0\f$ in box
 * support.)
 */
void define_POLY_CFLP(INSTANCE inp, int fType, IloModel & model, IloCplex & cplex,
                      int support)
{

    // select support type
    // here we get W and h, depending on the type of support
    switch (support)
    {
        case 1 :
            define_box_support(inp);
            break;
        case 2 : 
            define_budget_support(inp, false);
            break;
        default :
            cout << "ERROR : Support type not defined.\n" << endl;
            exit(123);
    }

    char varName[100];
    IloEnv env = model.getEnv();

    // location variables
    y_ilo = IloNumVarArray(env, inp.nF, 0, 1, ILOINT);

    // allocation variables (Note: Change here to switch between MS and SS)
    x_ilo = TwoD(env, inp.nF);
    for (int i = 0; i < inp.nF; i++)
        x_ilo[i] = IloNumVarArray(env, inp.nC, 0.0, 1.0, ILOFLOAT); // MS
        // x_ilo[i] = IloNumVarArray(env, inp.nC, 0.0, 1.0, ILOINT); // SS 

    // set var names
    for (int i = 0; i < inp.nF; i++)
    {
        sprintf(varName, "y.%d", (int)i);
        y_ilo[i].setName(varName);
        for (int j = 0; j < inp.nC; j++)
        {
            sprintf(varName, "x.%d.%d", (int)i, (int) j);
            x_ilo[i][j].setName(varName);
        }
    }

    // psi variables
    psi_ilo = TwoD(env, inp.nF);
    for (int i = 0; i < inp.nF; i++)
        psi_ilo[i] = IloNumVarArray(env, inp.nR, 0.0, IloInfinity, ILOFLOAT);

    // u variables
    u_ilo = IloNumVarArray(env, inp.nR, 0.0, IloInfinity, ILOFLOAT);
    
    // set vars names
    for (int i = 0; i < inp.nF; i++)
        for (int t = 0; t < inp.nR; t++)
        {
            sprintf(varName, "psi.%d.%d", (int)i, (int) t);
            psi_ilo[i][t].setName(varName);
        }
    for (int t = 0; t < inp.nR; t++)
    {
        sprintf(varName, "u.%d", (int)t);
        u_ilo[t].setName(varName);
    }

    // delta variable (for the objective function)
    delta_ilo = IloNumVarArray(env, 1, 0.0, IloInfinity, ILOFLOAT);
    delta_ilo[0].setName("delta");

    // customers demand
    for (int j = 0; j < inp.nC; j++)
    {
        IloExpr sum(env);
        for (int i = 0; i < inp.nF; i++)
            sum += x_ilo[i][j];
        model.add(sum == 1.0);
    }

    // "robust" demand - constr. W*psi >= x
    for (int i = 0; i < inp.nF; i++)
        for (int j = 0; j < inp.nC; j++)
        {
            IloExpr sum(env);
            for (int l = inp.start[j]; l < inp.start[j+1]; l++)
            {
                int t = inp.index[l];
                sum += inp.W[l]*psi_ilo[i][t];
            }
            sum -= x_ilo[i][j];

            model.add(sum >= 0.0);
        }

    // robust capacity h*psi <= s*y
    for (int i = 0; i < inp.nF; i++)
    {
        IloExpr sum(env);
        for (int t = 0; t < inp.nR; t++)
            sum += inp.h[t]*psi_ilo[i][t];

        sum -= inp.s[i]*y_ilo[i];
        model.add(sum <= 0.0);
    }

    // "robust" objective function: h*u <= delta
    IloExpr sum(env);
    for (int t = 0; t < inp.nR; t++)
        sum += inp.h[t]*u_ilo[t];

    sum -= delta_ilo[0];
    model.add(sum <= 0.0);


    // second robust obj function: W*u >= c*x
    for (int j = 0; j < inp.nC; j++)
    {
        IloExpr sum(env);
        for (int l = inp.start[j]; l < inp.start[j+1]; l++)
        {
            int t = inp.index[l];
            sum += inp.W[l]*u_ilo[t];
        }
        for (int i = 0; i < inp.nF; i++)
            sum -= inp.c[i][j]*x_ilo[i][j];

        model.add(sum >= 0.0);
    }

    // objective function: min f*y + delta
    IloExpr totCost(env);
    for (int i = 0; i < inp.nF; i++)
        totCost += y_ilo[i]*inp.f[i];
    totCost += delta_ilo[0];

    model.add(IloMinimize(env,totCost));
}



void define_new_POLY_CFLP(INSTANCE inp, int fType, IloModel & model, IloCplex & cplex,
                      int support)
{

    // select support type
    // here we get W and h, depending on the type of support
    switch (support)
    {
        case 1 :
            define_box_support(inp);
            break;
        case 2 : 
            // define_budget_support(inp, false);
            define_new_budget_support(inp, false);
            break;
        default :
            cout << "ERROR : Support type not defined.\n" << endl;
            exit(123);
    }

    char varName[100];
    char constrName[100];
    IloEnv env = model.getEnv();

    // location variables
    y_ilo = IloNumVarArray(env, inp.nF, 0, 1, ILOINT);

    // allocation variables (Note: Change here to switch between MS and SS)
    x_ilo = TwoD(env, inp.nF);
    for (int i = 0; i < inp.nF; i++)
        x_ilo[i] = IloNumVarArray(env, inp.nC, 0.0, 1.0, ILOFLOAT); // MS
        // x_ilo[i] = IloNumVarArray(env, inp.nC, 0.0, 1.0, ILOINT); // SS 

    // set var names
    for (int i = 0; i < inp.nF; i++)
    {
        sprintf(varName, "y.%d", (int)i);
        y_ilo[i].setName(varName);
        for (int j = 0; j < inp.nC; j++)
        {
            sprintf(varName, "x.%d.%d", (int)i, (int) j);
            x_ilo[i][j].setName(varName);
        }
    }

    // reconstruct the mapping
    int * setB = new int[inp.nC];
    for (int j = 0; j < inp.nC; j++)
        setB[j] = 0;
    for (int j = 0; j < inp.listB.size(); j++)
        setB[inp.listB[j]] = 1;

    cout << "Mapping (list of vars in budget constraints :: ";
    for (int j =0; j < inp.nC; j++)
        cout << " " << setB[j];
    cout << endl;
    cout << "Nr. of rows of matrix W = " << inp.nR << endl;

    // psi variables
    psi_ilo = TwoD(env, inp.nF);
    for (int i = 0; i < inp.nF; i++)
        psi_ilo[i] = IloNumVarArray(env, inp.nR, 0.0, IloInfinity, ILOFLOAT);

    // u variables
    u_ilo = IloNumVarArray(env, inp.nR, 0.0, IloInfinity, ILOFLOAT);
    
    // set vars names
    for (int i = 0; i < inp.nF; i++)
        for (int t = 0; t < inp.nR; t++)
        {
            sprintf(varName, "psi.%d.%d", (int)i, (int) t);
            psi_ilo[i][t].setName(varName);
        }
    for (int t = 0; t < inp.nR; t++)
    {
        sprintf(varName, "u.%d", (int)t);
        u_ilo[t].setName(varName);
    }

    // delta variable (for the objective function)
    delta_ilo = IloNumVarArray(env, 1, 0.0, IloInfinity, ILOFLOAT);
    delta_ilo[0].setName("delta");

    // customers demand
    for (int j = 0; j < inp.nC; j++)
    {
        IloExpr sum(env);
        for (int i = 0; i < inp.nF; i++)
            sum += x_ilo[i][j];
        sprintf(constrName, "demand.%d", (int)j);
        model.add(IloRange(env, 1.0, sum, 1.0, constrName));
    }

    // "robust" demand - constr. W*psi >= x
    int t;
    for (int i = 0; i < inp.nF; i++)
        for (int j = 0; j < inp.listB.size(); j++)
        {
            IloExpr sum(env);
            for (int l = inp.start[j]; l < inp.start[j+1]; l++)
            {
                t = inp.index[l];
                sum += inp.W[l]*psi_ilo[i][t];
            }
            sum -= x_ilo[i][inp.listB[j]];

            sprintf(constrName, "r.demand.%d.%d", (int)i, (int)inp.listB[j]);
            model.add(IloRange(env, 0.0, sum, IloInfinity, constrName));
        }


    // robust capacity h*psi <= s*y
    for (int i = 0; i < inp.nF; i++)
    {
        IloExpr sum(env);
        for (int t = 0; t < inp.nR; t++)
            sum += inp.h[t]*psi_ilo[i][t];

        for (int j = 0; j < inp.nC; j++)
            if (setB[j] > 0)
                sum += inp.d[j]*(1.0-_epsilon)*x_ilo[i][j];
            else
                sum += inp.d[j]*(1.0+_epsilon)*x_ilo[i][j];

        sum -= inp.s[i]*y_ilo[i];

        sprintf(constrName, "r.capacity.%d", (int)i);
        model.add(IloRange(env, -IloInfinity, sum, 0.0, constrName));
    }

    // "robust" objective function: h*u <= delta
    IloExpr sum(env);
    for (int t = 0; t < inp.nR; t++)
        sum += inp.h[t]*u_ilo[t];

    sum -= delta_ilo[0];
    sprintf(constrName, "r.delta");
    model.add(IloRange(env, -IloInfinity, sum, 0.0, constrName));


    // second robust obj function: W*u >= c*x
    for (int j = 0; j < inp.listB.size(); j++)
    {
        int el = inp.listB[j];  // element in budget constraints
        IloExpr sum(env);
        for (int l = inp.start[j]; l < inp.start[j+1]; l++)
        {
            t    = inp.index[l];
            sum += inp.W[l]*u_ilo[t];
        }
        for (int i = 0; i < inp.nF; i++)
            sum -= inp.c[i][el]*x_ilo[i][el];

        sprintf(constrName, "r.objective.%d", (int)el);
        model.add(IloRange(env, 0.0, sum, IloInfinity, constrName));
    }
    
    // thightening the model (does not seem to be beneficial)
    for (int i = 0; i < inp.nF; i++)
        for (int j = 0; j < inp.nC; j++)
            model.add(x_ilo[i][j] - y_ilo[i] <= 0.0);


    // objective function: min f*y + delta
    IloExpr totCost(env);
    for (int i = 0; i < inp.nF; i++)
        totCost += y_ilo[i]*inp.f[i];
    totCost += delta_ilo[0];


    // new terms (deterministic part of d)
    for (int j = 0; j < inp.nC; j++)
        {
            IloExpr coeff(env); 
            for (int i = 0; i < inp.nF; i++)
                coeff += inp.c[i][j]*x_ilo[i][j];

            if (setB[j] == 0)
                coeff *= inp.d[j]*(1+_epsilon);
            else
                coeff *= inp.d[j]*(1-_epsilon);

            totCost += coeff;
        }

    model.add(IloMinimize(env,totCost));
}

void define_new2_POLY_CFLP(INSTANCE inp, int fType, IloModel & model, IloCplex & cplex,
                      int support)
{

    // select support type
    // here we get W and h, depending on the type of support
    switch (support)
    {
        case 1 :
            define_box_support(inp);
            break;
        case 2 : 
            // define_budget_support(inp, false);
            define_new2_budget_support(inp, false);
            break;
        default :
            cout << "ERROR : Support type not defined.\n" << endl;
            exit(123);
    }
    cout << "EPSILON HERE IS " << _epsilon << endl;

    char varName[100];
    char constrName[100];
    IloEnv env = model.getEnv();

    // location variables
    y_ilo = IloNumVarArray(env, inp.nF, 0, 1, ILOINT);

    // allocation variables (Note: Change here to switch between MS and SS)
    x_ilo = TwoD(env, inp.nF);
    for (int i = 0; i < inp.nF; i++)
        x_ilo[i] = IloNumVarArray(env, inp.nC, 0.0, 1.0, ILOFLOAT); // MS
        // x_ilo[i] = IloNumVarArray(env, inp.nC, 0.0, 1.0, ILOINT); // SS 

    // set var names
    for (int i = 0; i < inp.nF; i++)
    {
        sprintf(varName, "y.%d", (int)i);
        y_ilo[i].setName(varName);
        for (int j = 0; j < inp.nC; j++)
        {
            sprintf(varName, "x.%d.%d", (int)i, (int) j);
            x_ilo[i][j].setName(varName);
        }
    }

    // reconstruct the mapping
    int * setB = new int[inp.nC];
    for (int j = 0; j < inp.nC; j++)
        setB[j] = 0;
    for (int j = 0; j < inp.listB.size(); j++)
        setB[inp.listB[j]] = 1;

    // psi variables
    psi_ilo = TwoD(env, inp.nF);
    for (int i = 0; i < inp.nF; i++)
        psi_ilo[i] = IloNumVarArray(env, inp.nR, 0.0, IloInfinity, ILOFLOAT);

    // u variables
    u_ilo = IloNumVarArray(env, inp.nR, 0.0, IloInfinity, ILOFLOAT);
    
    // set vars names
    for (int i = 0; i < inp.nF; i++)
        for (int t = 0; t < inp.nR; t++)
        {
            sprintf(varName, "psi.%d.%d", (int)i, (int) t);
            psi_ilo[i][t].setName(varName);
        }
    for (int t = 0; t < inp.nR; t++)
    {
        sprintf(varName, "u.%d", (int)t);
        u_ilo[t].setName(varName);
    }

    // delta variable (for the objective function)
    delta_ilo = IloNumVarArray(env, 1, 0.0, IloInfinity, ILOFLOAT);
    delta_ilo[0].setName("delta");

    // customers demand
    for (int j = 0; j < inp.nC; j++)
    {
        IloExpr sum(env);
        for (int i = 0; i < inp.nF; i++)
            sum += x_ilo[i][j];
        sprintf(constrName, "demand.%d", (int)j);
        model.add(IloRange(env, 1.0, sum, 1.0, constrName));
    }

    // "robust" demand - constr. W*psi >= x
    int t;
    for (int i = 0; i < inp.nF; i++)
        for (int j = 0; j < inp.nC; j++)
        {
            IloExpr sum(env);
            for (int l = inp.start[j]; l < inp.start[j+1]; l++)
            {
                t    = inp.index[l];
                sum += inp.W[l]*psi_ilo[i][t];
            }
            sum -= x_ilo[i][j];

            sprintf(constrName, "r.demand.%d.%d", (int)i, (int)j);
            model.add(IloRange(env, 0.0, sum, IloInfinity, constrName));
        }


    // robust capacity h*psi <= s*y
    for (int i = 0; i < inp.nF; i++)
    {
        IloExpr sum(env);
        for (int t = 0; t < inp.nR; t++)
            sum += inp.h[t]*psi_ilo[i][t];

        for (int j = 0; j < inp.nC; j++)
            // if (setB[j] >= 0)
                sum += inp.d[j]*(1.0-_epsilon)*x_ilo[i][j];

        sum -= inp.s[i]*y_ilo[i];

        sprintf(constrName, "r.capacity.%d", (int)i);
        model.add(IloRange(env, -IloInfinity, sum, 0.0, constrName));
    }

    // "robust" objective function: h*u <= delta
    IloExpr sum(env);
    for (int j = 0; j < inp.nC; j++)
        if (setB[j] == 1) // skip terms not in budget constraints
            sum += inp.h[j]*u_ilo[j];
    for (int t = inp.nC; t < inp.nR; t++)
            sum += inp.h[t]*u_ilo[t];

    sum -= delta_ilo[0];
    sprintf(constrName, "r.delta");
    model.add(IloRange(env, -IloInfinity, sum, 0.0, constrName));


    // second robust obj function: W*u >= c*x
    for (int j = 0; j < inp.nC; j++)
    {
        IloExpr sum(env);
        for (int l = inp.start[j]; l < inp.start[j+1]; l++)
        {
            int t = inp.index[l];
            sum += inp.W[l]*u_ilo[t];
        }
        for (int i = 0; i < inp.nF; i++)
            sum -= inp.c[i][j]*x_ilo[i][j];

        sprintf(constrName, "r.objective.%d", (int)j);
        model.add(IloRange(env, 0.0, sum, IloInfinity, constrName));
    }

    // objective function: min f*y + delta
    IloExpr totCost(env);
    for (int i = 0; i < inp.nF; i++)
        totCost += y_ilo[i]*inp.f[i];
    totCost += delta_ilo[0];
    
    // new terms (deterministic part of d)
    for (int j = 0; j < inp.nC; j++)
        {
            IloExpr coeff(env); 
            for (int i = 0; i < inp.nF; i++)
                coeff += inp.c[i][j]*x_ilo[i][j];

            if (setB[j] == 0)
                coeff *= inp.d[j]*(1+_epsilon);
            else
                coeff *= inp.d[j]*(1-_epsilon);

            totCost += coeff;
        }

    model.add(IloMinimize(env,totCost));
}


/// Set cplex parameters and solve the optimization problem
int solveCplexProblem(IloModel model, IloCplex cplex, INSTANCE inp, int solLimit, int timeLimit, int displayLimit)
{
    try
    {
        IloEnv env = model.getEnv();
        cplex.setOut(env.getNullStream());
        cplex.setParam(IloCplex::ClockType, 2); // 1 --> Cpu Time; 2 --> Wall-clock 
        cplex.setParam(IloCplex::ClockType, 2); // 1 --> Cpu Time; 2 --> Wall clock
        cplex.setParam(IloCplex::MIPInterval, 5000);
        cplex.setParam(IloCplex::MIPDisplay, displayLimit);
        cplex.setParam(IloCplex::IntSolLim, solLimit);
        cplex.setParam(IloCplex::TiLim, timeLimit);

        cplex.exportModel("cflp.lp");

        if (!cplex.solve())
        {
            env.error() << "Failed to Optimize MIP " << endl;
            throw(-1);
        }
        return 1;
    }
    catch (...)
    {
        cout << "Exception caught ... " << endl;
        return -1;
    }
}


/// Define \f$W\f$ and \f$\mathbf{h}\f$, given the value of `_epsilon`.
/**
 * We need to define \f$W\mathbf{d} \leq \mathbf{h}\f$, where:
 * * \f$W\f$ is of size [r x n]
 * * \f$\mathbf{h}\f$ is of size [r x 1]
 * where:
 * * \f$r = 2n\f$ (two elements for each column \f$j\f$)
 * * \f$n =\f$ number of customers
 * 
 * Since matrix \f$W\f$ is sparse, we use a column major format as follows. 
 * Consider the toy problem, for which matrix \f$W\f$ and vector \f$\mathbf{h}\f$
 * are: 
    \f[
    W = \left[\begin{array}{rrr}
    -1 & 0 & 0 \\
    0  & -1 & 0 \\
    0  & 0 & -1 \\
    1 & 0 & 0 \\
    0 & 1 & 0 \\
    0 & 0 & 1 \\
    \end{array}\right] \quad
    \textbf{d} = \left[\begin{array}{ccc}
    d_1 \\
    d_2 \\
    d_3
    \end{array}\right] \quad
    \textbf{h} =\left[
    \begin{array}{r}
    -d_0(1-\epsilon) \\ -d_1(1-\epsilon) \\ -d_2(1-\epsilon) \\ d_0(1+\epsilon) \\ d_1(1+\epsilon) \\ d_2(1+\epsilon) 
    \end{array}\right]
    \f]
 *
 * Vector h is kept as is it, while matrix W uses a column major format:
 *
 * > W     = [-1 1  -1 1  -1 1]
 *
 * > start = [0 2 4 6]
 *
 * > index = [0 3  1 4  2 5]
 * 
 * where elements of column \c j in \c W and \c index are found in positions 
 * going from \c start[j] to \c start[j+1] (note that vector \c start contains 
 * a final extra element to close the cycle). Thus, e.g.,  the elements of the 
 * second column of `W (j = 1)` are obtained as:
 * \code{.cpp}
 * for (int l = start[j]; l < start[j+1]; l++)
 *      W[l] is the element in row index[l] of the matrix
 * \endcode
 */
void define_box_support(INSTANCE & inp)
{
   read_parameters_box();
   inp.nR = 2*inp.nC;
   inp.h  = new double[inp.nR];
   for (int j = 0; j < inp.nC; j++)
   {
       inp.h[j]        = -inp.d[j]*(1.0-_epsilon);
       inp.h[j+inp.nC] =  inp.d[j]*(1.0+_epsilon);
       cout << "left and write " << inp.h[j] << " and " << inp.h[j+inp.nC] << endl;
   }

   // define matrix W in column major format
   inp.W     = new int[inp.nC*2];
   inp.index = new int[inp.nC*2];
   inp.start = new int[inp.nC+1]; // one extra element in last position

    int pos = 0;
    for (int j = 0; j < inp.nC; j++)
    {
        inp.start[j]   = pos;
        inp.W[pos]     = -1;
        inp.index[pos] = j;
        pos++;
        inp.W[pos]     = 1;
        inp.index[pos] = j+inp.nC;
        pos++;
    }
    inp.start[inp.nC] = pos;
    assert(pos == inp.nC*2);
}


/// Read parameters to define the Budget support.
/**
 * We allow for two options here:
 * * __case 1__: Create the budget support from the nominal values of the instance.
 * In this case, we read the parameters from the disk file (see below) and we
 * generate the polyhedron. This is a stochastic process, since the set of
 * columns included in each budget constrain is randomly generated. Therefore,
 * to ensure reproducibility, we save these values in a disk file, within the
 * folder "support."
 * 
 * * __case 2__: Read the budget support from the disk. This allows to reproduce the
 * results of an instance, by recreating the same budget support set.
 *
 * The file `paramsBudget.txt` has the following format:
 *      - first row : `_epsilon`, as above
 *      - second row: `_delta`, i.e., how the \f$b_l\f$ value (the r.h.s. value of each
 *                    budget constraint) is defined. We define \f$b_l\f$ as a
 *                    percentage of the total demand of customers included in
 *                    the budget constraint, i.e.:
 *                    \f[ b_l = \delta* (\sum_{j \in B_l} d_j), \quad 0 \leq
 *                    \delta \leq 1 \f]
 *      - third row : `_gamma`, i.e., percentage of columns included in each
 *                    budget constraint:
                      \f[|B_l| = \lfloor \_gamma n \rfloor, \quad 0 \leq \_gamma \leq 1\f]
 *      - forth row : \f$L \geq 1\f$, i.e., number of budget constraints.
 *
 *  The structure of \f$W\f$ and \f$\mathbf{h}\f$ and the column-major format
 *  used is simular to the one employed in define_box_support(). The first
 *  \f$2n\f$ rows of \f$W\f$ are identical to the box support. Next, we define
 *  \f$L\f$ extra rows, to include the budget constraints. Consider budget
 *  constraint \f$l\f$, associated to set \f$B_l\f$. This row is defined as
 *  follows:
 *  \f[
 *  w_{lj} =  \left\{
 *  \begin{array}{ll}
 *  1, & j \in B_l \\
 *  0, & \mbox{otherwise}
 *  \end{array}
 *  \right.
 *  \f]
 *  For example, if we want to add the following budget constraint (following
 *  the example from the paper):
 *  \f[
 *  d_1 + d_2 + d_3 \leq 4
 *  \f]
 *  where we have \f$L=1\f$, \f$B_l =\left\{1,2,3\right\}\f$, we obtain:
    \f[
    W = \left[\begin{array}{rrr}
    -1 & 0 & 0 \\
    0  & -1 & 0 \\
    0  & 0 & -1 \\
    1 & 0 & 0 \\
    0 & 1 & 0 \\
    0 & 0 & 1 \\
    1 & 1 & 1 \\
    \end{array}\right] \quad
    \textbf{d} = \left[\begin{array}{ccc}
    d_1 \\
    d_2 \\
    d_3
    \end{array}\right] \quad
    \textbf{h} =\left[
    \begin{array}{r}
    -d_0(1-\epsilon) \\ -d_1(1-\epsilon) \\ -d_2(1-\epsilon) \\ d_0(1+\epsilon) \\ d_1(1+\epsilon) \\ d_2(1+\epsilon) \\ 4 \\
    \end{array}\right]
    \f]
 *
 *  Consequently, in terms of column-major format, we get:

 * > W     = [-1 1 1 -1 1 1  -1 1 1]
 *
 * > start = [0 3 6 9]
 *
 * > index = [0 3 6 1 4 6  2 5 6]
 *
 */
void define_budget_support(INSTANCE & inp, bool fromDisk)
{

    int    nBl     = 0;
    int  **Bl;
    double *budget;

    read_parameters_budget();

    nBl    = floor(_gamma*(double)inp.nC); // cardinality of each B_l
    Bl     = new int*[L];
    budget = new double[L];
    for (int l = 0; l < L; l++)
        Bl[l] = new int[nBl];
    
    // cardinality of sets B_l
    cout << "[** |B_l| = " << nBl << "]\n" << endl;
    // nBl = 2;
    // cout << "HARD CODED HERE !!! " << endl;

    // randomly generate sets B_l and compute budget b_l
    std::vector<int> shuffled;
    for (int i = 0; i < inp.nC; ++i) 
        shuffled.push_back(i); // 0 2 3 ... nC-1

    if (readFromDisk==false)
    {
        // initialize random generator (uniform distribution)
        uniform_int_distribution<> d(0,inp.nC-1);
        for (int l = 0; l < L; l++)
        {
            // using built-in random generator:
            std::random_shuffle(shuffled.begin(), shuffled.end());
            budget[l] = 0.0;
            for (int k = 0; k < nBl; k++)
            {
                // int el = d(gen);
                int el = shuffled[k];
                Bl[l][k] = el;
                budget[l] += inp.d[el];
            }
            cout << endl;
        }
        // adjust b_l values
        for (int l = 0; l < L; l++)
            budget[l] = floor(_delta*budget[l]);

        // do we want to save the sets B_l (and the parameters?)
        bool save2Disk = true;
        if (save2Disk==true)
            save_instance_2_disk(_epsilon, _delta, _gamma, L, nBl, Bl, budget);
    }
    else
    {
        string  s1      = string(_FILENAME);
        s1              = s1.substr(s1.find_last_of("\\/"), 100);
        string filename = "support" + s1 + ".budget"; 

        ifstream fReader(filename, ios::in);
        if (!fReader)
        {
            cout << "Cannot open file '" << filename << "'." << endl;
            exit(111);
        }
        for (int l = 0; l < L; l++)
        {
            int temp = nBl;
            fReader >> nBl;
            cout << "t vs nBl " << temp << " " << nBl << endl;
            assert(temp == nBl);
            Bl[l] = new int[nBl];

            for (int k = 0; k < nBl; k++)
                fReader >> Bl[l][k];
            fReader >> budget[l];
        }

        fReader.close();
        cout << "[** Instance read from disk. File '" << filename << "']" << endl;
        cout << "[** Uncertainty Set Parameters :: _epsilon = " << _epsilon 
             << "; _delta = " << _delta << "; _gamma = " << _gamma << "; L = " << L 
             << "]\n" << endl;
    }

    for (int l = 0; l < L; l++)
    {
        cout <<"Budget constraint # " << l << ":: ";
        for (int k = 0; k < nBl; k++)
            cout << " " << Bl[l][k];
        cout << " :: Total Budget = " << budget[l] << endl;
    }

    // define mapping: list of budget constraints including column j 
    MyVect * mapping = new MyVect[inp.nC];
    for (int l = 0; l < L; l++)
        for (int k = 0; k < nBl; k++)
            mapping[Bl[l][k]].push_back(l);

   // total number of rows of W and h
   inp.nR = 2*inp.nC + L;
    // define vector h
    // a. box part
   inp.h  = new double[inp.nR];
   for (int j = 0; j < inp.nC; j++)
   {
       inp.h[j]        = -inp.d[j]*(1.0-_epsilon);
       // inp.h[j]        = 0.0;
       inp.h[j+inp.nC] =  inp.d[j]*(1.0+_epsilon);
   }


    
    // b. budget part
    for (int l = 0; l < L; l++)
        inp.h[2*inp.nC+l] = budget[l];

   // define matrix W in column major format
   int nEls  = 2*inp.nC + L*nBl;
   inp.W     = new int[nEls];
   inp.index = new int[nEls];
   inp.start = new int[inp.nC+1]; // one extra element in last position

    int pos = 0;
    for (int j = 0; j < inp.nC; j++)
    {
        inp.start[j]     = pos;
        inp.W[pos]       = -1;
        inp.index[pos++] = j;
        inp.W[pos]       = 1;
        inp.index[pos++] = j+inp.nC;
        if (mapping[j].size() > 0)
        {
            for (MyVect::iterator it = mapping[j].begin(); it != mapping[j].end(); it++)
            {
                inp.W[pos]       = 1;
                inp.index[pos++] = 2*inp.nC + *it;
            }
        }
    }
    inp.start[inp.nC] = pos;
    // NOTE : Remove ASSERT from final version
    assert(pos == nEls);

}


void define_new_budget_support(INSTANCE & inp, bool fromDisk)
{

    int    nBl     = 0;
    int  **Bl;
    double *budget;

    read_parameters_budget();

    nBl    = floor(_gamma*(double)inp.nC); // cardinality of each B_l
    Bl     = new int*[L];
    budget = new double[L];
    for (int l = 0; l < L; l++)
        Bl[l] = new int[nBl];
    
    // cardinality of sets B_l
    cout << "[** |B_l| = " << nBl << "]\n" << endl;
    /* nBl = 2;
     * cout << "HARD CODED HERE !!!!!!!!!" << endl; */

    // randomly generate sets B_l and compute budget b_l
    std::vector<int> shuffled;
    for (int i = 0; i < inp.nC; ++i) 
        shuffled.push_back(i); // 0 2 3 ... nC-1

    if (readFromDisk==false)
    {
        // initialize random generator (uniform distribution)
        uniform_int_distribution<> d(0,inp.nC-1);
        for (int l = 0; l < L; l++)
        {
            // using built-in random generator:
            std::random_shuffle(shuffled.begin(), shuffled.end());
            budget[l] = 0.0;
            for (int k = 0; k < nBl; k++)
            {
                // int el = d(gen);
                int el = shuffled[k];
                Bl[l][k] = el;
                budget[l] += inp.d[el];
            }
            cout << endl;
        }
        // adjust b_l values
        for (int l = 0; l < L; l++)
            budget[l] = floor(_delta*budget[l]);

        // do we want to save the sets B_l (and the parameters?)
        bool save2Disk = true;
        if (save2Disk==true)
            save_instance_2_disk(_epsilon, _delta, _gamma, L, nBl, Bl, budget);
    }
    else
    {
        string  s1      = string(_FILENAME);
        s1              = s1.substr(s1.find_last_of("\\/"), 100);
        string filename = "support" + s1 + ".budget"; 

        ifstream fReader(filename, ios::in);
        if (!fReader)
        {
            cout << "Cannot open file '" << filename << "'." << endl;
            exit(111);
        }
        for (int l = 0; l < L; l++)
        {
            int temp = nBl;
            fReader >> nBl;
            cout << "temp vs nBl " << temp << " " << nBl << endl;
            assert(temp == nBl);
            Bl[l] = new int[nBl];

            for (int k = 0; k < nBl; k++)
                fReader >> Bl[l][k];
            fReader >> budget[l];
        }

        fReader.close();
        cout << "[** Instance read from disk. File '" << filename << "']" << endl;
        cout << "[** Uncertainty Set Parameters :: _epsilon = " << _epsilon 
             << "; _delta = " << _delta << "; _gamma = " << _gamma << "; L = " << L 
             << "]\n" << endl;
    }

    for (int l = 0; l < L; l++)
    {
        cout <<"Budget constraint # " << l << ":: ";
        for (int k = 0; k < nBl; k++)
            cout << " " << Bl[l][k];
        cout << " :: Total Budget = " << budget[l] << endl;
    }
    cout << "Corrected Budget :: " << endl;
    for (int l = 0; l < L; l++)
    {
        for (int k = 0; k < nBl; k++)
            budget[l] -= inp.d[Bl[l][k]]*(1.0-_epsilon);
        cout << " :: Total Corrected Budget = " << budget[l] << endl;
    }


    // define mapping: list of budget constraints including column j 
    MyVect * mapping = new MyVect[inp.nC];
    for (int l = 0; l < L; l++)
        for (int k = 0; k < nBl; k++)
            mapping[Bl[l][k]].push_back(l);

    // count cardinality of budget constraints overall: how many vars are used
    // in all the budget constraints
    std::vector <int> listB;
    for (int j = 0; j < inp.nC; j++)
        if (mapping[j].size() > 0)
            inp.listB.push_back(j);

    cout << inp.listB.size() << " elements in constraints :: ";
    for (int j = 0; j < inp.listB.size(); j++)
        cout << " " << inp.listB[j];
    cout << endl;

   // total number of rows of W and h
   int cardB = inp.listB.size();
   inp.nR = cardB + L;
    
    // define vector h
    // a. box part
   inp.h  = new double[inp.nR];
   for (int j = 0; j < cardB; j++)
       inp.h[j] =  inp.d[inp.listB[j]]*(1.0+_epsilon);
    
    // b. budget part
    for (int l = 0; l < L; l++)
        inp.h[cardB+l] = budget[l];

   cout << "H VECTOR " << endl;
   for (int j = 0; j < cardB+L; j++)
       cout << "h[" << j << "] = " << inp.h[j] << endl;

   // define matrix W in column major format
   int nEls  = cardB + L*nBl;
   inp.W     = new int[nEls];
   inp.index = new int[nEls];
   inp.start = new int[cardB+1]; // one extra element in last position

    int pos = 0;
    for (int j = 0; j < cardB; j++)
    {
        inp.start[j]     = pos;
        inp.W[pos]       = 1;
        inp.index[pos++] = j;
        int el           = inp.listB[j]; // actual column index
        if (mapping[el].size() > 0)
        {
            for (MyVect::iterator it = mapping[el].begin(); it != mapping[el].end(); it++)
            {
                inp.W[pos]       = 1;
                inp.index[pos++] = cardB + *it;
            }
        }
    }
    inp.start[cardB] = pos;
    // NOTE : Remove ASSERT from final version
    assert(pos == nEls);

    cout << "Column Major Format :: " << endl;
    for (int j = 0; j < inp.listB.size(); j++)
    {
        cout << "  var " << j << " from " << inp.start[j] << " to " << inp.start[j+1] << " :: " << endl;
        for (int k = inp.start[j]; k < inp.start[j+1]; k++)
            cout << "... el " << inp.W[k] << " in row " << inp.index[k] << endl;
    }



}

void define_new2_budget_support(INSTANCE & inp, bool fromDisk)
{

    int    nBl     = 0;
    int  **Bl;
    double *budget;

    read_parameters_budget();

    nBl    = floor(_gamma*(double)inp.nC); // cardinality of each B_l
    Bl     = new int*[L];
    budget = new double[L];
    for (int l = 0; l < L; l++)
        Bl[l] = new int[nBl];
    
    // cardinality of sets B_l
    cout << "[** |B_l| = " << nBl << "]\n" << endl;
    // nBl = 2;
    // cout << "HARD CODED HERE !!!!!!!!!" << endl;

    // randomly generate sets B_l and compute budget b_l
    std::vector<int> shuffled;
    for (int i = 0; i < inp.nC; ++i) 
        shuffled.push_back(i); // 0 2 3 ... nC-1

    if (readFromDisk==false)
    {
        // initialize random generator (uniform distribution)
        uniform_int_distribution<> d(0,inp.nC-1);
        for (int l = 0; l < L; l++)
        {
            // using built-in random generator:
            std::random_shuffle(shuffled.begin(), shuffled.end());
            budget[l] = 0.0;
            for (int k = 0; k < nBl; k++)
            {
                // int el = d(gen);
                int el = shuffled[k];
                Bl[l][k] = el;
                budget[l] += inp.d[el];
            }
            cout << endl;
        }
        // adjust b_l values
        for (int l = 0; l < L; l++)
            budget[l] = floor(_delta*budget[l]);

        // do we want to save the sets B_l (and the parameters?)
        bool save2Disk = true;
        if (save2Disk==true)
            save_instance_2_disk(_epsilon, _delta, _gamma, L, nBl, Bl, budget);
    }
    else
    {
        string  s1      = string(_FILENAME);
        s1              = s1.substr(s1.find_last_of("\\/"), 100);
        string filename = "support" + s1 + ".budget"; 

        ifstream fReader(filename, ios::in);
        if (!fReader)
        {
            cout << "Cannot open file '" << filename << "'." << endl;
            exit(111);
        }
        for (int l = 0; l < L; l++)
        {
            int temp = nBl;
            fReader >> nBl;
            cout << "temp vs nBl " << temp << " " << nBl << endl;
            assert(temp == nBl);
            Bl[l] = new int[nBl];

            for (int k = 0; k < nBl; k++)
                fReader >> Bl[l][k];
            fReader >> budget[l];
        }

        fReader.close();
        cout << "[** Instance read from disk. File '" << filename << "']" << endl;
        cout << "[** Uncertainty Set Parameters :: _epsilon = " << _epsilon 
             << "; _delta = " << _delta << "; _gamma = " << _gamma << "; L = " << L 
             << "]\n" << endl;
    }

    for (int l = 0; l < L; l++)
    {
        cout <<"Budget constraint # " << l << ":: ";
        for (int k = 0; k < nBl; k++)
            cout << " " << Bl[l][k];
        cout << " :: Total Budget = " << budget[l] << endl;
    }
    cout << "Corrected Budget :: " << endl;
    for (int l = 0; l < L; l++)
    {
        for (int k = 0; k < nBl; k++)
            budget[l] -= inp.d[Bl[l][k]]*(1.0-_epsilon);
        cout << " :: Total Corrected Budget = " << budget[l] << endl;
    }


    // define mapping: list of budget constraints including column j 
    MyVect * mapping = new MyVect[inp.nC];
    for (int l = 0; l < L; l++)
        for (int k = 0; k < nBl; k++)
            mapping[Bl[l][k]].push_back(l);

    // count cardinality of budget constraints overall: how many vars are used
    // in all the budget constraints
    std::vector <int> listB;
    for (int j = 0; j < inp.nC; j++)
        if (mapping[j].size() > 0)
            inp.listB.push_back(j);

   // total number of rows of W and h
   inp.nR = inp.nC + L;
    
    // define vector h
    // a. box part
   inp.h  = new double[inp.nR];
   for (int j = 0; j < inp.nC; j++)
       inp.h[j] =  inp.d[j]*(2.0*_epsilon);
       // inp.h[j] =  inp.d[j]*(1.0+_epsilon);
    
    // b. budget part
    for (int l = 0; l < L; l++)
        inp.h[inp.nC+l] = budget[l];

   /* cout << "H VECTOR " << endl;
    * for (int j = 0; j < inp.nC+L; j++)
    *     cout << "h[" << j << "] = " << inp.h[j] << endl; */

   // define matrix W in column major format
   int nEls  = inp.nC + L*nBl;
   inp.W     = new int[nEls];
   inp.index = new int[nEls];
   inp.start = new int[inp.nC+1]; // one extra element in last position

    int pos = 0;
    for (int j = 0; j <inp.nC; j++)
    {
        inp.start[j]     = pos;
        inp.W[pos]       = 1;
        inp.index[pos++] = j;
        if (mapping[j].size() > 0)
        {
            for (MyVect::iterator it = mapping[j].begin(); it != mapping[j].end(); it++)
            {
                inp.W[pos]       = 1;
                inp.index[pos++] = inp.nC + *it;
            }
        }
    }
    inp.start[inp.nC] = pos;
    // NOTE : Remove ASSERT from final version
    assert(pos == nEls);

    /* cout << "Column Major Format :: " << endl;
     * for (int j = 0; j < inp.nC; j++)
     * {
     *     cout << "  var " << j << " from " << inp.start[j] << " to " << inp.start[j+1] << " :: " << endl;
     *     for (int k = inp.start[j]; k < inp.start[j+1]; k++)
     *         cout << "... el " << inp.W[k] << " in row " << inp.index[k] << endl;
     * } */
}

/// Save Budget Uncertainty Set Info on disk.
/** We save the parameters and the data needed to recreate the budget instance.
 *  This is done to ensure reproducibility of a robust instance. If we decide
 *  to save multiple versions of robust instances starting from the same
 *  nominal instances, we will have to come up with a different name coding.
 */
void save_instance_2_disk(double _epsilon, double _delta, double _gamma, int L, 
                          int nBl, int ** Bl, double * budget)
{

    string  s1      = string(_FILENAME);
    s1              = s1.substr(s1.find_last_of("\\/"), 100);
    string filename = "support" + s1 + ".budget"; 
    ofstream fWriter(filename, ios::out);

    for (int l = 0; l < L; l++)
    {
        fWriter << nBl;
        for (int k = 0; k < nBl; k++)
            fWriter << " " << Bl[l][k];
        fWriter << " " << budget[l] << endl;
    }

    fWriter.close();
    cout << "[** Instance saved on disk. File '" << filename << "']" << endl;
}

/// Read a robust instance based on budget uncertainty set.
/**
 *  This ensures reproducibility. We read the set \f$B_l\f$ and the budget
 *  values \f$b_l\f$. Thus, the robust instance can fully be reconstructed.
 * */
void read_instance_from_disk(double & _epsilon, double & _delta, double & _gamma, 
                             int & L, int & nBl, int ** Bl, double * budget)
{
    string  s1      = string(_FILENAME);
    s1              = s1.substr(s1.find_last_of("\\/"), 100);
    string filename = "support" + s1 + ".budget"; 
    /* string filename = "support" + s1 + "-budget-" + to_string(_epsilon) +"-" +
     *                   to_string(_delta) + "-" + to_string(_gamma) + "-" + to_string(L); */

    ifstream fReader(filename, ios::in);
    if (!fReader)
    {
        cout << "Cannot open file '" << filename << "'." << endl;
        exit(111);
    }
    fReader >> _epsilon >> _delta >> _gamma >> L;
    Bl     = new int*[L];
    budget = new double[L];
    for (int l = 0; l < L; l++)
    {
        fReader >> nBl;
        Bl[l] = new int[nBl];

        for (int k = 0; k < nBl; k++)
            fReader >> Bl[l][k];
        fReader >> budget[l];
    }

    fReader.close();
    cout << "[** Instance read from disk. File '" << filename << "']" << endl;

}

/// Benders decomposition algorithm
/** We use the Benders algorithm provided by cplex, with the following
 * strategy:
 * - `IloCplex::BendersUser`
 * In reality, the default strategy `BendersFull`, in which cplex chooses how
 * to assign variables to the master, would also work.
 * 
 */
void define_benders(IloModel & model, IloCplex & cplex, INSTANCE inp)
{
    IloEnv env = model.getEnv();

    // by default, variables are assigned to the subproblem (value 1)
    IloCplex::LongAnnotation benders = cplex.newLongAnnotation("cpxBendersPartition",1);

	 // Set benders strategy 
    cplex.setParam(IloCplex::Param::Benders::Strategy, IloCplex::BendersUser);
	// cplex.setParam(IloCplex::Param::Benders::Strategy, IloCplex::BendersFull);

    // assign location variables to the master
    for (int i = 0; i < inp.nF; i++)
        cplex.setAnnotation(benders, y_ilo[i], 0);

    // save annotation to disk
    cplex.writeBendersAnnotation("benders.ann");


}
