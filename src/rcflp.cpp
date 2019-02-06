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

  Description here

  \author 
  \version v. 1.0.0
  \date Begins: 25.05.18
  \date Ends:

  The model is:

  \f[ 
  \f]


  \file rcflp.cpp
  \brief General Implementation of the compact formulations for the (R)-CFLP.

  We define here the following variants of the CFLP, controlled via command 
  line with the flag -v:
  -v 1 :: Single Source
  -v 2 :: Multi Source

  In addition, we need to define the type of instance used. We are currently
  working with the following instances:
  -t 1 :: OR Library instances
  -t 2 :: Avella instances (Type 1, Type A, Type B)
  Note that these instances define the costs in different ways and, therefore
  the way in which the x_ij variables are defined changes. More precisely:
  - for the OR Library instances, the c_ij values are the cost of delivering
  the entire demand to the customer. Thus, we define x_ij as fractional and
  multiply c_ij*x_ij to get the real cost
  - for the Avella instances, c_ij is the cost of transporting one unit. Thus,
  x_ij is no longer the fraction of demand, but the actual number of units
  transported.
  Obviously, this different treatment of the x_ij variables has an effect on 
  the definition of the demand and capacity constraints as well.

  Therefore, to make the treatment of the two types of instances uniform, the
  cost c_ij of the OR Library is divided by the total demand:
    c_ij = c_ij / d_j
  To summarize, for both instance types (OR Library and Avella), we have:
  - x_ij \in [0,1] : percentage of demand
  - c_ij : cost of delivering one unit of demand
  - d_j  : total demand of customer j

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
    d_jx_{ij}c_{ij}
  \f]
*/

#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

#include <limits> 

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

/* #include "timer.h" */
#include "options.h"

using namespace std;

double INFTY = std::numeric_limits<double>::infinity();
const long _MAXRANDOM  = 2147483647;
const double ZERO      = 0.0e0;
const double EPSI      = 0.00001;
// Initialization of Random Number Generator
const int seed = 27;
// 64-bit Mersenne Twister by Matsumoto and Nishimura, 2000
mt19937_64 gen(seed);
/* mt19937_64 gen(seed()); */

/****************** VARIABLES DECLARATION ***************************/
char * _FILENAME;		//!< Instance name file
int fType;              //!< instance type (1-4)
int version;            //!< 1-SS; 2-MS; 3-SOCP
int support;            //!< 1-Box; 2-Budget
string instanceType;
string versionType;
string supportType;

typedef std::vector <int> MyVect;

/// Structure used to define the instance data
// NOTE: Change the same structure in the file inout.cpp !!!
struct INSTANCE {
    int nF;        // number of facilities
    int nC;        // number of customers
    double  *f;    // fixed costs
    double  *s;    // capacity
    double  *d;    // demand
    double **c;    // allocation costs
    double   totS; // total supply
    double   totD; // total demand

    int     nR;    // number of constraints polyhedron uncertainty set
    double  *h;    // rhs of polyhedron definining support
    int     *W;    // matrix W in column major format
    int *index;    // index of column major format for w
    int *start;    // starting position for elements of column j

};
INSTANCE inp; // instance data

// optimal solution and obj function value
struct SOLUTION {
    int  *ySol;
    double **xSol;
    double zStar;
    IloAlgorithm::Status zStatus;
};
SOLUTION opt; // solution data structure


/**** CPLEX DEFINITION ****/
typedef IloArray <IloNumVarArray> TwoD;
IloEnv env;
IloModel model(env, "cflp");
IloCplex cplex(model);
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
double Omega     = 1.0;
double epsi      = 0.1;

/****************** FUNCTIONS DECLARATION ***************************/
int readProblemData(char * _FILENAME, int fType, INSTANCE & inp);
void printOptions(char * _FILENAME, INSTANCE inp, int timeLimit);
void define_MS_CFLP(INSTANCE inp, int fType, IloModel & model, IloCplex & cplex);
void define_SS_CFLP(INSTANCE inp, int fType, IloModel & model, IloCplex & cplex);
void define_SOCP_CFLP(INSTANCE inp, int fType, IloModel & model, IloCplex & cplex);
void define_POLY_CFLP(INSTANCE inp, int fType, IloModel & model, IloCplex & cplex, int support);
int solveCplexProblem(IloModel model, IloCplex cplex, INSTANCE inp, int solLimit, int timeLimit, int displayLimit);
void getCplexSol(INSTANCE inp, IloCplex cplex, SOLUTION & opt);
void printSolution(char * _FILENAME, INSTANCE inp, SOLUTION opt, bool toDisk, int fullOutput);
void read_parameters_box(double & delta);
void define_box_support(INSTANCE & inp);
void define_budget_support(INSTANCE & inp, bool fromDisk);
void save_instance_2_disk(double epsilon, double delta, double gamma, int L, 
                          int nBl, int ** Bl, double * budget);
/****************** FUNCTIONS DECLARATION ***************************/

/************************ main program ******************************/
/// main program
/************************ main program ******************************/
int main(int argc, char *argv[])
{
    int err = parseOptions(argc, argv);
    if (err != 0) exit(1);

    readProblemData(_FILENAME, fType, inp);
    printOptions(_FILENAME, inp, timeLimit);



    if (version == 1) // single source nominal
        define_SS_CFLP(inp, fType, model, cplex);
    else if (version == 2) // multi source nominal
        define_MS_CFLP(inp, fType, model, cplex);
    else if (version == 3) // multi source ellipsoidal
        define_SOCP_CFLP(inp, fType, model, cplex);
    else if (version == 4) // robust polyhedral uncertainty set (both SS and MS)
        define_POLY_CFLP(inp, fType, model, cplex, support);

    solveCplexProblem(model, cplex, inp, solLimit, timeLimit, displayLimit);

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
    opt.ySol = new int[inp.nF];
    opt.xSol = new double*[inp.nF];
    for (int i = 0; i < inp.nF; i++)
        opt.xSol[i] = new double[inp.nC];

    opt.zStar = cplex.getObjValue();
    opt.zStatus = cplex.getStatus();

    for (int i = 0; i < inp.nF; i++)
        if (cplex.getValue(y_ilo[i]) >= 1.0-EPSI)
            opt.ySol[i] = 1;
        else
            opt.ySol[i] = 0;

    for (int i = 0; i < inp.nF; i++)
        for (int j = 0; j < inp.nC; j++)
            opt.xSol[i][j] = cplex.getValue(x_ilo[i][j]);
}


/// Print solution to screen (and disk, if required)
void printSolution(char * _FILENAME, INSTANCE inp, SOLUTION opt, bool toDisk, 
int fullOutput)
{
    cout << endl << "** SOLUTION **" << endl;
    cout << " ..z*     = " << setprecision(15) << opt.zStar << endl;
    cout << " ..status = " << opt.zStatus << endl;

    {
        ofstream fWriter("solution.txt", ios::out);
        fWriter << _FILENAME << "\t" << instanceType << "\t"
                << versionType << "\t" << supportType << "\t" 
                << setprecision(15) << opt.zStar << "\t" 
                << opt.zStatus << endl;
        fWriter.close();
    }

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

    char varName[100];
    IloEnv env = model.getEnv();

    // location variables 
    y_ilo = IloNumVarArray(env, inp.nF, 0, 1, ILOINT);

    // allocation variables
    x_ilo = TwoD(env, inp.nF);
    for (int i = 0; i < inp.nF; i++)
        x_ilo[i] = IloNumVarArray(env, inp.nC, 0.0, 1.0, ILOFLOAT);

    // set vars name
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
            sum += x_ilo[i][j]*x_ilo[i][j]*inp.c[i][j]*epsi*inp.c[i][j]*epsi;
    model.add(sum <= 0.0);

    // Q conic constraints
    for (int i = 0; i < inp.nF; i++)
    {
        IloExpr sum(env);
        // sum = -q_ilo[i];
        sum = -q_ilo[i]*q_ilo[i];
        for (int j = 0; j < inp.nC; j++)
            sum += x_ilo[i][j]*x_ilo[i][j]*epsi*epsi;
        model.add(sum <= 0.0);
    }

    // capacity constraints QUADRATIC
    for (int i = 0; i < inp.nF; i++)
    {
        IloExpr sum(env);
        for (int j = 0; j < inp.nC; j++)
            sum += x_ilo[i][j]*inp.d[j];

        sum += Omega*q_ilo[i];
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
    totCost += Omega*w_ilo;

    model.add(IloMinimize(env,totCost));
}

/// Define robust model based on polyhedral uncertainty set
/**
 *
 * We implement here both the SINGLE SOURCE and the MULTI SOURCE versions. Just
 * change the way in which the allocation variables are defined, to switch
 * between the SS and the MS versions.
 *
 * To test it, use the toy problem, whose instance values are defined in the
 * file 'toy.txt'. These values are the one reported in the paper, under the
 * example section. The budget constraint, though, has not been implemented
 * yet. The following tests can be executed:
 * 1. nominal version: Set lower_epsilon and upper_epsilon = 0
 * 2. multi-source: Define allocation variables as float (and give some slack
 * to the demand assigning positive values to lower and upper epsilon.
 * 3. single-source: Define allocation variables are integer.
 */
void define_POLY_CFLP(INSTANCE inp, int fType, IloModel & model, IloCplex & cplex,
                      int support)
{

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

    // allocation variables
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

/// Set cplex parameters and solve the optimization problem
int solveCplexProblem(IloModel model, IloCplex cplex, INSTANCE inp, int solLimit, int timeLimit, int displayLimit)
{
    try
    {
        IloEnv env = model.getEnv();
        /* cplex.setOut(env.getNullStream()); */
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

/// Define W and h, given the value of epsilon
/**
 * We need to define Wd <= h, where:
 *      W is of size [r x n]
 *      h is of size [r x 1]
 * where:
 * r = 2*n (two elements for each column j)
 * n = number of customers
 * 
 * Since matrix W is sparse, we use a column major format as follows. Consider
 * the toy problem, for which matrix W and vector h are: 
 *  
 *    |-1  0  0 |           |-d_0*(1-epsilon)|
 *    | 0 -1  0 |           |-d_1*(1-epsilon)|
 * W =| 0  0  0 |           |-d_2*(1-epsilon)|
 *    | 1  0  0 |       h = | d_0*(1+epsilon)|
 *    | 0  1  0 |           | d_1*(1+epsilon)|
 *    | 0  0  1 |           | d_2*(1+epsilon)|
 *
 * Vector h is kept as is it, while matrix W uses a column major format:
 * W     = [-1 1  -1 1  -1 1]
 * start = [0 2 4 6]
 * index = [0 3  1 4  2 5]
 * 
 * where elements of column j in W and index are found in positions going from
 * start[j] to start[j+1] (note that vector "start" contains a final extra
 * element to close the cycle). Thus, e.g.,  the elements of the second column
 * of W (j = 1) are obtained as:
 * for l = start[j], ..., start[j+1]
 *      W[l] is the element in row index[l] of the matrix
 */
void define_box_support(INSTANCE & inp)
{
   double epsilon = 0.0;
   read_parameters_box(epsilon);

   inp.nR = 2*inp.nC;
   inp.h  = new double[inp.nR];
   for (int j = 0; j < inp.nC; j++)
   {
       inp.h[j]        = -inp.d[j]*(1.0-epsilon);
       inp.h[j+inp.nC] =  inp.d[j]*(1.0+epsilon);
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


void define_budget_support(INSTANCE & inp, bool fromDisk)
{

    double epsilon = 0.0;
    double delta   = 0.0;
    double gamma   = 0.0;
    int    L       = 0;
    if (fromDisk==false)
    {
        read_parameters_budget(epsilon, delta, gamma, L);

        uniform_int_distribution<> d(0,inp.nC);

        // cardinality of sets B_l
        int nBl = floor(gamma*(double)inp.nC);
        cout << "[** |B_l| = " << nBl << "]" << endl;
        int ** Bl = new int*[L];
        double  * budget = new double[L];
        for (int l = 0; l < L; l++)
            Bl[l] = new int[nBl];

        for (int l = 0; l < L; l++)
        {
            budget[l] = 0.0;
            for (int k = 0; k < nBl; k++)
            {
                int el = d(gen);
                Bl[l][k] = el;
                budget[l] += inp.d[el];
            }
            cout << endl;
        }
        // adjust b_l values
        for (int l = 0; l < L; l++)
        {
            budget[l] = floor(delta*budget[l]);
            cout << "budget = " << budget[l] << endl;
        }
        for (int l = 0; l < L; l++)
        {
            for (int k = 0; k < nBl; k++)
                cout << " " << Bl[l][k];
            cout << endl;
        }

        // define mapping: list of budget constraints including column j 
        MyVect * mapping = new MyVect[inp.nC];
        for (int l = 0; l < L; l++)
            for (int k = 0; k < nBl; k++)
                mapping[Bl[l][k]].push_back(l);


       // total number of rows of W and h
       inp.nR = 2*inp.nC + L;
        // define vector h
       inp.h  = new double[inp.nR];
       for (int j = 0; j < inp.nC; j++)
       {
           inp.h[j]        = -inp.d[j]*(1.0-epsilon);
           inp.h[j+inp.nC] =  inp.d[j]*(1.0+epsilon);
       }

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

        // do we want to save the sets B_l (and the parameters?)
        bool save2Disk = false;
        if (save2Disk==true)
            save_instance_2_disk(epsilon, delta, gamma, L, nBl, Bl, budget);
    }
}

/// Save Budget Uncertainty Set Info on disk
/** We save the parameters and the data needed to recreate the budget instance.
 */
void save_instance_2_disk(double epsilon, double delta, double gamma, int L, 
                          int nBl, int ** Bl, double * budget)
{

    string  s1      = string(_FILENAME);
    s1              = s1.substr(s1.find_last_of("\\/"), 100);
    string filename = "support" + s1 + ".budget";

    ofstream fWriter(filename, ios::out);

    fWriter << epsilon << endl;
    fWriter << delta << endl;
    fWriter << gamma << endl;
    fWriter << L << endl;
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
