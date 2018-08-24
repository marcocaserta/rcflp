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
#include <cmath>
#include <climits>
#include <cfloat>
#include <cassert>

/* #include "timer.h" */
#include "options.h"

using namespace std;

double INFTY = std::numeric_limits<double>::infinity();
const long _MAXRANDOM  = 2147483647;
const double ZERO      = 0.0e0;
const double EPSI      = 0.00001;


/****************** VARIABLES DECLARATION ***************************/
char * _FILENAME;		//!< Instance name file
int fType;              //!< instance type (1-4)
int version;            //!< 1-SS; 2-MS; 3-SOCP
string instanceType;
string versionType;

/// Structure used to define the instance data
struct INSTANCE {
    int nF;        // number of facilities
    int nC;        // number of customers
    int nR;        // number of constraints polyhedron uncertainty set
    double  *f;    // fixed costs
    double  *s;    // capacity
    double  *d;    // demand
    double **c;    // allocation costs
    double   totS; // total supply
    double   totD; // total demand
};
INSTANCE inp; // instance data

// optimal solution and obj function value
struct SOLUTION {
    int  *ySol;
    double **xSol;
    double zStar;
    IloAlgorithm::Status zStatus;
};
SOLUTION opt;

typedef IloArray <IloNumVarArray> TwoD;

/**** CPLEX DEFINITION ****/
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
int timeLimit;
double Omega     = 1.0;
double epsi      = 0.1;
int solLimit     = 9999;
int displayLimit = 4;
/****************** VARIABLES DECLARATION ***************************/

/****************** FUNCTIONS DECLARATION ***************************/
int readProblemData(char * _FILENAME, int fType, INSTANCE & inp);
void printOptions(char * _FILENAME, int fType, int version, INSTANCE inp, int timeLimit);
void define_MS_CFLP(INSTANCE inp, int fType, IloModel & model, IloCplex & cplex);
void define_SS_CFLP(INSTANCE inp, int fType, IloModel & model, IloCplex & cplex);
void define_SOCP_CFLP(INSTANCE inp, int fType, IloModel & model, IloCplex & cplex);
void define_POLY_CFLP(INSTANCE inp, int fType, IloModel & model, IloCplex & cplex);
int solveCplexProblem(IloModel model, IloCplex cplex, INSTANCE inp, int solLimit, int timeLimit, int displayLimit);
void getCplexSol(INSTANCE inp, IloCplex cplex, SOLUTION & opt);
void printSolution(char * _FILENAME, INSTANCE inp, SOLUTION opt, bool toDisk, bool fullOutput);
/****************** FUNCTIONS DECLARATION ***************************/

/************************ main program ******************************/
/// main program
/************************ main program ******************************/
int main(int argc, char *argv[])
{
    int err = parseOptions(argc, argv);
    if (err != 0) exit(1);

    readProblemData(_FILENAME, fType, inp);
    printOptions(_FILENAME, fType, version, inp, timeLimit);


    if (version == 1) // single source nominal
        define_SS_CFLP(inp, fType, model, cplex);
    else if (version == 2) // multi source nominal
        define_MS_CFLP(inp, fType, model, cplex);
    else if (version == 3) // multi source ellipsoidal
        define_SOCP_CFLP(inp, fType, model, cplex);
    else if (version == 4)
        define_POLY_CFLP(inp, fType, model, cplex);

    solveCplexProblem(model, cplex, inp, solLimit, timeLimit, displayLimit);

    getCplexSol(inp, cplex, opt);
    printSolution(_FILENAME, inp, opt, true, true);

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
bool fullOutput)
{
    cout << endl << "** SOLUTION **" << endl;
    cout << " ..z*     = " << setprecision(15) << opt.zStar << endl;
    cout << " ..status = " << opt.zStatus << endl;

    {
        ofstream fWriter("solution.txt", ios::out);
        fWriter << _FILENAME << "\t" << instanceType << "\t"
                << versionType << "\t" << setprecision(15) << opt.zStar
                << "\t" << opt.zStatus << endl;
        fWriter.close();
    }

    if (fullOutput)
    {
        cout << "Open Facilities: ";
        for (int i = 0; i < inp.nF; i++)
            if (opt.ySol[i] == 1)
                cout << setw(4) << i;
        cout << endl;
        for (int i = 0; i < inp.nF; i++)
            for (int j = 0; j < inp.nC; j++)
                if (opt.xSol[i][j] > 0.0)
                    cout << "x(" << i << "," << j << ") = " << setprecision(3) << opt.xSol[i][j] << endl;
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

    // facility capacity (depends on instance type)
    for (int i = 0; i < inp.nF; i++)
    {
        IloExpr sum(env);
        for (int j = 0; j < inp.nC; j++)
            sum += x_ilo[i][j]*inp.d[j];
        sum -= y_ilo[i]*inp.s[i];

        model.add(sum <= 0.0);
    }

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

void define_POLY_CFLP(INSTANCE inp, int fType, IloModel & model, IloCplex & cplex)
{
    inp.nR = 2*inp.nC; // this is equivalent to a BOX UNCERTAINTY set
    int lower_epsilon = 0;   // box width around nominal value
    int upper_epsilon = 1;   // box width around nominal value

    char varName[100];
    IloEnv env = model.getEnv();

    // location variables
    y_ilo = IloNumVarArray(env, inp.nF, 0, 1, ILOINT);

    // allocation variables
    x_ilo = TwoD(env, inp.nF);
    for (int i = 0; i < inp.nF; i++)
        x_ilo[i] = IloNumVarArray(env, inp.nC, 0.0, 1.0, ILOFLOAT);
        /* x_ilo[i] = IloNumVarArray(env, inp.nC, 0.0, 1.0, ILOINT); */

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

    // delta variable
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

    // "robust" capacity demand
    for (int i = 0; i < inp.nF; i++)
        for (int j = 0; j < inp.nC; j++)
            model.add(-1.0*psi_ilo[i][j] + 1.0*psi_ilo[i][j+inp.nC] - x_ilo[i][j] >= 0.0);
    for (int i = 0; i < inp.nF; i++)
    {
        IloExpr sum(env);
        for (int j = 0; j < inp.nC; j++)
        {
            sum += -(inp.d[j]-lower_epsilon)*psi_ilo[i][j];
            sum +=  (inp.d[j]+upper_epsilon)*psi_ilo[i][j+inp.nC];
        }
        sum -= inp.s[i]*y_ilo[i];
        model.add(sum <= 0.0);
    }

    // "robust" objective function
    IloExpr sum(env);
    for (int j = 0; j < inp.nC; j++)
    {
        sum += -(inp.d[j]-lower_epsilon)*u_ilo[j];
        sum +=  (inp.d[j]+upper_epsilon)*u_ilo[j+inp.nC];
    }
    sum -= delta_ilo[0];
    model.add(sum <= 0.0);
    for (int j = 0; j < inp.nC; j++)
    {
        IloExpr sum(env);
        for (int i = 0; i < inp.nF; i++)
            sum += inp.c[i][j]*x_ilo[i][j];

        model.add(-1.0*u_ilo[j] + 1.0*u_ilo[j+inp.nC] - sum >= 0.0);
    }


    // objective function
    IloExpr totCost(env);
    for (int i = 0; i < inp.nF; i++)
        totCost += y_ilo[i]*inp.f[i];
    totCost += delta_ilo[0];

    model.add(IloMinimize(env,totCost));
    

}
