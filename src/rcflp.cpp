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

  We define here the following variants of the CFLP, controlled via coommand 
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
int version;            //!< 1-SS; 2-MS
string instanceType;
string versionType;

/// Structure used to define the instance data
struct INSTANCE {
    int nF;
    int nC;
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
int timeLimit;
double Omega     = 1.0;
double epsi      = 0.001;
int solLimit     = 9999;
int displayLimit = 4;
/****************** VARIABLES DECLARATION ***************************/

/****************** FUNCTIONS DECLARATION ***************************/
/* void readProblemData(char * _FILENAME, int fType); */
int readProblemData(char * _FILENAME, int fType, INSTANCE & inp);
void printOptions(char * _FILENAME, int fType, int version, INSTANCE inp, int timeLimit);
void define_MS_CFLP(INSTANCE inp, int fType, IloModel & model, IloCplex & cplex);
void define_SS_CFLP(INSTANCE inp, int fType, IloModel & model, IloCplex & cplex);
void define_SOCP_CFLP(INSTANCE inp, int fType, IloModel & model, IloCplex & cplex);
int solveCplexProblem(IloModel model, IloCplex cplex, INSTANCE inp, int solLimit, int timeLimit, int displayLimit);
void getCplexSol(INSTANCE inp, IloCplex cplex, SOLUTION & opt);
void printSolution(char * _FILENAME, INSTANCE inp, SOLUTION opt, bool toDisk);
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

    solveCplexProblem(model, cplex, inp, solLimit, timeLimit, displayLimit);

    getCplexSol(inp, cplex, opt);
    printSolution(_FILENAME, inp, opt, true);

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
void printSolution(char * _FILENAME, INSTANCE inp, SOLUTION opt, bool toDisk)
{
    cout << endl << "** SOLUTION **" << endl;
    cout << " ..z*     = " << setprecision(15) << opt.zStar << endl;
    cout << " ..status = " << opt.zStatus << endl;
    if (toDisk)
    {
        ofstream fWriter("solution.txt", ios::out);
        fWriter << _FILENAME << "\t" << instanceType << "\t"
                << versionType << "\t" << setprecision(15) << opt.zStar
                << "\t" << opt.zStatus << endl;
        fWriter.close();
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
    if (fType == 1)
        for (int i = 0; i < inp.nF; i++)
            x_ilo[i] = IloNumVarArray(env, inp.nC, 0.0, 1.0, ILOFLOAT);
    else if (fType == 2)
        for (int i = 0; i < inp.nF; i++)
            x_ilo[i] = IloNumVarArray(env, inp.nC, 0.0, inp.s[i], ILOFLOAT);

    // set var names
    for (int i = 0; i < inp.nF; i++)
    {
        sprintf(varName, "y.%d", (int)i);
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
        if (fType == 1)
        {
            for (int i = 0; i < inp.nF; i++)
                sum += x_ilo[i][j];
            model.add(sum == 1.0);
        }
        else if (fType == 2)
        {
            for (int i = 0; i < inp.nF; i++)
                sum += x_ilo[i][j];
            model.add(sum == inp.d[j]);
        }
    }

    // facility capacity (depends on instance type)
    for (int i = 0; i < inp.nF; i++)
    {
        IloExpr sum(env);
        if (fType == 1)
        {
            for (int j = 0; j < inp.nC; j++)
                sum += x_ilo[i][j]*inp.d[j];
        }
        else if (fType == 2)
        {
            for (int j = 0; j < inp.nC; j++)
                sum += x_ilo[i][j];
        }
        sum -= y_ilo[i]*inp.s[i];

        model.add(sum <= 0.0);
    }

    // does it tighten the formulation?
    /*     double coeff;
     *     for (int i = 0; i < inp.nF; i++)
     *         for (int j = 0; j < inp.nC; j++)
     *         {
     *             if (fType == 1 || version == 1)
     *                 coeff = 1.0;
     *             else
     *                 // only in case fTYpe = 2 (Avella) and version = 2 (Single-source)
     *                 coeff = inp.s[i];
     *
     *             model.add(x_ilo[i][j] - coeff*y_ilo[i] <= 0.0);
     *         } */


    // objective function
    IloExpr totCost(env);
    for (int i = 0; i < inp.nF; i++)
    {
        totCost += y_ilo[i]*inp.f[i];
        for (int j = 0; j < inp.nC; j++)
            totCost += x_ilo[i][j]*inp.c[i][j];
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

    // does it tighten the formulation?
    for (int i = 0; i < inp.nF; i++)
        for (int j = 0; j < inp.nC; j++)
            model.add(x_ilo[i][j] - y_ilo[i] <= 0.0);

    // objective function
    IloExpr totCost(env);
    double coeff = 0.0;
    for (int i = 0; i < inp.nF; i++)
    {
        totCost += y_ilo[i]*inp.f[i];
        for (int j = 0; j < inp.nC; j++)
        {
            if (fType == 1)
                coeff = inp.c[i][j];
            else if (fType == 2)
                coeff = inp.c[i][j]*inp.d[j];
            totCost += x_ilo[i][j]*coeff;
        }
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
    // single source case
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

    // second order cone W (for OR Library)
    IloExpr sum(env);
    sum = -w_ilo;
    // sum = -w_ilo*w_ilo;
    for (int i = 0; i < inp.nF; i++)
        for (int j = 0; j < inp.nC; j++)
            sum += x_ilo[i][j]*x_ilo[i][j]*inp.c[i][j]*epsi*inp.c[i][j]*epsi;
    model.add(sum <= 0.0);

    // Q conic constraints
    for (int i = 0; i < inp.nF; i++)
    {
        IloExpr sum(env);
        sum = -q_ilo[i];
        // sum = -q_ilo[i]*q_ilo[i];
        for (int j = 0; j < inp.nC; j++)
            sum += x_ilo[i][j]*x_ilo[i][j]*inp.c[i][j]*epsi*inp.c[i][j]*epsi;
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
            totCost += x_ilo[i][j]*inp.c[i][j];
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

