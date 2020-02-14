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

/*! \file binPacking.cpp 
  \brief Manage Bin Packing problem

 * We define, implement, and solve the Bin Packing problem.

*/

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <fstream>
#include <vector>
#include "binPacking.h"

using namespace std;


extern const double EPSI;
extern double _epsilon;
extern double _Omega;
extern double _delta;
extern double _gamma;
extern int    L;      

extern string instanceType;
extern string versionType;
extern string supportType;
extern int support;

extern IloNumVarArray y_ilo;
typedef IloArray <IloNumVarArray> TwoD;
extern TwoD x_ilo;


void define_box_support(InstanceBin & inp)
{
   inp.nR = 2*inp.nI;
   inp.h  = new double[inp.nR];
   for (int j = 0; j < inp.nI; j++)
   {
       inp.h[j]        = -inp.d[j]*(1.0-_epsilon);
       inp.h[j+inp.nI] =  inp.d[j]*(1.0+_epsilon);
   }

   // define matrix W in column major format
   inp.W     = new int[inp.nI*2];
   inp.index = new int[inp.nI*2];
   inp.start = new int[inp.nI+1]; // one extra element in last position

    int pos = 0;
    for (int j = 0; j < inp.nI; j++)
    {
        inp.start[j]   = pos;
        inp.W[pos]     = -1;
        inp.index[pos] = j;
        pos++;
        inp.W[pos]     = 1;
        inp.index[pos] = j+inp.nI;
        pos++;
    }
    inp.start[inp.nI] = pos;
    assert(pos == inp.nC*2);
}


void define_POLY_BINPACKING(InstanceBin & inpBin, IloModel & model, IloCplex & cplex, int support)
{
    // select support type
    // here we get W and h, depending on the type of support
    switch (support) {
        case 1 :
            define_box_support(inpBin);
            break;
        case 2 :
            // define_budget_support(inpBin, false);
            cout << "BUDGET SUPPORT not defined here!!!\n";
            break;
        default :
            cout << "ERROR : Support type not defined.\n" << endl;
            exit(123);
    }

    char varName[100];
    IloEnv env = model.getEnv();

    y_ilo = IloNumVarArray(env, inpBin.nI, 0, 1, ILOINT);
    x_ilo = TwoD(env, inpBin.nI);
    for (int i = 0; i < inpBin.nI; i++)
        x_ilo[i] = IloNumVarArray(env, inpBin.nI, 0, 1, ILOINT);

    // set vars names
    for (int i = 0; i < inpBin.nI; i++) {
        sprintf(varName, "y.%d", (int)i);
        y_ilo[i].setName(varName);
        for (int j = 0; j < inpBin.nI; j++) {
            sprintf(varName, "x.%d.%d", (int)i, (int)j);
            x_ilo[i][j].setName(varName);
        }
    }

    // packing constraints (each item in a bin)
    for (int i = 0; i < inpBin.nI; i++) {
        IloExpr sum(env);
        for (int j = 0; j < inpBin.nI; j++)
            sum += x_ilo[i][j];
        model.add(sum == 1.0);
    }

    // capacity constraint
    for (int j = 0; j < inpBin.nI; j++) {
        IloExpr sum(env);
        for (int i = 0; i < inpBin.nI; i++)
            sum += (1.+_epsilon)*inpBin.d[i]*x_ilo[i][j]; // update this in robust formulation
        sum -= y_ilo[j]*inpBin.totS;

        model.add(sum <= 0.0);
    }


    // objective function: min sum y_j
    IloExpr totCost(env);
    for (int j = 0; j < inpBin.nI; j++)
        totCost += y_ilo[j];

    model.add(IloMinimize(env,totCost));
}
void getCplexSolBIN(InstanceBin inpBin, IloCplex cplex, SOLUTION & opt)
{
    opt.nOpen = 0;
    opt.ySol = new int[inpBin.nI];
    opt.xSol = new double*[inpBin.nI];
    for (int i = 0; i < inpBin.nI; i++)
        opt.xSol[i] = new double[inpBin.nI];

    opt.zStar = cplex.getObjValue();
    opt.zStatus = cplex.getStatus();

    for (int i = 0; i < inpBin.nI; i++)
        if (cplex.getValue(y_ilo[i]) >= 1.0-EPSI)
        {
            opt.ySol[i] = 1;
            opt.nOpen++;
        }
        else
            opt.ySol[i] = 0;

    for (int i = 0; i < inpBin.nI; i++)
        for (int j = 0; j < inpBin.nI; j++)
            opt.xSol[i][j] = cplex.getValue(x_ilo[i][j]);

    if (support==1)
        versionType = "box";
    else if (support== 2)
        versionType = "budget";
    else {
        cout << "ERROR in WRITING SOLUTION: Version type not defined.\n" << endl;
        exit(123);
    }

    string  s1      = string(_FILENAME);
    s1              = s1.substr(s1.find_last_of("\\/"), 100);
    ostringstream obj;
    obj << _epsilon << "-" <<  _delta << "-" << _gamma << "-" << L;
    string filename = "solutions" + s1 + "-" + versionType + "-" + obj.str();

    ofstream fWriter(filename, ios::out);
    fWriter << inpBin.nI << endl;
    fWriter << setprecision(15) << opt.zStar << endl;
    fWriter << opt.zStatus << endl;
    fWriter << opt.nOpen << endl;
    for (int i = 0; i < inpBin.nI; i++)
        if (opt.ySol[i] == 1)
            fWriter << " " << i;
    fWriter << endl;
    for (int i = 0; i < inpBin.nI; i++)
        for (int j = 0; j < inpBin.nI; j++)
            if (opt.xSol[i][j] > 0.0)
                fWriter << i << " " << j << " " << opt.xSol[i][j] << endl;

    fWriter.close();
    cout << "Solution written to disk. ('" << filename <<"')" << endl;

    // TEST: check that used capacity does not exceed capacity of each bin
    for (int i = 0; i < inpBin.nI; i++) {
        if (opt.ySol[i] == 1) {
            float tot = 0.0;
            for (int j = 0; j < inpBin.nI; j++)
                if (opt.xSol[i][j] > 0.0)
                    tot += opt.xSol[i][j]*inpBin.d[j];
            // cout << "Facility " << i << " with capacity = " << inp.s[i] << " and consumption = " << tot << endl;
            assert(inpBin.totS >= tot);
        }
    }
}


