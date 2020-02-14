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

/*! \file mmcf.cpp 
  \brief Manage Bin Packing problem

 * We define, implement, and solve the Bin Packing problem.

*/

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <fstream>
#include <vector>
#include "mmcf.h"

using namespace std;



void define_POLY_MMCF(InstanceMMCF & inpMMCF, IloModel & model, IloCplex & cplex, int support)
{
    // select support type
    // here we get W and h, depending on the type of support
    switch (support) {
        case 1 :
            // define_box_support(inpBin);
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


    y_ilo = IloNumVarArray(env, inpMMCF.nArcs, 0, 1, ILOINT);
    x_ilo = TwoD(env, inpMMCF.nArcs);
    for (int i = 0; i < inpMMCF.nArcs; i++)
        x_ilo[i] = IloNumVarArray(env, inpMMCF.nK, 0, IloInfinity, ILOFLOAT);

    // set vars names
    for (int i = 0; i < inpMMCF.nArcs; i++) {
        sprintf(varName, "y.%d", (int)i);
        y_ilo[i].setName(varName);
        for (int j = 0; j < inpMMCF.nK; j++) {
            sprintf(varName, "x.%d.%d", (int)i, (int)j);
            x_ilo[i][j].setName(varName);
        }
    }


    // flow constraints
    // flow at node n of commodity k
    for (int k = 0; k < inpMMCF.nK; k++) {

        std::vector<double> d(inpMMCF.nNodes, 0.0);

        d[inpMMCF.source[k]] = -inpMMCF.d[k]*(1.0+_epsilon);
        d[inpMMCF.sink[k]]   = inpMMCF.d[k]*(1.0+_epsilon);
        /* cout << " \nComm " << k << ":: ";
         * for (int n = 0; n < inpMMCF.nNodes; n++)
         *     cout << " " << d[n]; */

        for (int n = 0; n < inpMMCF.nNodes; n++) {
            IloExpr sum(env);
            // inbound elements
            for (int i = 0; i < inpMMCF.inbound[n].size(); i++) {
                int el = inpMMCF.inbound[n][i];
                sum += x_ilo[el][k];
            }
            
            // outbound elements
            for (int i = 0; i < inpMMCF.outbound[n].size(); i++) {
                int el = inpMMCF.outbound[n][i];
                sum -= x_ilo[el][k];
            }

            model.add(sum >= d[n]);
        }
    }

    // logical constraints
    for (int i = 0; i < inpMMCF.nArcs; i++) {
        IloExpr sum(env);
        for (int k = 0; k < inpMMCF.nK; k++)
            sum += x_ilo[i][k];

        sum -= y_ilo[i]*inpMMCF.cap[i];

        model.add(sum <= 0.0);
    }

    IloExpr totCost(env);
    for (int i = 0; i < inpMMCF.nArcs; i++) {
        totCost += y_ilo[i]*inpMMCF.fixed[i];
        for (int k = 0; k < inpMMCF.nK; k++)
            totCost += x_ilo[i][k]*inpMMCF.c[i];
    }
    model.add(IloMinimize(env, totCost));
}


void getCplexSolMMCF(InstanceMMCF inpMMCF, IloCplex cplex, SOLUTION & opt)
{
    opt.nOpen = 0;
    opt.ySol = new int[inpMMCF.nArcs];
    opt.xSol = new double*[inpMMCF.nArcs];
    for (int i = 0; i < inpMMCF.nArcs; i++)
        opt.xSol[i] = new double[inpMMCF.nK];

    opt.zStar = cplex.getObjValue();
    opt.zStatus = cplex.getStatus();

    for (int i = 0; i < inpMMCF.nArcs; i++)
        if (cplex.getValue(y_ilo[i]) >= 1.0-EPSI)
        {
            opt.ySol[i] = 1;
            opt.nOpen++;
        }
        else
            opt.ySol[i] = 0;

    for (int i = 0; i < inpMMCF.nArcs; i++)
        for (int j = 0; j < inpMMCF.nK; j++)
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
    fWriter << inpMMCF.nNodes << " " << inpMMCF.nArcs << " " << inpMMCF.nK << endl;
    fWriter << setprecision(15) << opt.zStar << endl;
    fWriter << opt.zStatus << endl;
    fWriter << opt.nOpen << endl;
    for (int i = 0; i < inpMMCF.nArcs; i++)
        if (opt.ySol[i] == 1)
            fWriter << " " << i;
    fWriter << endl;
    for (int i = 0; i < inpMMCF.nArcs; i++)
        for (int j = 0; j < inpMMCF.nK; j++)
            if (opt.xSol[i][j] > 0.0)
                fWriter << i << " " << j << " " << opt.xSol[i][j] << endl;

    fWriter.close();
    cout << "Solution written to disk. ('" << filename <<"')" << endl;

    // TEST: check that used capacity does not exceed capacity of each bin
}

