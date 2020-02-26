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

extern int version;
extern string versionLabel[6];
extern string supportLabel[2];

extern IloNumVarArray y_ilo;
typedef IloArray <IloNumVarArray> TwoD;
extern TwoD x_ilo;
extern TwoD psi_ilo;
extern IloNumVarArray u_ilo;
extern IloNumVarArray delta_ilo;

typedef std::vector <int> MyVect;

void read_parameters_budgetMMCF()
{
    ifstream fReader("parameters/paramsBudget.txt", ios::in);
    if (!fReader)
    {
        cout << "Cannot open file 'paramsBudget.txt'." << endl;
        exit(111);
    }
    cout << "[** Reading parameters from file 'paramsBudget.txt'.]" << endl;
    string line;

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
void save_instance_2_diskMMCF(double _epsilon, double _delta, double _gamma, int L, 
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

void define_box_support(InstanceMMCF & inp)
{
   // define number of constraints in Wd <= h
   inp.nR = 2*inp.nK;
   inp.h  = new double[inp.nR];
   for (int k = 0; k < inp.nK; k++)
   {
       inp.h[k]        = -inp.d[k]*(1.0-_epsilon);
       inp.h[k+inp.nK] =  inp.d[k]*(1.0+_epsilon);
   }

   // define matrix W in column major format
   inp.W     = new int[inp.nK*2];
   inp.index = new int[inp.nK*2];
   inp.start = new int[inp.nK+1]; // one extra element in last position

    int pos = 0;
    for (int k = 0; k < inp.nK; k++)
    {
        inp.start[k]   = pos;
        inp.W[pos]     = -1;
        inp.index[pos] = k;
        pos++;
        inp.W[pos]     = 1;
        inp.index[pos] = k+inp.nK;
        pos++;
    }
    inp.start[inp.nK] = pos;
    assert(pos == inp.nK*2);

    cout << "Printing vector h " << endl;
    for (int k = 0; k < inp.nK; k++)
        cout << inp.h[k] << " " << inp.d[k] << " " << inp.h[k+inp.nK] << endl;
}

void define_budget_support(InstanceMMCF & inp, bool fromDisk)
{
    
    int    nBl     = 0;
    int  **Bl;
    double *budget;

    read_parameters_budgetMMCF();

    nBl    = floor(_gamma*(double)inp.nK); // cardinality of each B_l
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
    for (int i = 0; i < inp.nK; ++i) 
        shuffled.push_back(i); // 0 2 3 ... nC-1

    if (readFromDisk==false)
    {
        // initialize random generator (uniform distribution)
        uniform_int_distribution<> d(0,inp.nNodes-1);
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
        }
        // adjust b_l values
        for (int l = 0; l < L; l++)
            budget[l] = floor(_delta*budget[l]);

        // do we want to save the sets B_l (and the parameters?)
        bool save2Disk = true;
        if (save2Disk==true)
            save_instance_2_diskMMCF(_epsilon, _delta, _gamma, L, nBl, Bl, budget);
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
    MyVect * mapping = new MyVect[inp.nK];
    for (int l = 0; l < L; l++)
        for (int k = 0; k < nBl; k++)
            mapping[Bl[l][k]].push_back(l);

   // total number of rows of W and h
   inp.nR = 2*inp.nK + L;
    // define vector h
    // a. box part
   inp.h  = new double[inp.nR];
   for (int j = 0; j < inp.nK; j++)
   {
       inp.h[j]        = -inp.d[j]*(1.0-_epsilon);
       // inp.h[j]        = 0.0;
       inp.h[j+inp.nK] =  inp.d[j]*(1.0+_epsilon);
   }
    
    // b. budget part
    for (int l = 0; l < L; l++)
        inp.h[2*inp.nK+l] = budget[l];

   // define matrix W in column major format
   int nEls  = 2*inp.nK + L*nBl;
   inp.W     = new int[nEls];
   inp.index = new int[nEls];
   inp.start = new int[inp.nK+1]; // one extra element in last position

    int pos = 0;
    for (int j = 0; j < inp.nK; j++)
    {
        inp.start[j]     = pos;
        inp.W[pos]       = -1;
        inp.index[pos++] = j;
        inp.W[pos]       = 1;
        inp.index[pos++] = j+inp.nK;
        if (mapping[j].size() > 0)
        {
            for (MyVect::iterator it = mapping[j].begin(); it != mapping[j].end(); it++)
            {
                inp.W[pos]       = 1;
                inp.index[pos++] = 2*inp.nK + *it;
            }
        }
    }
    inp.start[inp.nK] = pos;
    // NOTE : Remove ASSERT from final version
    assert(pos == nEls);

}

void define_POLY_MMCF(InstanceMMCF & inpMMCF, IloModel & model, IloCplex & cplex, int support)
{
    // select support type
    // here we get W and h, depending on the type of support
    switch (support) {
        case 1 :
            define_box_support(inpMMCF);
            break;
        case 2 :
            define_budget_support(inpMMCF, false);
            break;
        default :
            cout << "ERROR : Support type not defined.\n" << endl;
            exit(123);
    }

    char varName[100];
    char constrName[100];
    IloEnv env = model.getEnv();


    y_ilo = IloNumVarArray(env, inpMMCF.nArcs, 0, 1, ILOINT);
    x_ilo = TwoD(env, inpMMCF.nArcs);
    for (int i = 0; i < inpMMCF.nArcs; i++)
        x_ilo[i] = IloNumVarArray(env, inpMMCF.nK, 0, 1.0, ILOFLOAT);

    // set vars names
    for (int i = 0; i < inpMMCF.nArcs; i++) {
        sprintf(varName, "y.%d", (int)i);
        y_ilo[i].setName(varName);
        for (int j = 0; j < inpMMCF.nK; j++) {
            sprintf(varName, "x.%d.%d", (int)i, (int)j);
            x_ilo[i][j].setName(varName);
        }
    }

    // psi variables
    psi_ilo = TwoD(env, inpMMCF.nArcs);
    for (int i = 0; i < inpMMCF.nArcs; i++)
        psi_ilo[i] = IloNumVarArray(env, inpMMCF.nR, 0.0, IloInfinity, ILOFLOAT);

    // u variables
    u_ilo = IloNumVarArray(env, inpMMCF.nR, 0.0, IloInfinity, ILOFLOAT);
    
    // set vars names
    for (int i = 0; i < inpMMCF.nArcs; i++)
        for (int t = 0; t < inpMMCF.nR; t++)
        {
            sprintf(varName, "psi.%d.%d", (int)i, (int) t);
            psi_ilo[i][t].setName(varName);
        }
    for (int t = 0; t < inpMMCF.nR; t++)
    {
        sprintf(varName, "u.%d", (int)t);
        u_ilo[t].setName(varName);
    }
    //
    // delta variable (for the objective function)
    delta_ilo = IloNumVarArray(env, 1, 0.0, IloInfinity, ILOFLOAT);
    delta_ilo[0].setName("delta.0");

    // flow constraints
    // flow at node n of commodity k: For each commodity k, the demand is zero
    // for all nodes, except for source and sink. Note that demands here are
    // for each commodity and each node, i.e., d_n^k.
    for (int k = 0; k < inpMMCF.nK; k++) {

        std::vector<double> d(inpMMCF.nNodes, 0.0);

        d[inpMMCF.source[k]] = -1.0;
        d[inpMMCF.sink[k]]   =  1.0;
        /* cout << " \nComm " << k << ":: ";
         * for (int n = 0; n < inpMMCF.nNodes; n++)
         *     cout << " " << d[n]; */

        for (int n = 0; n < inpMMCF.nNodes; n++) {
            IloExpr sum(env);
            // inbound elements
            for (int i = 0; i < inpMMCF.inbound[n].size(); i++) {
                int arc = inpMMCF.inbound[n][i];
                sum += x_ilo[arc][k];
            }
            
            // outbound elements
            for (int i = 0; i < inpMMCF.outbound[n].size(); i++) {
                int arc = inpMMCF.outbound[n][i];
                sum -= x_ilo[arc][k];
            }

            sprintf(constrName, "flow.%d.%d", k, n);
            model.add(IloRange(env, d[n], sum, IloInfinity, constrName));
        }
    }

    // robust demand W*psi >= x_ij^k
    for (int i = 0; i < inpMMCF.nArcs; i++) {
        for (int k = 0; k < inpMMCF.nK; k++) {
            IloExpr sum(env);

            for (int l = inpMMCF.start[k]; l < inpMMCF.start[k+1]; l++)
            {
                int t = inpMMCF.index[l];
                sum += inpMMCF.W[l]*psi_ilo[i][t];
            }
            sum -= x_ilo[i][k];

            sprintf(constrName, "robustDemand.%d.%d", (int)i, (int)k);
            model.add(IloRange(env, 0.0, sum, IloInfinity, constrName));
        }
    }


    // robust capacity h*psi <= s*y
    for (int i = 0; i < inpMMCF.nArcs; i++)
    {
        IloExpr sum(env);
        for (int t = 0; t < inpMMCF.nR; t++)
            sum += inpMMCF.h[t]*psi_ilo[i][t];

        sum -= inpMMCF.cap[i]*y_ilo[i];

        sprintf(constrName, "robustCapacity.%d", (int)i);
        model.add(IloRange(env, -IloInfinity, sum, 0.0, constrName));
    }

/*     // logical constraints
 *     for (int i = 0; i < inpMMCF.nArcs; i++) {
 *         IloExpr sum(env);
 *         for (int k = 0; k < inpMMCF.nK; k++)
 *             sum += x_ilo[i][k]*inpMMCF.d[k]*(1.0+_epsilon);
 *
 *         sum -= y_ilo[i]*inpMMCF.cap[i];
 *
 *         sprintf(constrName, "logical.%d", i);
 *         model.add(IloRange(env, -IloInfinity, sum, 0.0, constrName));
 *     } */
    
    // "robust" objective function: h*u <= delta
    IloExpr sum(env);
    for (int t = 0; t < inpMMCF.nR; t++)
        sum += inpMMCF.h[t]*u_ilo[t];

    sum -= delta_ilo[0];

    sprintf(constrName, "robustObj1");
    model.add(IloRange(env, -IloInfinity, sum, 0.0, constrName));
    
    // second robust obj function: W*u >= c*x
    for (int k = 0; k < inpMMCF.nK; k++)
    {
        IloExpr sum(env);
        for (int l = inpMMCF.start[k]; l < inpMMCF.start[k+1]; l++)
        {
            int t = inpMMCF.index[l];
            sum += inpMMCF.W[l]*u_ilo[t];
        }
        for (int i = 0; i < inpMMCF.nArcs; i++)
            sum -= inpMMCF.c[i]*x_ilo[i][k];

        sprintf(constrName, "robustObj2.%d", (int)k);
        model.add(IloRange(env, 0.0, sum, IloInfinity, constrName));
    }

    IloExpr totCost(env);
    for (int i = 0; i < inpMMCF.nArcs; i++) {
        totCost += y_ilo[i]*inpMMCF.fixed[i];
            // totCost += x_ilo[i][k]*inpMMCF.c[i]*inpMMCF.d[k]*(1.0+_epsilon);
    }
    totCost += delta_ilo[0];
    model.add(IloMinimize(env, totCost));
}


void define_POLY_MMCF_old(InstanceMMCF & inpMMCF, IloModel & model, IloCplex & cplex, int support)
{
    // select support type
    // here we get W and h, depending on the type of support
    switch (support) {
        case 1 :
            define_box_support(inpMMCF);
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
    char constrName[100];
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
    // flow at node n of commodity k: For each commodity k, the demand is zero
    // for all nodes, except for source and sink. Note that demands here are
    // for each commodity and each node, i.e., d_n^k.
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

            sprintf(constrName, "flow.%d%d", k, n);
            model.add(IloRange(env, d[n], sum, IloInfinity, constrName));
        }
    }

    // logical constraints
    for (int i = 0; i < inpMMCF.nArcs; i++) {
        IloExpr sum(env);
        for (int k = 0; k < inpMMCF.nK; k++)
            sum += x_ilo[i][k];

        sum -= y_ilo[i]*inpMMCF.cap[i];

        sprintf(constrName, "logical.%d", i);
        model.add(IloRange(env, -IloInfinity, sum, 0.0, constrName));
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

    string versionType = versionLabel[version-1];
    if (version > 3) 
        versionType += "-" + supportLabel[support-1];

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

