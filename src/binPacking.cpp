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
extern int readFromDisk;       //!< 0-No; (Generate a new Budget set B_l); 1-Yes

extern string versionLabel[6];
extern string supportLabel[2];
extern string instanceType;
// extern string versionType;
// extern string supportType;
extern int support;
extern int version;

extern IloNumVarArray y_ilo;
typedef IloArray <IloNumVarArray> TwoD;
extern TwoD x_ilo;
extern TwoD psi_ilo;
typedef std::vector <int> MyVect;


void read_parameters_budgetBIN()
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
    assert(pos == inp.nI*2);
}

void save_instance_2_diskBIN(double _epsilon, double _delta, double _gamma, int L, 
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

void define_budget_support(InstanceBin & inp, bool fromDisk)
{
   int    nBl = 0;
   int **Bl;
   double *budget;

   read_parameters_budgetBIN();

    nBl    = floor(_gamma*(double)inp.nI); // cardinality of each B_l
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
    for (int i = 0; i < inp.nI; ++i) 
        shuffled.push_back(i); // 0 2 3 ... nC-1


    if (readFromDisk==false)
    {
        // initialize random generator (uniform distribution)
        uniform_int_distribution<> d(0,inp.nI-1);
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
            save_instance_2_diskBIN(_epsilon, _delta, _gamma, L, nBl, Bl, budget);
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
    MyVect * mapping = new MyVect[inp.nI];
    for (int l = 0; l < L; l++)
        for (int k = 0; k < nBl; k++)
            mapping[Bl[l][k]].push_back(l);

   // total number of rows of W and h
   inp.nR = 2*inp.nI + L;
    // define vector h
    // a. box part
   inp.h  = new double[inp.nR];
   for (int j = 0; j < inp.nI; j++)
   {
       inp.h[j]        = -inp.d[j]*(1.0-_epsilon);
       // inp.h[j]        = 0.0;
       inp.h[j+inp.nI] =  inp.d[j]*(1.0+_epsilon);
   }

    
    // b. budget part
    for (int l = 0; l < L; l++)
        inp.h[2*inp.nI+l] = budget[l];

   // define matrix W in column major format
   int nEls  = 2*inp.nI + L*nBl;
   inp.W     = new int[nEls];
   inp.index = new int[nEls];
   inp.start = new int[inp.nI+1]; // one extra element in last position

    int pos = 0;
    for (int j = 0; j < inp.nI; j++)
    {
        inp.start[j]     = pos;
        inp.W[pos]       = -1;
        inp.index[pos++] = j;
        inp.W[pos]       = 1;
        inp.index[pos++] = j+inp.nI;
        if (mapping[j].size() > 0)
        {
            for (MyVect::iterator it = mapping[j].begin(); it != mapping[j].end(); it++)
            {
                inp.W[pos]       = 1;
                inp.index[pos++] = 2*inp.nI + *it;
            }
        }
    }
    inp.start[inp.nI] = pos;
    // NOTE : Remove ASSERT from final version
    assert(pos == nEls);
}


void define_POLY_BINPACKING(InstanceBin & inpBin, IloModel & model, IloCplex & cplex, int support)
{
    // select support type
    // here we get W and h, depending on the type of support
    switch (support) {
        case 1 :
            define_box_support(inpBin);
            cout << "BOX SUPPORT defined here!!!\n";
            break;
        case 2 :
            define_budget_support(inpBin, false);
            cout << "BUDGET SUPPORT not defined here!!!\n";
            break;
        default :
            cout << "ERROR : Support type not defined.\n" << endl;
            exit(123);
    }

    char varName[100];
    char constrName[100];
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

    // psi variables
    psi_ilo = TwoD(env, inpBin.nI);
    for (int j = 0; j < inpBin.nI; j++)
        psi_ilo[j] = IloNumVarArray(env, inpBin.nR, 0.0, IloInfinity, ILOFLOAT);

    // set vars names
    for (int j = 0; j < inpBin.nI; j++) 
        for (int t = 0; t < inpBin.nR; t++) {
        sprintf(varName, "psi.%d.%d", (int)j, (int)t);
        psi_ilo[j][t].setName(varName);
    }

    // packing constraints (each item in a bin)
    for (int i = 0; i < inpBin.nI; i++) {
        IloExpr sum(env);
        for (int j = 0; j < inpBin.nI; j++)
            sum += x_ilo[i][j];

        sprintf(constrName, "packing.%d", i);
        model.add(IloRange(env, 1.0, sum, 1.0, constrName));
    }

    // capacity constraint
    bool standard =false;
    if (standard) {
        for (int j = 0; j < inpBin.nI; j++) {
            IloExpr sum(env);
            for (int i = 0; i < inpBin.nI; i++)
                sum += (1.0 + _epsilon)*inpBin.d[i]*x_ilo[i][j]; // update this in robust formulation
            sum -= y_ilo[j]*inpBin.totS;

            sprintf(constrName, "capacity.%d", j);
            model.add(IloRange(env, -IloInfinity, sum, 0.0, constrName));
        }
    }
    else {
        
        // "robust" demand - constr. W*psi >= x
        // Rem. W is in column major format
        for (int j = 0; j < inpBin.nI; j++)
            for (int i = 0; i < inpBin.nI; i++)
            {
                IloExpr sum(env);
                for (int l = inpBin.start[i]; l < inpBin.start[i+1]; l++)
                {
                    int t = inpBin.index[l];
                    sum += inpBin.W[l]*psi_ilo[j][t];
                }
                sum -= x_ilo[i][j];

                sprintf(constrName, "robustDemand.%d%d", (int)j, (int)i);
                model.add(IloRange(env, 0.0, sum, IloInfinity, constrName));
            }

        // robust capacity h*psi <= s*y
        for (int j = 0; j < inpBin.nI; j++)
        {
            IloExpr sum(env);
            for (int t = 0; t < inpBin.nR; t++)
                sum += inpBin.h[t]*psi_ilo[j][t];

            sum -= inpBin.totS*y_ilo[j];

            sprintf(constrName, "robustCapacity.%d", (int)j);
            model.add(IloRange(env, -IloInfinity, sum, 0.0, constrName));
        }
    }


    // objective function: min sum y_j
    IloExpr totCost(env);
    for (int j = 0; j < inpBin.nI; j++)
        totCost += y_ilo[j];

    model.add(IloMinimize(env,totCost));
}

void define_POLY_BINPACKING_DETERM(InstanceBin & inpBin, IloModel & model, IloCplex & cplex, int support)
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
    char constrName[100];
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

        sprintf(constrName, "packing.%d", i);
        model.add(IloRange(env, 1.0, sum, 1.0, constrName));
    }

    // capacity constraint
    for (int j = 0; j < inpBin.nI; j++) {
        IloExpr sum(env);
        for (int i = 0; i < inpBin.nI; i++)
            sum += (1.+_epsilon)*inpBin.d[i]*x_ilo[i][j]; // update this in robust formulation
        sum -= y_ilo[j]*inpBin.totS;

        sprintf(constrName, "capacity.%d", j);
        model.add(IloRange(env, -IloInfinity, sum, 0.0, constrName));
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

    for (int j = 0; j < inpBin.nI; j++)
        if (cplex.getValue(y_ilo[j]) >= 1.0-EPSI)
        {
            opt.ySol[j] = 1;
            opt.nOpen++;
        }
        else
            opt.ySol[j] = 0;

    for (int i = 0; i < inpBin.nI; i++)
        for (int j = 0; j < inpBin.nI; j++)
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
    fWriter << inpBin.nI << endl;
    fWriter << setprecision(15) << opt.zStar << endl;
    fWriter << opt.zStatus << endl;
    fWriter << opt.nOpen << endl;
    for (int j = 0; j < inpBin.nI; j++)
        if (opt.ySol[j] == 1)
            fWriter << " " << j;
    fWriter << endl;
    for (int i = 0; i < inpBin.nI; i++)
        for (int j = 0; j < inpBin.nI; j++)
            if (opt.xSol[i][j] > 0.0)
                fWriter << i << " " << j << " " << opt.xSol[i][j] << endl;

    fWriter.close();
    cout << "Solution written to disk. ('" << filename <<"')" << endl;

    // TEST: check that used capacity does not exceed capacity of each bin
    for (int j = 0; j < inpBin.nI; j++)
        if (opt.ySol[j] == 1) {
            cout << "Assigned to bin " << j << " :: ";
            float tot = 0.0;
            for (int i = 0; i < inpBin.nI; i++) {
                if (opt.xSol[i][j] > 0.0)
                    cout << " " << i;
                    tot += opt.xSol[i][j]*(1.0+_epsilon)*inpBin.d[i];
            }
            cout << "\nBin " << j << " with capacity = " << inpBin.totS << " and consumption = " << tot << endl;
            assert(inpBin.totS >= tot);
        }
}


