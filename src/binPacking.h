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

/*! \file options.h
\brief Header file of options.cpp

*/
#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

extern char* _FILENAME;
extern int time_limit;

struct InstanceBin {

    double totS; // total capacity
    int nI;   // number of items
    int opt;  // optimal value 
    std::vector<double> d; // weight of each item (treated as demand here)

    int     nR;    //!< Number of constraints polyhedron uncertainty set
    double  *h;    //!< Rhs of polyhedron definining support
    int     *W;    //!< Matrix W in column major format
    int *index;    //!< Index of column major format for w
    int *start;    //!< Starting position for elements of column j
};
struct SOLUTION {
    int nOpen;
    int  *ySol;
    double **xSol;
    double zStar;
    IloAlgorithm::Status zStatus;
    IloNum startTime;
    IloNum cpuTime;
};

int parseOptions(int argc, char* argv[]);

void define_POLY_BINPACKING(InstanceBin & inpBin, IloModel & model, IloCplex & cplex, int support);
void getCplexSolBIN(InstanceBin inpBin, IloCplex cplex, SOLUTION & opt);
void define_box_support(InstanceBin & inp);


