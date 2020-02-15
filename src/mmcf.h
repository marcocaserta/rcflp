
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

struct InstanceMMCF {
    int nNodes;
    int nArcs;
    int nK;
    std::vector<int> origin;
    std::vector<int> destin;
    std::vector<double> c;
    std::vector<double> d;
    std::vector<double> fixed;
    std::vector<double> cap;
    std::vector<int> source;
    std::vector<int> sink;
    std::vector<std::vector<int>> outbound;
    std::vector<std::vector<int>> inbound;
    // int ** outbound;
    // int ** inbound;
    int     nR;    //!< Number of constraints polyhedron uncertainty set
    double  *h;    //!< Rhs of polyhedron definining support
    int     *W;    //!< Matrix W in column major format
    int *index;    //!< Index of column major format for w
    int *start;    //!< Starting position for elements of column j
};

void define_POLY_MMCF(InstanceMMCF & inpMMCF, IloModel & model, IloCplex & cplex, int support);
void getCplexSolMMCF(InstanceMMCF inpMMCF, IloCplex cplex, SOLUTION & opt);
void define_box_support(InstanceBin & inp);
