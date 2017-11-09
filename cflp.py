"""
/***************************************************************************
 *   copyright (C) 2015 by Marco Caserta                                   *
 *   marco dot caserta at ie dot edu                                       *
 *                                                                         *
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

 Algorithm for the Multi-source Robust Capacitated Facility Location Problem
 (MS-CFLP-R).

 Author: Marco Caserta (marco dot caserta at ie dot edu)
 Started : 26.04.2015
 Ended   :

 Command line options (see parseCommandLine):
-i inputfile

-t instance type:
1 -> orLibrary
2 -> Avella
3 -> SS-orLibrary
4 -> Robust orLibrary

"""

import sys, getopt
import cplex
from cplex.callbacks import HeuristicCallback
import bisect
from time import time

_INFTY = sys.float_info.max

inputfile = ""
ftype     = ""

x_var = []
y_var = []
q_var = []
w_var = -1


initT = 0

omega = 1.0
epsi  = 0.001

_EPSI = sys.float_info.epsilon

# Parse command line
def parseCommandLine(argv):
    global inputfile
    global ftype
    
    try:
        opts, args = getopt.getopt(argv, "ht:i:", ["help","type=", "ifile="])
    except getopt.GetoptError:
        print ("Command Line Erorr. Usage : python cflp.py -t <type> -i\
        <inputfile> ")

        sys.exit(2)

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print ("Usage : python cflp.py -t <type> -i <inputfile> ")
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-t", "--type"):
            ftype = arg
            
class Instance:

    # initialize instance reading from "inputfile"
    def __init__(self, inputfile, ftype):
        
        self.f = []
        self.s = []
        self.c = []
        self.d = []

        with open(inputfile) as ff:

            if ftype == "1" or ftype == "3" or ftype == "4":  # orLibrary
                data = ff.readline()
                # read nr facilities and nr customers
                self.m, self.n = [int(v) for v in data.split()]


                for i in range(self.m):
                    data = ff.readline()
                    cap, cost = [float(v) for v in data.split()]
                    self.f.append(cost)
                    self.s.append(cap)
                    
                self.read_demands(ff)
                self.read_costs(ff)

            elif ftype == "2": 
                data = [float(v) for v in open(inputfile).read().split()]
                self.n = int(data[0]) # customers
                self.m = int(data[1]) # facilities
                progr = 2
                #print("m and n are ", self.m, self.n)
                
                # read demands
                self.d = [data[i] for i in range(progr,progr + self.n)]
                progr += self.n
                
                # read capacities
                self.s = [data[i] for i in range(progr, progr + self.m)]
                progr += self.m
                
                # read fixed costs
                self.f = [data[i] for i in range(progr, progr + self.m)]
                progr += self.m

                for i in range(self.m):
                    dd = [data[j] for j in range(progr, progr + self.n)]
                    progr += self.n
                    self.c.append(dd)
            else:
                print("Problem type not defined (-t option). Use '-h' for help.")
                exit()

            self.totS = sum(self.s)
            self.totD = sum(self.d)
            self.surplus = self.totS - self.totD
            assert(self.surplus > 0.0)
            
    def read_demands(self, ff):
        # read demands
        data = ff.readline()
        self.d = [float(v) for v in data.split()]

    def read_costs(self, ff):
        # read costs (from each facility to customers)
        for i in range(self.m):
            data = ff.readline()
            ci = [float(v) for v in data.split()]
            self.c.append(ci)

"""
Multi-source Capacitated Facility Location Problem.

This is the DETERMINISTIC version of the MS-CFLP.
"""
def defineModelCFLP(inp, cpx, ftype):

    '''
    This function slightly changes depending on whether the instance is from the
    OR Library (ftype = 1) or Avella (ftype = 2):
    - For the OR Library, x_ij is the % of the demand of customer j satisfied by
      plant i. This implies 0 <= x_ij <= 1
    - For Avella instances, x_ij is the number of units of demands satisfied.
      Therefore, we have: 0 <= x_ij <= d_j
    '''
    
    global x_var
    global y_var
    
    cpx.objective.set_sense(cpx.objective.sense.minimize)

    # define binary variables (y_var)
    cpx.variables.add(obj   = inp.f,
                      lb    = [0]*inp.m,
                      ub    = [1]*inp.m,
                      types = ["B"]*inp.m)
    
    # define continuous variables (x_var)
    if ftype == "1":
        for i in range(inp.m):
            cpx.variables.add(obj   = inp.c[i],
                              lb    = [0.0]*inp.n,
                              ub    = [1.0]*inp.n,
                              types = ["C"]*inp.n)
    elif ftype == "2":
        for i in range(inp.m):
            cpx.variables.add(obj   = inp.c[i],
                              lb    = [0.0]*inp.n,
                              ub    = inp.d,
                              types = ["C"]*inp.n)

    # create indexes
    y_var = range(inp.m)
    x_var = [range(inp.n*i+inp.m, inp.n*i+inp.n+inp.m) for i in range(inp.m)]
            
    # capacity constraints
    for i in range(inp.m):
        index = [x_var[i][j] for j in range(inp.n)]
        
        if ftype == "1":
            value = [inp.d[j] for j in range(inp.n)]
        elif ftype == "2":
            value = [1.0 for j in range(inp.n)]
            
        index.append(y_var[i])
        value.append(-inp.s[i])
        capacity_constraint = cplex.SparsePair(ind=index, val=value)
        cpx.linear_constraints.add(lin_expr = [capacity_constraint],
                                   senses = ["L"],
                                   rhs = [0.0])
                       
    # demand constraints
    if ftype == "1":
        for j in range(inp.n):
            index = [x_var[i][j] for i in range(inp.m)]
            value = [1.0]*inp.m
            demand_constraint = cplex.SparsePair(ind=index, val=value)
            cpx.linear_constraints.add(lin_expr = [demand_constraint],
                                       senses   = ["E"],
                                       rhs      = [1.0])
    elif ftype == "2":
        for j in range(inp.n):
            index = [x_var[i][j] for i in range(inp.m)]
            value = [1.0]*inp.m
            demand_constraint = cplex.SparsePair(ind=index, val=value)
            cpx.linear_constraints.add(lin_expr = [demand_constraint],
                                       senses   = ["E"],
                                       rhs      = [inp.d[j]] )
    #  # link between x and y
    #  for i in range(inp.m):
    #      for j in range(inp.n):
    #          index = [x_var[i][j], y_var[i]]
    #          value = [1.0, -1.0]
    #          link_constraint = cplex.SparsePair(ind=index, val=value)
    #          cpx.linear_constraints.add(lin_expr = [link_constraint],
    #                                     senses   = ["L"],
    #                                     rhs      = [0.0])

    if ftype == "2":
        # variable upper bound (B) - v2
        for i in range(inp.m):
            for j in range(inp.n):
                index = [x_var[i][j]]
                value = [1.0]
                B_constraint = cplex.SparsePair(ind=index, val=value)
                cpx.linear_constraints.add(lin_expr = [B_constraint],
                                           senses   = ["L"],
                                           rhs      = [inp.d[j]])

    

"""
Robust Formulation for the MS-CFLP.
"""
def defineModelRCFLP(inp, cpx, ftype):
    '''
    Define model for robust version. The key idea is taken from Baron et al. (2010) and
    adapted for this version of the CFLP.
    
    Note: for the OR Library problems, the quadratic constraint W >= ... 
    must be modified w.r.t. the version presented in the paper. The problem is due to
    the fact that the c_ij values provided in the OR Library are not expressed in
    dollars/unit, as one might expect given the formulation but, rather, in dollars per
    percentage point.
    
    This is easily seen observing the objective function. Since x_ij measures the
    percentage of demand of j satisfied by i, the objective function term that
    accounts for the transportation cost should be c_ij*x_ij*d_j. However, we use
    c_ij*x_ij (and we get the right solution values, as published in the literature).
    This proves that c_ij is the cost of delivering the entire demand of j from i.
    
    Following the same reasoning, the constraint:
    W >= sqrt(sum_i,j c_ij^2*x_ij^2*epsilon^2*d_j^2)
    must be rewritten as:
    W >= sqrt(sum_i,j c_ij^2*x_ij^2*epsilon^2)
    since c_ij = c_ij*d_j
    
    Note that the same does not apply to the Q_i constraints. Since the cost does not
    appear in the constraints, the d_j term must be included.
    '''
    
    global x_var
    global y_var
    global q_var
    global w_var
    
    cpx.objective.set_sense(cpx.objective.sense.minimize)

    # define binary variables (y_var)
    cpx.variables.add(obj   = inp.f,
                      lb    = [0]*inp.m,
                      ub    = [1]*inp.m,
                      types = ["B"]*inp.m)
    
    # define continuous variables (x_var)
    for i in range(inp.m):
        cpx.variables.add(obj   = inp.c[i],
                          lb    = [0.0]*inp.n,
                          ub    = [1.0]*inp.n,
                          types = ["C"]*inp.n)

    # define binary variables (q_var)
    cpx.variables.add(obj   = [0.0]*inp.m,
                      lb    = [0.0]*inp.m,
                      types = ["C"]*inp.m)

    # define binary variables (w_var)
    cpx.variables.add(obj   = [omega],
                      lb    = [0.0],
                      types = ["C"])

    # create indexes
    y_var = range(inp.m)
    x_var = [range(inp.n*i+inp.m, inp.n*i+inp.n+inp.m) for i in range(inp.m)]
    q_var = range(inp.m + inp.n*inp.m, inp.m + inp.n*inp.m + inp.m)
    w_var = inp.m + inp.n*inp.m + inp.m


    # second order cone W (modified for the OR Library instances)
    index = [x_var[i][j] for j in range(inp.n) for i in range(inp.m)]
    #value = [(inp.c[i][j]*epsi*inp.d[j])**2 for j in range(inp.n) for i in range(inp.m)]
    value = [(inp.c[i][j]*epsi)**2 for j in range(inp.n) for i in range(inp.m)]
    index.append(w_var)
    value.append(-1.0)
    robust_cost = cplex.SparseTriple(ind1 = index, ind2 = index, val = value)
    cpx.quadratic_constraints.add(quad_expr = robust_cost,
                                  sense     = "L",
                                  rhs       = 0.0)

    # Q conic constraint
    for i in range(inp.m):
        index = [x_var[i][j] for j in range(inp.n)]
        value = [(epsi*inp.d[j])**2 for j in range(inp.n)]
        #value = [(epsi)**2 for j in range(inp.n)]
        index.append(q_var[i])
        value.append(-1.0)

        robust_capacity = cplex.SparseTriple(ind1= index, ind2 = index, val = value)
        cpx.quadratic_constraints.add(quad_expr = robust_capacity,
                                      sense     = "L",
                                      rhs       = 0.0)
    
    # capacity constraints QUADRATIC
    for i in range(inp.m):
        index = [x_var[i][j] for j in range(inp.n)]
        value = [inp.d[j]    for j in range(inp.n)]
            
        index.append(q_var[i])
        value.append(omega)

        index.append(y_var[i])
        value.append(-inp.s[i])
        
        capacity_constraint = cplex.SparsePair(ind=index, val=value)
        cpx.linear_constraints.add(lin_expr = [capacity_constraint],
                                   senses   = ["L"],
                                   rhs      = [0.0])
                       
    # demand constraints
    for j in range(inp.n):
        index = [x_var[i][j] for i in range(inp.m)]
        value = [1.0]*inp.m
        demand_constraint = cplex.SparsePair(ind=index, val=value)
        cpx.linear_constraints.add(lin_expr = [demand_constraint],
                                   senses   = ["E"],
                                   rhs      = [1.0])

    


def solveModel(inp, cpx, timeLimit=1000, solLimit=9999, withPrinting=True, display=2):
    try:

        cpx.parameters.mip.interval.set(2500) # how often to print info
        cpx.parameters.timelimit.set(timeLimit)
        cpx.parameters.mip.limits.solutions.set(solLimit)
        cpx.parameters.mip.display.set(display)
        #  cpx.parameters.mip.tolerances.mipgap.set(0.000000001)
        #  cpx.parameters.mip.tolerances.absmipgap.set(0.000000001)
        #cpx.parameters.mip.tolerances.integrality.set(0.01)

        
        cpx.solve()
        
    except CplexSolverError as e:
        print("Exception raised during solve: " + e)
    else:
        # get solution
        solution = cpx.solution
        if withPrinting:
            print ("\t ub[{0:4d}] = {1:15.4f} with solution status = {2:20s}".\
                   format(0,solution.get_objective_value(), solution.status[solution.get_status()]))
        return solution
        
def printSummary(inputfile, inp):
    print("\n\nmarco caserta (c) 2015 ")
    print("==========================================")
    print("Namefile       : {0:>25s}".format(inputfile))
    print("Nr. Facilities : {0:25d}".format(inp.m))
    print("Nr. Customers  : {0:25d}".format(inp.n))

    print("Instance Type  : {0:>25s}".format(ftype))
    print("==========================================")
    
def main(argv):
    '''
    Entry point.
    '''

    parseCommandLine(argv)
    
    inp = Instance(inputfile, ftype)
    printSummary(inputfile, inp)
        
    cpx = cplex.Cplex()
    if ftype == "1" or ftype == "2":
        defineModelCFLP(inp, cpx, ftype)
    #  elif ftype == "3": # removed from the code (consider Multi-source only)
    #      defineModelSSCFLP(inp, cpx, ftype)
    elif ftype == "4":
        defineModelRCFLP(inp, cpx, ftype)

    
    # solve full model using CPLEX
    initT = time()

    sol  = solveModel(inp, cpx, timeLimit=10000)
    zOpt = sol.get_objective_value()
    stat = sol.status[sol.get_status()]
    lb   = sol.MIP.get_best_objective()
    gap  = sol.MIP.get_mip_relative_gap()
    zTime = time() - initT
    
    #  printing summary statistics and saving solution
    print("Summary Instance :: {0:20s} {1:20.5f} {2:25s} {3:20.5f} {4:20.7f} {5:20.7f}\n".
          format(inputfile, zOpt, stat, lb, gap, zTime))
    with open("cplexResults.txt", "a") as outfile:
        outfile.write("{0:20s} {1:20.5f} {2:25s} {3:20.5f} {4:20.7f} {5:20.7f}\n".
                      format(inputfile, zOpt, stat, lb, gap, zTime))

    if ftype == "4":
        with open("optSol.txt", "a") as outfile2:
            ySol = [i for i in range(inp.m) if sol.get_values(y_var[i]) > _EPSI]
            
            for i in ySol:
                units = [sol.get_values(x_var[i][j]) for j in range(inp.n) if sol.get_values(x_var[i][j]) > _EPSI]
                aux   = [j for j in range(inp.n) if sol.get_values(x_var[i][j]) > _EPSI]
                
                outfile2.write("f({0:3d}) :: {1} -> {2}\n".format(i, aux, len(aux)))
                outfile2.write("f({0:3d}) :: {1}\n".format(i, units))

            for j in range(inp.n):
                yy = [i for i in range(inp.m) if sol.get_values(x_var[i][j]) > 0.0000001]
                outfile2.write("c({0:3d}) = {1}".format(j,yy))
                tot = 0.0
                for i in yy:
                    tot += sol.get_values(x_var[i][j])
                aux = [sol.get_values(x_var[i][j]) for i in yy]
                outfile2.write("\t {0} \t :: {1:10.2f} vs {2:10.2f}\n".format(aux, tot, inp.d[j]))

            # recomputing objective function value
            xSol  = []
            xSolR = []    
            for i in range(inp.m):     
                aux = [sol.get_values(x_var[i][j]) for j in range(inp.n)]
                auxR = [round(el, 5) for el in aux]
                xSol.append(aux)
                xSolR.append(auxR)

            fixed  = sum([inp.f[i] for i in ySol])
            trans  = sum([inp.c[i][j]*xSol[i][j] for i in range(inp.m) for j in range(inp.n)])
            transR = sum([inp.c[i][j]*xSolR[i][j] for i in range(inp.m) for j in range(inp.n)])
            cone   = omega*sol.get_values(w_var)
            print("Total COST = {0} + {1} + {2} = {3}".format(fixed, trans, cone, fixed + trans + cone))
            print("Total COST = {0} + {1} + {2} = {3}".format(fixed, transR, cone, fixed + transR + cone))

                
            
    print("Solution written to disk file * results.txt * ")
    
if __name__ == '__main__':
    main(sys.argv[1:])
    #unittest.main()
    

    
