"""
/***************************************************************************
 *   copyright (C) 2019 by Marco Caserta                                   *
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

 Analysis of the results on the robust CFLP.

 this code carries out the following steps:
 - read the output files of each evaluation
 - import them into a pandas dataframe and modify the dataframe to extract the
   parameters values
- create charts. The following charts can be created:
  * -t e : A chart which explores the effect of epsilon and scenario_epsi onto
           the infeasibility levels (either the counting, the sum, or the max)
  * -t g : A plot to explore the effects of the gamma parameter (the percentage
           of customers included in the budget constraint.)

  Other parameters of the command line are:
  * -i 'folder' : The folder containing the files which are the result of the
                  evaluation part
  * -s 'namefile' : The name of the chart file. By default, the chart is called
                    'plot.png'. Change this if we want to create multiple
                    charts.

 Author: Marco Caserta (marco dot caserta at ie dot edu)
 Started : 26.04.2015
 Ended   :

"""

import numpy as np
import pandas as pd
import os, sys, getopt
import glob
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

base = '/home/marco/gdrive/research/rcflp/analysis'
inputdir = ''
savedFile = 'plot.png'
chartType = ''

# Parse command line
def parseCommandLine(argv):
    """
    Note that two parameters need to be mandatorily defined, i.e., -i and -t.
    """
    global inputdir
    global savedFile
    global chartType
    
    try:
        opts, args = getopt.getopt(argv, "hs:t:i:", ["help", "sFile", "type", "idir="])
    except getopt.GetoptError:
        print ("Command Line Erorr. Usage : python cflp.py -i <inputdir> -t <chartType> -s <savedFile>")

        sys.exit(2)

    setFolder = False 
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print ("Usage : python cflp.py -i <inputdir> -t <chartType> -s <savedFile>")
            print("Chart Type :: ")
            print("\t 'e' : multiple lines over scenario_epsi")
            print("\t 'g' : multiple lines over gamma values")
            sys.exit()
        elif opt in ("-i", "--idir"):
            inputdir = arg
            setFolder = True
        elif opt in ("-t", "--type"):
            chartType = arg
            setChart = True
        elif opt in ("-s", "--sFile"):
            savedFile = arg

    if not setFolder:
        print("Folder not defined. Usage : python cflp.py -i <inputdir>. Use -h for help.")
        sys.exit(2)
    if not setChart:
        print("Chart type not defined. Usage : python cflp.py -t <chartType>. Use -h for help.")
        sys.exit(2)


def createPlotsGridEpsi(dfE):
    """
    This module creates a plot for the scenario_epsi parameter. Each line on the
    chart will give the infeasibility (y axis) while changing epsilon (x axis)
    for a given value of scenario_epsi (in 0.1, ...., 0.9).
    """

    L_vals = dfE['L'].unique()
    L_vals = np.sort(L_vals)

    plt.figure(figsize = (20,10))
    gs1 = gridspec.GridSpec(2, 3)
    gs1.update(wspace=0.25, hspace=0.25)

    ax = []

    i = 0
    j = 0
    count = 0

    for l in L_vals:
        ax.append(plt.subplot(gs1[i,j]))
        j += 1
        if j % 3 == 0:
            j = 0
            i += 1

        support       = 'budget'
        delta         = 1.1
        gamma         = 0.3
        L             = l

        ddfE = dfE[ (dfE['support_type']  == support) & 
                    (dfE['delta'] == delta) &
                    (dfE['gamma'] == gamma) &
                    (dfE['L']     == L)        
                ]
        
        
        #ddf['sum'] = ddf['sum']*(-1)

        # this is used to plot both the deterministic and the robust obj
        # function
        #  ddfE['scenario_epsi'] = ['900-D']*len(ddfE)
        #  ddf['scenario_epsi'] = ['900-R']*len(ddf)
        
        ddfE.groupby(by=['epsi', 'scenario_epsi'])['infeasible'].sum().unstack().plot(ax=ax[count])

        title = "Budget L = " + str(L)
        plt.title(title);
        count += 1

    plt.savefig(savedFile)
    print("Plot saved to disk - '{0}'".format(savedFile))


def createPlotsGridG(dfE):
    """
    This module creates a plot for the gamma parameter. Each line on the
    chart will give the infeasibility (y axis) while changing epsilon (x axis)
    for a given value of gamma (in 0.1, ...., 0.5).
    """

    L_vals = dfE['L'].unique()
    L_vals = np.sort(L_vals)

    plt.figure(figsize = (20,10))
    gs1 = gridspec.GridSpec(2, 3)
    gs1.update(wspace=0.25, hspace=0.25)

    ax = []

    i = 0
    j = 0
    count = 0

    for l in L_vals:
        ax.append(plt.subplot(gs1[i,j]))
        j += 1
        if j % 3 == 0:
            j = 0
            i += 1

        support       = 'budget'
        delta         = 1.1
        L             = l

        ddfE = dfE[ (dfE['support_type']  == support) & 
                    (dfE['delta'] == delta) &
                    (dfE['L']     == L)        
                ]
        
        ddfE.groupby(by=['epsi', 'gamma'])['infeasible'].sum().unstack().plot( ax=ax[count])

        title = "Budget L = " + str(L)
        plt.title(title);
        count += 1

    plt.savefig(savedFile)
    print("Plot saved to disk - '{0}'".format(savedFile))


def createDataFrame():
    """
    Read the files contained in the target directory. All the observations are 
    stored inside a dataframe.
    """

    listOfFiles = glob.glob(os.path.join(base,inputdir, '*.*'))
    i = 0
    totFiles = len(listOfFiles)
    for file in listOfFiles:
        aux = pd.read_csv(file, header=None, sep=';')
        if i%1000 == 0:
            print("[{0:6d}/{1:6d}] '{2}' \t size :: {3}".format(i+1,
            totFiles, file, aux.shape))
        if i == 0:
            dfE = aux
        else:
            dfE = dfE.append(aux,ignore_index=True, sort=False)
       
        i += 1

    print("\n")
    print("-"*80)
    print("Total Nr. of Files :: ", i )
    print("DataFrame Size     :: ", dfE.shape)
    print("-"*80)

    return dfE


def transformDataFrame(dfE):
    """
    Add new columns and extract parameter values from the name of the file.
    NOTE: This part must be modified if the 'm' part of the name (discrete or
    robust objective functions) is or is not present.
    """

    # add columns names
    dfE.columns=['scenario', 'solution', 'fixed', 'variable', 'sum', 'max']

    # create new columns
    dfE['infeasible'] = np.where(dfE['max'] > 0.001, 1, 0)
    dfE['z']          = dfE['fixed'] + dfE['variable']

    # extract parameters
    temp = dfE['solution'].str.split('-', n = 7, expand = True)
    dfE['support_type']  = temp[2]
    dfE['epsi']          = temp[3]
    dfE['delta']         = temp[4]
    dfE['gamma']         = temp[5]
    dfE['L']             = temp[6]
    dfE["epsi"]          = dfE.epsi.astype(float)
    dfE["delta"]         = dfE.delta.astype(float)
    dfE["gamma"]         = dfE.gamma.astype(float)
    dfE["L"]             = dfE.L.astype(int)

    temp                 = dfE['scenario'].str.split('_', n = 6, expand = True)
    dfE['scenario_epsi'] = temp[1]
    dfE['scenario_epsi'] = dfE.scenario_epsi.astype(int)

    return dfE

def main(argv):
    '''
    Entry point.
    '''

    parseCommandLine(argv)

    print("*"*50)
    print("* Parameters ")
    print("*"*50)
    print("* Folder with files : ", inputdir)
    print("* Type of chart     : ", chartType)
    print("* Plot name         : ", savedFile)
    print("*"*50)
    print("\n")

    dfE = createDataFrame()
    dfE = transformDataFrame(dfE)

    if chartType == 'e':
        createPlotsGridEpsi(dfE)
    elif chartType == 'g':
        createPlotsGridG(dfE)


    
if __name__ == '__main__':
    main(sys.argv[1:])
    #unittest.main()

