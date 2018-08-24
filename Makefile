#####################################################
#
# This makefile is used to compile the 
# sources (c++) that implements the
# algorithm for the Robust Facility Location Problem
# 
#####################################################
SRCDIR        = src
BINDIR        = bin
OBJDIR        = obj

# ---------------------------------------------------------------------
# ILOG CPLEX stuff
SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic
# ---------------------------------------------------------------------
CONCERTDIR    = /home/marco/opt/ILOG/cplex127/concert
CPLEXDIR      = /home/marco/opt/ILOG/cplex127/cplex
CPLEXBINDIR   = $(CPLEXDIR)/bin/$(BINDIST)
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include

# BE CAREFULL: All warnings have been disables!!!!
CCOPT     = -m64 -O -w -fPIC -fexceptions  -DIL_STD ## -DOPTIMAL -DNDEBUG
CCLNFLAGS = -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert -lm -pthread 
CCFLAGS   = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) 


# set default name for the executable
EXEC ?= rcflp

# set compiler
CC = g++

# set debug options
DEBUG = -ggdb

#set optimization level
OPTLEVEL = -O -DEBIAN_BUILDARCH=pentium

#set flags
FLAGS =  -fomit-frame-pointer -pipe -Wimplicit -Wparentheses -Wreturn-type -Wcast-qual -Wpointer-arith -Wwrite-strings
#-Wconversion

default: $(OBJDIR)/rcflp.o $(OBJDIR)/options.o $(OBJDIR)/inout.o 
	$(CC) $(CCFLAGS) $(OBJDIR)/options.o $(OBJDIR)/inout.o $(OBJDIR)/rcflp.o -o $(BINDIR)/$(EXEC) $(CCLNFLAGS) 
$(OBJDIR)/rcflp.o: $(SRCDIR)/rcflp.cpp $(SRCDIR)/inout.cpp 
	$(CC) -c -Wignored-attributes $(CCFLAGS) $(SRCDIR)/rcflp.cpp -o $(OBJDIR)/rcflp.o
$(OBJDIR)/options.o: $(SRCDIR)/options.cpp
	$(CC) -c $(FLAGS) $(SRCDIR)/options.cpp -o $(OBJDIR)/options.o
$(OBJDIR)/inout.o: $(SRCDIR)/inout.cpp
	$(CC) -c $(FLAGS) $(SRCDIR)/inout.cpp -o $(OBJDIR)/inout.o
# $(OBJDIR)/timer.o: $(SRCDIR)/timer.cpp
#     $(CC) -c $(FLAGS) $(SRCDIR)/timer.cpp -o $(OBJDIR)/timer.o

# clean backup files
clean:
	rm *.*~ 
	rm *.o

# create doxygen documentation using "doxygen.conf" file
# the documentation is put into the directory Doc	
doc: $(SRCDIR)/rcflp.cpp doxygen.conf
	doxygen doxygen.conf


#	
