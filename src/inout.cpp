#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <fstream>

using namespace std;

struct INSTANCE {
    int nF;
    int nC;
    int nR;
    double  *f;
    double  *s;
    double  *d;
    double **c;
    double   totS;
    double   totD;
};

extern string instanceType;
extern string versionType;

/// Read benchmark instances
/**
 * Currently, two types of instances can be imported:
 * type 1: OR Library
 * type 2: Avella (Test Bed 1, Test Bed A. Test Bed B)
 */
int readProblemData(char * _FILENAME, int fType, INSTANCE & inp)
{
    inp.totS = 0.0;
    inp.totD = 0.0;
    ifstream fReader(_FILENAME, ios::in);
    if (!fReader)
    {
        cout << "cannot open file " << _FILENAME << endl;
        exit(1);
    }
    // read OR Library instances
    if (fType == 1)
    {
        
        fReader >> inp.nF >> inp.nC;
        inp.s = new double[inp.nF];
        inp.f = new double[inp.nF];
        inp.d = new double[inp.nC];
        inp.c = new double*[inp.nF];
        for (int i = 0; i < inp.nF; i++)
            inp.c[i] = new double[inp.nC];

        for (int i = 0; i < inp.nF; i++)
        {
            fReader >> inp.s[i] >> inp.f[i];
            inp.totS += inp.s[i];
        }

        for (int j = 0; j < inp.nC; j++)
        {
            fReader >> inp.d[j];
            inp.totD += inp.d[j];
        }

        for (int i = 0; i < inp.nF; i++)
            for (int j = 0; j < inp.nC; j++)
            {
                fReader >> inp.c[i][j];
                inp.c[i][j] /= inp.d[j];
            }
    }
    // read Avella instances
    else if (fType == 2)
    {
        fReader >> inp.nC >> inp.nF;
        inp.f = new double[inp.nF];
        inp.s = new double[inp.nF];
        inp.d = new double[inp.nC];
        inp.c = new double*[inp.nF];
        for (int i = 0; i < inp.nF; i++)
            inp.c[i] = new double[inp.nC];
        for (int j = 0; j < inp.nC; j++)
        {
            fReader >> inp.d[j];
            inp.totD += inp.d[j];
        }
        for (int i = 0; i < inp.nF; i++)
        {
            fReader >> inp.s[i];
            inp.totS += inp.s[i];
        }
        for (int i = 0; i < inp.nF; i++)
            fReader >> inp.f[i];

        for (int i = 0; i < inp.nF; i++)
            for (int j = 0; j < inp.nC; j++)
                fReader >> inp.c[i][j];
    }
    else
        cout << "Problem type not defined (-t option). Use '-h' for help. " << endl;

    fReader.close();

    return 1;
}

/// Print instance info and algorithmic parameters.
void printOptions(char * _FILENAME, int fType, int version, INSTANCE inp, int timeLimit)
{
   cout << "-------------------------------------" << endl;
   cout << "- OPTIONS : " << endl;
   cout << "-------------------------------------" << endl;
   cout << "  DATA FILE      = " << _FILENAME        << endl;
   cout << "  Instance type  = " << instanceType << endl;
   cout << "  Version        = " << versionType << endl;
   cout << "  Time Limit     = " << timeLimit << endl;
   cout << "  Nr. Facilities = " << inp.nF << endl;
   cout << "  Nr. Customers  = " << inp.nC << endl;
   cout << "-------------------------------------" <<  endl << endl;   
}

