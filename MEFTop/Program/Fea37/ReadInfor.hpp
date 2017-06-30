//===================================================================================
//
// Read information from files
//
// AUTHORS: daiwentao
// E-mail: daiwentao@moon.ibp.ac.cn
// Date: 2012.01.04
//              
//
//===================================================================================


#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <list>
#include <iostream>
#include <cstring>
#include <vector>
#include <fstream>
#include "GenSSE.hpp"

using namespace std;

typedef struct _Quepre
{
	string seq;
	string ss;
	string rs;
	vector< vector<double> > CM8a;
	vector< vector<double> > CM12a;
	double radius;
	vector<ELE> sse;
} QUE;


typedef struct _Mod
{
	vector<vector< vector<double> > > coord; //Coord of main chain
	vector<vector< vector<double> > > acoord; //Coord of all atoms
	string seq;
	string ss;
	string rs;
	vector<int> Acc; 
	vector<vector<double> > CAmatrix;
	double radius;
	vector<ELE> sse;
}MOD;

int ReadID(vector<string>& allid, string IDfile);
void ReadSequence(const char* fastafile,string& sequence);
int ReadCM(string& cmfile,string& ss,string& rs, vector< vector<double> >& CMA, string& seq);
int Readsssa(string& cmfile,string& ss,string& rs);
int Readpdbcoord(vector<vector< vector<double> > >& coord, string name);
int ReadMCcoord(vector<vector< vector<double> > >& coord, string name);
int ReadDSSP(string& second, vector<int>& acc, string secondfile);
int GenCAmatrix(vector<vector<double> >& matrix, vector<vector< vector<double> > > xyz);
int SurfRes(string res, vector<int>& acc, string& SIs );
void CleQUE(QUE& querry );
double RGvalue(vector<vector< vector<double> > > xyz);
double PredRG(string sequence);
double ReadTM(string name);




