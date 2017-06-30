//===================================================================================
//
// Calculate the NACACO score as feature for MS
//
// AUTHORS: daiwentao
// E-mail: daiwentao@moon.ibp.ac.cn
// Date: 2012.06.18
//              
//
//===================================================================================

#ifndef NCACO_H
#define NCACO_H


#include <vector>
#include <fstream>
#include <string>
#include <cstring>
#include <iostream>
#include <cmath>
#include <cstdio>
#include <map>
#include <ctime>

using namespace std;


typedef vector<vector<double > > DV2;
typedef vector<vector<vector<vector<double > > > > DV4;
typedef vector<vector<vector<vector<vector<vector<double > > > > > > DV6;

extern DV2 VDWRADIUS;
extern DV2 CONTABLE;
extern DV4 TRPTABLE;
extern DV2 SOLVETABLE;
extern DV4 CENTABLE;
extern DV4 StruEnergySet;
//extern double G_EXTRA;

struct residueInfor
{
	string resid;
	int index;
	vector<double> CAcoordinate;
};

typedef struct _ChainCoord	//to store chain's five-beads model information
{
	int resinum;
	int atomnum;

	double beginxyz[3];
	double endxyz[6];
	
	vector<vector<string > > atomname;
	vector<string> resiname;
	vector<int> resiserial;
	vector<vector<vector<double > > > xyz;

} ChainCoord;

typedef struct _LocalFrag
{
	
	string rname[3];
	double angle[4];
	
} LocalFrag;

void ReadAllPar(DV2& vdwradtable, const string vdwfile, DV2& contable, const string confile, DV4& trptable, const string trpfile, DV2& solvetable, const string solvefile, DV4& centable, const string centroidINFOfile, DV4&struenergyset, const string stru_energy_file);
void ReadVdwRadius(int& INTERRESI, DV2& VDWRADIUS, const string vdwfile);
int ReadConPar(DV2& contable, const string parfile);
int ReadTrpPar(DV4& trptable, const string parfile);
int ReadCentroid(DV4& centable, const string cenfile);
int ReadSolvPar(DV2 &solvetable, const string solvefile);
void ReadStruEnergyTable(const string energy_file,DV4 &StruEnergySet);

int Read4AtomPdb(ChainCoord& chain, char* pdbfile);
void AddCentroid(ChainCoord& chain, DV4& centable);
double EconOfFiveBead99_99(ChainCoord& chain, DV2& VDWRADIUS, DV2& contable, int &collidetnum);
double Etrp(ChainCoord& chain, DV4& trptable);
double Esolve(ChainCoord& chain,const vector<vector<double> >& solvetable);
double Estru(ChainCoord& chain, const vector<vector<vector<vector<double > > > > &StruEnergySet);

int Chain2LocalFrag(vector<LocalFrag>& localfrag, ChainCoord& chain);
int Get4Angle(double angle[4], double xyz[12][3]);
bool CalContact(const double ContactDisCutoff,const vector<vector<double> > &CAcoordinate1,
				const vector<vector<double> > & CAcoordinate2,double &dis);
bool internal2cartesian (double *c1, double *c2, double *c3, double *p, double *c4);
void randomShuffle(vector<string>& array);
double Etrp_n(ChainCoord& chain, DV4& trptable);
void GenCENmatrix(ChainCoord& chain, vector<vector<double> >& Cmatrix);

#endif
