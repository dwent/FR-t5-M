//===========================================================================================================
//Edited by Wu Aiping (Email: wuaiping@moon.ibp.ac.cn)
//Energy function based residue-five-beads model (N,CA,C,O, and centroid of side-chain)
//Energy Function: V = Econ + a*Etrp + b*Ehp
//	1)Econ : atom-pairs contact MIU-potential;
//	2)Etrp : local 3-residues' sequence-dependent conformation MIU-potential;
//	3)Ehp  : hydrogen-bonding potential, realized by MIU-poteintial or Quasi-chemical
//			approximation;
//	4)a, b and c in Ehp: three weigth factors to balance different energy sets; 
//===========================================================================================================

#ifndef FIVEBEADENERGY_H
#define FIVEBEADENERGY_H


#include <vector>
#include <fstream>
#include <string>
#include <iostream>
#include <cmath>
#include <cstdio>
#include "MathTool.hpp"
#include "IOprotein.hpp"

using namespace std;


//to store chain's five-beads model information
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


//VECTOR TYPEDEFINE
typedef vector<vector<double > > DV2;
typedef vector<vector<vector<vector<double > > > > DV4;
typedef vector<vector<vector<vector<vector<vector<double > > > > > > DV6;
//PARAMETER TABLES
extern DV2 VDWRADIUS;
extern DV2 CONTABLE;
extern DV4 TRPTABLE;
extern DV6 HBTABLE;
extern DV4 CENTABLE;

//MOST FUNCTIONS CAN BE FOUND IN THE CPP FILE
//AND THEIR FUNCTIONS ARE EXPLAINED THERE

void ReadAllPar(DV2& vdwradtable, const char* vdwfile, DV2& contable, const char* confile, DV4& trptable, const char* trpfile, DV6& hbtable, const char* hbfile, DV4& centable, const char* centroidINFOfile);
void ReadVdwRadius(int& INTERRESI, DV2& VDWRADIUS, const char* vdwfile);
int ReadConPar(DV2& contable, const char* parfile);
int ReadTrpPar(DV4& trptable, const char* parfile);
int ReadHbPar(DV6& hbtable, const char* parfile);
int ReadCentroid(DV4& centable, const char* cenfile);


//calculate Econ
int Read4AtomPdb(ChainCoord& chain, const char* pdbfile);
std::string ReadHoriz(char* ssfile);
double EconOfFiveBead99_99(ChainCoord& chain, DV2& VDWRADIUS, DV2& contable);

int ReadPdb(ChainCoord& chain, const char* pdbfile);
int Aname2Index(string rname, string aname);
double EconOfAllAtom84_84(ChainCoord& chain, int INTERRESI, DV2& VDWRADIUS, DV2& contable);


//calculate Etrp
typedef struct _LocalFrag
{
	
	string rname[3];
	double angle[4];
	
} LocalFrag;

int Get4Angle(double angle[4], double xyz[12][3]);
int Chain2LocalFrag(vector<LocalFrag>& localfrag, ChainCoord& chain);
double Etrp(ChainCoord& chain, DV4& trptable);


//calculate Ehp
typedef struct _Frag6Vector
{
	int rindex;
	
	double b1[2][3];
	double b2[2][3];
	double b3[2][3];
	
	double P1[2][3];
	double P2[2][3];
	double P3[2][3];
	
	double Oatom[3];
	double Hatom[3];
		
} Frag6Vector;

int Get6Angle(Frag6Vector& frag1, Frag6Vector& frag2, double angle[6]);
double VectorAngle(double v1[2][3], double v2[2][3]);
void Resi2Vector(double bv[2][3], double Pv[2][3], double resi[4][3]);
void OxyHyDis(double d[2], Frag6Vector& frag1, Frag6Vector& frag2);
void Get6Vector_2resi(Frag6Vector& frag6vector, double xyz[12][3], int rindex);
int Chain2FragVector_2resi(vector<Frag6Vector>& chain6v, ChainCoord& chain);
double Ehb_2resi(ChainCoord& chain, DV6& hbtable);


//add centroid atom to 4-atom's residue (N,CA,C,O)
void AddCentroid(ChainCoord& chain, DV2& tolerance, DV4& centable);
//bool internal2cartesian (double *c1, double *c2, double *c3, double *p, double *c4);


#endif
