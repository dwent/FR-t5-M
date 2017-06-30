				/*********************************************************
				Copyright	    	 : Edited by Wuaiping, 2006
				Program name	 	 : PTADEE
				Program version	 : 0.2
				Begin time		 : Dec. 28 2005
				Refined time 	 : Jun. 14 2007
				Email    		 : wuaiping@moon.ibp.ac.cn
				*********************************************************/
				
#ifndef MATHTOOL_H
#define MATHTOOL_H

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <string>
#include <iostream>
#include "BasicDef.hpp"


using namespace std;

/////////////////////changed from caoyang//////////////
typedef vector<vector<double> > mymatrix;
typedef vector<double> myvector;

inline void SetMatrix(mymatrix &sm, int m, int n);
inline void MatrixTimesMatrix(const mymatrix &a, const mymatrix &v, mymatrix &x, int m, int n, int l);
inline bool TransVectorTimesVector(const myvector &trans, const myvector &vctor, mymatrix &mtx);
inline bool MatrixTimesTransVector(const mymatrix &mtx, const myvector &tvt, myvector &vct);
inline void RealTimesMatrix(double at, const mymatrix &mx, mymatrix &mc);
inline bool MatrixAddMatrix(const mymatrix &ma, const mymatrix &mb, mymatrix &mc);
inline bool Norm2Vector(const myvector &c, myvector &cc);
inline void ShowMatrix(mymatrix &mtx);
inline bool RotationMatrix(const myvector &axis, double angle, mymatrix &romtx);
//bool CoordinateRotation(myvector &pointA, const myvector &axis, double angle, myvector &pointB);
bool CoordinateRotation(const myvector &pointA, const myvector &axisA, const myvector &axisB, double angle, myvector &pointB);


double Distance(double *p1, double *p2);
double Angle(double *p1, double *p2, double *p3);
double Dihedral(double *p1, double *p2, double *p3, double *p4);
bool internal2cartesian (double *c1, double *c2, double *c3, double *p, double *c4);


void crossproduct(double *c1, double *c2, double *cc);
double innerproduct(double *c1, double *c2);
void vectorsum(double *c1, double *c2, double *cc);
void subtract(double *c1, double *c2, double *cc);
void norm(double *c);
void multi(double coefficient, double *c);

double degtorad(double deg);
double rad2deg(double rad);
///////////////////////////////////////////////////////////////////

double CompareAngle(double angle1, double angle2);

////////////////come from panqing//////////////////
void VectorRotation (double v0[3], double v1[3], double mymatrix[9]);
void CrossProductNormalize (double[], double[], double[]);
void RotateMatrix_z (double, double, double, double[][3]);
void CrossProduct (double[], double[], double[]);
void DNormalizationOfVector (double[]);
void GaussInverse3 (double *, double, int *);
void MatrixTimesMatrix (double *, double *, double *, int, int, int);
double InnerProduct (double[], double[]);
void Exchangerowcolumn (double *, int, int, int, int);
////////////////////////////////////////////////////////////////

double MeanVal(vector<double> x);
double Std(vector<double> x);
double PearsonCoefficient(vector<double> x, vector<double> y);
void ReadProfPssm(std::ifstream &infile,IV2& fary,IV2& mary,std::string &seq);
#endif
