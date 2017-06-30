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

using namespace std;

/////////////////////changed from caoyang//////////////
typedef vector<vector<double> > matrix;
typedef vector<double> myvector;

inline void SetMatrix(matrix &sm, int m, int n);
inline void MatrixTimesMatrix(const matrix &a, const matrix &v, matrix &x, int m, int n, int l);
inline bool TransVectorTimesVector(const myvector &trans, const myvector &vctor, matrix &mtx);
inline bool MatrixTimesTransVector(const matrix &mtx, const myvector &tvt, myvector &vct);
inline void RealTimesMatrix(double at, const matrix &mx, matrix &mc);
inline bool MatrixAddMatrix(const matrix &ma, const matrix &mb, matrix &mc);
inline bool Norm2Vector(const myvector &c, myvector &cc);
inline void ShowMatrix(matrix &mtx);
inline bool RotationMatrix(const myvector &axis, double angle, matrix &romtx);
//bool CoordinateRotation(myvector &pointA, const myvector &axis, double angle, myvector &pointB);
bool CoordinateRotation(const myvector &pointA, const myvector &axisA, const myvector &axisB, double angle, myvector &pointB);


double Distance(double *p1, double *p2);
double MyDistance(const vector<double> &p1,const vector<double> &p2);
double Angle(double *p1, double *p2, double *p3);
double MyAngle(const vector<double> &p1, const vector<double> &p2, const vector<double> &p3);//write by Tian
double Dihedral(double *p1, double *p2, double *p3, double *p4);

//calculate dihedral between two plane(plane1: p1,p2,p3 plane2: p4,p5,p6), (-180,180), write by Tian
double MyDihedral(const vector<double> &p1, const vector<double> &p2, const vector<double> &p3, 
                  const vector<double> &p4, const vector<double> &p5, const vector<double> &p6);
				  
//calculate rotate angle between two vector(vector1: p1,p2 vector2: p3,p4), (-180,180), write by Tian
double MyRotateAngle(const vector<double> &p1, const vector<double> &p2, const vector<double> &p3, const vector<double> &p4);

void Mysubtract(const vector<double> &c1, const vector<double> &c2, double *cc);//write by Tian

void crossproduct(double *c1, double *c2, double *cc);
double innerproduct(double *c1, double *c2);
void vectorsum(double *c1, double *c2, double *cc);
void subtract(double *c1, double *c2, double *cc);
void norm(double *c);
void multi(double coefficient, double *c);

double deg2rad(double deg);
double rad2deg(double rad);
///////////////////////////////////////////////////////////////////

double CompareAngle(double angle1, double angle2);

////////////////come from panqing//////////////////
void VectorRotation (double v0[3], double v1[3], double matrix[9]);
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

#endif
