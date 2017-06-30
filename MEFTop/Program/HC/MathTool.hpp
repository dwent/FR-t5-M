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


using namespace std;

/////////////////////changed from caoyang//////////////
typedef vector<vector<double> > MATRIX;
typedef vector<double> myvector;

inline void SetMatrix(MATRIX &sm, int m, int n);
inline void MatrixTimesMatrix(const MATRIX &a, const MATRIX &v, MATRIX &x, int m, int n, int l);
inline bool TransVectorTimesVector(const myvector &trans, const myvector &vctor, MATRIX &mtx);
inline bool MatrixTimesTransVector(const MATRIX &mtx, const myvector &tvt, myvector &vct);
inline void RealTimesMatrix(double at, const MATRIX &mx, MATRIX &mc);
inline bool MatrixAddMatrix(const MATRIX &ma, const MATRIX &mb, MATRIX &mc);
inline bool Norm2Vector(const myvector &c, myvector &cc);
inline void ShowMatrix(MATRIX &mtx);
inline bool RotationMatrix(const myvector &axis, double angle, MATRIX &romtx);
//bool CoordinateRotation(myvector &pointA, const myvector &axis, double angle, myvector &pointB);
bool CoordinateRotation(const myvector &pointA, const myvector &axisA, const myvector &axisB, double angle, myvector &pointB);


double Distance(double *p1, double *p2);
double VDistance(vector<double> p1, vector<double> p2);
double Angle(double *p1, double *p2, double *p3);
double Dihedral(double *p1, double *p2, double *p3, double *p4);
bool internal2cartesian (double *c1, double *c2, double *c3, double *p, double *c4);


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
void VectorRotation (double v0[3], double v1[3], double MATRIX[9]);
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



////////// added by TianLiqing ////////////////////////
inline void CalAver(const vector<double> &v1,const vector<double> &v2,vector<double> &v3)
{
	int i;
	for(i=0;i<3;i++)
	{
		v3[i]=(v1[i]+v2[i])/2;
	}
}


inline double CalDistance(const vector<double> &p1, const vector<double> &p2)
{
	double result;

	result = sqrt ((p1[0] - p2[0]) * (p1[0] - p2[0]) +
		 		   (p1[1] - p2[1]) * (p1[1] - p2[1]) +
		 		   (p1[2] - p2[2]) * (p1[2] - p2[2]));
	
	return result;
}


//vector need define size out of the function
inline void Vsubtract(const vector<double> &c1, const vector<double> &c2, vector<double> &cc)//c1-c2
{
	cc[0] = c1[0] - c2[0];
	cc[1] = c1[1] - c2[1];
	cc[2] = c1[2] - c2[2];
}

inline void Vcrossproduct(const vector<double> &c1, const vector<double> &c2, vector<double> &cc)//c1*c2
{
   cc[0] = c1[1] * c2[2] - c1[2] * c2[1];
   cc[1] = c1[2] * c2[0] - c1[0] * c2[2];
   cc[2] = c1[0] * c2[1] - c2[0] * c1[1];
}

inline bool Vnorm(vector<double> &c)
{

	double len = sqrt(c[0]*c[0] + c[1]*c[1] +c[2]*c[2]);
	
//if(c[0]<1e-3 && c[1]<1e-3 && c[2]<1e-3 && c[0]>-1e-3 && c[1]>-1e-3 && c[2]>-1e-3) cout<<c[0]<<" "<<c[1]<<" "<<c[2]<<" "<<len<<endl;
	if(len<1e-1) return false;//special, note!!!
		
	c[0] /= len;c[1] /= len;c[2] /= len;
	
	return true;
}


#endif
