				/*********************************************************
				Copyright	    	 : Edited by Wuaiping, 2006
				Program name	 	 : PTADEE
				Program version	 : 0.2
				Begin time		 : Dec. 28 2005
				Refined time 	 : Jun. 14 2007
				Email    		 : wuaiping@moon.ibp.ac.cn
				*********************************************************/

#include "MathTool.hpp"
#include <iostream>

const double PI=3.1415926535897932384626433832795028841971693993751;
const double EXTRA=1.0e-10;

using namespace std;



//==============================================================//
//Method to initialize a MATRIX
//
//==============================================================//
inline void SetMatrix(MATRIX &sm, int m, int n)
{
   vector<double> tmp(n, 0);
   sm.assign(m, tmp);
   //tmp.~vector<double>();
}

inline void MatrixTimesMatrix(const MATRIX &a, const MATRIX &v, MATRIX &x, int m, int n, int l)
{
   int i=0, j=0, k=0;
   double sum=0;
   for(i=0; i<m; ++i)
   {
      for(k=0; k<l; ++k)
      {
         for(j=0, sum=0; j<n; ++j)
         {
            //cout<<a[i][j]<<" * "<<v[j][k]<<" + ";
            sum+=a[i][j]*v[j][k];
         }
         x[i][k]=sum;
      }
   }
}

inline bool TransVectorTimesVector(const myvector &trans, const myvector &vctor, MATRIX &mtx)
{
   if(trans.size()!=vctor.size())
   {
      cout<<"Error in TransVectorTimesVector()"<<endl;
      return false;
   }
   int i=0, j=0;
   SetMatrix(mtx, trans.size(), vctor.size());
   for(i=0; i<trans.size(); ++i)
      for(j=0; j<vctor.size(); ++j)
         mtx[i][j]=trans[i]*vctor[j];
   
   return true;
}

inline bool MatrixTimesTransVector(const MATRIX &mtx, const myvector &tvt, myvector &vct)
{
   if(mtx[0].size()!=tvt.size())
   {
      cout<<"Error in MatrixTimesTransVector()"<<endl;
      return false;
   }
   int i=0, j=0;
   vct.assign(mtx.size(), 0);
   for(i=0; i<mtx.size(); ++i)
      for(j=0; j<mtx[i].size(); ++j)
         vct[i]+=mtx[i][j]*tvt[j];
   
   return true;
}

inline void RealTimesMatrix(double at, const MATRIX &mx, MATRIX &mc)
{
   int i=0, j=0;
   for(i=0; i<mx.size(); ++i)
      for(j=0; j<mx[i].size(); ++j)
         mc[i][j]=at*mx[i][j];
}

inline bool MatrixAddMatrix(const MATRIX &ma, const MATRIX &mb, MATRIX &mc)
{
   if(ma.size()!=mb.size() || ma[0].size()!=mb[0].size()
   ||ma.size()!=mc.size() || ma[0].size()!=mc[0].size())
   {
      cout<<"Error size! MatrixAddMatrix()"<<endl;
      return false;
   }
   int i=0, j=0;
   for(i=0; i<ma.size(); ++i)
      for(j=0; j<ma[i].size(); ++j)
         mc[i][j]=ma[i][j]+mb[i][j];
    
   return true;     
}

inline bool Norm2Vector(const myvector &c, myvector &cc)
{
   if(c.size()!=cc.size())
   {
      cout<<"Error in Norm2Vector() "<<endl;
      return false;
   }
   int i=0; 
   double len=0;
   for(i=0; i<c.size(); ++i)
   {
      len+=c[i]*c[i];
   }
   len=sqrt(len);
   if(len<EXTRA) return false;
   for(i=0; i<c.size(); ++i)
      cc[i]=c[i]/len;
   
   return true;
   
}


inline void ShowMatrix(MATRIX &mtx)
{
   cout<<"ShowMatrix"<<endl;
   for(int i=0; i<mtx.size(); ++i)
   {
      for(int j=0; j<mtx[i].size(); ++j)
         cout<<"["<<mtx[i][j]<<"]	";
      cout<<endl;
   }
}


/************************Rotation Matrix****************************
* Give the Rotation Axis and the Rotation Angle
* Generate the Rotation Matrix
* Author: C.Y.
* Date 2006.11.30.
*************************End***************************************/
inline bool RotationMatrix(const myvector &axis, double angle, MATRIX &romtx)
{
   if(axis.size()!=3) 
   {
      cout<<"Error 1 in RotationMatrix()"<<endl;
      return false;
   }
   myvector ouc(3, 0);
   if(!Norm2Vector(axis, ouc)) 
   {
      cout<<"Error 2 in RotationMatrix()"<<endl;
      return false;
   }
   int i=0, j=0;
   MATRIX su, u_i, unit, tmp;
   
   SetMatrix(su, 3, 3);
   su[0][0]=0; su[0][1]=-ouc[2]; su[0][2]=ouc[1]; 
   su[1][0]=ouc[2]; su[1][1]=0; su[1][2]=-ouc[0]; 
   su[2][0]=-ouc[1]; su[2][1]=ouc[0]; su[2][2]=0;
   RealTimesMatrix(sin(angle), su, su);
   TransVectorTimesVector(ouc, ouc, u_i);
   
   SetMatrix(unit, 3, 3);
   SetMatrix(tmp, 3, 3);
   unit[0][0]=1; unit[1][1]=1; unit[2][2]=1;
   RealTimesMatrix(-1, u_i, tmp);
   MatrixAddMatrix(unit, tmp, tmp);
   RealTimesMatrix(cos(angle), tmp, tmp); 
   
   MatrixAddMatrix(u_i, tmp, tmp);
   MatrixAddMatrix(su, tmp, tmp);
   
   SetMatrix(romtx, 4, 4);
   for(i=0; i<tmp.size(); ++i)
      for(j=0; j<tmp[i].size(); ++j)
         romtx[i][j]=tmp[i][j];
   romtx[0][3]=0; romtx[1][3]=0; romtx[2][3]=0; 
   romtx[3][0]=0; romtx[3][1]=0; romtx[3][2]=0; 
   romtx[3][3]=1;
   
   return true;
}
/***************************CoordinateRotation********************
* Give the Rotation Axis and the Rotation Angle for a PointA
* Generate the coordinate of PointB after operation
* Author: C.Y.
* Date 2006.12.2.  
*****************************END**********************************/
/*
bool CoordinateRotation(myvector &pointA, const myvector &axis, double angle, myvector &pointB)
{
   if(pointA.size()!=3)
   {
      cout<<"Error in CoordinateRotation()"<<endl;
      return false;
   }
   MATRIX rotmtx;
   if(!RotationMatrix(axis, deg2rad(angle), rotmtx)) return false;
   
   pointA.push_back(1);
   if(!MatrixTimesTransVector(rotmtx, pointA, pointB))
   { 
      pointA.pop_back(); 
      return false;
   }
   pointA.pop_back(); pointB.pop_back();
   return true;
}
*/
/***************************CoordinateRotation********************
* For a PointA, Given the Coordinates of 2 points(axisA, axisB) 
* in line of Rotation Axis and the Rotation Angle 
* Generate the coordinate of PointB after operation
* Author: C.Y.
* Date 2006.12.2.  
*
*****************************END**********************************/
bool CoordinateRotation(const myvector &pointA, const myvector &axisA, const myvector &axisB, double angle, myvector &pointB)
{
   if(pointA.size()!=3 || axisA.size()!=3 || axisB.size()!=3)
   {
      cout<<"Error in CoordinateRotation()"<<endl;
      return false;
   }
   
   myvector axis(3, 0);
   axis[0]=axisB[0]-axisA[0]; 
   axis[1]=axisB[1]-axisA[1];
   axis[2]=axisB[2]-axisA[2]; 
   
   MATRIX rotmtx;
   if(!RotationMatrix(axis, deg2rad(angle), rotmtx)) return false;
   
   myvector point_A(3, 0);
   point_A[0]=pointA[0]-axisA[0];
   point_A[1]=pointA[1]-axisA[1];
   point_A[2]=pointA[2]-axisA[2];
   
   point_A.push_back(1);
   if(!MatrixTimesTransVector(rotmtx, point_A, pointB))
   { 
      point_A.pop_back(); 
      return false;
   }
   point_A.pop_back(); pointB.pop_back();
   
   pointB[0]+=axisA[0];
   pointB[1]+=axisA[1];
   pointB[2]+=axisA[2];
   
   return true;
}




/***************************************************************
//Calculate the distance between point p1 and p2
//p1[x, y, z], p2[x, y, z]
***************************************************************/
double Distance(double *p1, double *p2)
{
	double distance=0;
	
	distance += (p1[0]-p2[0])*(p1[0]-p2[0]);
	distance += (p1[1]-p2[1])*(p1[1]-p2[1]);
	distance += (p1[2]-p2[2])*(p1[2]-p2[2]);

	return sqrt(distance);
}



/***************************************************************
//Calculate the distance between point p1 and p2
//p1[x, y, z], p2[x, y, z]
***************************************************************/
double VDistance(vector<double> p1, vector<double> p2)
{
	double distance=0;
	
	distance += (p1[0]-p2[0])*(p1[0]-p2[0]);
	distance += (p1[1]-p2[1])*(p1[1]-p2[1]);
	distance += (p1[2]-p2[2])*(p1[2]-p2[2]);

	return sqrt(distance);
}


/***************************************************************
//Calculate the angle between p1-p2-p3
//p1[x, y, z], p2[x, y, z], p3[x, y, z]
***************************************************************/
double Angle(double *p1, double *p2, double *p3)
{
	double angle=0.0;
	
	double c1[3], c2[3];
	c1[0]=p1[0]-p2[0];
	c1[1]=p1[1]-p2[1];
	c1[2]=p1[2]-p2[2];
	
	c2[0]=p3[0]-p2[0];
	c2[1]=p3[1]-p2[1];
	c2[2]=p3[2]-p2[2];	
	
	double numerator=c1[0]*c2[0]+c1[1]*c2[1]+c1[2]*c2[2];
	double denominator1=sqrt(c1[0]*c1[0]+c1[1]*c1[1]+c1[2]*c1[2]);
	double denominator2=sqrt(c2[0]*c2[0]+c2[1]*c2[1]+c2[2]*c2[2]);
	
	if((denominator1<1.0e-10 && denominator1>-1.0e-10) || (denominator2<1.0e-10 && denominator2>-1.0e-10))
	{
		angle=0.0;
	}
	else
	{
		double tempval=numerator/(denominator1*denominator2);
		if(tempval>1.0 || tempval<-1.0)
		{
			cout << "\nAngle()--Error: illegal three points!" << endl;
			return 360;
		}
		else
		{
			angle=acos(numerator/(denominator1*denominator2));
			angle=(angle/PI)*180.0;
		}
	}
	
	return angle;
}


/***************************************************************************
//Calculate the dihedral between plane p1-p2-p3 and p2-p3-p4
//p1[x, y, z], p2[x, y, z], p3[x, y, z], p4[x, y, z]
***************************************************************************/
double Dihedral(double *p1, double *p2, double *p3, double *p4)
{
   double vector1[3];
   subtract(p1, p2, vector1);
   double vector2[3];
   subtract(p2, p3, vector2);
   double vector3[3];
   subtract(p3, p4, vector3);
   
   double v1[3];
   crossproduct(vector2, vector1, v1);
   double v2[3];
   crossproduct(vector3, vector2, v2);
   
   norm (v1);
   norm (v2);
   
   double dihedral = innerproduct(v1, v2);
   
   if (dihedral>1)
   {
      dihedral=1;
   }
   else if(dihedral<-1)
   {
      dihedral=-1;
   }
  // else if(dihedral>1+EXTRA || dihedral<-1-EXTRA)
   //{
  //    cout<<"Error, double Dihedral()\n";
  //   exit(-1);
   //}
   
   double v5[3];
   crossproduct(v2, v1, v5);
   double direction = innerproduct(v5, vector2);
   
 
   if (direction>0)
    {
       return  (acos(dihedral)/PI)*180;
    }  
    else
    {  
       return -(acos(dihedral)/PI)*180;
    }
   
}



//=====================================================================
//Translate Internal data into Cartesian data for 
//	the last of the four points
//====================================================================
bool internal2cartesian (double *c1, double *c2, double *c3, double *p, double *c4)
{
    double EXTRA=1.0e-20;
    
   double d1[3];
   double d2[3];
   double xp[3];

   subtract(c1, c2, d1);
   subtract(c3, c2, d2);
   crossproduct(d1, d2, xp);

   if( (xp[0]<EXTRA && xp[0]<-EXTRA) && (xp[1]<EXTRA && xp[1]>-EXTRA) && (xp[2]<EXTRA && xp[2]>-EXTRA) )
    {
        cout<<"Error! Points 1, 2, 3 are in one line & no plane!\n";

        return false;        
    }

   double d3[3];
   double yp[3];
   double r[3] ;
   double ypp[3];
   double tmp1[3];

   crossproduct(d2, xp, yp);
   double ang1 = deg2rad(p[2]);
   norm(xp);
   norm(yp);
   multi(cos(ang1), xp);
   multi(sin(ang1), yp);
   vectorsum( xp, yp, r);

   crossproduct(d2, r, ypp);
   double ang2 = deg2rad(p[1]);
   norm(d2 );
   norm(ypp);
   multi(-cos(ang2), d2);
   multi(sin(ang2), ypp);
   vectorsum( d2, ypp, d3 );
   multi(p[0], d3);
   vectorsum(c3, d3, tmp1);

   c4[0] = tmp1[0];
   c4[1] = tmp1[1];
   c4[2] = tmp1[2];


   return true;

}



/***************************************************
>> 
***************************************************/
void crossproduct(double *c1, double *c2, double *cc)
{
   cc[0] = c1[1] * c2[2] - c1[2] * c2[1];
   cc[1] = c1[2] * c2[0] - c1[0] * c2[2];
   cc[2] = c1[0] * c2[1] - c2[0] * c1[1];

}


/***************************************************
>> 
***************************************************/
double innerproduct(double *c1, double *c2)
{
   return c1[0] * c2[0] + c1[1] * c2[1] + c1[2] * c2[2];
}


/***************************************************
>> 
***************************************************/
void vectorsum(double *c1, double *c2, double *cc)
{

   cc[0] = c1[0] + c2[0];
   cc[1] = c1[1] + c2[1];
   cc[2] = c1[2] + c2[2];

}



/***************************************************
>> 
***************************************************/
void subtract(double *c1, double *c2, double *cc)
{

   cc[0] = c1[0] - c2[0];
   cc[1] = c1[1] - c2[1];
   cc[2] = c1[2] - c2[2];

}



/***************************************************
>> 
***************************************************/
void norm(double *c)
{

   double len = sqrt(c[0]*c[0] + c[1]*c[1] +c[2]*c[2]);
   c[0] /= len;
   c[1] /= len;
   c[2] /= len;

}



/***************************************************
>> 
***************************************************/
void multi(double coefficient, double *c)
{
   c[0] *= coefficient;
   c[1] *= coefficient;
   c[2] *= coefficient;

}


/***************************************************
>> 
***************************************************/
double deg2rad(double deg)
{
   return deg*PI/180;
}



/***************************************************
>> 
***************************************************/
double rad2deg(double rad)
{
   return rad*180/PI;
}


/*******************************************************************
>> Calculate the difference between two angles 1 and 2
>> The output value must be positive and between 0-180
********************************************************************/
double CompareAngle(double angle1, double angle2)
{
	double outval=0.0;
	int i=0, j=0;
	
	double temp=angle1*angle2;
	if(temp>=0)
	{
		outval=abs(angle1-angle2);
	}
	else
	{
		outval=abs(angle1-angle2);
		if(outval>180)
		{
			outval = 360-outval;
		}
	}
	
	return outval;
}


/*-------------------------------------------------------------------------------------
 * this function is rotation: vector v1 rotate to vector v0 needs the MATRIX
 * 
 * [v1x,v1y,v1z]^T = M * [v0x,v0y,v0z]^T
 	v0 is the target vector, v1 is the vector to be rotated.
 * ----------------------------------------------------------------------------------*/
void
VectorRotation (double v0[3], double v1[3], double MATRIX[9])
{
  int i, j;
  double MATRIX0[3][3], MATRIX1[3][3], tmp0[9], tmp1[9], tmp[9];


  RotateMatrix_z (v0[0], v0[1], v0[2], MATRIX0);
  RotateMatrix_z (v1[0], v1[1], v1[2], MATRIX1);

  for (i = 0; i < 3; i++)
    {
      for (j = 0; j < 3; j++)
	{
	  tmp0[i * 3 + j] = MATRIX0[i][j];
	  tmp1[i * 3 + j] = MATRIX1[i][j];
	}
    }


  i = 0;
  GaussInverse3 (tmp1, 10e-14, &i);
  
/*  if (i==0)
  {
  printf("the MATRIX is singular\n");
  exit(0);
  
  }
  
 */
  
  MatrixTimesMatrix (tmp0, tmp1, tmp, 3, 3, 3);

  for (i = 0; i < 3; i++)
    {
      for (j = 0; j < 3; j++)
	{
	  MATRIX[i * 3 + j] = tmp[j * 3 + i];
	}
    }


}


/* --------------------------------------------------------------------------*/
void
CrossProduct (double u[], double v[], double w[])
{
  w[0] = u[1] * v[2] - u[2] * v[1];
  w[1] = u[2] * v[0] - v[2] * u[0];
  w[2] = u[0] * v[1] - v[0] * u[1];
}

/*--------------------------------------------------------------------------*/
void
CrossProductNormalize (double U[], double V[], double W[])
{
  double Normal;

  /* cross product               */
  W[0] = U[1] * V[2] - U[2] * V[1];
  W[1] = U[2] * V[0] - U[0] * V[2];
  W[2] = U[0] * V[1] - U[1] * V[0];
  /* normal                      */
  Normal = W[0] * W[0] + W[1] * W[1] + W[2] * W[2];
  Normal = sqrt (Normal);

  W[0] = W[0] / Normal;
  W[1] = W[1] / Normal;
  W[2] = W[2] / Normal;
}

/*------------------------------------------------------------------------*/
double
InnerProduct (double u[], double v[])
{
  double w;

  w = u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
  return (w);
}

/*---------------------------------------------------------------------------
 * Distance_Two_points -- Distance of two 3D points
 * ----------------------------------------------------------------------------*/
double
Distance_Two_points3D (double p1[], double p2[])
{
  double result;

  result = sqrt ((p1[0] - p2[0]) * (p1[0] - p2[0]) +
		 (p1[1] - p2[1]) * (p1[1] - p2[1]) +
		 (p1[2] - p2[2]) * (p1[2] - p2[2]));

  return (result);
}


/*-----------------------------------------------------------------------
*  Rotate the vector (nx,ny,nz)^T to be z-axis        
*                 That is, MATRIX^T*(nx,ny,nz)^T = [0,0, |n|]^T          
--------------------------------------------------------------------------*/
void
RotateMatrix_z (double nx, double ny, double nz, double MATRIX[3][3])
{
  double c1, c2, s1, s2, normal, normalz;

  normal = sqrt (nx * nx + ny * ny + nz * nz);
  normalz = sqrt (nx * nx + ny * ny);
  c1 = nz / normal;
  c2 = -ny / normalz;
  s1 = normalz / normal;
  s2 = nx / normalz;

  if (normalz < 0.001)
    {
      MATRIX[0][0] = 1.0;
      MATRIX[0][1] = 0.0;
      MATRIX[0][2] = 0.0;
      MATRIX[1][0] = 0.0;
      MATRIX[1][1] = 1.0;
      MATRIX[1][2] = 0.0;
      MATRIX[2][0] = 0.0;
      MATRIX[2][1] = 0.0;
      MATRIX[2][2] = 1.0;
    }
  if (normalz >= 0.001)
    {
      MATRIX[0][0] = c2;
      MATRIX[0][1] = -c1 * s2;
      MATRIX[0][2] = s1 * s2;
      MATRIX[1][0] = s2;
      MATRIX[1][1] = c1 * c2;
      MATRIX[1][2] = -s1 * c2;
      MATRIX[2][0] = 0.0;
      MATRIX[2][1] = s1;
      MATRIX[2][2] = c1;
    }
}

/*--------------------------------------------------------------------------*/
void
DNormalizationOfVector (double vector[])
{
  double length;

  length = sqrt (vector[0] * vector[0] + vector[1] * vector[1] +
		 vector[2] * vector[2]);
  if (length == 0.0)
    {
      printf ("the length of the vector is zero\n");
      return;
    }

  vector[0] = vector[0] / length;
  vector[1] = vector[1] / length;
  vector[2] = vector[2] / length;
}

/*----------------------------------------------------------------------------
 gaussInverse3                                                          
 - Inverse 3 * 3  MATRIX A by Gauss elimination methods                   
*----------------------------------------------------------------------------*/
void
GaussInverse3 (double *a, double eps, int *message)
{
  double max;
  int n, k, ik, jk, i, j, z[6];

  n = 3;
  *message = 1;
  for (k = 0; k < n; k++)
    {
      max = 0.0;
      for (i = k; i < n; i++)
	for (j = k; j < n; j++)
	  if (fabs (*(a + i * n + j)) > max)
	    {
	      ik = i;
	      jk = j;
	      max = fabs (*(a + i * n + j));
	    }
      if (max < eps || max == 0.0)
	{
	  *message = 0;
	  printf ("The MATRIX in gaussinverse3 is singular, %lf\n", max);
	  return;
	}
      max = 1.0 / *(a + ik * n + jk);
      *(a + ik * n + jk) = 1.0;
      z[2 * k] = ik;
      z[2 * k + 1] = jk;
      Exchangerowcolumn (a, n, k, ik, jk);
      for (j = 0; j < n; j++)
	*(a + k * n + j) = *(a + k * n + j) * max;
      for (i = 0; i < n; i++)
	if (i != k)
	  {
	    max = *(a + i * n + k);
	    *(a + i * n + k) = 0.0;
	    for (j = 0; j < n; j++)
	      *(a + i * n + j) = *(a + i * n + j) - max * *(a + k * n + j);
	  }
    };
  for (k = n - 2; k > -1; k--)
    {
      ik = z[2 * k + 1];
      jk = z[2 * k];
      Exchangerowcolumn (a, n, k, ik, jk);
    }
}

/*---------------------------------------------------------------------------------
* Exchangerowcolumn                                                     
* - Exchange two rows and two columns of the MATRIX A ,
*   exchange the (k,ik) rows
*                (k,ik) column
----------------------------------------------------------------------------------*/
void Exchangerowcolumn (double *a, int n, int k, int ik, int jk)
{
  double b;
  int j;

  /*  exchange the (k,ik) rows */
  if (ik != k)
    {
      for (j = 0; j < n; j++)
	{
	  b = *(a + ik * n + j);
	  *(a + ik * n + j) = *(a + k * n + j);
	  *(a + k * n + j) = b;
	}
    }
  /*  exchange the (k,ik) column  */
  if (jk != k)
    {
      for (j = 0; j < n; j++)
	{
	  b = *(a + j * n + jk);
	  *(a + j * n + jk) = *(a + j * n + k);
	  *(a + j * n + k) = b;
	}
    }
}

/*------------------------------------------------------------------------
* MatrixTimesMatrix                                                     
* times a MATRIX a with MATRIX v, the result is x. i.e., x = a*v        
*   a;       pointe to a (m,n) MATRIX               
*   v;       pointe to a (n,l) MATRIX                 
*   x;      pointe to a (m,l) MATRIX         
-------------------------------------------------------------------------*/
void
MatrixTimesMatrix (double *a, double *v, double *x, int m, int n, int l)
{

  int i, j, k;
  double sum;

  for (k = 0; k < l; k++)
    for (i = 0; i < m; i++)
      {
	sum = 0.0;
	for (j = 0; j < n; j++)
	  sum = sum + *(a + i * n + j) * *(v + j * l + k);
	*(x + i * l + k) = sum;
      }
}



//========================================================================//
//Calculate the mean value of the input vector
//
//========================================================================//
double MeanVal(vector<double> x)
{
	double mean=0.0;
	
	for(int i=0; i<x.size(); i++)
	{
		mean += x[i];
	}
	mean = mean/x.size();
	
	return mean;
}



//========================================================================//
//Calculate the standard square deviation of the input vector
//
//========================================================================//
double Std(vector<double> x)
{
	double std=0.0;
	
	double mean=MeanVal(x);
	int N=x.size();
	for(int i=0; i<N; i++)
	{
		std += (x[i]-mean)*(x[i]-mean);
	}
	
	std /= N;
	std = sqrt(std);
	
	
	return std;
}




//========================================================================//
//Calculate the Pearson Coefficient of the two input vectors
//
//========================================================================//
double PearsonCoefficient(vector<double> x, vector<double> y)
{
	double coe=0.0;
	int N=x.size();
	
	double meanx=MeanVal(x);
	double meany=MeanVal(y);
	double stdx=Std(x);
	double stdy=Std(y);
	for(int i=0; i<N; i++)
	{
		coe += ((x[i]-meanx)/stdx)*((y[i]-meany)/stdy);
	}
	coe /= N;
	
	
	return coe;
}


