//=======================================================
//
// Description:
//	Generate center-axis of alpha-helix.
//		(Changed from Tianliqing)
//
// Contact: wuaiping@moon.ibp.ac.cn, 2010_05_11
//
//=======================================================


#include "HelixAxis.hpp"


using namespace std;



//================================================================
// Calculate the angular bisector of the input three points.
//================================================================
void CalBiosection(const vector<double> v1,
			const vector<double> &v2,
			const vector<double> &v3,
			vector<double> &vdtemp)
{
	int i;
	
	vector<double> d1(3);
	vector<double> d2(3);
	vector<double> d3(3);
	
	Vsubtract(v3,v1,d1);
	Vsubtract(v2,v1,d2);
	Vsubtract(v3,v2,d3);
	
	vector<double> xp(3);
	
	Vcrossproduct(d2,d3,xp);
	Vcrossproduct(xp,d1,vdtemp);
	Vnorm(vdtemp);
}



//=======================================================
// Calculate the common perpendicular of two vectors
//=======================================================
void CalPerpendicular(const vector<double> &point1,const vector<double> &biosection1,
			const vector<double> &point2,const vector<double> &biosection2,
			pair<vector<double>,vector<double> > &pair_temp)
{
	int i,j;
	
	double x1,y1,z1,x2,y2,z2;
	x1=point1[0];y1=point1[1];z1=point1[2];
	x2=point2[0];y2=point2[1];z2=point2[2];
	
	double A1,B1,C1,A2,B2,C2;
	A1=biosection1[0];B1=biosection1[1];C1=biosection1[2];
	A2=biosection2[0];B2=biosection2[1];C2=biosection2[2];
	
	
	//Perpendicular vector
	vector<double> vtemp(3);
	Vcrossproduct(biosection1,biosection2,vtemp);
	Vnorm(vtemp);
	
	double a,b,c;
	a=vtemp[0];b=vtemp[1];c=vtemp[2];
	
	
	//
	double k;
	k=(a*x2+b*y2+c*z2-a*x1-b*y1-c*z1)/(a*a+b*b+c*c);
	
	//point1 mirror in plane(point2 biosection2 and Perpendicular vector(a b c) )
	double x3,y3,z3;
	x3=x1+a*k;
	y3=y1+b*k;
	z3=z1+c*k;


	//
	double m,t;
	//m=(B1*x3-B1*x2-A1*y3+A1*y2)/(B1*A2-A1*B2);
	t=(x2*B2-x3*B2-y2*A2+y3*A2)/(A1*B2-B1*A2);
	
	//Perpendicular of point2
	double x4,y4,z4;
	x4=x3+A1*t;
	y4=y3+B1*t;
	z4=z3+C1*t;


	//
	double p,q;
	p=(x1*B1-x4*B1-y1*A1+y4*A1)/(a*B1-b*A1);
	
	//Perpendicular of point1
	double x5,y5,z5;
	x5=x4+a*p;
	y5=y4+b*p;
	z5=z4+c*p;


	//set
	vector<double> O1(3),O2(3);
	O1[0]=x5,O1[1]=y5,O1[2]=z5;
	O2[0]=x4,O2[1]=y4,O2[2]=z4;
	
	pair_temp.first=O1;
	pair_temp.second=O2;
}



//=========================================================
// Calculate the center-axis of the input coordinates,
//	especially used to alpha-helix element.
//=========================================================
void HelixAxis(vector<vector<double > >& axis,
		const vector<vector<double> > &coor)
{
	int i,j;
	
	
	//calculate bisection
	vector<vector<double> > bisection;
	vector<vector<double> > CA;//erase first and last residue of the helix, CA.size()==bisection.size()
	
	for(i=1;i<coor.size()-1;i++)
	{
		CA.push_back(coor[i]);
		vector<double> vdtemp(3);
		CalBiosection(coor[i-1],coor[i],coor[i+1],vdtemp);
		bisection.push_back(vdtemp);
	}
	
	
	//calculate two point OnOn+1 which is perpendicular to both ln and ln+1, ln and ln+1 are biosection which pass CAn and CAn+1
	vector<pair<vector<double>,vector<double> > > perpendicular;
	
	for(i=0;i<CA.size()-1;i++)
	{
		pair<vector<double>,vector<double> > pair_temp;
		CalPerpendicular(CA[i],bisection[i],CA[i+1],bisection[i+1],pair_temp);
		
		perpendicular.push_back(pair_temp);
	}
	
	
	//calculate helix axis
	axis.push_back(perpendicular[0].first);
	for(i=0;i<perpendicular.size()-1;i++)
	{
		vector<double> vdtemp(3);
		CalAver(perpendicular[i].second,perpendicular[i+1].first,vdtemp);
		axis.push_back(vdtemp);
	}
	axis.push_back(perpendicular[perpendicular.size()-1].second);
	
}




