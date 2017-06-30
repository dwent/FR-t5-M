//=======================================================
//
// Description:
//	Generate center-axis of alpha-helix.
//		(Changed from Tianliqing)
//
// Contact: wuaiping@moon.ibp.ac.cn, 2010_05_11
//
//=======================================================


#ifndef HELIXAXIS_H
#define HELIXAXIS_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <map>
#include "MathTool.hpp"


using namespace std;


void CalBiosection(const vector<double> v1,
			const vector<double> &v2,
			const vector<double> &v3,
			vector<double> &vdtemp);

void CalPerpendicular(const vector<double> &point1,const vector<double> &biosection1,
			const vector<double> &point2,const vector<double> &biosection2,
			pair<vector<double>,vector<double> > &pair_temp);

void HelixAxis(vector<vector<double > >& axis,
		const vector<vector<double> > &coor);


#endif

