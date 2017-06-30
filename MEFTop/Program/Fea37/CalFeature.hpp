//===================================================================================
//
// Caculate the feature score based on the querry and model
//
// AUTHORS: daiwentao
// E-mail: daiwentao@moon.ibp.ac.cn
// Date: 2012.03.12
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
#include "ReadInfor.hpp"

using namespace std;

double resSS(string qss, string mss);
double resBE(string qrs, string mrs);
double resSS1(string qss, string mss);
double resBE1(string qrs, string mrs);
void pairSSA(vector<double>& score, QUE querry, MOD model);
double cosine(vector<double> ratio);
double correl(vector<double> ratio);
double expdist(vector<double> ratio);
double dotprod(vector<double> ratio);
void cmFeature(vector<double>& score, vector< vector<double> > cm, 
vector< vector<double> > ca, double cutoff, int nu);
void slenFeature(vector<double>& score, QUE querry, MOD model, int nu);
double ssCon(QUE querry, MOD model, vector<int> smindex); 
void sconFeature(vector<double>& score, QUE querry, MOD model, int nu);
void radFeature(vector<double>& score, double rad1, double rad2);
void cmFeatureIn(vector<int> si, vector<double>& score, vector< vector<double> > cm, 
vector< vector<double> > ca, double cutoff, int nu);
void sconFeatureIn(vector<int> resi, vector<double>& score, QUE querry, MOD model, int nu);
void slenFeatureIn(vector<int> ri, vector<double>& score, QUE querry, MOD model, int nu);
int ReadInf(string& cseq, string& chf, string& cmark, int& num, double& HCrad, const char* inffile);
void HCFeature(vector<double>& score, MOD model, vector<vector<double> > hc);


