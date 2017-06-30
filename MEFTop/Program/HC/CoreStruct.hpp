//=======================================================
//
// Description:
//	Seperate hydrophobic-core from protein and
//		output its structure.
//
// Contact: wuaiping@moon.ibp.ac.cn, 2010_05_19
//
//=======================================================


#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cstdlib>
#include "HydroFace.hpp"
#include "HydroCluster.hpp"


using namespace std;



typedef struct _HcoreStr
{
	vector<int> hcmark;		//hcmark[i]=1 or 0, "1" meas hydrophobic-core positon, "0" means not.
	
	string gapseq;			//sequence with gap mark "-", means chain with break.
	string gap2nd;			//second structure with gap mark "-", means chain with break.
	
	vector<char> rname;		//residue names, 1-letter.
	vector<int> rserial;		//residue serials, not begin from zero.
	
	vector<vector<vector<double > > > xyz;	//coordinates of each residue: N, CA, C, O, four atoms.
	
	int Facenum;
	double radius;
} HcoreStr;



void GetHcoreStr(HcoreStr& hc, 
		  const ChainCoord& chain, 
		  const vector<int>& hcid, 
		  const vector<char>& second, 
		  const vector<ELE>& allele,
		  int Fnum);


void OutputHCore(string outname, const HcoreStr& hc);

void AAname1to3(char name1, string& name3);


