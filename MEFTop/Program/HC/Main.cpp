//=======================================================
//
// Description:
//	Define and recognize hydrophobic-core.
//
// Contact: wuaiping@moon.ibp.ac.cn, 2010_05_13
//
//=======================================================



#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cstdlib>
#include "CoreStruct.hpp"


using namespace std;


int ReadID(vector<string>& allid, string IDfile);
int HCRadius(HcoreStr& hc);


//=================================================================
// Compile:  g++ -O2 -o coredistribut Main.cpp HydroCluster.cpp HydroFace.cpp IOprotein.cpp MathTool.cpp HelixAxis.cpp
// Usage  : ./coredistribut  ../../data/domain40.ID  /tjjiang/wuaiping/DATABASE/SCOP-v1.73/coordfile/allpdb/  /tjjiang/wuaiping/DATABASE/SCOP-v1.73/dssp/seq40/
//=================================================================
//=================================================================
// 进行所有蛋白结构内疏水内核分布信息的统计；
//=================================================================
int main(int argc, char* argv[])
{
	int i=0, j=0, k=0, p=0;

	string prtid=argv[1];
	//string idfile="/tjjiang/daiwentao/Database/HydroCore/ssSCOP/Decoy/test.id";




	vector<int> Ncore, Nface, Nresi;
	string pdbpath =argv[2];
	string dssppath=argv[3];
	string outpath =argv[4];
	//string pdbpath="/tjjiang/daiwentao/Database/HydroCore/ssSCOP/Decoy/test/";
	//string dssppath="/tjjiang/daiwentao/Database/HydroCore/ssSCOP/Decoy/testdssp/";
	//string outpath="/tjjiang/daiwentao/Database/HydroCore/ssSCOP/Decoy/";
	//string pdbpath="/tjjiang/wuaiping/DATABASE/SCOP-v1.73/coordfile/allpdb/";
	//string dssppath="/tjjiang/wuaiping/DATABASE/SCOP-v1.73/dssp/seq40/";



	//== read pdb and add the side-chain centroid for each residue;
	//string temppdb=pdbpath+ID[p]+".ent";
	//string temppdb=pdbpath+ID[p]+".A1.f01.full.Lf.pdb";
	//string temppdb=pdbpath+ID[p]+".cut";
	string temppdb=pdbpath+prtid+".pdb"; 
	string tempdssp=dssppath+prtid+".dssp";		

	char pdbfile[500], dsspfile[500];
	strcpy(pdbfile, temppdb.c_str());
	strcpy(dsspfile, tempdssp.c_str());
	ChainCoord chain;
	Read4AtomPdb(chain, pdbfile);
	bool tf=true;
	for(i=0; i<chain.xyz.size(); i++)
	{
		if(chain.xyz[i].size()!=4) {tf=false; break;}
	}
	if(!tf) 
		return -1;
	//continue;

	vector<char> resi, second;
	vector<double> acc;
	ReadDSSP(resi, second, acc, dsspfile);
	if(chain.resinum!=resi.size() || chain.resinum>1000) 
		return -1;
	//continue;		


	//if(second.size()!=chain.xyz.size()) continue;
	char cenfile[500]="./bin/centroidrotamer";
	DV4 centable;
	centable.resize(19);
	for(i=0; i<19; i++)
	{
		centable[i].resize(36);
		for(j=0; j<36; j++)
		{
			centable[i][j].resize(36);
			for(k=0; k<36; k++)
			{
				centable[i][j][k].resize(6);
			}
		}
	}
	ReadCentroid(centable, cenfile);
	AddCentroid(chain, centable);


	//== calculate the distance matrix among all Centroid atoms;
	vector<char> rname;
	for(i=0; i<chain.resiname.size(); i++)
	{
		char ch;
		AAname3to1(ch, chain.resiname[i]);
		rname.push_back(ch);
	}
	vector<char> HPseq;
	Seq2HP(HPseq, rname);

	vector<int> hi;
	for(i=0; i<HPseq.size(); i++)
	{
		if(HPseq[i]=='H') hi.push_back(i);
	}


	DV2 Dcent;	// the distance between each Cc-Cc pair
	DV3 NVcent;	// normalized vector of each Cc-Cc pair
	CentroidMatrix(Dcent, NVcent, chain, hi);


	//== cluster all hydrophobic residues into classes
	vector<vector<int > > map;
	ClusterMap(map, Dcent, NVcent, chain, hi);


	//== hydro-face-based hydrophobic-cores' recognization
	vector<ELE> allele;
	InitialSSE(allele, second, chain);
	for(j=0; j<allele.size(); j++)
	{
		string tseq="", hp="";
		for(k=allele[j].pos; k<allele[j].pos+allele[j].len; k++)
		{
			tseq.push_back(resi[k]);
			hp.push_back(HPseq[k]);
		}
		allele[j].seq=tseq;
		allele[j].HPseq=hp;

		HydroFace(allele[j]);
	}

	vector<vector<int > > FC;
	FaceCluster(FC, map, Dcent, hi, allele);

	/*
	for(i=0; i<FC.size(); i++)
	{
	for(j=0; j<FC[i].size(); j++)
	{
	cout << FC[i][j]+chain.resiserial[0] << ' ';
	}cout << endl;
	}cout << endl;	*/


	//== the number of hydrophobic cores in each structure.
	int num=0;
	for(i=0; i<FC.size(); i++)
	{
		if(FC[i].size()>=6) //At least 6 residues in one HydroCore
		{	
			int fnum=FaceInCluster(allele, FC[i]);
			if(fnum>1) //At least 2 HydroFaces in one HydroCore
			{
				num++;
				cout << "Core " << num << ": ";
				for(j=0; j<FC[i].size(); j++) cout << chain.resiserial[FC[i][j]] << ' ';
				cout << endl;
				Nface.push_back(fnum);
				Nresi.push_back(FC[i].size());
			}

		}
	}
	cout << p+1 << '\t' << prtid << '\t' << chain.resinum << '\t' << hi.size() << '\t' << num << endl;
	Ncore.push_back(num);

	//== BE CAREFUL: all core-structs were outputed, used FC[i].size()>=6 and HFace>1 constraint.
	string core="_core";
	string suffix=".pdb";
	for(i=0; i<FC.size(); i++)
	{
		int fnum=FaceInCluster(allele, FC[i]);
		if(FC[i].size()>=6 && fnum>1)
		{								
			HcoreStr hc;
			for(j=0; j<FC[i].size(); j++)
			{
				GetHcoreStr(hc, chain, FC[i], second, allele, fnum);
			}


			//===Statistic the probilty of Hydrophobic Core
			HCRadius(hc);

			//===Output the Hydro Core
			char chindex[10]="";
			sprintf(chindex, "%d", i+1);		
			string outname=outpath+prtid+core+chindex+suffix;
			OutputHCore(outname, hc);
		}
	}



	/*	
	vector<int> bin;
	bin.resize(20);
	for(i=0; i<Ncore.size(); i++)
	{
	if(Ncore[i]>19) bin[19]++;
	else bin[Ncore[i]]++;
	}
	cout << "\nThe distribution of hydro-core number:" << endl;
	for(i=0; i<20; i++)
	{
	cout << i << '\t' << bin[i]*1.0/Ncore.size() << endl;
	}

	vector<int> fbin;
	fbin.resize(20);
	for(i=0; i<Nface.size(); i++)
	{
	if(Nface[i]>19) fbin[19]++;
	else fbin[Nface[i]]++;
	}
	cout << "\nThe distribution of hydro-faces in hydro-core:" << endl;
	for(i=0; i<20; i++)
	{
	cout << i << '\t' << fbin[i]*1.0/Nface.size() << endl;
	}

	vector<int> rbin;
	rbin.resize(100);
	for(i=0; i<Nresi.size(); i++)
	{
	if(Nresi[i]>99) rbin[99]++;
	else rbin[Nresi[i]]++;
	}
	cout << "\nThe distribution of residues in hydro-core:" << endl;
	for(i=0; i<100; i++)
	{
	cout << i << '\t' << rbin[i]*1.0/Nresi.size() << endl;
	}
	*/	


	cout << "\nProgram exit normally!" << endl << endl;

	return 0;
}




/*
//=================================================================
// Compile:  g++ -O2 -o hydrocore Main.cpp HydroCluster.cpp HydroFace.cpp IOprotein.cpp MathTool.cpp HelixAxis.cpp
// Usage  : ./hydrocore Input/1A0FA.pdb Input/1A0FA.dssp | tee log
//=================================================================
//=================================================================
// 进行蛋白结构内疏水残基的聚类，疏水内核的识别；
//=================================================================
int main(int argc, char* argv[])
{
int i=0, j=0, k=0;


//== read pdb and add the side-chain centroid for each residue;
//char pdbfile[500]="Input/1A0FA.pdb";
char pdbfile[500];
strcpy(pdbfile, argv[1]);
ChainCoord chain;
Read4AtomPdb(chain, pdbfile);

char cenfile[500]="/tjjiang/wuaiping/HydroCore/data/centroidrotamer";
DV4 centable;
centable.resize(19);
for(i=0; i<19; i++)
{
centable[i].resize(36);
for(j=0; j<36; j++)
{
centable[i][j].resize(36);
for(k=0; k<36; k++)
{
centable[i][j][k].resize(6);
}
}
}
ReadCentroid(centable, cenfile);
AddCentroid(chain, centable);


//== calculate the distance matrix among all Centroid atoms;
vector<char> rname;
for(i=0; i<chain.resiname.size(); i++)
{
char ch;
AAname3to1(ch, chain.resiname[i]);
rname.push_back(ch);
}
vector<char> HPseq;
Seq2HP(HPseq, rname);

vector<int> hi;
for(i=0; i<HPseq.size(); i++)
{
if(HPseq[i]=='H') hi.push_back(i);
}


DV2 Dcent;	// the distance between each Cc-Cc pair
DV3 NVcent;	// normalized vector of each Cc-Cc pair
CentroidMatrix(Dcent, NVcent, chain, hi);


//== cluster all hydrophobic residues into classes
vector<vector<int > > map;
ClusterMap(map, Dcent, NVcent, chain, hi);

//for(i=0; i<map.size(); i++)
//{ {for(j=0; j<map[i].size(); j++)
//   cout << map[i][j];
// }cout << endl;
//}cout << endl;


//vector<vector<int > > cluster;
//HydroCluster(cluster, map, Dcent, hi);

//for(i=0; i<cluster.size(); i++)
//{
//	for(j=0; j<cluster[i].size(); j++)
//	{
//		cout << cluster[i][j]+chain.resiserial[0] << ' ';
//	}cout << endl;
//}cout << endl;


//== hydro-face-based hydrophobic-cores' recognization
vector<char> resi, second;
vector<double> acc;
//char dsspfile[500]="Input/1A0FA.dssp";
char dsspfile[500];
strcpy(dsspfile, argv[2]);
ReadDSSP(resi, second, acc, dsspfile);
vector<ELE> allele;
InitialSSE(allele, second, chain);
for(j=0; j<allele.size(); j++)
{
string tseq="", hp="";
for(k=allele[j].pos; k<allele[j].pos+allele[j].len; k++)
{
tseq.push_back(resi[k]);
hp.push_back(HPseq[k]);
}
allele[j].seq=tseq;
allele[j].HPseq=hp;

HydroFace(allele[j]);
}

vector<vector<int > > FC;
FaceCluster(FC, map, Dcent, hi, allele);

for(i=0; i<FC.size(); i++)
{
for(j=0; j<FC[i].size(); j++)
{
cout << FC[i][j]+chain.resiserial[0] << ' ';
}cout << endl;
}cout << endl;



cout << "\nProgram exit normally!" << endl << endl;

return 0;
}
*/



//===================================================================
// Read all IDs from the ID-list file
//===================================================================
int ReadID(vector<string>& allid, string IDfile)
{
	ifstream fcid;
	fcid.open(IDfile.c_str());
	if( !fcid.is_open() )
	{
		cout << "Can not open " << IDfile << endl;
		exit(-1); 
	}


	string stemp="";
	while( !fcid.eof() )
	{
		stemp="";
		fcid >> stemp;
		if(stemp!="") allid.push_back(stemp);
	}
	fcid.close();

	return 0;	
}

//===================================================================
// Read all IDs from the ID-list file
//===================================================================
int HCRadius(HcoreStr& hc)
{
	int hclen=hc.rserial.size();
	double Center[3];
	double SumX =0.0, SumY =0.0, SumZ =0.0;
	double (*CAcoor)[3]=new double [hclen][3];	

	for(int i=0; i<hclen; i++)
	{
		CAcoor[i][0]=hc.xyz[i][1][0];
		CAcoor[i][1]=hc.xyz[i][1][1];
		CAcoor[i][2]=hc.xyz[i][1][2];

		SumX +=CAcoor[i][0];
		SumY +=CAcoor[i][1];
		SumZ +=CAcoor[i][2];
	}
	Center[0] =SumX/hclen;
	Center[1] =SumY/hclen;
	Center[2] =SumZ/hclen;

	SumX =0.0;SumY =0.0;SumZ =0.0;

	for(int i=0; i<hclen; i++)
	{
		SumX +=(CAcoor[i][0]-Center[0])*(CAcoor[i][0]-Center[0]);
		SumY +=(CAcoor[i][1]-Center[1])*(CAcoor[i][1]-Center[1]);
		SumZ +=(CAcoor[i][2]-Center[2])*(CAcoor[i][2]-Center[2]);
	}
	double Sum = (SumX + SumY +SumZ)/hclen;
	hc.radius=pow(Sum, 0.5); 

	if(CAcoor!=NULL)
		delete [] CAcoor;
	return 0;
}




