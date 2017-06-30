//===================================================================================
//
// Generate features from sequence and model
//
//Modify the function slenFeature() and slenFeatureIn() as model.sse.size=1
//
// AUTHORS: daiwentao
// E-mail: daiwentao@moon.ibp.ac.cn
// Date: 2012.06.18
//      
//	Add the HC modual, give a signal of hl file       
//
//===================================================================================


#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <list>
#include <iostream>
#include <ctime>
#include <cstring>
#include <vector>
#include <fstream>


#include "CalFeature.hpp"
#include "NCACO.hpp"

using namespace std;

DV2 VDWRADIUS;
DV2 CONTABLE;
DV4 TRPTABLE;
DV2 SOLVETABLE;
DV4 CENTABLE;
DV4 StruEnergySet;

char AAname3to1(string name3);
void ModelSequence(vector<int>& ind, vector<int>& resn, string& seq, string name);
void OrderSequence(vector<int>& seqin, vector<int>& resin, QUE& que, MOD& mod);

//=============================================================================
//
//
//=============================================================================
int main(int argc, const char *argv[])
{
	cout<<"!***********************!"<<endl;
	cout<< "Generate Feature 0.1!" <<endl;
	/****************************************************/
	string idfile=argv[1];
	//string idfile="/tjjiang/daiwentao/Database/CASP5_8/SubModel/TrainCasp7/Simple.id";
	vector<string> ID;
	ReadID(ID, idfile);
	QUE que;
	
	//== Read all parameters of Scoring-Function
	const string vdwfile = "parameter/vdwmatrix99_99";
	const string confile="parameter/EconTable99_99";
	const string trpfile="parameter/EtrpTable";
	const string stru_energy_file="parameter/EbetasheetTable";
	const string solvefile="parameter/EsolvTable";
	const string centroINFOfile="parameter/centroidrotamer";
	
	ReadAllPar(VDWRADIUS, vdwfile, CONTABLE, confile, TRPTABLE, trpfile, SOLVETABLE, solvefile, CENTABLE, centroINFOfile, StruEnergySet, stru_energy_file);
	
	for(int i=0; i<ID.size(); i++)
	{	
	MOD mod;
	double tmscore=1.0;
	//===Read information from files
	string seqpath=argv[2];
	string cmpath8=argv[3];
	string cmpath12=argv[4];
	string sssapath=argv[5];	
	string pdbpath=argv[6]+ID[i]+".pdb";
	//string pdbpath=argv[6]+ID[i]+".pdb"; //Some targets of casp_good
	string dssppath=argv[7]+ID[i]+".dssp";
	string HCpath=argv[8];
	string HCpath1=HCpath+ID[i]+"/"+ID[i]+".hl";
	string outpath=argv[9];
	//string tmpath1=argv[10];
	//string tmpath=tmpath1+"/"+ID[i]+".txt";
	//===Querry prediction
	if(i==0 || (ID[i].substr(0,5) != ID[i-1].substr(0,5) ) )
	{
	CleQUE(que);
	ReadSequence(seqpath.c_str(), que.seq);
	ReadCM(cmpath8, que.ss, que.rs, que.CM8a, que.seq);
	ReadCM(cmpath12, que.ss, que.rs, que.CM12a, que.seq);
	Readsssa(sssapath, que.ss, que.rs);
	que.radius=PredRG(que.seq);
	InitialSSE(que.sse, que.ss);
	}
	//===Model files
	Readpdbcoord(mod.acoord, pdbpath);
	ReadMCcoord(mod.coord, pdbpath);
	ReadDSSP(mod.ss, mod.Acc, dssppath);	
	//tmscore=ReadTM(tmpath);
	
	//==HC information
	vector<vector<double> > hcinfor;
	ifstream hfile;
	hfile.open(HCpath1.c_str());
	if( !hfile.is_open() )
	{
		//cout<<"Here"<<endl;
		vector<double> linf;
		linf.resize(3, 0.0);
		hcinfor.push_back(linf);
	}
	else
	{
		vector<string> HID;
		ReadID(HID, HCpath1);
		//cout<<HID.size()<<endl;	
		//hcinfor[0]=double(HID.size());
		for(int j=0; j<HID.size() ; j++)
		{
			//cout<<HID[j]<<endl;
			vector<double> linf;
			linf.resize(3, 0.0);
			string cseq, cmark, csse;
			int HFnum=0;
			double HCrad=0.0;
			string hcfile=HCpath+ID[i]+"/"+HID[j]+".inf";
			char inffile[500];
			strcpy(inffile, hcfile.c_str());				
			ReadInf(cseq, csse, cmark, HFnum, HCrad, inffile);
			//===Hydrophobic residues proportion
			int allnum=0;		
			for(int k=0; k<cmark.size(); k++)
			{
			if(cmark[k]=='1' || cmark[k]=='0')
			{
				allnum++;
				//if(cmark[j]=='1') HRnum++;
			}
			}
			linf[0]=HCrad;
			linf[1]=double(HFnum);
			linf[2]=double(allnum);
			hcinfor.push_back(linf);
			linf.clear();
		}
	}
	hfile.close();
	
	//cout<<hcinfor.size()<<endl;
	//for(int m=0; m<hcinfor.size(); m++)
	//{
		//cout<<"HC radius: "<<hcinfor[m][0]<<endl;
		//cout<<"HF number: "<<hcinfor[m][1]<<endl;
		//cout<<"HC seq: "<<hcinfor[m][2]<<endl;

	//}
	
	
	
	//===Contact matrix and surface residues
	SurfRes(que.seq, mod.Acc, mod.rs);
	//GenCAmatrix(mod.CAmatrix, mod.coord);
	mod.radius=RGvalue(mod.acoord);
	//===Secondary element states
	InitialSSE(mod.sse, mod.ss);
	
	//cout<<ID[i]<<endl;
	//cout<<que.seq<<'\n'<<que.ss<<'\n'<<que.rs<<endl;
	//cout<<mod.ss<<'\n'<<mod.rs<<endl;
	//cout<<mod.CAmatrix.size()<<'\n'<<que.CM8a.size()<<endl;
	//cout<<que.radius<<" "<<mod.radius<<endl;
	//cout<<que.sse.size()<<" "<<mod.sse.size()<<endl;
	//cout<<que.sse[1].type<<" "<<mod.sse[1].type<<endl;
	//cout<<que.sse[1].len<<" "<<mod.sse[1].len<<endl;
	
	//===Calculate the NCACO score
	ChainCoord chain;
	char pdbfile[500]; 
	strcpy(pdbfile, pdbpath.c_str()); 
	Read4AtomPdb(chain, pdbfile); 
	AddCentroid(chain, CENTABLE); 
	int colliedtnum=0;
	double conscore=0.0, trpscore=0.0, solvescore=0.0,trpscore1=0.0,trpscore2=0.0;
	conscore=EconOfFiveBead99_99(chain, VDWRADIUS, CONTABLE, colliedtnum);
	double sl=double(mod.ss.size());
	double conscore1=0.00114289*sl*sl-3.98118*sl+77.1005;
	conscore1=exp(log(2)*(-(2*(conscore-conscore1)/(conscore1+conscore))*(2*(conscore-conscore1)/(conscore1))));
	trpscore=Etrp(chain, TRPTABLE);
	/*
	solvescore=Esolve(chain, SOLVETABLE);
	for(int j=0; j<mod.ss.size(); j++)
	{
		
		randomShuffle(chain.resiname); 
		double lscore=Etrp(chain, TRPTABLE); 
		trpscore2 = trpscore2+lscore; 
	} 
	trpscore2 =trpscore2/mod.ss.size(); 
	trpscore2=exp(log(2)*(-(2*(trpscore-trpscore2)/(trpscore2+trpscore))*(2*(trpscore-trpscore2)/(trpscore2)))); 
	trpscore1=Etrp_n(chain, TRPTABLE);
	cout<<conscore<<" "<<trpscore<<" "<<solvescore<<" "<<colliedtnum<<endl;
	cout<<tmscore<<" "<<conscore<<" "<<conscore1<<" "<<trpscore<<" "<<trpscore1<<" "<<trpscore2<<" "<<mod.ss.size()<<endl;
	cout<<conscore<<" "<<mod.ss.size()<<endl;
	*/
	//vector<vector<double> > CENmatrix; 
	GenCENmatrix(chain, mod.CAmatrix);
	
	
	if(que.ss.size()==mod.ss.size())
	{
		//cout<<resSS1(que.ss, mod.ss)<<" "<<resBE1(que.rs, mod.rs)<<endl;
		
		vector<double> passa;
		
		passa.resize(14, 0.0);
		pairSSA(passa, que, mod);
		passa.resize(19, 0.0);
		cmFeature(passa, que.CM8a, mod.CAmatrix, 8.0, 14);
		passa.resize(24, 0.0);
		cmFeature(passa, que.CM12a, mod.CAmatrix, 12.0, 19);
		passa.resize(27, 0.0);
		slenFeature(passa, que, mod, 24);
		passa.resize(30, 0.0);
		sconFeature(passa, que, mod, 27);
		passa.resize(32, 0.0);
		radFeature(passa, que.radius, mod.radius);
		
		passa.resize(35, 0.0);
		HCFeature(passa, mod, hcinfor);
		passa.resize(37, 0.0);
		passa[35]=conscore1;
		passa[36]=trpscore;
		
		
		//cout<<que.sse.size()<<" "<<mod.sse.size()<<endl;
		ofstream outfile(outpath.c_str(), ios::app);
		if(!outfile)
		{
		cerr << "\nError: unable to open output file!\n";
		return -1;
		}
		outfile <<tmscore<<" ";
		for(int j=0; j<passa.size(); j++) outfile<<j+1<<":"<<passa[j]<<" ";
		outfile <<endl;
		outfile.close();
		//cout<<tmscore<<" ";
		//for(int j=0; j<passa.size(); j++) cout<<j+1<<":"<<passa[j]<<" ";
		//cout<<endl;
	}
	else 
	{
		cout<<"Error: "<<ID[i]<<endl;
		cout<<que.ss.size()<<endl;
		cout<<mod.ss.size()<<endl;
		
		vector<int> seqin;
		vector<int> resin;
		seqin.resize(que.seq.size(), -1); 
		//cout<<seqin.size()<<endl;
		ModelSequence(seqin, resin, mod.seq, pdbpath);
		//cout<<resin.size()<<endl;
		for(int k=0; k<resin.size(); k++)
		{
			seqin[resin[k]-1]=k;
		}
		
		OrderSequence(seqin, resin, que, mod);
		for(int k=0; k<seqin.size(); k++)
		{
			if(seqin[k] == -1)
			{
				cout<<"-";				
			}
			else cout<<que.seq[k];
		}cout<<endl;
		
		vector<double> passa;
		
		passa.resize(14, 0.0);
		pairSSA(passa, que, mod);
		passa.resize(19, 0.0);
		cmFeatureIn(seqin, passa, que.CM8a, mod.CAmatrix, 8.0, 14);
		passa.resize(24, 0.0);
		cmFeatureIn(seqin, passa, que.CM12a, mod.CAmatrix, 12.0, 19);
		passa.resize(27, 0.0);
		slenFeatureIn(resin, passa, que, mod, 24); 
		passa.resize(30, 0.0);
		sconFeatureIn(resin, passa, que, mod, 27); //cout<<"OK!"<<endl;
		passa.resize(32, 0.0);
		radFeature(passa, que.radius, mod.radius); 
		
		passa.resize(35, 0.0);
		HCFeature(passa, mod, hcinfor);
		passa.resize(37, 0.0);
		passa[35]=conscore1;
		passa[36]=trpscore;
		
		ofstream outfile(outpath.c_str(), ios::app);
		if(!outfile)
		{
		cerr << "\nError: unable to open output file!\n";
		return -1;
		}
		outfile <<tmscore<<" ";
		for(int j=0; j<passa.size(); j++) outfile<<j+1<<":"<<passa[j]<<" ";
		outfile <<endl;
		outfile.close();
	}
	
	}
	/*****************************************************/
	cout<< "\nProgram exit normally!\n\n";
      
	return 0;
}

/**********************************************************************
>> Translate residue names from 3-lette to 1-letter form
>> 
**********************************************************************/
char AAname3to1(string name3)
{
	char name1;
	if( name3=="ALA" )	name1='A';
	else if( name3=="ARG" ) name1='R';
	else if( name3=="ASN" ) name1='N';
	else if( (name3=="ASP")||(name3=="ASX") ) name1='D';
	else if( name3=="CYS" ) name1='C';
	else if( (name3=="GLN")||(name3=="GLX") ) name1='Q';
	else if( name3=="GLU" ) name1='E';
	else if( name3=="GLY" ) name1='G';
	else if( name3=="HSD" ) name1='H';
	else if( (name3=="HIS")||(name3=="HSE")||(name3=="HSD")||(name3=="HSP") ) name1='H';
	else if( name3=="ILE" ) name1='I';
	else if( name3=="LEU" ) name1='L';
	else if( name3=="LYS" ) name1='K';
	else if( name3=="MET" ) name1='M';
	else if( name3=="PHE" ) name1='F';
	else if( name3=="PRO" ) name1='P';
	else if( name3=="SER" ) name1='S';
	else if( name3=="THR" ) name1='T';
	else if( name3=="TRP" ) name1='W';
	else if( name3=="TYR" ) name1='Y';
	else if( name3=="VAL" ) name1='V';
	else
	{
		name1='?';
		cout << "\nAAname3to1()--Error: no matched residue name.\t" << name3 << endl;
	}
	
	return name1;
}


//==================================================================
//Get the residue string from Model
//==================================================================
void ModelSequence(vector<int>& ind, vector<int>& resn, string& seq, string name)
{
	
	ifstream fin(name.c_str());
	if( !fin.is_open() )
	{
		cout << "Error: can not open " << name << endl;
		//return -2;
	}

	char buf[200];
	bool MultiModel=false;
	while( !fin.eof() )
	{
		fin.getline(buf, 200);
		
		string sbuf(buf);
		
		if((sbuf.substr(0,5)=="MODEL")&&MultiModel) break;//MODEL 2 start
		else if(sbuf.substr(0,5)=="MODEL") MultiModel=true;
		
		if(strncmp(buf, "ATOM", 4)==0)
		{
			if(sbuf.substr(12,4)!=" CA ") continue;
		
			if((buf[16]!=' ')&&(buf[16]!='A')) continue;//retain A atom if the atom can not been determined
				
			string resid=sbuf.substr(17,3);
			
			seq.push_back(AAname3to1(resid));
			
			string num="";
			for(int i=22; i<26; i++)
			{
				if(buf[i] != ' ') num.push_back(buf[i]);
			}
			resn.push_back(atoi(num.c_str()) );
		}
	}
	fin.close();
	
	
}

//==================================================================
//Generate index between sequence and model
//==================================================================
void OrderSequence(vector<int>& seqin, vector<int>& resin, QUE& que, MOD& mod)
{
	int mark=0;
		int de1=1;
		int delta=0;
		if(que.seq.size()>mod.seq.size()) delta=seqin.size()-resin.size();
		else delta=resin.size()-seqin.size();
		for(int k=0; k<seqin.size(); k++)
		{		
				
			if((k != resin[mark]-de1) && (mark == 0) )
			{
				delta--;
				continue;
			}
			else
			{
				for(int l=mark; l<mark+1+delta; l++)
				{	
					int flag=0;			
					if( (resin[l]-de1)==k && que.seq[k]==mod.seq[l])
					{
						seqin[k]=l;
						mark=l+1;
						break;
					}
					else
					{
						flag++;
						de1=de1+flag;
						if( (resin[l]-de1)==k && que.seq[k]==mod.seq[l])
						{
							seqin[k]=l;
							mark=l+1;
							delta--;							
							break;
						}
					}
				}
			}
		}
		//test index
		//for(int k=0; k<seqin.size(); k++)
		//{
			//if(seqin[k] == -1)
			//{
				//cout<<"-";				
			//}
			//else cout<<que.seq[k];
		//}cout<<endl;
}




