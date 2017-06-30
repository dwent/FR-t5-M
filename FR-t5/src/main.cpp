//===================================================================================
//
// Implementation of Threading Sequence Alignment
//
// AUTHORS:     Huyun
//
// DESCIRPTION: Training and Testing Datasets
//               LOCAL/GlOBAL alignment. (Feb, 2009)
//===================================================================================
				/**********************************************************
				Copyright(c)2012, IBP, CHINESE ACADEMY OF SCIENCE
				All rights reserved.
				NAME:		main.cpp
				ABSTRACT:	THREADING MAIN
				VERSION:	0.01
				AUTHOR:	Miao Zhichao
				DATE:		2012.03.19 Mon
				CONTACT:	chichaumiau@gmail.com

				NOTICE: This is free software and the source code is freely
				available. You are free to redistribute or modify under the
				conditions that (1) this notice is not removed or modified
				in any way and (2) any modified versions of the program are
				also available for free.
						** Absolutely no Warranty **
				***********************************************************/

#include <cstdlib>
#include <cstdio>
#include <list>
#include <cstring>

#include "Benchmark.hpp"
#include "Smith-Waterman.hpp"
//#include "weight.h"
#include "FiveBeadEnergy.hpp"
#define EXCLUDE
//HUNYUN PARAMETERS
//THE DIR WHERE PUT THE PARAMETER FILES
string pararoot="/defusion/data/taijiao/jtj03/miaozhichao/casp9/FR-t5/parameter/";
//THE ID FILE
const char *IDfile ="/defusion/data/taijiao/jtj03/miaozhichao/DATABASE/RATHREADER_package/ID.lib";
//TEMPLATE PDB DIR
const char pdbdir[200]="/defusion/data/taijiao/jtj03/miaozhichao/DATABASE/RATHREADER_package/ent/";
//TEMPLATE SEQ PROFILE DIR
const char* profdir="/defusion/data/taijiao/jtj03/miaozhichao/DATABASE/RATHREADER_package/temp/psi/";	
//TEMPLATE DSSP DIR
const char* dsspdir="/defusion/data/taijiao/jtj03/miaozhichao/DATABASE/RATHREADER_package/temp/dssp/";
//ERROR ID IN F0125 CALCULATION
const char *ID_f0125_erro ="/defusion/data/taijiao/jtj01/Database/PDBLib/LibCheck/f0125.wrong.956.id";
//LESS THEN 40 AA TEMPLATE ID
const char *ID_less40AA ="/defusion/data/taijiao/jtj01/Database/PDBLib/less40AA.id";

DV2 VDWRADIUS;//HUNYUN PARAMETERS
DV2 CONTABLE;//HUNYUN PARAMETERS
DV4 TRPTABLE;//HUNYUN PARAMETERS
DV6 HBTABLE;//HUNYUN PARAMETERS
DV4 CENTABLE;//HUNYUN PARAMETERS

//=============================================================================
int Read_Pdb_Lib(string Terms,vector<string> allID,vector<SeqStr> &alltargstr);//HUYUN FUNCTION
void Sub_full_Align(string fullseq,vector<int> sub, AlignResult& alignresult);//HUYUN FUNCTION
void GetSubSite(string sub,vector<int> &Subsite);//HUYUN FUNCTION
void RemovePdb_f0125_erro(vector<string> &allID);//HUYUN FUNCTION
void RemovePdb_less40AA(vector<string> &allID);//HUYUN FUNCTION
//
//=============================================================================
#ifdef HUYUN
//HUYUN MAIN FUNCTION OF THREADING
/*
 * USAGE: BIN [QUERY ID] [QUERY PATH] 
 */ 
int main(int argc, const char *argv[])
{
	//clock_t s;double d;s=clock();
	//cout << "usage:" <<endl;
	//DEFINE ALGORITHM AND SCORING TERMS
	string Agorithm ="A1";
	string Terms ="f0125";
	string PredID =argv[1];
	string QueryPath=argv[2];
	string QL ="full";
	//DEFINE THE PARAMETER FILES
	string vdwfile=pararoot+"vdwmatrix99_99_0";
	string confile=pararoot+"EconTable99_99";
	string trpfile=pararoot+"EtrpTable";
	string hbfile=pararoot+"EhbTable_2resi";
	string centroidINFOfile=pararoot+"centroidrotamer";
	//READ PARAMETERS INTO MULTI-DIM VECTORS
	ReadAllPar(VDWRADIUS, vdwfile.c_str(), 
	           CONTABLE, confile.c_str(), 
	           TRPTABLE, trpfile.c_str(), 
	           HBTABLE, hbfile.c_str(), 
	           CENTABLE, centroidINFOfile.c_str());
	//FINISH
	cout << "Read AllPair In!\n\n";
	
	//READ ALL TEMPLATE IDS
	vector<string> allID;
      
	ReadTemplateID(allID, IDfile); 
	int templatenum=allID.size();
      cout << "the number of the template orginal: " << templatenum << endl;
	//READ ID EXCLUSION FOR LESS THEN 40 AA 
	RemovePdb_less40AA(allID);	
	templatenum=allID.size();
      cout << "the number of the template which more than 40 AA: " << templatenum << endl;
 	//READ ID EXCLUSION FOR ERROR ID IN CALCULATING SCORING TERM 0125 
 	RemovePdb_f0125_erro(allID);
	templatenum=allID.size();
      cout << "the number of the template after checking erro: " << templatenum << endl;

///////////////////////////////////////////////////////////////////////////////

	//SET THE OUTPUT FILES: ALIGNMENT FILE
	string OutPath=QueryPath;
	string outname1 =OutPath +PredID +".align";
	ofstream outfile1(outname1.c_str(), ios::out | ios::trunc);
	if(!outfile1)
	{
		cerr << "\nError: unable to open output file1!\n";
	}
	//SET THE OUTPUT FILES: FULL RANK FILE
	string  outname2 =OutPath +PredID+".Lf.rank";
	ofstream outfile2(outname2.c_str(), ios::out | ios::trunc);
	if(!outfile2)
	{
		cerr << "\nError: unable to open output file2!\n";
	}
////////////////////////////////////////////////////////////////////////////////
	//SET THE OUTPUT FILES: PART RANK FILE
	string  outname3 =OutPath +PredID+".Lp.rank";
	ofstream outfile3(outname3.c_str(), ios::out | ios::trunc);
	if(!outfile3)
	{
		cerr << "\nError: unable to open output file3!\n";
	}
////////////////////////////////////////////////////////////////////////////////


	//INPUT FILES IN THE QUERY PATH
	//FASTA FILE
	//ASK: WHY NEED A FASTA FILE??????BY MIAO
	string fastapath= QueryPath+PredID+".fasta";
	//PSI FILE
	string pssmpath= QueryPath+PredID+".psi";
	//PSIPRED FILE
	string secondpath= QueryPath+PredID+".horiz";

	char sequence[2000]="";	
	string fullseq ="";
	//== Read Profile and 2nd-structure infomation of query sequence
	SeqStr queueseq;
	vector<int> Subsite;
	if(QL =="full")
	{
		//ASK: WHY NEED A FASTA FILE??????BY MIAO
		ReadNormFasta(fastapath.c_str(), sequence);
		PSSMofSeqStr(queueseq, pssmpath.c_str());
		PSIPREDofSeqStr(queueseq, secondpath.c_str());
	}
	else
	{

		char Fullsequence[2000]="";	
		ReadNormFasta(fastapath.c_str(),Fullsequence);
//////////////////////////////////////
		for(int i=0; i<strlen(Fullsequence); i++)
		{
			fullseq +=Fullsequence[i];
		}

////////////////////////////////////////////
		GetSubSite(QL,Subsite);
		SubPSSMofSeqStr(Subsite,queueseq, pssmpath.c_str());
		SubPSIPREDofSeqStr(Subsite,queueseq, secondpath.c_str());

		int len = Subsite.size() ;
		for(int i=0; i<len; i++)
		{
			int index =Subsite[i];	
			sequence[i] = Fullsequence[index];
		}
		sequence[len]='\0';;

		//PSSMofSeqStr(queueseq, pssmpath.c_str());
		//PSIPREDofSeqStr(queueseq, secondpath.c_str());
	}
	//TRANSFORM AA TYPES TO INDICES
	TransAA2Index(queueseq.resiindex, queueseq.sequence);


	/*for(int i=0; i<strlen(sequence); i++)
	{
		cout << queueseq.second[i];
	}
	cout << endl;*/

      vector<AlignResult> AllResult;
      vector<double>allscore1;      vector<double>allscore2;

	vector<SeqStr> alltargstr;
	//READ TEMPLATE LIBRARY
	Read_Pdb_Lib(Terms,allID,alltargstr);
	//for(int i=0; i<50; i++)
	//THE MAIN THREADING PART
	for(int i=0; i<templatenum; i++)
	{
		AlignResult alignresult;
           	cout <<"Threading_1v1:" << i << "\t" << alltargstr[i].ID << endl; 
		
		char a_str[5000];
		char b_str[5000];
		strcpy(a_str, sequence);
		strcpy(b_str, alltargstr[i].sequence.c_str());
		//CHOOSE THE DP ALGORITHM FOR THREADING
		Threading_1v1(Agorithm,Terms,queueseq,alltargstr[i],alignresult,a_str,b_str);
		//Threading_1v1_zhang_f0125(queueseq,alltargstr[i],alignresult,a_str,b_str);

		AllResult.push_back(alignresult);	

		allscore1.push_back(alignresult.ZhNormScore1);
		allscore2.push_back(alignresult.ZhNormScore2);
	}

       vector<double>allZscore1; vector<double>allZscore2;
       CalNormScoreVector(allscore1,allZscore1); CalNormScoreVector(allscore2,allZscore2);
       allscore1.clear();  allscore2.clear();

	//SORT THE RESULTS TO RANK TOP ONES
	for(int i=0; i<templatenum; i++)
	{
		outfile1 << alltargstr[i].ID  << " Rs=" <<AllResult[i].rawscore<<" ZS1=" << allZscore1[i]<< " ZS2=" << allZscore2[i] <<endl;
		//outfile1  <<" q.L=" << AllResult[i].qseq.size() <<" gapq.L=" << AllResult[i].gapqseq.size()<< " Na=" <<AllResult[i].Nalign ;
		//outfile1  << " Rs=" <<AllResult[i].rawscore<< " Vs=" <<AllResult[i].Vs << " Vr=" <<AllResult[i].Vr;
		//outfile1  << " Ql=" <<AllResult[i].qseq.size() << " Tl=" <<AllResult[i].tseq.size() << endl;
		//outfile1  << " L_f=" <<AllResult[i].L_full << " L_p=" <<AllResult[i].L_partitial << " L_p2=" <<AllResult[i].L_partitial2 << endl;
		if(QL =="full")
		{	
			outfile1 <<"Q::" << AllResult[i].gapqseq << endl;
			outfile1 <<"T::" << AllResult[i].gaptseq << endl;
		}
		else
		{
			Sub_full_Align(fullseq,Subsite,AllResult[i]);
			outfile1 <<"Q::" << AllResult[i].gapqseq << endl;
			outfile1 <<"T::" << AllResult[i].gaptseq << endl;
		}
		outfile1 << endl;
	}
	//FIND THE TOP ALIGNMENTS
	vector<int> TopIndex1;
	vector<int> TopIndex2;
	int TopN =50;
     	TopNIndex(TopIndex1,allZscore1,TopN);
     	TopNIndex(TopIndex2,allZscore2,TopN);
	//OUTPUT THE RESULTS OF TOP ALIGNMENTS
	for(int i=0; i<TopN; i++)
	{
		int ti1 =TopIndex1[i] ;	int ti2 =TopIndex2[i] ;
		outfile2 << "Top " << i+1 << ":" << ti1 << '\t' << alltargstr[ti1].ID<< '\t' << allZscore1[ti1] << endl;
		outfile3 << "Top " << i+1 << ":" << ti2 << '\t' << alltargstr[ti2].ID<< '\t' << allZscore2[ti2] << endl;		
	}



	//d=(double)(clock()-s)/CLOCKS_PER_SEC;printf("TIME_DURATION_MTHREAD: [%7.2f] sec\n",d);
	//FINISH ALL
	cout << "\nProgram exit normally!\n\n";
	
	return 0;

}
#endif
//== Read all information of structure templates in Folds-Library 
//THIS IS TIME-CONSUMING
int Read_Pdb_Lib(string Terms,vector<string> allID,vector<SeqStr> &alltargstr)
{
	//== Read all information of structure templates in Folds-Library 
   
	int templatenum=allID.size();
      cout << "the number of the template: " << templatenum << endl;

	char pdbfile[500];
	///DEFINE SUFFICES
	const char pdbsuffix[10]=".ent";
	const char* profsuffix=".psi";
	const char* dsspsuffix=".dssp";

	char targpssmfile[500];
	char targsecondfile[500];
   	for(int i=0; i<templatenum; i++)
	{
            cout << i << "\t" << allID[i] << endl; 
		SeqStr targstr;
		targstr.ID =allID[i];
		
		//strcpy(targpssmfile, rootpath);
		strcpy(targpssmfile, profdir);
		strcat(targpssmfile, allID[i].c_str());
		strcat(targpssmfile, profsuffix);
		//THE PSI FILES
		PSSMofSeqStr(targstr, targpssmfile);
		//strcpy(targsecondfile, rootpath);
		strcpy(targsecondfile, dsspdir);
		strcat(targsecondfile, allID[i].c_str());
		strcat(targsecondfile, dsspsuffix);    
		//THE DSSP FILES         
		DSSPofSeqStr(targstr, targsecondfile);
        //TRANSFORM AA TYPES INTO INDICES
		TransAA2Index(targstr.resiindex, targstr.sequence);
		//CHECK THE SCORING TERMS
 		if(Terms !="f0" && Terms !="f01"&& Terms !="f1" )
		{             
			ChainCoord chain;
			//strcpy(pdbfile, rootpath);
			strcpy(pdbfile, pdbdir);
			strcat(pdbfile, allID[i].c_str());
			strcat(pdbfile, pdbsuffix);
			//READ THE PDB FILES
			Read4AtomPdb(chain, pdbfile);
             //TRIPLETS FOR THE PDB FILE
			TRPofSeqStr(targstr, chain);   
		}         
		//STORE ALL THE INFORMATION    
		alltargstr.push_back(targstr);
	}


	return 0;
}
//HUYUN FUNCTION
//AS I KNOW, IT IS USED FOR DELETE THE TWO TERMINALS OF ALIGNMENT FOR THE QUERY SEQ
void Sub_full_Align(string fullseq,vector<int> sub,AlignResult& alignresult)
{
	//cout << "qseq ="<< alignresult.qseq<< endl;
	//cout << "fqseq="<<fullseq << endl;
	//cout << "tseq ="<< alignresult.tseq<< endl;

	//cout << "max_i ="<< max_i<< "\tmax_j ="<< max_j<< endl;
	string gapqseq =alignresult.gapqseq;
	string gaptseq =alignresult.gaptseq;	
	//cout << "gapqseq ="<< gapqseq<< endl;
	//cout << "gaptseq ="<< gaptseq<< endl;

	string QB ="";string QE ="";
	string TB ="";string TE ="";

	vector<vector<int > > b;vector<int> a;
	vector<int> c;	vector<int> state ;
	int len=fullseq.size();
	for(int i=0; i<len; i++)
	{
		bool flag =false;
		for(int j=0; j<sub.size(); j++)
		{
			if( i==sub[j]) {flag =true; break;}
		}		
		if(flag) {	a.push_back(1); }
		else{	a.push_back(0);}
	}
	vector<int> Vira;
	for(int i=0; i<=len; i++)
	{
		if(i<len) {Vira.push_back(a[i]);}
		else {Vira.push_back(a[len-1]);}
	}
	for(int i=0; i<len; i++)
	{
		cout << Vira[i] ;
		if(Vira[i+1] ==Vira[i] && i!=len-1){c.push_back(i);}
		if(Vira[i+1] ==Vira[i] && i==len-1){c.push_back(i);b.push_back(c);c.clear(); state.push_back(Vira[i]);}
		if(Vira[i+1] !=Vira[i] )
		{	
			{b.push_back(c);c.clear(); state.push_back(Vira[i]);}
		}
	}
		
	//cout <<"b.size()= "<< b.size()<< endl;
	for(int i=0; i<b.size(); i++)
	{
		cout << i << "\t" << state[i] << endl;
	}
//////下面只处理了sub.a1_b1的情况，而没有处理sub.a1_b1.a2_b2的情况

	for(int i=0; i<len; i++)
	{
		if(Vira[i] ==0)
		{	
			if(i<sub[0])  // sub没有包括N 端
			{
			 	QB += fullseq[i];
				TB += '-';
			}
			else if(i> sub[sub.size()-1]) // sub没有包括C 端
			{
				QE +=fullseq[i];
				TE +='-';
			}
		}
	}
	gapqseq = QB+gapqseq+QE;
	gaptseq = TB+gaptseq+TE;
	//cout << "gapqseq2 ="<< gapqseq<< endl;
	//cout << "gaptseq2 ="<< gaptseq<< endl;

	int num =0;	
	for(int i=0; i<gapqseq.size(); i++)
	{
		if(gapqseq[i] !='-') num++;
	}
	cout <<" num: "<< num << "\tlen: " << len << endl;
	if(num !=len)
	{
		cout << "In main() Sub_full_Align len !=num" <<endl;
		exit(-1);
	}
	alignresult.gapqseq = QB +alignresult.gapqseq +QE;
    	alignresult.gaptseq = TB +alignresult.gaptseq +TE;
}
/*
void Sub_full_Align1(vector<int> sub, AlignResult& alignresult)
{
	//cout << "qseq ="<< alignresult.qseq<< endl;
	//cout << "tseq ="<< alignresult.tseq<< endl;

	//cout << "max_i ="<< max_i<< "\tmax_j ="<< max_j<< endl;
	//cout << "gapqseq ="<< alignresult.gapqseq<< endl;
	//cout << "gaptseq ="<< alignresult.gaptseq<< endl;

	const char* chstr1=alignresult.gapqseq.c_str();
	const char* chstr2=alignresult.gaptseq.c_str();

	int L_gapqseq = alignresult.gapqseq.size();
      int L_gaptseq = alignresult.gaptseq.size();
      //cout << "L_gapqseq  " << L_gapqseq <<";L_gaptseq  " << L_gaptseq << endl;
    
      //caculate the gap numbers of the gapqseq and gaptseq
      int G_gapqseq = 0,G_gaptseq = 0;
      for(int j=0; j<L_gapqseq; j++)
	{
		if(chstr1[j]=='-') G_gapqseq++;
		if(chstr2[j]=='-') G_gaptseq++;    
        }
      //cout << "G_gapqseq  " <<  G_gapqseq <<";G_gaptseq  " << G_gaptseq << endl;

      //caculate the start rserial[] of the local aligment
      int gapqseq_start ,gaptseq_start ;
       //gapqseq_start = max_i-(L_gapqseq-G_gapqseq)+1;
       //gaptseq_start = max_j-(L_gaptseq-G_gaptseq)+1;
      gapqseq_start = max_i-(L_gapqseq-G_gapqseq);
      gaptseq_start = max_j-(L_gaptseq-G_gaptseq);
    	//cout << "max_i  " << max_i<<"  ,max_j  " << max_j <<endl; 
    	//cout << "gapqseq_start  " << gapqseq_start <<"  ,gaptseq_start  " << gaptseq_start <<endl; 

	string QB ="";string QE ="";
	string TB ="";string TE ="";
	for(int j=0; j< gaptseq_start; j++)
	{
		QB +='-';
		TB +=alignresult.tseq[j];
	}
	for(int j=0; j<gapqseq_start; j++)
	{
		QB +=alignresult.qseq[j];
		TB +='-';
	}
	int Q_L =alignresult.qseq.size();
	int T_L =alignresult.tseq.size();

	for(int j=max_i; j<Q_L; j++)
	{
		QE +=alignresult.qseq[j];
		TE +='-';
	}
	for(int j=max_j; j<T_L; j++)
	{
		QE +='-';
		TE +=alignresult.tseq[j];
	}


	alignresult.gapqseq = QB +alignresult.gapqseq +QE;
    	alignresult.gaptseq = TB +alignresult.gaptseq +TE;

}
*/
//HUYUN FUNCTION
//AS I KNOW, IT IS USED TO GET THE INDICES OF THE SUB STRING
void GetSubSite(string sub,vector<int> &Subsite)
{
	//sub.a1_b1(.a2_b2)
	vector<int> a;vector<int> b;
	vector<int> s;vector<int> e;
	for(int i=0; i<sub.size(); i++)
	{
		if(sub[i] == '.') a.push_back(i);
		if(sub[i] == '_') b.push_back(i);
	}
	a.push_back(sub.size());
	for(int i=0; i<b.size(); i++)
	{	
		string stemp1 ="";
		for(int j=a[i]+1; j<b[i]; j++)
		{
			stemp1 +=sub[j];
		}
		s.push_back(atoi(stemp1.c_str()));
	}	
	for(int i=0; i<b.size(); i++)
	{	
		string stemp2 ="";
		for(int j=b[i]+1; j<a[i+1]; j++)
		{
			stemp2 +=sub[j];
		}
		e.push_back(atoi(stemp2.c_str()));
	}

	for(int i=0; i<s.size(); i++)
	{
		for(int j=s[i]; j<=e[i]; j++)
		{
			Subsite.push_back(j);
			//cout << j << endl;
		}
	}
}
//FR-T5 THREADING CAN NOT CARRY OUT CALCULATION FOR ALL PDB STRUCTURES BECAUSE OF BAD CODING
//SO THIS FUNCTION IS USED TO ELIMINATE THE PDB IDS WHICH CAN NOT BE MEASURED BY FR-T5
void RemovePdb_f0125_erro(vector<string> &allID)
{

	
	vector<string> ID2;
	ReadTemplateID(ID2, ID_f0125_erro);
	
	vector<string> ID3;	
	cout << "in all 146726 pdbs f0125 erro id num=" << ID2.size() <<endl;
	long int num =0;
	for(long int i=0;i<allID.size();i++)
	{
		bool flag= true;
		for(long int j=0;j<ID2.size();j++)
		{
			if(allID[i].substr(0,4) ==ID2[j].substr(0,4) ) 
			{
				flag= false;
				break;

			}
		}
		if(flag)
		{
			ID3.push_back(allID[i]);
		}
	}
	allID.clear();
	for(long int i=0;i<ID3.size();i++)
	{
		allID.push_back(ID3[i]);
	}

}
//LESS THEN 40 AA TEMPALTES ARE TOO SHORT FOR THREADING
//SO EXCLUDE THESE IDS
void RemovePdb_less40AA(vector<string> &allID)
{
	
	vector<string> ID2;
	ReadTemplateID(ID2, ID_less40AA);

	vector<string> ID3;	
	cout << "in all 146726 pdbs less than 40 AA  num=" << ID2.size() <<endl;
	long int num =0;
	for(long int i=0;i<allID.size();i++)
	{
		bool flag= true;
		for(long int j=0;j<ID2.size();j++)
		{
			if(allID[i] ==ID2[j] ) 
			{
				flag= false;
				break;

			}
		}
		if(flag)
		{
			ID3.push_back(allID[i]);
		}
	}

	allID.clear();
	for(long int i=0;i<ID3.size();i++)
	{
		allID.push_back(ID3[i]);
	}

}

#ifdef MIAOZHICHAO
//MIAOZHICHAO MAIN FUNCTION OF THREADING
/*
 * EXECUTE [DSSPFILE] [SEQUENCE PROFILE] [PARAMETERFILE] [OUTPUT NUMBER] [OUTPUT FILE]
 */

int main(int argc, const char *argv[])
{
	if(argc<4)
	{
		printf("Usage: EXECUTE [DSSPFILE] [SEQUENCE PROFILE] [OUTPUT NUMBER] [OUTPUT FILE]\n");exit(0);
	}
	//READ PARAMETERS FORM THE PARAMETER FILE
	cost_t wt1,wt2,wt3,gp1,gp2;
	cost_t wtc0,wtc1,wtc2,wtc3,gpc1,gpc2;
	cost_t wt4[15];
	std::string progName(argv[0]);
	std::string buf,rootdir;
	wt1=1.59688;wt2=0.629216;wt3=1.29089;gp1=-8.05998;gp2=-0.366436;
	wtc0=0.0401365;wtc1=1.98704;wtc2=2.17441;gpc1=1.64107;gpc2=1.90402;wtc3=1.72;
	wt4[0]=-0.126808;wt4[1]=-0.0618381;wt4[2]=-0.0772241;wt4[3]=-0.0299057;
	wt4[4]=-0.0684414;wt4[5]=-0.321477;wt4[6]=0.93563;wt4[7]=0.0365641;
	wt4[8]=-0.172593;wt4[9]=-0.0249622;wt4[10]=-1.02738;wt4[11]=0.0134908;
	wt4[12]=-0.0866414;wt4[13]=-0.081371;wt4[14]=-0.00192084;

	int i,j,k,l,m,n;
	std::string file1,file2,file3,file4;
	buf=progName.substr(0,progName.rfind('/'));
	rootdir=buf.substr(0,buf.rfind('/'));
	if(buf==rootdir)
		rootdir="../";
	else
		rootdir+="/";
	file1=(std::string)rootdir+"data/ID.lib";
	file2=(std::string)rootdir+"data/dssp.lib";
	file3=(std::string)rootdir+"data/psi.lib";
	file4=(std::string)rootdir+"data/uff.lib";
#ifdef EXCLUDE
	std::string file5=(std::string)rootdir+"data/exception";//FOR TEST
	SV1 exID;ReadLibID((char*)file5.c_str(),exID);bool exflag;//FOR TEST
#endif
	SV1 allID;
	ReadLibID((char*)file1.c_str(),allID);
	int Tnum=allID.size();
	std::ifstream ifDsspLib(file2.c_str()),ifProfLib(file3.c_str()),ifUffLib(file4.c_str());
	if(!ifDsspLib){std::cerr<<"\nCan not open: "<<file2<<'\n';exit(0);}
	if(!ifProfLib){std::cerr<<"\nCan not open: "<<file3<<'\n';exit(0);}
	if(!ifUffLib){std::cerr<<"\nCan not open: "<<file4<<'\n';exit(0);}
	
	//READ THE QUERY (DSSP PROFILE PSIPRED)
	std::string seq1,seq2;
	std::string seq1_al,seq2_al;
	IV2 fary1,fary2,mary1,mary2;
	char* Pf=(char*)argv[2];
	ReadProfPssm(Pf,fary1,mary1,seq1);
	char* Df=(char*)argv[1];
	std::string dsse1=ReadHoriz(Df),dsse2;
	int ls1=seq1.length(),ls2;
	//READ ANOTHER PARAMETER
	RV2 uffmat;
	RV1 ufftab;ReadUffTab(ufftab);
	//SET VALUES FOR THE MATRIX
	std::map<char,int> mapaa;
	for(i=0;i<20;i++)
		mapaa[aatype[i]]=i;
	IV1 ss1,ss2;
	for(i=0;i<ls1;i++)
		ss1.push_back(mapaa[seq1[i]]);
	cost_t matrix1[20][20],matrix[20][20];
	for(i=0;i<20;i++)matrix1[i][i]=BLOSUM62[i][i]*1.0;
	for(i=0;i<20;i++)
	for(j=0;j<20;j++)
	{
		matrix1[i][j]=BLOSUM62[i][j]*1.0;
		matrix[i][j]=matrix1[i][j]/sqrt(matrix1[i][i]*matrix1[j][j]);
	}
	//SET THE CONSERVATION VALUES
	bool consflag=false;
	RV1 cons1(ls1,0),cons2;
	RV1 con1(ls1,0),con2;
	int lc=2;
	int ib,ie;
	if(findhomo(mary1))
	{
		consflag=true;
		for(j=0;j<20;j++)
		for(k=j;k<20;k++)
		{
			for(i=0;i<ls1;i++)
			cons1[i]+=mary1[i][j]*mary1[i][k]*matrix[j][k]/10000.0;
		}
		for(i=0;i<ls1;i++)
		{
			ib=-lc;if(i+ib<0)ib=-i;
			ie=lc;if(i+ie>=ls1)ie=ls1-i-1;
			for(j=ib;j<=ie;j++)con1[i]+=cons1[i+j];
			con1[i]/=(cost_t)(ie-ib+1);
		}
	}
	
	cost_t sc=0;
	RV1 rtmp;
	RV2 score;
	cost_t hcascore=0,uffscore=0,conscore=0;
	Rank atmp(i,sc,sc,seq1,seq2,sc);
	std::vector<Rank> result;
	std::string seqquery;
	cost_t sda,meana=0,meana2=0;
	cost_t sdb,meanb=0,meanb2=0;
	//THE THREADING PART
	for(l=0;l<Tnum;l++)
	{
		getline(ifDsspLib,dsse2);
		ReadProfPssm(ifProfLib,fary2,mary2,seq2);
		ReadUff(ifUffLib,uffmat,ss2);
#ifdef EXCLUDE
		exflag=false;
		for(i=0;i<(int)exID.size();i++)
		if(exID[i]==allID[l])
		{
			exflag=true;
			break;
		}
#endif
		ls2=seq2.length();
		cons2.assign(ls2,0);con2.assign(ls2,0);
		if(consflag&&findhomo(mary2))
		{
			for(j=0;j<20;j++)
			for(k=j;k<20;k++)
			for(i=0;i<ls2;i++)
				cons2[i]+=mary2[i][j]*mary2[i][k]*matrix[j][k]/10000.0;
			for(i=0;i<ls2;i++)
			{
				ib=-lc;if(i+ib<0)ib=-i;
				ie=lc;if(i+ie>=ls2)ie=ls2-i-1;
				for(j=ib;j<=ie;j++)con2[i]+=cons2[i+j];
				con2[i]/=(cost_t)(ie-ib+1);
			}
		}

		rtmp.assign(ls2,wt3);
		score.assign(ls1,rtmp);
		for(i=0;i<ls1;i++)
		for(j=0;j<ls2;j++)
		{
			//SCORING TERMS
			if(hca(seq1[i])&&hca(seq2[j]))hcascore=1;
			else if(seq1[i]=='P'&&seq2[j]=='P')hcascore=1;
			else if(seq1[i]==seq2[j]&&!hca(seq1[i]))hcascore=0.7;
			else hcascore=0;
			conscore=0;
			ib=-lc;if(i+ib<0)ib=-i;if(j+ib<0)ib=-j;
			ie=lc;if(i+ie>=ls1)ie=ls1-i-1;if(j+ie>=ls2)ie=ls2-j-1;
			for(k=ib;k<=ie;k++)conscore+=cons1[i+k]*cons2[j+k];
			conscore/=(cost_t)(ie-ib+1);
			for(k=0;k<20;k++)score[i][j]+=mary1[i][k]*fary2[j][k]*0.01*(1+wtc0*conscore);
			score[i][j]+=(dsse1[i]==dsse2[j]?1.0:-1.0)*wt1*(1+wtc1*conscore)
				+hcascore*wt2*(1+wtc2*conscore);	
			for(k=0;k<5;k++)
			{
				m=ss2[j]*1800+k*360+ss1[i]*18;
				uffscore=0;
				for(n=0;n<18;n++)
				uffscore+=uffmat[j*5+k][n]*ufftab[m+n];
				score[i][j]+=wt4[ss2[j]*5+k]*uffscore*(1+wtc3*conscore);
			}
		}
#ifdef EXCLUDE
		if(exflag)sc=-100;
		else
#endif
		//GET THE SCORE AND PRESERVE
		sc=NwAlign3(score,gp1,gp2,seq1,seq2,gpc1,gpc2,dsse1,dsse2,seq1_al,seq2_al,con1,con2);
		j=seq1_al.length();
		for(i=0;i<j;i++)if(seq1_al[i]!='-')break;
		j-=i;
		for(i=seq1_al.length()-1;i>0;i--)if(seq1_al[i]!='-')break;
		j-=seq1_al.length()-1-i;
		atmp.ix=l;
		sda=sc;
		sdb=sc/(cost_t)j;
		atmp.s1=sda;
		atmp.s2=sdb;
		meana+=sda;meana2+=sda*sda;
		meanb+=sdb;meanb2+=sdb*sdb;
		atmp.ss1=seq1_al;
		atmp.ss2=seq2_al;
		result.push_back(atmp);
	}
	ifDsspLib.close();ifProfLib.close();ifUffLib.close();
	fflush(stdout);
	printf("Finish!\n");
	//MEASURE THE Z-SCORE
	meana/=(cost_t)Tnum;meana2/=(cost_t)Tnum;
	meanb/=(cost_t)Tnum;meanb2/=(cost_t)Tnum;
	sda=sqrt(meana2-meana*meana);
	sdb=sqrt(meanb2-meanb*meanb);
	j=k=-1;uffscore=conscore=-1;
	for(i=0;i<Tnum;i++)
	{
		if(result[i].s1>uffscore)
		{
			j=i;
			uffscore=result[i].s1;
		}
		if(result[i].s2>conscore)
		{
			k=i;
			conscore=result[i].s2;
		}
	}
	if(seqid(result[j])>=seqid(result[k]))
	{
		for(i=0;i<Tnum;i++)result[i].z=(result[i].s1-meana)/sda;
	}
	else
	{
		for(i=0;i<Tnum;i++)result[i].z=(result[i].s2-meanb)/sdb;
	}
	if(Tnum>1)sort(result.begin(),result.end(),std::less<Rank>());
	//OUTPUT THE RESULTS
	std::ofstream fout(argv[4]);
	if(!fout){std::cerr<<"\nCan not open: "<<argv[4]<<'\n';exit(0);}
	for(i=0;i<atoi(argv[3]);i++)
	{
		fout<<'>'<<allID[result[i].ix]<<' '<<result[i].z<<endl;
		fout<<result[i].ss2<<endl<<result[i].ss1<<endl;
	}
	fout.close();
	return 0;
}
#endif
