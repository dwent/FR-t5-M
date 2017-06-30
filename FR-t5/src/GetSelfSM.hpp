//===================================================================================
//
// Get N*M score matrix defined by oneself
//
// AUTHORS:     Wu Aiping
//
// DESCIRPTION: This score matrix is asymmetrically, N is the length
//		   of sequence 1 and M is the length of sequence 2.
//		   (Aug, 2007)
////===================================================================================
//
// Implementation of Threading Sequence Alignment
//
// AUTHORS:     Huyun
//
// DESCIRPTION: Training and Testing Datasets
//               LOCAL/GlOBAL alignment. (Feb, 2009)
//===================================================================================
//===================================================================================


#ifndef GETSELFSM_H
#define GETSELFSM_H


#include <vector>
#include <string>
#include <iostream>
#include "FiveBeadEnergy.hpp"


using namespace std;



/*
int ReadTrpPar(DV4& trptable, const char* parfile);


//=======================================================
// Struct ChainCoord, generated for recording
//	each template pdb structures
//=======================================================

typedef struct _ChainCoord
{
	int resinum;
	int atomnum;

	double beginxyz[3];
	double endxyz[6];
	
	vector<vector<string > > atomname;
	vector<string> resiname;
	vector<int> resiserial;
	vector<vector<vector<double > > > xyz;

} ChainCoord;


int Read4AtomPdb(ChainCoord& chain, const char* pdbfile);
*/


//========================================================
// Struct SeqStr, generated for queue sequence
//	and each template structure
//========================================================
typedef struct SeqStr_
{
	int resinum;
	
	string sequence;
	string ID;
	
	vector<vector<double > > prof;
	vector<vector<double > > freq;
	vector<char> second;
	vector<vector<double > > strfreq;
	
	vector<int> resiindex;
	vector<vector<int > > angleindex;
	
} SeqStr;



//========================================================
// Struct AlignResult, generated for storing
//	the result of alignment.
//========================================================
//========================================================
// Struct AlignResult, generated for storing
//	the result of alignment.
//========================================================
/*typedef struct AlignResult_  //Global
{
	string qseq;
	string tseq;
	
	string gapqseq;
	string gaptseq;
	
	double maxscore;
	
} AlignResult;


typedef struct AlignResult_  //Local
{
	string qseq;
	string tseq;
	
	string gapqseq;
	string gaptseq;
	
	int Nalign;
	vector<double> Vscore;
	vector<double> Vref;
	
	double normscore;
	
} AlignResult;
*/


typedef struct AlignResult_  
{
	string qseq;
	string tseq;
	
	string gapqseq;
	string gaptseq;
	
	int Nalign;
	vector<double> Vscore;
	vector<double> Vref;

	double Vs;
	double Vr;
	
	double normscore;
	double maxscore;//used only in Local align

	double rawscore;//like Zhang Muster
	int L_full;//like Zhang Muster
	int L_partitial;//like Zhang Muster
	int L_partitial2;// as following

//GAPVPVDENDEGLQRALQFAMAEYNRA--------------------S---------------NDKYSSRVVRVISAKRQLVSGIKYILQVEIGRTTCPKSSGDLQSCEFHDEPEMAKYTTCTFVVYSIPWLNQIKLLESK--------CQ
//-------GDKPIWEQIGSSFIQHYYQLFDNDRTQLGAIYIDASCLTWEGQQFQGKAAIVEKLSSLPFQKIQHSITAQDHQPTPDSCIISMVVGQLKA-----------------DEDPIMGFHQMFLLKNINDAWVCTNDMFRLALHNF--
//-----CQ  转化为
//GAPVPVDENDEGLQRALQFAMAEYNRA--------------------S---------------NDKYSSRVVRVISAKRQLVSGIKYILQVEIGRTTCPKSSGDLQSCEFHDEPEMAKYTTCTFVVYSIPWLNQIKLLESKCQ--------
//-------GDKPIWEQIGSSFIQHYYQLFDNDRTQLGAIYIDASCLTWEGQQFQGKAAIVEKLSSLPFQKIQHSITAQDHQPTPDSCIISMVVGQLKA-----------------DEDPIMGFHQMFLLKNINDAWVCTNDM--FRLALHNF
	double ZhNormScore1;//like Zhang Muster
	double ZhNormScore2;//like Zhang Muster
} AlignResult;

typedef struct Block_    //store the trim information
{
	vector<int> initpos;
	vector<int> endpos;
	vector<int> length;
	vector<int> state;		// state[i]=0 means aligned-segment, state[i]=-1 means unaligned-segment.eg. loop region;
	
} Block;

void PSSMofSeqStr(SeqStr& seqstr, const char* pssmfile);
void PSIPREDofSeqStr(SeqStr& seqstr, const char* secondfile);

void SubPSSMofSeqStr(vector<int> Subsite,SeqStr& seqstr, const char* pssmfile);
void ReadLibID(char* file,SV1 &ss);
void SubPSIPREDofSeqStr(vector<int> Subsite,SeqStr& seqstr, const char* secondfile);

void DSSPofSeqStr(SeqStr& seqstr, const char* secondfile);
void STRFREQofSeqStr(SeqStr& seqstr, const char* strfreqfile);
void TRPofSeqStr(SeqStr& seqstr, ChainCoord& chain);

void ReadNormFasta(const char* fastafile, char* sequence);
int ReadPSSM(string& sequence, vector<vector<double > >& freq, vector<vector<double > >& prof, const char* pssmfile);
int ReadPSIPRED(vector<char>& second, const char* secondfile);
int ReadDSSP(vector<char>& second, const char* secondfile);
int TransAA2Index(vector<int>& resiindex, string sequence);
int ReadStrFreq(vector<vector<double > >& strfreq, const char* strfreqfile);

int Get4AngleIndex(vector<vector<int > >& angleindex, ChainCoord& chain);
//int Get4Angle(double angle[4], double xyz[12][3]);

int GetSelfSM(vector<vector<double > >& selfsm, SeqStr& queueseq, SeqStr& targstr, double factor[8]);


double Eseqseq(vector<double>& freq, vector<double>& prof);
double E2nd(char origsecond, char targsecond);
double Etrpbin(int qpos, int tpos, vector<int>& resiindex, vector<vector<int > >& angleindex);
double Eseqstr(vector<double>& prof, vector<double>& strfreq);
double HCA(char origseq, char targseq);



typedef struct PosPair_
{
	int id1, id2;
	char rname1, rname2;
	
} PosPair;


int ReadProSupFile(vector<PosPair>& stdalign, const char* prosupfile);
int ReadProSupFilePro(vector<PosPair>& stdalign, const char* prosupfile);
int CompAlign_Local(AlignResult& alignresult, vector<PosPair>& stdalign, int initpos1, int initpos2,int max_i,int max_j,
		vector<int>& rserial1, vector<int>& rserial2);

int CompAlign_Global(AlignResult& alignresult, vector<PosPair>& stdalign, int initpos1, int initpos2,
		vector<int>& rserial1, vector<int>& rserial2);













#endif
