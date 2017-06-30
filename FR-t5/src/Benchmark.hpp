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
				NAME:		BENCHMARK.Hpp
				ABSTRACT:	THREADING BENCHMARK
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
//HEAD
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <list>
#include <iostream>
#include <vector>
#include <fstream>
#include <cstring>
#include "Smith-Waterman.hpp"
//NAMESPACE
//ASK: WHY USE THIS IN ALL FILES???????WHY USE STD ALL THE TIME???BY MIAO
using namespace std;
//READ ALL THE INFORMATION OF PROSUP DATASET
int ProsupInfor
	(	
		vector<string>& ID1,
		vector<string>& ID2,
		vector<int>& POS1,
		vector<int>&POS2,
		vector<vector<PosPair > >& allstdalign
	
	);
//READ THE PAIR INFORMATION OF PROSUP DATASET
void ProsupPairPDB(vector<string> ID1,vector<string> ID2,vector<ChainCoord>& allchain1,vector<ChainCoord>& allchain2);
//READ SEQUENCE AND STRUCTURE INFORMATION OF PROSUP DATASET
void ProsupSeq_Str_Info(vector<string> ID1,vector<string> ID2,vector<ChainCoord> allchain2,vector<SeqStr>& allqueueseq,vector<SeqStr>& alltargstr);
//GLOBAL TRAINING FUNCTION FOR THE PROSUP DATASET
void GlobalTrainning(vector<string> ID1,vector<string> ID2,vector<ChainCoord> allchain1,vector<ChainCoord> allchain2,vector<int> POS1,
	vector<int> POS2,vector<vector<PosPair > >allstdalign,vector<SeqStr> allqueueseq,vector<SeqStr> alltargstr);
//LOCAL TRAINING FUNCTION FOR THE PROSUP DATASE
void LocalTrainning(vector<string> ID1,vector<string> ID2,vector<ChainCoord> allchain1,vector<ChainCoord> allchain2,vector<int> POS1,
	vector<int> POS2,vector<vector<PosPair > >allstdalign,vector<SeqStr> allqueueseq,vector<SeqStr> alltargstr);
//THE MAIN FUNCTION FOR TRAINING
int Prosup_Trainning_Stage();

//READ THE PAIR INFORMATION OF LINDAHL DATASET
void Read_LindahlPair_Info(const string file, string ID, vector<string> &info);
//READ LINDAHL INFORMATION OF FAMILY FOLD AND SUPERFAMILY
int LindahlInfor(string ID);
//== Read all information of structure templates in Folds-Library  //976 proteins
int Read_Lindahl_Sets(vector<SeqStr> &alltargstr);
//TEST THE LINDAHL SET
int Lindahl_Testing(string PredID,string OutPath);

//FUNCTIONS BELOW ARE USED FOR DYNAMIC PROGRAMMING
//ZHANG FOR ZHANGYANG'S ALGORITHM
//NW FOR NEEDLEMAN-WUNSCH
//SW FOR SMITH-WATERMAN
//GL FOR GLOBAL-LOCAL
//LG FOR LOCAL-GLOBAL
//F0 F1 F01 F025 F0125 FOR DIFFERENT SCORING TERMS
//ASK F0125 ONLY CONSISTS 4 SCORING TERMS, WHERE IS THE FIFTH TERM??? BY MIAO
void Threading_1v1(string Agorithm,string Terms,SeqStr queueseq,SeqStr targstr,AlignResult& alignresult,const char* a_str,const char* b_str);
void Threading_1v1_sw(SeqStr queueseq,SeqStr targstr,AlignResult& alignresult,const char* a_str,const char* b_str);
void Threading_1v1_nw(SeqStr queueseq,SeqStr targstr,AlignResult& alignresult,const char* a_str,const char* b_str);
void Threading_1v1_zhang(SeqStr queueseq,SeqStr targstr,AlignResult& alignresult,const char* a_str,const char* b_str);
//SW FOR SMITH-WATERMAN
void Threading_1v1_sw_f0(SeqStr queueseq,SeqStr targstr,AlignResult& alignresult,const char* a_str,const char* b_str);
void Threading_1v1_sw_f1(SeqStr queueseq,SeqStr targstr,AlignResult& alignresult,const char* a_str,const char* b_str);
void Threading_1v1_sw_f01(SeqStr queueseq,SeqStr targstr,AlignResult& alignresult,const char* a_str,const char* b_str);
void Threading_1v1_sw_f025(SeqStr queueseq,SeqStr targstr,AlignResult& alignresult,const char* a_str,const char* b_str);
void Threading_1v1_sw_f0125(SeqStr queueseq,SeqStr targstr,AlignResult& alignresult,const char* a_str,const char* b_str);
//ZHANG FOR ZHANGYANG'S ALGORITHM
void Threading_1v1_zhang_f0(SeqStr queueseq,SeqStr targstr,AlignResult& alignresult,const char* a_str,const char* b_str);
void Threading_1v1_zhang_f1(SeqStr queueseq,SeqStr targstr,AlignResult& alignresult,const char* a_str,const char* b_str);
void Threading_1v1_zhang_f01(SeqStr queueseq,SeqStr targstr,AlignResult& alignresult,const char* a_str,const char* b_str);
void Threading_1v1_zhang_f012(SeqStr queueseq,SeqStr targstr,AlignResult& alignresult,const char* a_str,const char* b_str);
void Threading_1v1_zhang_f025(SeqStr queueseq,SeqStr targstr,AlignResult& alignresult,const char* a_str,const char* b_str);
void Threading_1v1_zhang_f0125(SeqStr queueseq,SeqStr targstr,AlignResult& alignresult,const char* a_str,const char* b_str);
//NW FOR NEEDLEMAN-WUNSCH
void Threading_1v1_nw_f0(SeqStr queueseq,SeqStr targstr,AlignResult& alignresult,const char* a_str,const char* b_str);
void Threading_1v1_nw_f1(SeqStr queueseq,SeqStr targstr,AlignResult& alignresult,const char* a_str,const char* b_str);
void Threading_1v1_nw_f01(SeqStr queueseq,SeqStr targstr,AlignResult& alignresult,const char* a_str,const char* b_str);
void Threading_1v1_nw_f025(SeqStr queueseq,SeqStr targstr,AlignResult& alignresult,const char* a_str,const char* b_str);
void Threading_1v1_nw_f0125(SeqStr queueseq,SeqStr targstr,AlignResult& alignresult,const char* a_str,const char* b_str);
//GL FOR GLOBAL-LOCAL
void Threading_1v1_gl_f0(SeqStr queueseq,SeqStr targstr,AlignResult& alignresult,const char* a_str,const char* b_str);
void Threading_1v1_gl_f1(SeqStr queueseq,SeqStr targstr,AlignResult& alignresult,const char* a_str,const char* b_str);
void Threading_1v1_gl_f01(SeqStr queueseq,SeqStr targstr,AlignResult& alignresult,const char* a_str,const char* b_str);
void Threading_1v1_gl_f025(SeqStr queueseq,SeqStr targstr,AlignResult& alignresult,const char* a_str,const char* b_str);
void Threading_1v1_gl_f0125(SeqStr queueseq,SeqStr targstr,AlignResult& alignresult,const char* a_str,const char* b_str);
//LG FOR LOCAL-GLOBAL
void Threading_1v1_lg_f0(SeqStr queueseq,SeqStr targstr,AlignResult& alignresult,const char* a_str,const char* b_str);
void Threading_1v1_lg_f1(SeqStr queueseq,SeqStr targstr,AlignResult& alignresult,const char* a_str,const char* b_str);
void Threading_1v1_lg_f01(SeqStr queueseq,SeqStr targstr,AlignResult& alignresult,const char* a_str,const char* b_str);
void Threading_1v1_lg_f025(SeqStr queueseq,SeqStr targstr,AlignResult& alignresult,const char* a_str,const char* b_str);
void Threading_1v1_lg_f0125(SeqStr queueseq,SeqStr targstr,AlignResult& alignresult,const char* a_str,const char* b_str);
//READ THE UFF PARAMETER TABLE 5400 REAL NUMBERS
void ReadUffTab(RV1 &mat);
