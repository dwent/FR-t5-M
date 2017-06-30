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
				NAME:		SMITHWATERMAN.Hpp
				ABSTRACT:	THREADING THE DYNAMIC PROGRAMMING FUNCTIONS
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
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <list>
#include <iostream>
#include <vector>
#include <fstream>
#include <cstring>


#include "Matrices.hpp"
#include "Trace.hpp"
#include "TraceArrow.hpp"
#include "Options.hpp"
#include "GetSelfSM.hpp"
#include "FiveBeadEnergy.hpp"
#include "MathTool.hpp"
//#include "weight.h"

using namespace std;

//LOCAL DYNAMIC PROGRAMMING: SMITH-WATERMAN
int Align_Local(	AlignResult& alignresult,
		SeqStr& queueseq,
		SeqStr& targstr,
            int& max_i,
            int& max_j,
		double factor[8],
		const char* a_str,
		const char* b_str,
		const char* sm_filename,
		cost_t match,
		cost_t mismatch,
		cost_t indel_var,
		cost_t indel_const,
		bool logarithmic,
		bool local,
		bool quiet,
		bool sm_opt,
		bool self_sm_opt
	     );
//GLOBAL DYNAMIC PROGRAMMING: NEEDLEMAN-WUNSCH
int Align_Global(	AlignResult& alignresult,
		SeqStr& queueseq,
		SeqStr& targstr,
		double factor[8],
		const char* a_str,
		const char* b_str,
		const char* sm_filename,
		cost_t match,
		cost_t mismatch,
		cost_t indel_var,
		cost_t indel_const,
		bool logarithmic,
		bool local,
		bool quiet,
		bool sm_opt,
		bool self_sm_opt
	     );
// f0:1 f1:1.2 f2:-1.19 f3:0 f4:-2.3 f5:0.8 con:-9 val:-2 score:	0.658474
// f0:1 f1:0 f2:0 f3:0 f4:0.7 f5:0 con:-6.7 val:-0.65 score:	0.615288
// f0:1 f1:1 f2:0 f3:0 f4:-2.5 f5:0 con:-7.25 val:-2 score:	0.649101
//ZHANGYANG GLOBAL DYNAMIC PROGRAMMING NEEDLEMAN-WUNSCH
int Align_Global_Zhang(	AlignResult& alignresult,
		SeqStr& queueseq,
		SeqStr& targstr,
		double factor[8],
		const char* a_str,
		const char* b_str,
		const char* sm_filename,
		cost_t match,
		cost_t mismatch,
		cost_t indel_var,
		cost_t indel_const,
		bool logarithmic,
		bool local,
		bool quiet,
		bool sm_opt,
		bool self_sm_opt
	     );
//ANOTHER TYPE OF NEEDLEMAN-WUNSCH
cost_t NwAlign3(RV2 &score,cost_t gap_open,cost_t gap_extn,
                        std::string seq1,std::string seq2,
                        cost_t gpc1,cost_t gpc2,
                        std::string dsse1,std::string dsse2,
                        std::string &seq1_al,std::string &seq2_al,
                        RV1 &con1,RV1 &con2);
//GLOBAL-LOCAL ALIGNMENT
int Align_Global_Local(	AlignResult& alignresult,
		SeqStr& queueseq,
		SeqStr& targstr,
		double factor[8],
		const char* a_str,
		const char* b_str,
		const char* sm_filename,
		cost_t match,
		cost_t mismatch,
		cost_t indel_var,
		cost_t indel_const,
		bool logarithmic,
		bool local,
		bool quiet,
		bool sm_opt,
		bool self_sm_opt
	     );
//LOCAL-GLOBAL ALIGNMENT
int Align_Local_Global(	AlignResult& alignresult,
		SeqStr& queueseq,
		SeqStr& targstr,
		double factor[8],
		const char* a_str,
		const char* b_str,
		const char* sm_filename,
		cost_t match,
		cost_t mismatch,
		cost_t indel_var,
		cost_t indel_const,
		bool logarithmic,
		bool local,
		bool quiet,
		bool sm_opt,
		bool self_sm_opt
	     );

// reference: assessing the performance of fold recognition methods by means of a compehensive benchmark
//The "global-local" alignment algorithm doesnot penalize unmatched N- or C-termini segments in the probe sequence(as in the local alignment),
// but does penalize any gaps in the target strucutre(as in the global aligment with ends penalization)

// Attention , i did not consider the "local-global" algorithm , because it penalize any aligned amino acids from the sequences. thus,its 
// application is limited to special cases.
//READ IDS OF THE TEMPLATE LIBRARY
int ReadTemplateID(vector<string>& allID, const char* IDfile);
//INITIALIZE THE BLOCK STRUCT ACCORDING TO THE SEQUENCE
void InitialBlock(string gapseq,Block & b);
//================================================================================
// Record the position of all matched pairs in the query sequence
//	and the template sequence;
//================================================================================
void GetAlignedPair(vector<vector<int > >& alignedpair, string gapqseq, string gaptseq);
//=========================================================================
// Get the aligned score from the score-matrix according to
// 	the index information in alignedpair[][]
//=========================================================================
double GetAlignedScore(vector<double>& alignedscore, vector<vector<int > >& alignedpair, vector<vector<double > >& selfsm);
//PART LENGTH OF ALIGNMENT USED FOR NORMALIZATION
void L_partitial( string gapqseq,string gaptseq,AlignResult& alignresult);
//NORMALIZE THE RAW SCORE BOTH BY LFULL AND LPART
void SW_L_partitial_full( string gapqseq,string gaptseq,AlignResult& alignresult);
//MAKE THE FULL ALIGNED SEQUENCE PAIR
void SW_full_Align(int max_i,int max_j, AlignResult& alignresult);
//get 	double rawscore;//like Zhang Muster
//	int L_full;//like Zhang Muster
//	int L_partitial;//like Zhang Muster
//=========================================================================
// Calculate all the Z-scores of query-template alignments
// Z(V) = (Z-<Z>)/(sqrt(<Z^2>-<Z>^2));
//=========================================================================
void CalNormScore(vector<AlignResult>& allalignresult);
//=========================================================================
// Calculate all the Z-scores of query-template alignments
// Z(V) = (Z-<Z>)/(sqrt(<Z^2>-<Z>^2));
//=========================================================================
void CalNormScoreVector(vector<double>allscore,vector<double>& allnormscore);
//=========================================================================
// Select the top selectN number of query-template alignments
// 	according to the normscore value;
//=========================================================================
void TopNIndex(vector<int>& TopIndex,vector<double>allnormscore,int selectN);
//=========================================================================
// Select the top selectN number of query-template alignments
// 	according to the allalignresult[i].normscore value;
//=========================================================================
void OptTemplate(vector<int>& TopIndex, vector<AlignResult>& TopAlignResult, vector<AlignResult>& allalignresult, int selectN);
