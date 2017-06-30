/*
  Copyright (C) 1999 by 
      Rolf Backofen, Peter Clote, Sebastian Will and Wiley & Sons, Inc.  
  All Rights Reserved.
  
  Permission to use, copy, modify, and distribute this
  software and its documentation for NON-COMMERCIAL purposes
  and without fee is hereby granted provided that this
  copyright notice appears in all copies.
  
  THE AUTHORS AND PUBLISHER MAKE NO REPRESENTATIONS OR
  WARRANTIES ABOUT THE SUITABILITY OF THE SOFTWARE, EITHER
  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
  IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
  PARTICULAR PURPOSE, OR NON-INFRINGEMENT. THE AUTHORS
  AND PUBLISHER SHALL NOT BE LIABLE FOR ANY DAMAGES SUFFERED
  BY LICENSEE AS A RESULT OF USING, MODIFYING OR DISTRIBUTING
  THIS SOFTWARE OR ITS DERIVATIVES.
 */

//============================================================================
//
// Changed by Wu Aiping (wuaiping@moon.ibp.ac.cn), Aug, 2007
// CHANGED BY MIAO ZHICHAO (chichaumiau@gmail.com) MAR  2012
//============================================================================
				/**********************************************************
				Copyright(c)2012, IBP, CHINESE ACADEMY OF SCIENCE
				All rights reserved.
				NAME:		BACISDEF.h
				ABSTRACT:	THREADING BACIS DEFINITIONS
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

#ifndef BASICDEF_H
#define BASICDEF_H
//HEAD
#include <cmath>
#include <fstream>
#include <sstream>
#include <map>
#include <algorithm>
#include <ctime>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#define NDEBUG

typedef double cost_t;//A REAL TYPE

// don't compare two floating point values by ==,
// this may almost always lead to errors.
// instead use the function approx_equal
#define APPROX_EQUAL 0.00001

//bool approx_equal(cost_t x, const cost_t &y, cost_t ae=APPROX_EQUAL);
bool approx_equal(cost_t x, const cost_t &y, cost_t ae);

const int PRECISION=6;    // printlength of DistanceMatrix values
const int ROWFILL=  2;    // min distance between two successive row-values
const int NOFILL=   '0';  // this fillcharacter means that no-fill is used.

//4-DIM REAL TYPES (VECTORS)
typedef std::vector<cost_t> RV1;
typedef std::vector<RV1> RV2;
typedef std::vector<RV2> RV3;
typedef std::vector<RV3> RV4;
#ifndef TIME
const std::string aatype="ARNDCQEGHILKMFPSTWYV"; //20 AA TYPES
#endif

//2-DIM STRING TYPES
typedef std::vector<std::string> SV1;
typedef std::vector<SV1> SV2;
typedef std::vector<bool> BV1;

using std::cout;
using std::endl;

const cost_t INFTY = 10000000.0;  //A VERY TINY REAL NUMBER
const cost_t deg2rad=0.01745329252; //A REAL NUMBER FOR ANGLE CONVERSION
const cost_t MINUSINFTY = -INFTY; //A VERY TINY -REAL NUMBER

//A PAIR CLASS USED IN DYNAMIC PROGRAMMING
class Pair {
 public:
  Pair(int i, int j);
  int i,j;
};

//4-DIM INTEGER TYPES
typedef std::vector<int> IV1;
typedef std::vector<IV1> IV2;
typedef std::vector<IV2> IV3;
typedef std::vector<IV3> IV4;
//=====================================================================
//DEFINE IF A SEQUENCE HAS FOUND SOME HOMOLOGOUS SEQUENCES
//INPUT THE 2-D MATRIX OF SEQUENCE PROFILE AND MEASURE 
//RETURN TRUE WHEN FIND SOME HOMOLOGOUS SEQUENCES FOR MORE THAN 40% SITES, ELSE FALSE
//=====================================================================
bool findhomo(IV2 &mary);

#endif

