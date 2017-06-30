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
//
//============================================================================


#ifndef PRINTARRAY_H
#define PRINTARRAY_H

#include <iostream>
#include "BasicDef.hpp"
#include "Matrix.hpp"

using namespace std;


// =======================================================================
//
// CLASS:        PrintArray
//
// DESCRIPTION:  PrintArray contains an array of char (which is
//               graphical object for printing). One can set chars,
//               and let the object print itself out onto cout.
//
//
// AUTHOR:       Rolf Backofen <backofen@informatik.uni-muenchen.de>
//
// =======================================================================
				/**********************************************************
				Copyright(c)2012, IBP, CHINESE ACADEMY OF SCIENCE
				All rights reserved.
				NAME:		PRINTARRAY.Hpp
				ABSTRACT:	THREADING PRINTARRAY
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
//A CLASS USED TO SORT ACCORDING TO Z
class Rank
{
public:
	Rank(int &m_ix,cost_t &m_s1,cost_t &m_s2,std::string &m_ss1,std::string m_ss2,cost_t m_z):ix(m_ix),s1(m_s1),s2(m_s2),ss1(m_ss1),ss2(m_ss2),z(m_z){}
	int ix;
	cost_t s1,s2;
	std::string ss1,ss2;
	cost_t z;
	bool operator < (const Rank &m)const 
	{return z>m.z;}
};
//========================================================================
//MEASURE THE SEQUENCE IDENTITY PERCENTAGE FROM THE RANK STRUCT
//========================================================================
cost_t seqid(Rank &a);
//AN CLASS USED TO PRINTARRAY DATA
class PrintArray: matrix<char> {
public:
  PrintArray(int cols, int rows);
  void print();
  void set(int x, int y, char c) {
      (*this)[x][y]=c;
  }
  char get(int x,int y) { return (*this)[x][y];}
};

#endif
