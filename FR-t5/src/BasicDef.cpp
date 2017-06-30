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


#include "BasicDef.hpp"

// =======================================================================
//
// CLASS:        Pair
//
// AUTHOR:       Rolf Backofen <backofen@informatik.uni-muenchen.de>
//               Sebastian Will <wills@informatik.uni-muenchen.de>
//
// =======================================================================
				/**********************************************************
				Copyright(c)2012, IBP, CHINESE ACADEMY OF SCIENCE
				All rights reserved.
				NAME:		BACISDEF.cpp
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
//UNKNOW FUNCTION USED BY HUYUN
bool approx_equal(cost_t x, const cost_t &y, cost_t ae=APPROX_EQUAL)
{
  return fabs( x-y ) < ae;
}

//=====================================================================
//DEFINE IF A SEQUENCE HAS FOUND SOME HOMOLOGOUS SEQUENCES
//INPUT THE 2-D MATRIX OF SEQUENCE PROFILE AND MEASURE 
//RETURN TRUE WHEN FIND SOME HOMOLOGOUS SEQUENCES FOR MORE THAN 40% SITES, ELSE FALSE
//=====================================================================
bool findhomo(IV2 &mary)
{
	int i,j,n=0,L=mary.size();
	for(i=0;i<L;i++)
	{
		for(j=0;j<20;j++)
		if(mary[i][j]>99)
		{
			n++;
			break;
		}
	}
	if((cost_t)n/(cost_t)L>0.4)return false;
	return true;
}
