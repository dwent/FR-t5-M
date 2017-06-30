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

 
#include "PrintArray.hpp"

using namespace std;


// =======================================================================
//
// CLASS:        PrintArray
//
// AUTHOR:       Rolf Backofen <backofen@informatik.uni-muenchen.de>
//
// =======================================================================

				/**********************************************************
				Copyright(c)2012, IBP, CHINESE ACADEMY OF SCIENCE
				All rights reserved.
				NAME:		PRINTARRAY.cpp
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

PrintArray::PrintArray(int prows, int pcols)
  : matrix<char>(pcols,prows)	//changed by wuaiping from ": matrix(pcols,prows)"
{
  for (int x=0;x<rows();++x)
    for (int y=0;y<cols();++y)
      set(x,y,' ');
}

void PrintArray::print()
{
  for (int y=0;y<cols();++y, cout<<endl)
    for (int x=0;x<rows(); ++x)
        cout << get(x,y);
}
//========================================================================
//MEASURE THE SEQUENCE IDENTITY PERCENTAGE FROM THE RANK STRUCT
//========================================================================
cost_t seqid(Rank &a)
{
	int i,j=0,k=0;
	for(i=0;i<(int)a.ss1.length();i++)
	if(a.ss1[i]!='-')
	{
		j++;
		if(a.ss2[i]==a.ss1[i])k++;
	}
	return (cost_t)k/(cost_t)j;
}
