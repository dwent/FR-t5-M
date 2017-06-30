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

#ifndef MATRICES_H
#define MATRICES_H

#include <iostream>
#include <assert.h>
#include <fstream>
#include <cstdio>
#include <cstdlib>

#include "BasicDef.hpp"
#include "PrintArray.hpp"
#include "Matrix.hpp"
#include "Sequence.hpp"

using namespace std;

// =======================================================================
//
// CLASS:        DistanceMatrix
//
// AUTHOR:       Rolf Backofen <backofen@informatik.uni-muenchen.de>
//
// DESCRIPTION:  Creates and stores a distance or similarity matrix
//               for an alignment algorithm. Basically, it is a matrix of
//               entries of type cost_t. Basic
//               operations for DistanceMatrix are
//                    M(a,b)                  construction (a, b strings).
//                    draw(pr)                draw array into PrintArray
//                    << M                    print array only
//
// =======================================================================
//THE DistanceMatrix CLASS
class DistanceMatrix : public matrix<cost_t> {
private:
  const Sequence &a_;
  const Sequence &b_;
public:
  DistanceMatrix(const Sequence &a, const Sequence &b);

  void draw(PrintArray &pr);

  friend ostream& operator << (ostream &,DistanceMatrix &);
};
//READ THE UFF LIBRARY
void ReadUff(std::ifstream &infile,RV2 &uffmat,IV1 &ss);
//READ A PSI FILE
void ReadProfPssm(char* Psifile,IV2& fary,IV2& mary,std::string &seq);

#endif
