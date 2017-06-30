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

 
#ifndef TRACE_H
#define TRACE_H

#define STACKLEN 10000

#include <string>
#include <list>
#include <iostream>
#include <stack>
#include <assert.h>
#include <fstream>
#include "Matrix.hpp"
#include "BasicDef.hpp"

#include "PrintArray.hpp"
#include "Sequence.hpp"


using namespace std;


// =================================================================
// AUTHORS:      Rolf Backofen <backofen@informatik.uni-muenchen.de>
//               Sebastian Will <wills@informatik.uni-muenchen.de>
// FILEDESC:     This file contains the classes implementing trace matrices
//               (which consist of TraceMatrix and TMatrixEntryElem),
//               and classes implementing the backtracing (which
//               consists of BackTracer and TraceElem).
// =================================================================


class TraceArrow;
class DistanceMatrix;



// =================================================================
// CLASS:        TMatrixEntryElem
//
// DESCRIPTION:  Element of an entry (which is a set) in TraceMatrix
//

class TMatrixEntryElem {
 public:
  TMatrixEntryElem(const  TraceArrow &, int mult=1);
  
  //TraceArrow const & TMatrixEntryElem::get_arrow() const;
  TraceArrow const & get_arrow() const;	//changed by wuaiping
  
  int get_mult() const;
  void draw(PrintArray &pr) const;
  // calls the draw method of the tracearrow
  // rowsstring = true means we print out the string associated with the
  // rows, false means the string associated with columns
  void print_alignment_string
    (int i, int j, bool rowstring, const Sequence &s) const; 

  void store_alignment_string (string& alignedS,
				     int i, int j, 
                               bool rowstring, 
                               const Sequence &s) const;
  // call the correspoding method with mult_
  void print_alignment_context(int i, int j,
                               bool rowcontext, bool prefix, bool gaps,
                               const Sequence &a,
                               const Sequence &b) const;
 private:
  const TraceArrow& ta_;
  int mult_;
};

// =================================================================
// CLASS:        TraceMatrix
//
// DESCRIPTION:  Matrix to store trace information 
//
typedef list<TMatrixEntryElem> TMatrixEntry;

class TraceMatrix: public matrix<TMatrixEntry> {
 private: 
  const DistanceMatrix &Dmatr;
 public:
  TraceMatrix(int n, int m, const DistanceMatrix &pDmatr);
  const DistanceMatrix &get_dmatr() const;
};






// =================================================================
// CLASS:        TraceElem
//
// DESCRIPTION:  Entry in a backtrace
//
class TraceElem {
 private:
  // endposition of the arrow
  const int i_;
  const int j_;
  // the traceentryelement itself
  const TMatrixEntryElem &tee_;
 public:
  TraceElem(int i, int j, const TMatrixEntryElem &tee);
  // calls the draw method of the stored traceentryelement with the stored
  // endposition (i,j)
  void draw(PrintArray &pr) const;
  // calls the print_alignment_string method of the stored traceentryelement
  // with the stored endposition (i,j)
  void print_alignment_string(bool rowstring, const Sequence &s) const;
  
  // calls the store_alignment_string method of the stored traceentryelement
  // with the stored endposition (i,j)
  void store_alignment_string(string& alignedS, bool rowstring, const Sequence &s) const;
  
  // call the corresponding methods
  void print_alignment_context(bool rowcontext, bool prefix, bool gaps,
                               const Sequence &a,
                               const Sequence &b) const;
};

// =================================================================
// CLASSES:      BackTracer 
//  
// DESCRIPTION:  This class has to purpose: First to store a trace 
//               matrix, second to calculate a backtrace (and store it), 
//               and third to print out a backtrace. 

class BackTracer {
 public:
  BackTracer();
  
  ~BackTracer();
  
  // calculate the traceback, starting in TraceMatrix tm at position i,j
  // end_type determines where to end
  void calc_traceback(const TraceMatrix & tm, int i, int j, int end_type);
  
  // draw the traceback to PrintArray pr
  void draw_traceback(PrintArray &pr);

  // print out alignment
  void print_alignment(const Sequence &a, const Sequence &b);
  
  //store alignment
  void store_alignment(string& alignedA, string& alignedB, const Sequence& a, const Sequence& b);
  
  // print out alignment, where the context of alignment is shown in a
  // (which is useful for local alignment)
  void print_alignment_show_contexta(const Sequence &a, const Sequence &b) ;
  void print_alignment_show_contextb(const Sequence &a, const Sequence &b) ;
    
private:
  list<TraceElem> trace;
};

#endif
