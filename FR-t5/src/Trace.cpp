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
				/**********************************************************
				Copyright(c)2012, IBP, CHINESE ACADEMY OF SCIENCE
				All rights reserved.
				NAME:		TRACE.cpp
				ABSTRACT:	THREADING TRACE (NO CHANGE FROM THE PREVIEWS VERSION)
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
 


#include "Trace.hpp"
#include "TraceArrow.hpp"
#include "Matrices.hpp"
#include "BasicDef.hpp"

using namespace std;


// =================================================================
// AUTHORS:      Rolf Backofen <backofen@informatik.uni-muenchen.de>
//               Sebastian Will <wills@informatik.uni-muenchen.de>
//
// =================================================================



// =================================================================
// CLASS:        TMatrixEntryElem
//
TMatrixEntryElem::TMatrixEntryElem(const TraceArrow &ta, int mult)
  : ta_(ta), mult_(mult)
{
}

TraceArrow const & TMatrixEntryElem::get_arrow() const {
  return ta_;
}


int TMatrixEntryElem::get_mult() const {
  return mult_;
}
void TMatrixEntryElem::print_alignment_string (int i, int j, 
                                               bool rowstring, 
                                               const Sequence &s) const
{
  ta_.print_alignment_string(i,j,mult_,rowstring,s);
}

void TMatrixEntryElem::store_alignment_string (string& alignedS,
							  int i, int j, 
                                               bool rowstring, 
                                               const Sequence &s) const
{
  ta_.store_alignment_string(alignedS,i,j,mult_,rowstring,s);
}

void TMatrixEntryElem::print_alignment_context (int i, int j, 
                                                bool rowcontext, bool prefix,
                                                bool gaps,
                                                const Sequence &a,
                                                const Sequence &b
                                                ) const
{
  ta_.print_alignment_context(i,j,mult_,rowcontext,prefix,gaps,a,b);
}

// =================================================================
// CLASS:        TraceMatrix
//
TraceMatrix::TraceMatrix(int n, int m, const DistanceMatrix &pDmatr)
  : matrix<TMatrixEntry>(n+1,m+1),Dmatr(pDmatr) {}

const DistanceMatrix &TraceMatrix::get_dmatr() const
{
  return Dmatr;
}



// =================================================================
// CLASS:        TraceElem
//
TraceElem::TraceElem(int i,int j, const TMatrixEntryElem &tee)
  : i_(i), j_(j), tee_(tee) {}

// calls the draw method of the stored traceentryelement with the stored
// endposition (i,j)
void TraceElem::draw(PrintArray &pr) const
{
  const TraceArrow & ta = tee_.get_arrow();
  ta.draw( pr , i_, j_,tee_.get_mult() );
}

// calls the print_alignment_string method of the stored traceentryelement
// with the stored endposition (i,j)
void TraceElem::print_alignment_string(bool rowstring, const Sequence &s) const
{
  tee_.print_alignment_string(i_,j_,rowstring,s);
}

// calls the store_alignment_string method of the stored traceentryelement
// with the stored endposition (i,j)
void TraceElem::store_alignment_string(string& alignedS, bool rowstring, const Sequence &s) const
{
	tee_.store_alignment_string(alignedS, i_, j_, rowstring, s);
}

// if rowscontext and prefix is true, then s is printed up to i_; 
// otherwise, s is printed up to j_; 
// similar holds for the other combinations
void TraceElem::print_alignment_context(bool rowcontext,
                                        bool prefix,
                                        bool gaps,
                                        const Sequence &a,
                                        const Sequence &b) const
{
  tee_.print_alignment_context(i_,j_,rowcontext,prefix,gaps,a,b);
}


// =================================================================
// CLASS:        BackTracer
//
BackTracer::BackTracer() {}
BackTracer::~BackTracer() {}
void BackTracer::calc_traceback(const TraceMatrix &tm,
                                int i, int j, int end_type) 

{
  // end-type 0 means stop at the cell (0,0)
  //
  assert(end_type!=0 || (i>=0 && j>=0));

  // get the traceset at (i,j) and the first first element
  // (which is a pointer to the TraceMatrix, and a TraceArrow
  //
  const TMatrixEntry &tme = tm[i][j];
  assert(!tme.empty());
  const TMatrixEntryElem &tee = tme.front();
  const TraceArrow & ta = tee.get_arrow();
  
  
#ifndef NDEBUG  
  // debug information
  cout << "Input BackTracePos: (" << i << "," << j << "), Out: (";
#endif

  // now get the end-position
  //
  pair<int,int> p = ta.end_position(i,j,tee.get_mult());
  int to_i = p.first;
  int to_j = p.second;
  
#ifndef NDEBUG  
  // debug information
  cout << to_i << "," << to_j << ")" << endl;
#endif
  
  const TraceMatrix &to_tmatr = ta.get_to_tmatr();
  
  // push that element on the backtrace stack
  //
  trace.push_front( TraceElem( to_i, to_j, tee ) );
  
  // check end condition and call recursive
  //     end_type == 0    ===>    stop at cell (0,0)
  //
  if ((end_type == 0) && ((to_i > 0) || (to_j > 0)))
    calc_traceback(  to_tmatr , to_i, to_j, end_type );
  else if ((end_type==1) && ((to_tmatr.get_dmatr())[to_i][to_j] > 0))
    calc_traceback(  to_tmatr , to_i, to_j, end_type );
}

// draw the trace into the printarray. Just call the draw function
// for every element of the backtrace
//
void BackTracer::draw_traceback(PrintArray &pr) {
  for(list<TraceElem>::iterator i=trace.begin() ; i != trace.end() ; ++i)
    i->draw(pr);
}

void BackTracer::print_alignment(const Sequence &a, const Sequence &b) {
  
  cout << "Alignment" << endl << "      a: ";

  // print "a"-line
  //
  for(list<TraceElem>::iterator e=trace.begin() ; e != trace.end() ; ++e)
    e->print_alignment_string(true,a);

  cout << endl << "      b: ";

  for(list<TraceElem>::iterator e=trace.begin() ; e != trace.end() ; ++e)
    e->print_alignment_string(false,b);
  
  cout << endl;
}


//Added by Wu Aiping, store alignment result into two strings
void BackTracer::store_alignment(string& alignedA,
					 string& alignedB,
					 const Sequence& a,
					 const Sequence& b)
{
  // store "a"-line
  //
  for(list<TraceElem>::iterator e=trace.begin() ; e != trace.end() ; ++e)
    e->store_alignment_string(alignedA, true, a);

  // store "b"-line
  for(list<TraceElem>::iterator e=trace.begin() ; e != trace.end() ; ++e)
    e->store_alignment_string(alignedB, false, b);

}



// print out alignment, where the context of alignment is shown in a
// (which is useful for local alignment)
void BackTracer::print_alignment_show_contexta(const Sequence &a,
                                               const Sequence &b) {
  
  cout << "Alignment" << endl << "      a: ";


  // print "a"-line
  //
  list<TraceElem>::iterator last_e = trace.begin();
  trace.begin()->print_alignment_context(true,true,false,a,b);
  for(list<TraceElem>::iterator e=trace.begin() ; e != trace.end() ; ++e) {
    last_e = e;
    e->print_alignment_string(true,a);
  }
  last_e->print_alignment_context(true,false,false,a,b);


  // print "b"-line
  //
  cout << endl << "      b: ";

  last_e = trace.begin();
  trace.begin()->print_alignment_context(true,true,true,a,b);
  for(list<TraceElem>::iterator e=trace.begin() ; e != trace.end() ; ++e) {
    last_e = e;
    e->print_alignment_string(false,b);
  }
  last_e->print_alignment_context(true,false,true,a,b);

  cout << endl;
}
// the same for b
void BackTracer::print_alignment_show_contextb(const Sequence &a,
                                               const Sequence &b) {
  
  cout << "Alignment" << endl << "      a: ";

  // print "a"-line
  //
  trace.begin()->print_alignment_context(false,true,false,a,b);
  for(list<TraceElem>::iterator e=trace.begin() ; e != trace.end() ; ++e)
    e->print_alignment_string(true,a);

  cout << endl << "      b: ";

  trace.begin()->print_alignment_context(false,true,true,a,b);
  for(list<TraceElem>::iterator e=trace.begin() ; e != trace.end() ; ++e)
    e->print_alignment_string(false,b);
  
  cout << endl;
}

