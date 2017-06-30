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

 
#ifndef TRACEARROW_H
#define TRACEARROW_H

#include <vector>
#include <iostream>
//#include "BasicDef.hpp"
//#include "PrintArray.hpp"
#include "Trace.hpp"
#include "ScoreMatrix.hpp"

using namespace std;


// ---------------------------------------------------------
//
// CLASS:        TraceArrow
//
// DESCRIPTION:  TraceArrows are basically vectors, pointing to a 
//               TraceMatrix (described by their index in an array of
//               TraceMatrices, which are stored in the Trace object).
//               TraceArrows are defined by the number of columns and
//               rows they span, which are positive integers. The
//               represented vector 
//                          (-rows,-cols).
//               TraceArrow has the following methods which are needed
//               for defining the recursion equations. They are of the
//               form
//                   f(int i, int j, int mult),
//               where (i,j) is the current cell, and mult is the
//               multiplicator used (e.g., mult = 1 for
//               Needleman/Wunsch, and possibly mult=k in
//               Smith/Waterman).
//               The functions are
//                   cost(int i, int j, int mult)
//                       returns cost for this arrow
//                       basic implementation: levenstein distance,
//                       i.e., mult*(rows + cols);
//                   new_position(int i, int j, int mult)
//                       returns (i-mult*rows,j-mult*cols) 
//               other functions:
//                   draw(PrintArray &pr, int end_i, int end_j, int mult)
//                       draw a connection between the cells 
//                              (end_i,end_j) 
//                       with
//                              (end_i+mult*rows,end_j+mult*cols) 
//                       in pr. IMPORTANT: Overwrite existing entries  
//


class TraceArrow {
protected:
  const int rows,cols;
  const TraceMatrix & to_tmatr;
  static const char rchar = '-';
  static const char dchar = '|';
  static const char diagchar = '\\';
public:
  TraceArrow(const TraceMatrix &pto_tmatr, int prows, int pcols);
  const TraceMatrix &get_to_tmatr() const;
  virtual cost_t cost(int i, int j, int mult=1) const;
  virtual cost_t cost(int i, int j, 
		      const Sequence &, const Sequence &, int mult=1) const;

  virtual pair<int,int> end_position(int i, int j, int mult=1) const;
  // the span of the arrow. Note that this has to be done for
  // wraparound arrows differently, but is submitted for the moment
  // virtual pair<int,int> span(int mult=1) const;
  virtual void draw(PrintArray &pr, int end_i, int end_j, int mult) const;
  virtual void print_alignment_string(int end_i, int end_j, int mult, 
                              bool rowstring, const Sequence &s) const;
  
  void store_alignment_string(string& alignedS,
			       int end_i, int end_j, int mult, 
                          bool rowstring,
                          const Sequence &s) const;
  
  // prints out the context of a possibly local alignment
  // if prefix is true, then the prefix context is print, otherwise the
  // suffix context
  void print_alignment_context(int i, int j, int mult,
                                      bool rowstring, bool prefix,  bool gaps,
                                      const Sequence &a,
                                      const Sequence &b) const;
};


// ---------------------------------------------------------
// SUBCLASSES of TRACEARROW for specific alignment algorithms
// ---------------------------------------------------------

// ---------------------------------------------------------
//
// CLASS:        SubstArrow
//
// DESCRIPTION:  Subclass of TraceArrow that models (possibly k-times)
//               substitution. The cost substitution is independent
//               of the substited arrows (only difference counst)
//
class SubstArrow : public TraceArrow {
 protected:
  const cost_t subst;
 public:
  SubstArrow(const TraceMatrix &pto_tmatr, int prows_cols, cost_t psubst);
  cost_t cost(int i, int j,
              const Sequence &a, const Sequence &b, int mult=1) const;
};
  


// ---------------------------------------------------------
//
// CLASS:        SimilarityArrow
//
// DESCRIPTION:  Subclass of TraceArrow that models (possibly k-times)
//               substitution. Calculates cost based on similarity given scores
//               for match and mismatch
//
class SimilarityArrow : public TraceArrow {
 protected:
  const cost_t match;
  const cost_t mismatch;
 public:
  SimilarityArrow(const TraceMatrix &to_tmatr, int
                  rows_cols, cost_t pmatch, cost_t pmismatch);
  cost_t cost(int i, int j,
              const Sequence &a, const Sequence &b, int mult=1) const;
};


// ---------------------------------------------------------
//
// CLASS:        ScoreMatrixArrow
//
// DESCRIPTION:  Subclass of TraceArrow that models (possibly k-times)
//               substitution. Calculates cost based on a given
//               a score matrix
//
class ScoreMatrixArrow : public TraceArrow {
 protected:
  const ScoreMatrix &sm_;
 public:
  ScoreMatrixArrow(const TraceMatrix &to_tmatr, int rows_cols, const ScoreMatrix &sm);
  cost_t cost(int i, int j,
              const Sequence &a, const Sequence &b, int mult=1) const;
};
  

// ---------------------------------------------------------
//
// CLASS:        SelfScoreMatrixArrow
//
// DESCRIPTION:  Subclass of TraceArrow that models (possibly k-times)
//               substitution. Calculates cost based on a given
//               a score matrix. This score matrix can be not symmetrical
//		    matrix such as BLOSUM and PAM, it can be defined as
//		    n*m asymmetrical format by oneself 
//
class SelfScoreMatrixArrow : public TraceArrow {
 protected:
   const vector<vector<double > >& selfsm_;
 public:
  SelfScoreMatrixArrow(const TraceMatrix &to_tmatr, int rows_cols, const vector<vector<double > >& selfsm);
  cost_t cost(int i, int j,
              const Sequence &a, const Sequence &b, int mult=1) const;
};




// ---------------------------------------------------------
//
// CLASS:        AffineIndelArrow
//
// DESCRIPTION:  Subclass of TraceArrow that models Indel-arrows
//               with affine gap penality. In the constructor, 
//               whether it is an up-arrow or left arrow is model
//               by setting either rows or cols to 0. One mustn't have
//               both rows and cols different from 0.
//
class AffineIndelArrow : public TraceArrow {
 protected:
  const cost_t gap_var;
  const cost_t gap_const;
 public:
  AffineIndelArrow(const TraceMatrix &pto_tmatr, int prows, int pcols,
                   cost_t pgap_var, cost_t pgap_const=0);
  cost_t cost(int i, int j, int mult=1) const;
};

// ============================================================
//
// CLASS:        LogarithmicIndelArrow
//
// DESCRIPTION:  Subclass of TraceArrow that models Indel-arrows
//               with logarithmic gap penality. In the constructor, 
//               whether it is an up-arrow or left arrow is model
//               by setting either rows or cols to 0. One mustn't have
//               both rows and cols different from 0.
//

class LogarithmicIndelArrow : public TraceArrow {
protected:
  const cost_t gap_var;
  const cost_t gap_const;
 public:
  LogarithmicIndelArrow(const TraceMatrix &pto_tmatr, int prows, int pcols,
                   cost_t pgap_var, cost_t pgap_const=0);
  virtual cost_t cost(int i, int j, int mult=1) const;
};


// ---------------------------------------------------------
//
// CLASS:        WrapAroundArrow
//
class WrapAroundArrow : public TraceArrow {
protected:
  static const char revdiagchar = '/';
  const int len_b;
  const cost_t subst_or_indel_cost;
public:
  WrapAroundArrow(const TraceMatrix &pto_tmatr,
                   int prows, int pcols, int plen_b, cost_t cost=1);
  virtual void draw(PrintArray &pr, int end_i, int end_j, int mult=1) const;
  
  virtual pair<int,int> end_position(int i, int j, int mult=1) const;
  virtual void print_alignment_string(int end_i, int end_j, int mult, 
                              bool rowstring, const Sequence &s) const;
  // TO BE DONE: 

  // we have a special cost function for   WrapAroundArrow,
  // since this will be called from main program.
  // hence, we know the speciality of wraparound
  virtual cost_t cost(int i, const Sequence &a, const Sequence &b) const;
};

#endif



