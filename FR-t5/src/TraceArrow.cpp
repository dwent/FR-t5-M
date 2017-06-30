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


#include <string>
#include <assert.h>
#include <fstream>
#include <cmath>

//#include "PrintArray.hpp"
#include "TraceArrow.hpp"
//#include "BasicDef.hpp"

using namespace std;


// ---------------------------------------------------------
//
// CLASS:        TraceArrow
//
// AUTHORS:      Rolf Backofen <backofen@informatik.uni-muenchen.de>
//

// ---------------------------------------------------------------------
// Constructor. Checks whether in row printing is possible for
// arrows that are only in one row.
// 
TraceArrow::TraceArrow(const TraceMatrix &pto_tmatr, int prows, int pcols)
  : rows(prows), cols(pcols), to_tmatr(pto_tmatr) {}

// -----------------------------------------------------------------
// get functions
const TraceMatrix &TraceArrow::get_to_tmatr() const
{
  return to_tmatr;
}


// ---------------------------------------------------------------------
// functions for defining the recursion equation:
// They are of the form
//     f(int i, int j, int mult),
// where (i,j) is the current cell, and mult is the
// multiplicator used (e.g., mult = 1 for
// Needleman/Wunsch, and possibly mult=k in
// Smith/Waterman).
// 

// ---------------------------------------------------------------------
// to be safe, test applicability first. 
//
pair<int,int> TraceArrow::end_position(int i, int j, int mult) const
{
  assert((i-mult*rows >= 0) && (j-mult*cols >= 0));
  return(make_pair(i-mult*rows,j-mult*cols));
}

// ---------------------------------------------------------------------
// prints-out the context of a possibly local alignment. If gaps is true,
// then ' ' is printed. Otherwise, if rowcontext is true, then  the string of a
// is printed. Otherwise, b is printed
//
void TraceArrow::print_alignment_context(int i, int j, int mult,
                                         bool rowcontext, bool prefix,
                                         bool gaps,
                                         const Sequence &a,
                                         const Sequence &b) const {

  int start, end, start_i, start_j, end_i, end_j;

  // // debug
  // cout << endl << "hier:" << i <<"|"<< 
  // j <<"|"<< mult <<"|"<< rowcontext <<"|"<< prefix <<"|"<< gaps << endl; 


  // region for printing
  //
  if (prefix) {
    start_i = 1;
    start_j = 1;
    end_i=i;
    end_j=j;
  }
  else {
    start_i = i+mult*rows+1;
    start_j = j+mult*cols+1;
    end_i=a.length();
    end_j=b.length();
  }

  // If row context is printed, then we use only _i;
  // otherwise, we use the _j
  //
  if (rowcontext) {
    start=start_i;
    end=end_i;
  }
  else {
    start=start_j;
    end=end_j;
  }
   
  // // debug
  // cout << "start:" << start << " end:" <<  end << endl;

  // now print the context
  //
  for (int k=start;k<=end;k++)
      if (gaps)
          cout << '.';
      else if (rowcontext) 
        cout << a[k];
      else
        cout << b[k];
}


// ---------------------------------------------------------------------
// the basic definition is the Levenshtein distance
//
cost_t TraceArrow::cost(int i, int j, int mult) const
{
  return(mult*(rows+cols));
}

cost_t TraceArrow::cost(int i, int j, 
			const Sequence &a, const Sequence &b, int mult) const {
  return(mult*(rows+cols));
}

// -----------------------------------------------------------------
// print the part of the alignment string for the row resp. column
// string.
//
void TraceArrow::print_alignment_string(int end_i, int end_j, int mult, 
                                        bool rowstring,
                                        const Sequence &s) const
{
  // if string is a rowstring, then stringspan is rows. Otherwise, it
  // is cols. The not_stringspan is how much the arrow is longer in the
  // dimension different of the dimension the string is written.
  // 
  int stringspan;
  int not_stringspan;
  int startprint;

  if (rowstring == true)
    {
      stringspan = mult*rows;
      not_stringspan = max(0,mult*(cols-rows));
      startprint = end_i;
    }
  else
    {
      stringspan = mult*cols;
      not_stringspan = max(0,mult*(rows-cols));
      startprint = end_j;
    }
  
  // print elements of from s[i] for the part 
  //
  for (int l=0; l < stringspan; ++l)
    cout << s[startprint+l+1];
  for (int l=0; l < not_stringspan; ++l)
    cout << '-';
}


// -----------------------------------------------------------------
// store the part of the alignment string for the row resp. column
// string.
//
void TraceArrow::store_alignment_string(string& alignedS,
						 int end_i, int end_j, int mult, 
                                        bool rowstring,
                                        const Sequence &s) const
{
  // if string is a rowstring, then stringspan is rows. Otherwise, it
  // is cols. The not_stringspan is how much the arrow is longer in the
  // dimension different of the dimension the string is written.
  // 
  int stringspan;
  int not_stringspan;
  int startprint;

  if (rowstring == true)
    {
      stringspan = mult*rows;
      not_stringspan = max(0,mult*(cols-rows));
      startprint = end_i;
    }
  else
    {
      stringspan = mult*cols;
      not_stringspan = max(0,mult*(rows-cols));
      startprint = end_j;
    }
  
  // print elements of from s[i] for the part 
  //
  for (int l=0; l < stringspan; ++l)
    alignedS.push_back( s[startprint+l+1] );
  for (int l=0; l < not_stringspan; ++l)
    alignedS.push_back( '-' );
}



// ---------------------------------------------------------------------
// Print arrow on printarray. Prints all the way from
// cell (end_i,end_j) to cell (end_i+mult*rows,end_j+mult*cols)
// Overwrites existing entries, i.e.,  (end_i,end_j) and
// (end_i+mult*rows,end_j+mult*cols) will be connected by a line.
//
void TraceArrow::draw(PrintArray &pr, int end_i, int end_j, int mult) const
{
  int current_i = end_i;
  int current_j = end_j;
  int x,y;

#ifndef NDEBUG
  cout << "Draw EndPos: (" << end_i << "," << end_j << ")" << endl;
#endif
  
  assert(mult!=0 /*Error: Arrow with 0-multiplication invalid*/);

  /* print diagonal */
  /*                */
  while ((current_i-end_i < mult*rows) && (current_j-end_j < mult*cols))
    {
      y = 1 + 2*current_i + 1;
      x = 2 + (current_j+1)*(PRECISION+ROWFILL);
      pr.set(x,y,diagchar);
      /* draw a continuing dashed line in next row (up to the final  */
      /*   cell)                                                     */
      /*                                                             */
      ++y;
      for (int l=1;l<PRECISION+ROWFILL;++l)
        {
          ++x;
          if (l==1)
            pr.set(x,y,diagchar);
          else
            pr.set(x,y,rchar);
        }
      current_i++;
      current_j++;
    }
  /* print down     */
  /*                */
  while (current_i-end_i <  mult*rows)
    {
      y = 1 + 2*current_i + 1;
      x = 2 + (current_j+1)*(PRECISION+ROWFILL) - 1;
      pr.set(x,y,dchar);
      current_i++;
    }
  /* print right    */
  /*                */
  while (current_j-end_j <  mult*cols)
    {
      y = 1 + 2*current_i;
      x = 2 + (current_j+1)*(PRECISION+ROWFILL);
      pr.set(x,y,rchar);
      for (int l=1;l<PRECISION+ROWFILL;++l)
        {
          ++x;
          pr.set(x,y,rchar);
        }
      current_j++;
    }
}



// ---------------------------------------------------------
//
// CLASS:        SubstArrow
//
SubstArrow::SubstArrow(const TraceMatrix &to_tmatr, int rows_cols,
                       cost_t psubst)
  : TraceArrow(to_tmatr,rows_cols,rows_cols),subst(psubst) {}
     
cost_t SubstArrow::cost(int i, int j, 
                        const Sequence &a, const Sequence &b, int mult) const
{
  cost_t c=0;
  for (int l=0; l < mult*rows; ++l) 
    if (a[i+l] != b[j+l])
        c = c+subst;
  return c;
}



// ---------------------------------------------------------
//
// CLASS:        SimilarityArrow
//
SimilarityArrow::SimilarityArrow(const TraceMatrix &to_tmatr, int
                                 rows_cols, cost_t pmatch, cost_t pmismatch) 
  : TraceArrow(to_tmatr,rows_cols,rows_cols),
  match(pmatch),mismatch(pmismatch) {}
     
cost_t SimilarityArrow::cost(int i, int j, 
                             const Sequence &a, const Sequence &b, int mult) const
{
  cost_t c=0;
  for (int l=0; l < mult*rows; ++l) 
    if (a[i+l] == b[j+l])
        c = c+match;
    else
      c = c+mismatch;
  return c;
}

// ---------------------------------------------------------
//
// CLASS:        ScoreMatrixArrow
//
ScoreMatrixArrow::ScoreMatrixArrow(const TraceMatrix &to_tmatr, int
                                 rows_cols, const ScoreMatrix &sm) 
  : TraceArrow(to_tmatr,rows_cols,rows_cols),
    sm_(sm) {}

cost_t ScoreMatrixArrow::cost(int i, int j, 
			      const Sequence &a, const Sequence &b, int mult) const
{
  cost_t c=0;
  for (int l=0; l < mult*rows; ++l) 
    c += sm_.get( a[i+l], b[j+l] );
  return c;
}


// ---------------------------------------------------------
//
// CLASS:        SelfScoreMatrixArrow
//
SelfScoreMatrixArrow::SelfScoreMatrixArrow(const TraceMatrix &to_tmatr, int
                                 rows_cols, const vector<vector<double > >& selfsm) 
  : TraceArrow(to_tmatr,rows_cols,rows_cols), selfsm_(selfsm){}

cost_t SelfScoreMatrixArrow::cost(int i, int j, 
			      const Sequence &a, const Sequence &b, int mult) const
{
  cost_t c=0;
  for (int l=0; l < mult*rows; ++l) 
    c += selfsm_[i+l-1][j+l-1];
  return c;
}



// ---------------------------------------------------------
//
// CLASS:        AffineIndelArrow
//

AffineIndelArrow::AffineIndelArrow(const TraceMatrix &pto_tmatr,
                                   int prows, int pcols,
                                   cost_t pgap_var, cost_t pgap_const)
  : TraceArrow(pto_tmatr,prows,pcols),gap_var(pgap_var),gap_const(pgap_const)
{
  assert((rows == 0) || (cols == 0));
}

cost_t AffineIndelArrow::cost(int i, int j, int mult) const
{
  assert(mult>0);
  return gap_const + gap_var*(mult*(rows+cols));
}

// ---------------------------------------------------------
//
// CLASS:        LogarithmicIndelArrow
//

LogarithmicIndelArrow::LogarithmicIndelArrow(const TraceMatrix &pto_tmatr,
                                   int prows, int pcols,
                                   cost_t pgap_var, cost_t pgap_const)
  : TraceArrow(pto_tmatr,prows,pcols),gap_var(pgap_var),gap_const(pgap_const)
{
  assert((rows == 0) || (cols == 0));
}

cost_t LogarithmicIndelArrow::cost(int i, int j, int mult) const
{
  assert(mult>0);
  return(gap_const + gap_var*((rows+cols)*log(mult*1.0)));
}


// ---------------------------------------------------------
//
// CLASS:        WrapAroundArrow
//
//
WrapAroundArrow::WrapAroundArrow(const TraceMatrix &pto_tmatr, 
                                 int prows, int pcols,
                                 int plen_b, cost_t cost)
  : TraceArrow(pto_tmatr,prows,pcols),len_b(plen_b),subst_or_indel_cost(cost)
{
  assert(((rows == 1) && (cols == 1))
         ||
         ((rows == 0) && (cols == 1)));
}

// -----------------------------------------------------------------
// cost function for wraparound arrow. Is always 1 insertion or deletion
// cost (since wraparound can have only length 1)
cost_t WrapAroundArrow::cost(int i, 
                             const Sequence &a, const Sequence &b) const
{
  if (rows == 0)
    return(subst_or_indel_cost);
  else {
    if (a[i] == b[1])
      return(0);
    else
      return(subst_or_indel_cost);
  }
}

// -----------------------------------------------------------------
// end-position: always len_b
//
pair<int,int> WrapAroundArrow::end_position(int i, int j, int mult) const
{
  assert((mult==1) && (j==1));
  return(make_pair(i-rows,len_b));
}

// -----------------------------------------------------------------
// print could be imagine as starting by column 0
//
void WrapAroundArrow::print_alignment_string(int end_i, int end_j, int mult, 
                                             bool rowstring,
                                             const Sequence &s) const
{
  assert((mult==1) && (end_j==len_b));
  TraceArrow::print_alignment_string(end_i, 0, 1, rowstring, s);
}

void WrapAroundArrow::draw(PrintArray &pr, int end_i, int end_j, int mult)
 const
{
  int y,start_x,end_x;

  // wraparound arrows cannot be multiplied, and it must connect
  // the first and last column
  //
  assert(mult == 1);
  assert(end_j == len_b);
  
  // calculate the y-value for the row before the printing,
  // and the x-value of the beginning and end of the line 
  // between the to cells.
  //
  start_x = 2 + 1*(PRECISION+ROWFILL);
  end_x = 2 + (len_b+1)*(PRECISION+ROWFILL);
  if (rows==1) 
    y = 1 + 2*end_i + 1;
  else 
    y = 1 + 2*(end_i-1) + 1;
  
  // print 
  //            \------------------------------------------------/
  //
  pr.set(start_x,y,diagchar);
  for (int l=1;l<end_x-start_x+1;++l)
      pr.set(start_x+l,y,rchar);

  if (rows==1) 
    pr.set(end_x+1,y,diagchar);
  else
    pr.set(end_x+1,y,revdiagchar);
  
  // print the connection to cell (end_i,end_j)
  //
  if (rows == 1)
    pr.set(end_x,y-1,diagchar);
  else
    pr.set(end_x,y+1,revdiagchar);

  // print the connection to the starting point of the arrow
  // (which is either (1,end_i), or (1,end_i+1)).
  //
  ++start_x;
  pr.set(start_x,y+1,diagchar);
  for (int l=2;l<PRECISION+ROWFILL;++l)
    {
      ++start_x;
      pr.set(start_x,y+1,rchar);
    }
}
