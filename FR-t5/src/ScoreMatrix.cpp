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

 
#include <fstream>
#include "ScoreMatrix.hpp"

using namespace std;



ScoreMatrix::ScoreMatrix(const char *filename)
{
  ifstream in(filename);
  if (!in.good()) {
    cerr<<"Can't open file "<<filename<<" for input"<<endl;
    exit(-1);
  }
  
  read_seqelems(in);
  read_matr(in);
}

ScoreMatrix::~ScoreMatrix() {
  delete matr;
}

// get an entry of the matrix, this is implemented efficiently. The
// translation of sequence element names to indices is done via the
// table "lookup".
cost_t ScoreMatrix::get(char x,char y) const {
  return (*matr)[lookup[x]][lookup[y]];  
}

// read one line from an input stream all characters in this line are
// regarded as sequence elements. A sequence element is a symbol for
// e.g a base or a amino acid.
// Further, initialize lookup for fast access to matrix entries by
// sequence elements
void ScoreMatrix::read_seqelems(istream &in) {
  static const int buf_siz=256; 
  char buf[buf_siz];
  
  if (in.eof()) {
    cerr << "ScoreMatrix: Can't read sequence elements"<<endl;
    exit(-1);
  }
  
  in.getline(buf,buf_siz);
  
  for (int i=0; buf[i]!=0; i++)
    if (buf[i]!=' ') seqelems.push_back(buf[i]);
  
  lookup.resize(256);
  
  for (size_t i=0; i<seqelems.size(); ++i)
    lookup[seqelems[i]] = i;
}

// read matrix from an input stream
// before the seqelems has to be initialized (e.g., by read_seqelems)
void ScoreMatrix::read_matr(istream &in) {
  
  matr = new matrix<cost_t>( seqelems.size(), seqelems.size() );
  
  for (size_t i=0; i<seqelems.size(); ++i) {
    string s;
    in >> s;
    // check that rows correspond to columns
    if ( s[0] != seqelems[i] ) {
      cerr << "ScoreMatrix: Row names are not equal to column names."<<endl;
      exit(-1);
    }
    
    for (size_t j=0; j<seqelems.size(); ++j) {
      if (in.eof()) {
	cerr << "ScoreMatrix: Can't read matrix, not enough entries"<<endl;
	exit(-1);
      }
      in >> (*matr)[i][j];
      
      // check symmetry
      if ( j<i )
	if ( ! approx_equal((*matr)[i][j], (*matr)[j][i], APPROX_EQUAL) ) {	//changed by wuaiping from "if ( ! approx_equal((*matr)[i][j], (*matr)[j][i]) ) {	"
	  cerr << "ScoreMatrix: Score-matrix is not symmetric"<<endl;
	  exit(-1);
	}
    }
  }
}
