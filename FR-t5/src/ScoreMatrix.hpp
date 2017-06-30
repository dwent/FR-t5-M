// ============================================================
//
// CLASS:       ScoreMatrix
//
// DESCRIPTION: class to read in a matrix of scores for each pair
//              of sequence elements. Store this matrix, and provide
//              access to the scores.
//
// AUTHOR:      Rolf Backofen, Sebastian Will

//============================================================================
//
// Changed by Wu Aiping (wuaiping@moon.ibp.ac.cn), Aug, 2007
//
//============================================================================


#ifndef SCOREMATRIX_H
#define SCOREMATRIX_H

//class istream;
#include <iostream>
#include <cstring>
#include <cstdlib>
#include <cstdio>

#include "Matrix.hpp"
#include "BasicDef.hpp"

using namespace std;

//A CLASS USED TO STORE SCORING MATRIX DATA FOR DYNAMIC PROGRAMMING
class ScoreMatrix {
public:
  // read matrix from filename
  // format-example:
  /*
       A    C    G    T
  A   2.0 -1.6 -1.2 -1.1 
  C  -1.6  1.5 -1.4 -1.0
  G  -1.2 -1.4  1.3 -2.0
  T  -1.1 -1.0 -2.0  2.0
  */
  // the matrix must be symmetric, names for cols and rows have to be equal
  ScoreMatrix(const char *filename);
  
  ~ScoreMatrix();
  
  // access an entry by sequence elements, i.e. characters
  cost_t get(char,char) const;
  
private:
  matrix<cost_t> *matr;
  
  vector<char> seqelems;
  
  vector<int> lookup;
  
  void read_matr(istream &in);
  
  void read_seqelems(istream &in);
};

#endif

