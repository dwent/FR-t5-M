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


#ifndef SEQUENCE_H
#define SEQUENCE_H

#include <assert.h>
#include <cstring>

using namespace std;

//
// CLASS:       Sequence
//
// DESCRIPTION: we do not want to use a string for a sequence, since
//              for the algorithm the elements should be accessed by
//              values 1..n, instead of 0..n-1
//              

class Sequence {
public:
  Sequence(const char *str) : length_(strlen(str)), str_(str-1) {}
  int length() const {return length_;}
  char operator [] (int i) const {
    assert (1<=i && i<=length());
    return str_[i];
  }
private:
  int length_;
  const char *str_;
};

#endif // SEQUENCE_HH
