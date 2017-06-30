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

 
#include "Matrices.hpp"

using namespace std;


// =======================================================================
//
// CLASS:        DistanceMatrix
//
// AUTHOR:       Rolf Backofen <backofen@informatik.uni-muenchen.de>
//
// =======================================================================
				/**********************************************************
				Copyright(c)2012, IBP, CHINESE ACADEMY OF SCIENCE
				All rights reserved.
				NAME:		MATRICES.cpp
				ABSTRACT:	THREADING MATRICES
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
// ----------------------------------------------------------------------
// DistanceMatrix construction. a and b are the strings to be aligned
// sets row, col, a,b and initializes darray.
// 
DistanceMatrix::DistanceMatrix(const Sequence &a, const Sequence &b) 
  : matrix<cost_t>(a.length()+1,b.length()+1),a_(a), b_(b)
{
#ifndef NDEBUG
  cout <<"create distance matrix "<<a.length()+1<<" x "<<b.length()+1<<endl;
#endif
}

// ------------------------------------------------------------
// draw to printarray
//
void DistanceMatrix::draw(PrintArray &pr)
{
  char floatstr[PRECISION+1];
  int x,y;
  
  int n = a_.length();
  int m = b_.length();

  /* print b-row, with prefix for a-column   */
  /*                                         */
  for (int j=1; j<=m; ++j)
    {
      x = 2 + j*(PRECISION+ROWFILL) + (PRECISION+ROWFILL-1);
      y = 0;
      pr.set(x,y,b_[j]);
    }

  /* print darray, with a-entry first.       */
  /*                                         */
  for (int i=0; i<=n; ++i) 
    {
      
      y = 1+2*i;

      // a-col
      if (i>0) {
	pr.set(1,y,a_[i]);
      }
      
      // rest
      for (int j=0; j<=m; ++j)
        {
          if ((*this)[i][j] == INFTY)
            {
             strcpy(floatstr, "   INF");
            //floatstr = "   INF";
            }
          else if (MINUSINFTY == (*this)[i][j])
            { 
          	strcpy(floatstr, "  -INF");
            //floatstr = "  -INF";
            }
          else
            sprintf(floatstr,"%6.2f",(double)((*this)[i][j]));
          // cout << "|" << floatstr << "|" << endl;
          for (int l=0;l<PRECISION;++l)
            {
              x = 2 + j*(PRECISION+ROWFILL) + (ROWFILL+l);
              if (floatstr[l] != ' ')
                pr.set(x,y,floatstr[l]);
            }
        }
    }
}


// ----------------------------------------------------------------------
// Printing. The darray with the a preceeding row for b_[] and a preceeding 
// column for a_[] is printed. Every value is printed with precision
// PRECISION, and ROWFILL times a ' ' is added between every value.
//
ostream& operator << (ostream &out,DistanceMatrix &mat) {

  char floatstr[PRECISION+1];

  int n = mat.rows()-1;
  int m = mat.cols()-1;
  
  /* print b-row, with prefix for a-column   */
  /*                                         */
  cout << "   ";for (int k=0; k<PRECISION+ROWFILL-1;++k) cout<<' ';
  for (int j=1; j<=m; ++j)
    {
      for (int k=0; k<PRECISION+ROWFILL-1;++k)
        cout << ' ';
      cout << mat.b_[j];
    }
  cout<<endl;

  /* print darray, with a-entry first.       */
  /*                                         */
  for (int i=0; i<=n; ++i, cout<<endl) 
    {
      // a-col
      if (i>0) cout << ' ' << mat.a_[i];
      else     cout << "  ";

      for (int j=0; j<=m; ++j)
        {
          for (int l=0;l<ROWFILL;++l)
            cout << ' ';
          if (mat[i][j] == INFTY)
            cout << "   INF";
          else if (MINUSINFTY == mat[i][j])
            cout <<  "  -INF";
          else {
            sprintf(floatstr,"%6.2f",(double)(mat[i][j]));
            cout << floatstr;
          }
        }
    }
        
  return out;
}
//=======================================================================
//READ A SEQUENCE PROFILE FILE. AND 2 MATRICES IN THE PROFILE
//ARE STORED INT HE FARY AND MARY MATRICES, THE PDB SEQUENCE IN THE 
//STRING SEQ
//=======================================================================
void ReadProfPssm(char* Psifile,IV2& fary,IV2& mary,std::string &seq)
{
	mary.clear();fary.clear();seq.clear();
	IV1 itmp;
	int i;
	std::string line;
	std::string buf;
	std::stringstream sline;	
	std::ifstream infile(Psifile);
	if(!infile){std::cerr<<"\nCan not open: "<<Psifile<<'\n';exit(0);}
	getline(infile,buf);getline(infile,buf);getline(infile,buf);
	while(!infile.eof())
	{
		line.clear();
		getline(infile,line);
		if(line.empty())break;
		seq.push_back(line[6]);
		sline<<line;
		sline>>buf;sline>>buf;
		for(i=0;i<20;i++)
		{
			sline>>buf;
			itmp.push_back(atoi(buf.c_str()));
		}
		fary.push_back(itmp);itmp.clear();
		for(i=0;i<20;i++)
		{
			sline>>buf;
			itmp.push_back(atoi(buf.c_str()));
		}
		mary.push_back(itmp);itmp.clear();
		sline.str("");
	};infile.close();
}
//=======================================================================
//READ THE LSPP LIBRARY DATA OF THE TEMPLATE LIBRARY
//=======================================================================
void ReadUff(std::ifstream &infile,RV2 &uffmat,IV1 &ss)
{
	uffmat.clear();ss.clear();
	int i;
	RV1 vac(18,0),rtmp(18,0);
	cost_t phi,psi;
	cost_t phit,psit;
	std::string line;
	while(!infile.eof())
	{
		line.clear();
		getline(infile,line);
		if(line.empty())continue;
		if(line[0]=='>')continue;
		if(line.substr(0,4)=="END ")break;
		ss.push_back(atoi(line.substr(5,2).c_str())-1);
		for(i=0;i<5;i++)
		{
			phi=atof(line.substr(7+i*14,7).c_str());
			psi=atof(line.substr(14+i*14,7).c_str());
			if(phi>200||psi>200)
			{
				uffmat.push_back(vac);continue;
			}
			phit=phi*deg2rad;
			psit=psi*deg2rad;
			//FOURIER EXPENSION
			rtmp.assign(18,0);
			rtmp[0]=cos(phit);
			rtmp[1]=sin(phit);
			rtmp[2]=cos(2.*phit);
			rtmp[3]=sin(2.*phit);
			rtmp[4]=cos(3.*phit);
			rtmp[5]=sin(3.*phit);
			rtmp[6]=cos(psit);
			rtmp[7]=sin(psit);
			rtmp[8]=cos(2.*psit);
			rtmp[9]=sin(2.*psit);
			rtmp[10]=cos(3.*psit);
			rtmp[11]=sin(3.*psit);
			rtmp[12]=cos(phit)*cos(psit);
			rtmp[13]=cos(2.*phit)*cos(2.*psit);
			rtmp[14]=cos(3.*phit)*cos(3.*psit);
			rtmp[15]=sin(phit)*sin(psit);
			rtmp[16]=sin(2.*phit)*sin(2.*psit);
			rtmp[17]=sin(3.*phit)*sin(3.*psit);
			uffmat.push_back(rtmp);
		}
	};
}
