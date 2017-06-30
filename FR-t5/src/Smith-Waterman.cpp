//===================================================================================
//
// Implementation of Threading Sequence Alignment
//
// AUTHORS:     Huyun
//
// DESCIRPTION: Training and Testing Datasets
//               LOCAL/GlOBAL alignment. (Feb, 2009)
//===================================================================================
				/**********************************************************
				Copyright(c)2012, IBP, CHINESE ACADEMY OF SCIENCE
				All rights reserved.
				NAME:		SMITHWATERMAN.cpp
				ABSTRACT:	THREADING THE DYNAMIC PROGRAMMING FUNCTIONS
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

#include "Smith-Waterman.hpp"
//#include "weight.h"

using namespace std;

int Align_Global_Zhang(	AlignResult& alignresult,
		SeqStr& queueseq,
		SeqStr& targstr,
		double factor[8],
		const char* a_str,
		const char* b_str,
		const char* sm_filename,
		cost_t match,
		cost_t mismatch,
		cost_t indel_var,
		cost_t indel_const,
		bool logarithmic,
		bool local,
		bool quiet,
		bool sm_opt,
		bool self_sm_opt
	     )
{
/*	// ---------------------------------------------------------  
	// options
	//
	char *a_str=0;
	char *b_str=0;
	char *sm_filename=0;
	cost_t match;
	cost_t mismatch;
	cost_t indel_var;
	cost_t indel_const;
	bool logarithmic;
	bool local;  
	bool quiet;
	bool sm_opt;


	option_def options[] = 
	{
		{0,0,0,O_ARG_STRING,&a_str,0,"sequence a","The first Sequence"},
		{0,0,0,O_ARG_STRING,&b_str,0,"sequence b","The second Sequence"},
		{"match",'m',0,O_ARG_DOUBLE,&match,"4.0","num","Cost of match"},
		{"mismatch",'i',0,O_ARG_DOUBLE,&mismatch,"-4.0","num","Cost of mismatch"},
		{"scorematrix",'s',&sm_opt,O_ARG_STRING,&sm_filename,0,"filename","File containing a score matrix"},
		{"gapconst",'c',0,O_ARG_DOUBLE,&indel_const,"-4.0","num","Cost to begin insertion or deletion "},
		{"gapvar",'v',0,O_ARG_DOUBLE,&indel_var,"-2.0","num","Cost to begin insertion or deletion "},
		{"logarithmical",'l',&logarithmic,O_NO_ARG,0,0,0,"Logarithmical gap costs"},
		{"local",'o',&local,O_NO_ARG,0,0,0,"Local alignment"},
		{"quiet",'q',&quiet,O_NO_ARG,0,0,0,"Quiet mode. Just print alignment"},
		// STOP
		{0,0,0,0,0,0,0,0}
	};

 
	if (!process_options(argc,argv,options))
	{
		print_help(argv[0],options);
		exit(-1);
	};
	
	print_options(options);
*/
	Sequence a(a_str);
	Sequence b(b_str);

	//Sequence a(b_str);
	//Sequence b(a_str);

	// ---------------------------------------------------------
	// auxiliary variables
	//
	int n = a.length();
	int m = b.length();	

	// ---------------------------------------------------------
	// define the distance matrix and trace matrix
	//
	DistanceMatrix D(a,b);
	TraceMatrix TraceD(n,m,D);


	// ---------------------------------------------------------
	// the different arrows
	//
	// define arrows
	//              
	TraceArrow *diag_arrow;
	TraceArrow *left_arrow;
	TraceArrow *up_arrow;
	ScoreMatrix *scorematrix=0;
	
	vector<vector<double > > selfsm;
	if (sm_opt)
	{
		if( self_sm_opt )
		{
			
			GetSelfSM(selfsm, queueseq, targstr, factor);
			
			//cout << "\nSelf Defined Score Matrix:\n";
		/*
			for(int i=0; i<selfsm.size(); i++)
			{
				for(int j=0; j<selfsm[i].size(); j++)
				{
					if(i>=149 && j>=25)
					 {cout <<" i=" << i << " " << a[i]<< " j=" <<j << " " << b[j] <<  "  ";
					cout << selfsm[i][j] <<endl;}
				}
			}
		*/	
			diag_arrow = new SelfScoreMatrixArrow(TraceD, 1, selfsm);
						
		}
		else
		{
			scorematrix = new ScoreMatrix( sm_filename );
			diag_arrow = new ScoreMatrixArrow(TraceD, 1, *scorematrix);  
		}
	}
	else
	{
		diag_arrow = new SimilarityArrow(TraceD, 1, match, mismatch);
	}

	if (logarithmic)
	{
		left_arrow = new LogarithmicIndelArrow(TraceD,0,1,indel_var,indel_const);
		up_arrow = new LogarithmicIndelArrow(TraceD,1,0,indel_var,indel_const);
	}
	else
	{
		left_arrow = new AffineIndelArrow(TraceD,0,1,indel_var,indel_const);
		up_arrow = new AffineIndelArrow(TraceD,1,0,indel_var,indel_const);
	}
    
	// ---------------------------------------------------------
	// Intialization
	//
	D[0][0]=0;//
	TraceArrow *left_arrow_startgap;
	TraceArrow *up_arrow_startgap;
	//left_arrow_startgap = new AffineIndelArrow(TraceD,0,1,0,0);
	//up_arrow_startgap = new AffineIndelArrow(TraceD,1,0,0,0);

	left_arrow_startgap = new AffineIndelArrow(TraceD,0,1,0.000,0.000);
	up_arrow_startgap = new AffineIndelArrow(TraceD,1,0,0.000,0.000);


	for (int i=1;i<=n;++i)
	{

		D[i][0] = up_arrow_startgap->cost(i,0,i);
		TraceD[i][0].push_back( TMatrixEntryElem(*up_arrow_startgap) ); //默认缺失时候，i=1
		//TraceD[i][0].push_back( TMatrixEntryElem(*up_arrow_startgap,i) );
	
    	}
    	
	for (int j=1;j<=m;++j)
	{
		D[0][j] = left_arrow_startgap->cost(0,j,j); 
		TraceD[0][j].push_back( TMatrixEntryElem(*left_arrow_startgap) );
		//TraceD[0][j].push_back( TMatrixEntryElem(*left_arrow_startgap,j) );

	}

	cost_t maxD=MINUSINFTY;
	int max_i=0, max_j=0;

	cost_t v_diag, v_up[n], v_left[m];
	
	// ---------------------------------------------------------
	// recursive definition of D-Matrix
	//
	for (int i=1;i<=n;++i)
	{
		for (int j=1;j<=m;++j)
		{
			// calculate D_ij
			//
			//  for debuging
			//cout << "hier i: " << i << "  j:" << j << endl;
        
			if (local) D[i][j] = 0;
			else D[i][j] = MINUSINFTY;

			v_diag= D.get(diag_arrow->end_position(i,j)) + diag_arrow->cost(i,j,a,b);
			D[i][j] = max(D[i][j],v_diag);
			//if(i>160 && j>140) cout <<" i=" << i << " " << a[i]<< " j=" <<j << " " << b[j]<< " D.get(diag_arrow->end_position(i,j))=" <<D.get(diag_arrow->end_position(i,j)) << "  diag_arrow->cost(i,j,a,b)=" <<diag_arrow->cost(i,j,a,b)<<  endl;
			
			bool flag = (targstr.second[j-1]==queueseq.second[i-1] && (targstr.second[j-1]=='H' || targstr.second[j-1]=='E'));
			//flag =false;		
			if(!flag)
			{
				for (int k=1; k<=i; k++)
				{ 	if( i !=n)
					{
						v_up[k] = D.get(up_arrow->end_position(i,j,k)) + up_arrow->cost(i,j,k);
					}
					else
					{
						v_up[k] = D.get(up_arrow->end_position(i,j,k))+0;
					}
					D[i][j] = max(D[i][j],v_up[k]);
					//if(i>=149 && j>=25) cout <<" i=" << i << " " << a[i]<< " j=" <<j << " " << b[j]<< "  k=" << k << " D.get(up_arrow->end_position(i,j,k))=" <<D.get(up_arrow->end_position(i,j,k))<< "  up_arrow->cost(i,j,k)=" <<up_arrow->cost(i,j,k)<<  endl;
				}

				for (int k=1; k<=j; k++)
				{
					if( j !=m)
					{
						v_left[k] = D.get(left_arrow->end_position(i,j,k)) + left_arrow->cost(i,j,k);
					}
					else
					{	
						v_left[k] = D.get(left_arrow->end_position(i,j,k))+0 ;
					}
					D[i][j] = max(D[i][j],v_left[k]);
				}

				// set trace entries
				if (v_diag == D[i][j])
				{
					TraceD[i][j].push_back(  TMatrixEntryElem(*diag_arrow) );
        			}
        		
				for (int k=1; k<=i; k++)
				{
					if (v_up[k] == D[i][j])
					{
						TraceD[i][j].push_back(  TMatrixEntryElem(*up_arrow,k) );
	    				}
	    			}
	    		
				for (int k=1; k<=j; k++)
				{
					if (v_left[k] == D[i][j])
					{
						TraceD[i][j].push_back(  TMatrixEntryElem(*left_arrow,k) );
      					}
      				}
			}
			else
			{

				TraceD[i][j].push_back(  TMatrixEntryElem(*diag_arrow) );

			}
			
			if (D[i][j] > maxD)
			{
				max_i = i;
				max_j = j;
				maxD = D[i][j];
			}

		}
	}

	#ifndef NDEBUG
		cout<<D;
	#endif


	// ------------------------------------------------------------
	// do backtrace
	BackTracer bt;
  
	// ---------------------------------------------------------
	// calc traceback
	//
	if (local)
		bt.calc_traceback(TraceD,max_i,max_j,1);   // 1 mean stop at 0;
	else
		bt.calc_traceback(TraceD,n,m,0);           // 0 means stop at cell (0,0)
  
	// ---------------------------------------------------------
	// print D and traceback
	// we have to calculate how many characters are printed maximally 
	// in lines and columns
	// 
	PrintArray pr(1+2*(n+1),2+(m+1)*(ROWFILL+PRECISION));

	if (!quiet)
	{
		// first, the trace is drawn on D (connecting all cells in
		// backtrace with
		//
		bt.draw_traceback(pr);
		
		// then draw the matrix overwriteing part of the trace
		//
		D.draw(pr);
    
		// now print the distance matrix and trace stored in pr
		//
		pr.print();
	}
	
	// finally, print out alignment
	//
	//bt.print_alignment(a,b);
  	//bt.print_alignment_show_contexta(a, b);
  	
  	// finally, store the alignment into a struct --- AlignResult and used later

  	alignresult.qseq = queueseq.sequence;
  	alignresult.tseq = targstr.sequence;
  	
  	string gapqseq="";
  	string gaptseq="";
  	
  	bt.store_alignment(gapqseq, gaptseq, a, b);
  	
  	alignresult.gapqseq = gapqseq;
  	alignresult.gaptseq = gaptseq;


	//cout <<"Q::" << gapqseq << endl;
	//cout <<"T::" << gaptseq << endl;
// Yang zhang defined ranking scheme 
	alignresult.rawscore =D[n][m];
	L_partitial(gapqseq,gaptseq,alignresult);
  	
	//== Get score-vector of all aligned positions and relate self-reference score-vector
  	/*vector<vector<int > > alignedpair;
  	GetAlignedPair(alignedpair, gapqseq, gaptseq);
  	vector<double> alignedscore;
  	alignresult.Vs =GetAlignedScore(alignedscore, alignedpair, selfsm);
	//cout <<"maxD " << maxD << endl; MaxD 记录的是最大值
  	
  	vector<double> refscore;
  	vector<vector<double > > refsm;
  	GetSelfSM(refsm, targstr, targstr, factor);
  	vector<vector<int > > refpair;
  	for(int i=0; i<alignedpair.size(); i++)
  	{
  		vector<int> vtemp;
  		vtemp.push_back(alignedpair[i][1]);
  		vtemp.push_back(alignedpair[i][1]);
  	//cout << alignedpair[i][0] << '-' << alignedpair[i][1] << ' ';	
  		refpair.push_back(vtemp);
  		vtemp.clear();
  	}
   	alignresult.Vr =GetAlignedScore(refscore, refpair, refsm);
   	
   	alignresult.Nalign=alignedpair.size();
   	alignresult.Vscore=alignedscore;
   	alignresult.Vref=refscore;*/
	// ---------------------------------------------------------
	// done
	
	
	selfsm.clear();
	delete diag_arrow;
	delete left_arrow;
	delete up_arrow;
	delete scorematrix;

	
	return 0;
}
//gapqseq:EEEEEECCEEEEEEECCCCEEEEEECCCEEEECCCCCHHHHHHHHCCCCHHHHHHHHCCCCCCCEEECCCCEEEEEEEEC
//gaptseq:EEEEEECCEEEEEEECCCC-EEEEECCC----CCEEHHHHHHHCCCCHHHHHHHHHHHCCCCCCEEECCCCCCCCEEEEC
//==============================================================================
//THE NEEDLEMAN-WUNSCH ALGORITHM, THIS IS THE THIRD VERSION BY MIAOZHICHAO
//BECAUSE THE PREVIOUS VERSIONS DO NOT PERFORM WELL.
//IN PUT THE GAP OPEN, GAP EXTENSION, SEQUENCE 1 AND SEQUENCE 2, 
//CONSERVATION GAP OPEN AND CONSERVATION GAP EXTENSION
//TWO STRINGS OF DSSP PREDICTION RESULTS
//AND TWO VECTORS OF CONSERVATION VALUES
//IT WILL GIVE THE ALIANMENT RESULTS
//==============================================================================
cost_t NwAlign3(RV2 &score,cost_t gap_open,cost_t gap_extn,
                        //GAP OPEN, GAP EXTENSION
                        std::string seq1,std::string seq2,
                        //SEQUENCE 1, SEQUENCE  2
                        cost_t gpc1,cost_t gpc2,
                        //CONSERVATION GAP OPEN AND CONSERVATION GAP EXTENSION
                        std::string dsse1,std::string dsse2,
                        //TWO STRINGS OF DSSP PREDICTION RESULTS
                        std::string &seq1_al,std::string &seq2_al,
                        //ALIANMENT RESULTS
                        RV1 &con1,RV1 &con2)
                        //TWO VECTORS OF CONSERVATION VALUES
{
	int i,j,k;
	seq1_al.clear();
	seq2_al.clear();
	int nseq1=seq1.length();
	int nseq2=seq2.length();
	RV1 rtmp(nseq2+1,0.0);
	RV2 val(nseq1+1,rtmp),preV,preH;preV=preH=val;
	IV1 itmp(nseq2+1,0);
	IV2 idir(nseq1+1,itmp),jpV,jpH;jpV=jpH=idir;
	
	cost_t D,V,H,val1,val2,score0;
	int it;
	IV1 j2i(nseq2+1,-1);
	//INITIALIZE
	val[0][0]=0.0;
	for(i=1;i<=nseq1;i++)
	{
		val[i][0]=gap_extn*i;
		preV[i][0]=val[i][0];
		idir[i][0]=0;
		jpV[i][0]=1;
		jpH[i][0]=i;
	}
	for(i=1;i<=nseq2;i++)
	{
		val[0][i]=gap_extn*i;
		preH[0][i]=val[0][i];
		idir[0][i]=0;
		jpV[0][i]=i;
		jpH[0][i]=1;
	}
    //ALIGN
    for(j=1;j<=nseq2;j++)
    {
		for(i=1;i<=nseq1;i++)
		{
			D=val[i-1][j-1]+score[i-1][j-1];
			jpH[i][j]=1;
            val1=val[i-1][j]+gap_open*(1+gpc1*con1[i-1]);
            val2=preH[i-1][j]+gap_extn*(1+gpc2*con1[i-1]);
            if(val1>val2) H=val1;
            else 
            {
				H=val2;
				if(i>1)jpH[i][j]=jpH[i-1][j]+1;
            }
            jpV[i][j]=1;
            val1=val[i][j-1]+gap_open*(1+gpc1*con2[j-1]);
            val2=preV[i][j-1]+gap_extn*(1+gpc2*con2[j-1]);
            if(val1>val2) V=val1;
            else
            {
				V=val2;
				if(j>1)jpV[i][j]=jpV[i][j-1]+1;
			}
            preH[i][j]=H;
            preV[i][j]=V;
            if(dsse1[i-1]==dsse2[j-1]&&(dsse1[i-1]=='H'||dsse1[i-1]=='E'))
            {
				idir[i][j]=1;
				val[i][j]=D;
			}
            else if(D>H&&D>V)
            {
				idir[i][j]=1;
				val[i][j]=D;
			}
            else if(H>V)
            {
				idir[i][j]=2;
				val[i][j]=H;
			}
            else
            {
				idir[i][j]=3;
				val[i][j]=V;
			}
		}
	}
	score0=val[nseq1][nseq2];
	for(j=1;j<=nseq2;j++)
	j2i[j]=-1;
	i=nseq1;
	j=nseq2;
	//TRACEBACK
	while(i>0&&j>0)
	{
		if(idir[i][j]==1)
		{
			j2i[j]=i;
			seq1_al+=seq1[i-1];
			seq2_al+=seq2[j-1];
			i--;
			j--;
			
		}
		else if(idir[i][j]==2)
		{
			it=jpH[i][j];
			for(k=1;k<=it;k++)
			{
				if(i>0)
				{
					seq1_al+=seq1[i-1];
					seq2_al+='-';
					i--;					
				}
			}
		}
		else
		{
			it=jpV[i][j];
			for(k=1;k<=it;k++)
			{
				if(j>0)
				{
					seq1_al+='-';
					seq2_al+=seq2[j-1];
					j--;					
				}
			}
		}
	};
	if(i>0)
	{
		for(k=i-1;k>=0;k--)
		{
			seq1_al+=seq1[k];
			seq2_al+='-';
		}
	}
	else if(j>0)
	{
		for(k=j-1;k>=0;k--)
		{
			seq1_al+='-';
			seq2_al+=seq2[k];
		}
	}
	reverse(seq1_al.begin(),seq1_al.end());
	reverse(seq2_al.begin(),seq2_al.end());

	return score0;
}
void InitialBlock(string gapseq,Block & b)
{
	int i,j,k;

	const char* chstr= gapseq.c_str();
	int len = gapseq.size();

	vector<int> Flag;

	for(i=0; i<len; i++)
	{
		Flag.push_back(0);
	}

	for(i=0; i<len; i++)
	{

		if(chstr[i]=='-' && Flag[i] == 0)
		{
			int number = 0;
			b.state.push_back(-1);
			b.initpos.push_back(i);
			for(j=0; j<len-i;j++)
			{
				if(chstr[i+j]=='-')
			       {	
					number ++;
					Flag[i+j] =1;
					if((i+j) == len-1)
                                	{
						b.length.push_back(number);
				       	b.endpos.push_back(i+j);
						break;
					}
				}
				else
                                {
					b.length.push_back(number);
				      b.endpos.push_back(i+j-1);
					break;
				}
			}				
		}
		if(chstr[i]!='-' && Flag[i] == 0)
		{			
			int number = 0;
			b.state.push_back(0);
			b.initpos.push_back(i);
			for(j=0; j<len-i;j++)
			{
				if(chstr[i+j]!='-')
			       {	
					number ++;
					Flag[i+j] =1;
					if((i+j) == len-1)
                                	{
						b.length.push_back(number);
				       	b.endpos.push_back(i+j);
						break;
					}
				}
				else
                                {
					b.length.push_back(number);
				       b.endpos.push_back(i+j-1);
					break;
				}
			}			
		}
	}


}
void SW_L_partitial_full( string gapqseq,string gaptseq,AlignResult& alignresult)
{
	
	int q1 = alignresult.qseq.length();
	int q2 = alignresult.tseq.length();
	int AL =gapqseq.size();

	int al1=0;int al2=0;
	int g1=0;int g2=0;

	for(int i=0; i<AL; i++)
	{
		if(gapqseq[i]!='-') al1 ++;
		if(gaptseq[i]!='-') al2 ++;
	}
	g1 =q1-al1;	g2 =q2-al2;

	alignresult.L_full = AL +g1 +g2;
	alignresult.L_partitial2 =alignresult.L_full-g1;
	alignresult.L_partitial =alignresult.L_full-g2;
	alignresult.ZhNormScore1 =(alignresult.rawscore)/(alignresult.L_full);
	alignresult.ZhNormScore2 =(alignresult.rawscore)/(alignresult.L_partitial);		
}
void SW_full_Align(int max_i,int max_j, AlignResult& alignresult)
{
	//cout << "qseq ="<< alignresult.qseq<< endl;
	//cout << "tseq ="<< alignresult.tseq<< endl;

	//cout << "max_i ="<< max_i<< "\tmax_j ="<< max_j<< endl;
	//cout << "gapqseq ="<< alignresult.gapqseq<< endl;
	//cout << "gaptseq ="<< alignresult.gaptseq<< endl;

	const char* chstr1=alignresult.gapqseq.c_str();
	const char* chstr2=alignresult.gaptseq.c_str();

	int L_gapqseq = alignresult.gapqseq.size();
      int L_gaptseq = alignresult.gaptseq.size();
      //cout << "L_gapqseq  " << L_gapqseq <<";L_gaptseq  " << L_gaptseq << endl;
    
      //caculate the gap numbers of the gapqseq and gaptseq
      int G_gapqseq = 0,G_gaptseq = 0;
      for(int j=0; j<L_gapqseq; j++)
	{
		if(chstr1[j]=='-') G_gapqseq++;
		if(chstr2[j]=='-') G_gaptseq++;    
        }
      //cout << "G_gapqseq  " <<  G_gapqseq <<";G_gaptseq  " << G_gaptseq << endl;

      //caculate the start rserial[] of the local aligment
      int gapqseq_start ,gaptseq_start ;
       //gapqseq_start = max_i-(L_gapqseq-G_gapqseq)+1;
       //gaptseq_start = max_j-(L_gaptseq-G_gaptseq)+1;
      gapqseq_start = max_i-(L_gapqseq-G_gapqseq);
      gaptseq_start = max_j-(L_gaptseq-G_gaptseq);
    	//cout << "max_i  " << max_i<<"  ,max_j  " << max_j <<endl; 
    	//cout << "gapqseq_start  " << gapqseq_start <<"  ,gaptseq_start  " << gaptseq_start <<endl; 

	string QB ="";string QE ="";
	string TB ="";string TE ="";
	for(int j=0; j< gaptseq_start; j++)
	{
		QB +='-';
		TB +=alignresult.tseq[j];
	}
	for(int j=0; j<gapqseq_start; j++)
	{
		QB +=alignresult.qseq[j];
		TB +='-';
	}
	int Q_L =alignresult.qseq.size();
	int T_L =alignresult.tseq.size();

	for(int j=max_i; j<Q_L; j++)
	{
		QE +=alignresult.qseq[j];
		TE +='-';
	}
	for(int j=max_j; j<T_L; j++)
	{
		QE +='-';
		TE +=alignresult.tseq[j];
	}


	alignresult.gapqseq = QB +alignresult.gapqseq +QE;
    	alignresult.gaptseq = TB +alignresult.gaptseq +TE;

}
void L_partitial( string gapqseq,string gaptseq,AlignResult& alignresult)
{
 	Block b1;Block b2;
	InitialBlock( gapqseq,b1);
	InitialBlock( gaptseq,b2);
	int query_end_gap =0;
	int query_begin_gap =0;
	int exchange_gap =0;
	
	
 //except	begining end and ending gap for the query
//Q::---------------------MFKVYGYDSNIHKCGPCDNAKRLLTV-KKQPFEFINIMPEKGVFDDEKIAELLTKLGRDTQIGLTMPQVFAPDGSHIGGFDQLREYFK----
//T::GKGDPKKPRGKMSSYAFFVQTSREEHKKKHPDASVNFSEFSKKCSERWKTMS-----------AKEKGKFEDMAKADKARYE-------------------REMKTYIPPKGE
	int len = b1.state.size();		
	int len2 = b2.state.size();
	if(b1.state[0] == -1)
	{
		query_begin_gap =b1.length[0];	
	}	
	if(b1.state[len-1] == -1)
	{
		query_end_gap =b1.length[len-1];	
	}

	bool flag1 = (len>1 && (b1.state[len-1] == 0) && (b1.state[len-2] == -1) );	
	bool flag2 =true;
	if(flag1 )
	{
		for(int j=b1.initpos[len-1];j<=b1.endpos[len-1];j++)
		{
			if(gaptseq[j] !='-')  
			{
				flag2 =false;break;
			}
		}
		if(flag2 ) { exchange_gap =b1.length[len-2];}
	}


//NDKYSSRVVRVISAKRQLVSGIKYILQVEIGRTTCPKSSGDLQSCEFHDEPEMAKYTTCTFVVYSIPWLNQIKLLESK--------CQ
//SLPFQKIQHSITAQDHQPTPDSCIISMVVGQLKA-----------------DEDPIMGFHQMFLLKNINDAWVCTNDMFRLALHNF--
//-----CQ  转化为
//NDKYSSRVVRVISAKRQLVSGIKYILQVEIGRTTCPKSSGDLQSCEFHDEPEMAKYTTCTFVVYSIPWLNQIKLLESKCQ--------
//SLPFQKIQHSITAQDHQPTPDSCIISMVVGQLKA-----------------DEDPIMGFHQMFLLKNINDAWVCTNDM--FRLALHNF
	alignresult.L_full =gapqseq.size();
	alignresult.L_partitial =alignresult.L_full-query_end_gap-query_begin_gap;
	alignresult.L_partitial2 =alignresult.L_partitial-exchange_gap;
	alignresult.ZhNormScore1 =(alignresult.rawscore)/(alignresult.L_full);
	alignresult.ZhNormScore2 =(alignresult.rawscore)/(alignresult.L_partitial);
	//cout << "rawscore ="<< alignresult.rawscore<< endl;
	//cout << "query_begin_gap ="<< query_begin_gap<< endl;
	//cout << "query_end_gap ="<< query_end_gap<< endl;
	//cout << "exchange_gap ="<< exchange_gap<< endl;
}
//================================================================================
// Record the position of all matched pairs in the query sequence
//	and the template sequence;
//================================================================================
void GetAlignedPair(vector<vector<int > >& alignedpair, string gapqseq, string gaptseq)
{
	int i1=-1, i2=-1;
	const char* chstr1=gapqseq.c_str();
	const char* chstr2=gaptseq.c_str();	
	
	for(int i=0; i<gapqseq.size(); i++)
	{
		if(chstr1[i]!='-') i1++;
		if(chstr2[i]!='-') i2++;
		
		if(chstr1[i]!='-' && chstr2[i]!='-')
		{
			vector<int> vtemp;
			vtemp.push_back(i1);
			vtemp.push_back(i2);
			
			alignedpair.push_back(vtemp);
			vtemp.clear();
		}
	}
	
}



//=========================================================================
// Get the aligned score from the score-matrix according to
// 	the index information in alignedpair[][]
//=========================================================================
double GetAlignedScore(vector<double>& alignedscore, vector<vector<int > >& alignedpair, vector<vector<double > >& selfsm)
{
	int i=0;
	int ti=-1, tj=-1;
	double all =0.0;
	for(i=0; i<alignedpair.size(); i++)
	{
		ti=alignedpair[i][0];
		tj=alignedpair[i][1];
		
		alignedscore.push_back(selfsm[ti][tj]);
		//cout << alignedscore[i] << ' ';
		all +=selfsm[ti][tj];
	}
	if(all <0.001) all=0.001;
	return all;
	//cout << "GetAlign Score =" << all <<endl;
	//cout << endl;
}

//=========================================================================
// Calculate all the Z-scores of query-template alignments
// Z(V) = (Z-<Z>)/(sqrt(<Z^2>-<Z>^2));
//=========================================================================
void CalNormScoreVector(vector<double>allscore,vector<double>& allnormscore)
{
       int num=allscore.size(); 
       
       double mean=0.0, mean2=0.0;
       int i;
	for(i=0; i<num; i++)
	{
		mean += allscore[i];
		mean2 += allscore[i]*allscore[i];
	}
	mean /= num;
	mean2 /= num;
	double delta = sqrt(mean2-mean*mean);
	
	for(i=0; i<num; i++)
	{
		allnormscore.push_back((allscore[i]-mean)/delta);
	}	

}

//=========================================================================
// Calculate all the Z-scores of query-template alignments
// Z(V) = (Z-<Z>)/(sqrt(<Z^2>-<Z>^2));
//=========================================================================
void CalNormScore(vector<AlignResult>& allalignresult)
{
	int Nquery=allalignresult[0].qseq.size();
	int i=0;
	int num=allalignresult.size();
	
	
	vector<double> allscore;
	double tempscore=0.0;
	for(i=0; i<num; i++)
	{
		tempscore=PearsonCoefficient(allalignresult[i].Vscore, allalignresult[i].Vref);
		tempscore=tempscore*allalignresult[i].Nalign/Nquery;
		
		allscore.push_back(tempscore);
	}
	
	double mean=0.0, mean2=0.0;
	for(i=0; i<num; i++)
	{
		mean += allscore[i];
		mean2 += allscore[i]*allscore[i];
	}
	mean /= num;
	mean2 /= num;
	double delta = sqrt(mean2-mean*mean);
	
	for(i=0; i<num; i++)
	{
		allalignresult[i].normscore = (allscore[i]-mean)/delta;
	}
}



//=========================================================================
// Select the top selectN number of query-template alignments
// 	according to the allalignresult[i].normscore value;
//=========================================================================
void OptTemplate(vector<int>& TopIndex, vector<AlignResult>& TopAlignResult, vector<AlignResult>& allalignresult, int selectN)
{
	int i=0, j=0;
	vector<double> allnormscore;
	int num=allalignresult.size();
	
	for(i=0; i<num; i++)
	{
		allnormscore.push_back(allalignresult[i].normscore);
	}
	
	
	int ti=-1;
	for(i=0; i<selectN; i++)
	{
		ti=-1;
		ti=(int) (max_element(allnormscore.begin(), allnormscore.end())-allnormscore.begin());
		
		TopIndex.push_back(ti);
		TopAlignResult.push_back(allalignresult[ti]);
		
		allnormscore[ti]=-1.0e20;
		
		if(allnormscore.empty())
			break;
	}
	allnormscore.clear();
	
}
//=========================================================================
// Select the top selectN number of query-template alignments
// 	according to the normscore value;
//=========================================================================

void TopNIndex(vector<int>& TopIndex,vector<double>allnormscore,int selectN)
{
       int i=0, j=0;	
	int num=allnormscore.size();

       int ti=-1;
	for(i=0; i<selectN; i++)
	{
		ti=-1;
		ti=(int) (max_element(allnormscore.begin(), allnormscore.end())-allnormscore.begin());
		
		TopIndex.push_back(ti);				
		allnormscore[ti]=-1.0e20;
		
		if(allnormscore.empty())
	       break;
	}
	
}


int ReadTemplateID(vector<string>& allID, const char* IDfile)
{
	string ID="";
	
	ifstream fin(IDfile);
	if( !fin.is_open() )
	{
		cout << "Can not open " << IDfile << endl;
		exit(-1);
	}
	
	while( !fin.eof() )
	{
		ID="";
		
		fin >> ID;
		
		if(ID!="")
		{
			allID.push_back(ID);
		}
	}
	fin.close();
	
	
	return 0;
}




