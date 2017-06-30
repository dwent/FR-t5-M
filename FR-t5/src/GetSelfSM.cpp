//===================================================================================
//
// Get N*M score matrix defined by oneself
//
// AUTHORS:     Wu Aiping
//
// DESCIRPTION: This score matrix is asymmetrically, N is the length
//		   of sequence 1 and M is the length of sequence 2.
//		   (Aug, 2007)
////===================================================================================
//
// Implementation of Threading Sequence Alignment
//
// AUTHORS:     Huyun
//
// DESCIRPTION: Training and Testing Datasets
//               LOCAL/GlOBAL alignment. (Feb, 2009)
//===================================================================================
//===================================================================================
				/**********************************************************
				Copyright(c)2012, IBP, CHINESE ACADEMY OF SCIENCE
				All rights reserved.
				NAME:		GETSELFSM.cpp
				ABSTRACT:	THREADING GETSELFSM
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


#include <fstream>
#include "GetSelfSM.hpp"

using namespace std;


/*
//=====================================================================================================
//Read parameter file of 20*20*20*(6*6*6*6)= 20*20*20*1296 bins of Etrp table;
//For a 3-residues local fragment, there are 20*20*20 kinds of amino acid
//	compositions, and cut Phi, Psi dihedrals into 6 bins by 60-degrees per-cell,
//	and cut two angles of b(i)-b(i+2) and P(i)-P(i+2) into 6 bins by 30-degrees
//	per-cell;
//=====================================================================================================
int ReadTrpPar(DV4& trptable, const char* parfile)
{
	int a=0, b=0, c=0;
	int i=0, j=0;
	
	ifstream fin(parfile);
	if( !fin.is_open() )
	{
		cout << "\nReadTrpPar()--Error: can not open " << parfile << endl;
		exit(-1);
	}
	
	
	string data[1296];
	for(i=0; i<8000; i++)
	{	
		a=i/400;
		b=(i-400*a)/20;
		c=i-400*a-20*b;
		
		for(j=0; j<1296; j++)
		{
			fin >> data[j];
			
			trptable[a][b][c][j] = atof(data[j].c_str());
		}
	}
	fin.close();
	
	
	return 0;
}




//=============================================================================================
//Read protein structure pdb file
//Substract the infomation of four-beads (N, CA, C, O) of each residue
//=============================================================================================
int Read4AtomPdb(ChainCoord& chain, const char* pdbfile)
{
	int i=0, j=0, k=0, l=0;
	
	ifstream fin(pdbfile);
	if( !fin.is_open() )
	{
		cout << "\nError: can not open " << pdbfile << endl;
		return -1;
	}
	
	int atomindex=0;
	char buf[100];
	while( !fin.eof() )
	{
		strcpy(buf, "");
		fin.getline(buf, 100);
		
		if(strlen(buf)>10)
		{
			atomindex++;
		}
	}
	fin.close();
	
	if(atomindex<4)
	{
		cout << "\nNo residue in " << pdbfile << endl;
		return -2;
	}
	
		
	PDB pdbobj;
	int atomnum=pdbobj.AtomNum(pdbfile);
	int resinum=pdbobj.ResiNum(pdbfile);
	
	chain.resinum=resinum;
	chain.atomnum=atomnum;
	chain.beginxyz[0]=0, chain.beginxyz[1]=0, chain.beginxyz[2]=0;
	chain.endxyz[0]=1, chain.endxyz[1]=1, chain.endxyz[2]=1;
	chain.endxyz[3]=2.3, chain.endxyz[4]=2.3, chain.endxyz[5]=2.3;

	
	
	double x[atomnum], y[atomnum], z[atomnum];
	string aname[atomnum], rname[atomnum];
	int rserial[atomnum];
	
	pdbobj.ReadCoord(pdbfile, x, "XCOORD");
	pdbobj.ReadCoord(pdbfile, y, "YCOORD");
	pdbobj.ReadCoord(pdbfile, z, "ZCOORD");
	pdbobj.ReadName(pdbfile, aname, "ATOMNAME");
	pdbobj.ReadName(pdbfile, rname, "RESINAME");
	pdbobj.ReadSerial(pdbfile, rserial, "RESISERIAL");

	
	vector<vector<double > > resi;
	vector<double> atom;
	int weight=1;
	double radius=0.0;
	double cc[3];
	for(i=0; i<atomnum; i++)
	{
		atom.clear();
		if(aname[i]=="N")
		{
			if(i>0)
			{
				chain.xyz.push_back(resi);
				resi.clear();
			}
			
			atom.clear();
			atom.push_back(x[i]);
			atom.push_back(y[i]);
			atom.push_back(z[i]);
			
			resi.push_back(atom);
						
			chain.resiname.push_back(rname[i]);
			chain.resiserial.push_back(rserial[i]);
			
		}
		else if(aname[i]=="CA" || aname[i]=="C" || aname[i]=="O")
		{
			atom.push_back(x[i]);
			atom.push_back(y[i]);
			atom.push_back(z[i]);
			
			resi.push_back(atom);
		}
		
		if(i==atomnum-1)
		{
			chain.xyz.push_back(resi);
			resi.clear();
		}
		
	}

	return 0;	
	
}

*/
//===================================================================//
//Read protein sequences that to be predicted
//Fasta file predefined by user in Config file
//File format: the standard Fasta sequence file
//===================================================================//
void ReadNormFasta(const char* fastafile, char* sequence)
{
	ifstream fin(fastafile);
	if( !fin.is_open() )
	{
		cout << "\nCan not open fasta file in ReadFasta()!!\n";
		exit(-1);	
	}
	
	char buf[2000];
	while( !fin.eof() )
	{
		strcpy(buf, "");
		fin.getline(buf, 2000);
		
		if(buf[0]!='>' && strlen(buf)>0)
		{
			strcat(sequence, buf);
		}
	}
	
	fin.close();
}

//================================================================================
// Fill the PSSM information of the SeqStr struct, 
// The PSSM information read from the output file of the PSI-blast
//================================================================================
void PSSMofSeqStr(SeqStr& seqstr, const char* pssmfile)
{
	string sequence;
	vector<vector<double > > freq;
	vector<vector<double > > prof;
	ReadPSSM(sequence, freq, prof, pssmfile);
	
	seqstr.resinum=sequence.size();
	seqstr.sequence=sequence;	
	seqstr.freq=freq;
	seqstr.prof=prof;
	
	freq.clear();
	prof.clear();
}
void SubPSSMofSeqStr(vector<int> Subsite,SeqStr& seqstr, const char* pssmfile)
{
	string sequence;
	vector<vector<double > > freq;
	vector<vector<double > > prof;
	ReadPSSM(sequence, freq, prof, pssmfile);
	
	string Subsequence;
	vector<vector<double > > Subfreq;
	vector<vector<double > > Subprof;
	for(int i=0; i<Subsite.size(); i++)
	{
		int index =Subsite[i];
		Subsequence +=sequence[index];
		Subfreq.push_back(freq[index]);
		Subprof.push_back(prof[index]);
	}
	seqstr.resinum =Subsite.size();
	seqstr.sequence=Subsequence;	
	seqstr.freq=Subfreq;
	seqstr.prof=Subprof;
	
	freq.clear();
	prof.clear();
	
	Subfreq.clear();
	Subprof.clear();
}


//================================================================================
// Fill the Second Structure (SS) information of the SeqStr struct, 
// The SS information read from the output file of the PsiPred
//================================================================================
void PSIPREDofSeqStr(SeqStr& seqstr, const char* secondfile)
{
	vector<char> second;
	ReadPSIPRED(second, secondfile);
	
	for(int i=0; i<second.size(); i++)
	{
		seqstr.second.push_back(second[i]);
	}	

	second.clear();
}

void SubPSIPREDofSeqStr(vector<int> Subsite,SeqStr& seqstr, const char* secondfile)
{
	vector<char> second;
	ReadPSIPRED(second, secondfile);
	
	for(int i=0; i<Subsite.size(); i++)
	{
		int index =Subsite[i];
		seqstr.second.push_back(second[index]);
	}	

	second.clear();
}

//================================================================================
// Fill the Second Structure (SS) information of the SeqStr struct, 
// The SS information read from the output file of the DSSP
//================================================================================
void DSSPofSeqStr(SeqStr& seqstr, const char* secondfile)
{
	vector<char> second;
	ReadDSSP(second, secondfile);
	
	for(int i=0; i<second.size(); i++)
	{
		seqstr.second.push_back(second[i]);
	}
	
	if(seqstr.second.size()==seqstr.resinum-1)	//ATTENTATION: the last residue has no secondary structure in DSSP file.
	{
		seqstr.second.push_back('C');
	}
	
	second.clear();
}



//================================================================================
// Fill the StrFreq information of the SeqStr struct, 
// The StrFreq information read from the pre-calculated strfreq file
//================================================================================
void STRFREQofSeqStr(SeqStr& seqstr, const char* strfreqfile)
{
	//get frequency-distributions based on local 9-residues' strucutres
	vector<vector<double > > strfreq;
	ReadStrFreq(strfreq, strfreqfile);
	
	seqstr.strfreq=strfreq;
	strfreq.clear();
}


//================================================================================
// Fill the 3-residue angle index information of the SeqStr struct, 
// These information transfered from the local 3-residues' structures;
//================================================================================
void TRPofSeqStr(SeqStr& seqstr, ChainCoord& chain)
{
	//get the 4-angles index reprenstation (0,1,2,3,4,5) in Etrp of each local 3-residues' structure
	vector<vector<int > > angleindex;
	Get4AngleIndex(angleindex, chain);
	
	seqstr.angleindex=angleindex;
	angleindex.clear();
	
}




//======================================================================
// Read PSSM profile that generated from PSI-BLAST;
//
//======================================================================
int ReadPSSM(string& sequence, vector<vector<double > >& freq, vector<vector<double > >& prof, const char* pssmfile)
{
	int i=0, j=0, k=0;
	
	char buf[200];
	
	ifstream fin(pssmfile);
	if( !fin.is_open() )
	{
		cout << "\nReadPSSM()--Error: can not open " << pssmfile << endl;
		exit(-1);
	}
	
	int index=0;
	vector<double> tempscore;
	vector<double> tempfreq;
	while( !fin.eof() )
	{
		fin.getline(buf, 200);
		
		if( index>=3 && strlen(buf)>=160 )
		{
			tempscore.clear();
			tempfreq.clear();
			
			//get residue name
			sequence.push_back( buf[6] );
			
			//get score
			i=9;
			for(j=0; j<20; j++)
			{
				i=9+3*j;
				while( buf[i]==' ' )
					i++;
				tempscore.push_back( atof(&buf[i]) );
				
			}
			prof.push_back(tempscore);
			
			//get freq
			i=70;
			for(j=0; j<20; j++)
			{
				i=70+4*j;
				while( buf[i]==' ' )
					i++;
				tempfreq.push_back( atof(&buf[i])/100.0 );
			}
			freq.push_back(tempfreq);
		}
		
		index++;
	}
	fin.close();

	tempscore.clear();
	tempfreq.clear();
		
	return 0;
}
//==========================================================================
// READ TEMPLATE LIBRARY ID FROM THE ID FILE
// 
//==========================================================================
void ReadLibID(char* file,SV1 &ss)
{
	ss.clear();
	std::string buf;
	std::ifstream infile;
	infile.open(file);
	if(!infile){std::cerr<<"\nCan not open: "<<file<<'\n';exit(0);}
	while(!infile.eof())
	{
		buf.clear();
		getline(infile,buf);
		if(buf.empty())break;
		ss.push_back(buf);
	};infile.close();
}

//==========================================================================
// Read Secondary Structure file that generated from PSIPRED;
//
//==========================================================================
int ReadPSIPRED(vector<char>& second, const char* secondfile)
{
	//Read the PSIPRED predicted *.horiz file --- secondfile
	ifstream fin(secondfile);
	if( !fin.is_open() )
	{
		cout << "\nReadPSIPRED()--Error: can not open " << secondfile << endl;
		exit(-1);
	}
	
	
	//suggest the largest residue num of the chain is 2000
	int i=0;
	char Pred[2000]="";
	string record[2];
	while( !fin.eof() )
	{
		record[0]="";
		record[1]="";
		fin >> record[0] >> record[1];
		if(record[0]=="Pred:")
		{
			strcat(Pred, record[1].c_str());
		}
	}
	fin.close();
	
	
	//
	i=0;
	while( Pred[i]!='\0' )
	{
		second.push_back(Pred[i]);
		i++;
	}
	
	
	return 0;
}



//==========================================================================
// Read Secondary Structure file that generated from DSSP;
// Changed from Tianliqing's program.
//==========================================================================
int ReadDSSP(vector<char>& second, const char* secondfile)
{
	string stemp;
	
	ifstream infileDssp(secondfile);
	if(!infileDssp)
	{
		cerr<<"ReadDSSP()---error: can not open dssp file " << secondfile << endl;
		exit(-1);
	}
		
	do
	{
		infileDssp>>stemp;
		if(stemp!="#") 
		{
			getline(infileDssp,stemp);
			continue;
		}
		else
		{
			getline(infileDssp,stemp);
			break;
		}
		
	}while(!infileDssp.eof());
	
	int index=0;			
	do
	{
		stemp.clear();
		getline(infileDssp,stemp);
		
		if(stemp.size()==0) break;

		//if(stemp[13]=='!') continue;
		if(stemp[10]!=' ') continue;
		
		//if(stemp[13]=='!' && index==0)
		//{
		//	second.push_back('C');
		//}
		//else if(stemp[13]=='!' && index>0)
		//{
		//	second.push_back(second[index-1]);
		//}
		if(stemp[13]=='!')
		{
			continue;
		}
		else if((stemp[16]=='H')||(stemp[16]=='G')||(stemp[16]=='I'))
		{
			second.push_back('H');
		}
		else if(stemp[16]=='E')
		{
			second.push_back('E');
		}
		else
		{
			second.push_back('C');
		}
		
		index++;
		
	}while(!infileDssp.eof());
	infileDssp.close();

	
	return 0;
}



//==========================================================================
// Translate amino acid type into index in all 20 types;
// 20 amino acid types sorted as: ARNDCQEGHILKMFPSTWYV
//==========================================================================
int TransAA2Index(vector<int>& resiindex, string sequence)
{
	//amino acides' serials: {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
	//		    		 "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL" };
	int i=0, j=0;
	const char* AAI="ARNDCQEGHILKMFPSTWYV";
	
	bool flag=false;
	const char* seqch=sequence.c_str();
	for(i=0; i<sequence.size(); i++)
	{
		flag=false;
		for(j=0; j<20; j++)
		{
			if(seqch[i]==AAI[j])
			{
				resiindex.push_back(j);
				flag=true;
				break;
			}
		}
		
		if( !flag )
		{
			cout << "\nTransAA2Index()--Error: can not recognize amino acid type " << seqch[i] << endl;
			return(-1);
		}
	}
	
	
	return 0;
}



//==========================================================================
// Read frequency-distribution file generated based on local
//	 9-residues' structures;
//
//==========================================================================
int ReadStrFreq(vector<vector<double > >& strfreq, const char* strfreqfile)
{
	int i=0, j=0, k=0;
	
	ifstream fin(strfreqfile);
	if( !fin.is_open() )
	{
		cout << "\nReadStrFreq()--Error: can not open " << strfreqfile << endl;
		exit(-1);
	}
	

	string record[21];
	vector<double> tempfreq;
	while( !fin.eof() )
	{
		record[0]="";
		
		for(i=0; i<21; i++)
		{
			fin >> record[i];
		}
		
		if(record[0]=="" || record[0]=="RES") continue;
		
		tempfreq.clear();
		for(i=1; i<21; i++)
		{
			tempfreq.push_back( atof(record[i].c_str()) );
		}
		
		strfreq.push_back(tempfreq);
	}
	fin.close();

	tempfreq.clear();
	
	
	return 0;
}



//==========================================================================
// According to the pdb structure of the template, get each
//	 local 3-residue fragment's 4 angle representation
// ATTENTATION: N-terminal and C-terminal assigned all index
//	 with -1, not used later;
//==========================================================================
int Get4AngleIndex(vector<vector<int > >& angleindex, ChainCoord& chain)
{
	int i=0, j=0, k=0;
	double angle[4];
	double xyz[12][3];
	
	//At the N-terminal and C-terminal of template structure,
	//	assign the angle[4] all zero, these values will not used
	
	int tempi=-1;
	vector<int> vangle;
	
	//i=0;
	vangle.push_back(-1);
	vangle.push_back(-1);
	vangle.push_back(-1);
	vangle.push_back(-1);
	angleindex.push_back(vangle);
	vangle.clear();
	
	for(i=1; i<chain.resinum-1; i++)
	{
		vangle.clear();
		
		if(chain.xyz[i-1].size()!=4 || chain.xyz[i].size()!=4 || chain.xyz[i+1].size()!=4 || chain.resiserial[i]!=chain.resiserial[i-1]+1 || chain.resiserial[i+1]!=chain.resiserial[i]+1)
		{
			vangle.push_back(-1);
			vangle.push_back(-1);
			vangle.push_back(-1);
			vangle.push_back(-1);
			angleindex.push_back(vangle);
			vangle.clear();
			
			continue;
		}
		
		for(j=0; j<12; j++)
		{//cout << i << ' ' << j << endl;
			xyz[j][0]=chain.xyz[j/4+i-1][j%4][0];
			xyz[j][1]=chain.xyz[j/4+i-1][j%4][1];
			xyz[j][2]=chain.xyz[j/4+i-1][j%4][2];
		}
			
		Get4Angle(angle, xyz);

		tempi = (int)floor(angle[0]+180)/60;
		vangle.push_back(tempi);
			
		tempi = (int)floor(angle[1]+180)/60;
		vangle.push_back(tempi);
			
		tempi = (int)floor(angle[2])/30;
		vangle.push_back(tempi);
			
		tempi = (int)floor(angle[3])/30;
		vangle.push_back(tempi);

		angleindex.push_back( vangle );
	}
	
	vangle.clear();

	//i=chain.resinum-1;
	vangle.push_back(-1);
	vangle.push_back(-1);
	vangle.push_back(-1);
	vangle.push_back(-1);
	angleindex.push_back(vangle);
	vangle.clear();

	
	return 0;
}


/*
//=========================================================================================
//Substract 2 dihedrals and 2 angles from a local 3-residues' fragment
//For fragment A(i)A(i+1)A(i+2), 2 dihedrals are Phi(i+1) and Psi(i+1),
//	2 angles are angle between P(i) and P(i+2), angle between b(i) and
//	b(i+2)
//=========================================================================================
int Get4Angle(double angle[4], double xyz[12][3])
{
	int i=0, j=0, k=0, m=0;
	
	double dih=0.0;
	
	//==calculate the first dihedral ( Phi-A(i+1) )==//
	dih = Dihedral(xyz[2], xyz[4], xyz[5], xyz[6]);
	angle[0]=dih;
	
	//==calculate the second dihedral ( Psi-A(i+1) )==//	
	dih = Dihedral(xyz[4], xyz[5], xyz[6], xyz[8]);
	angle[1]=dih;
	

	//==calculate the angle between P-A(i) and P-A(i+2), as the third angle==//
	double x1, y1, z1, x2, y2, z2, x3, y3, z3;
	x1=xyz[1][0], y1=xyz[1][1], z1=xyz[1][2];
	x2=xyz[0][0], y2=xyz[0][1], z2=xyz[0][2];
	x3=xyz[2][0], y3=xyz[2][1], z3=xyz[2][2];
	
	double A, B, C;
	A = (y2-y1)*(z3-z1)-(z2-z1)*(y3-y1);
	B = -( (x2-x1)*(z3-z1)-(z2-z1)*(x3-x1) );
	C = (x2-x1)*(y3-y1)-(y2-y1)*(x3-x1);
	
	double DummyP0[3], RealP0[3];
	DummyP0[0]=x1+1.0;
	DummyP0[1]=(B/A)*(DummyP0[0]-x1)+y1;
	DummyP0[2]=(C/A)*(DummyP0[0]-x1)+z1;
	RealP0[0]=x1, RealP0[1]=y1, RealP0[2]=z1;
	
	x1=xyz[9][0], y1=xyz[9][1], z1=xyz[9][2];
	x2=xyz[8][0], y2=xyz[8][1], z2=xyz[8][2];
	x3=xyz[10][0], y3=xyz[10][1], z3=xyz[10][2];

	A = (y2-y1)*(z3-z1)-(z2-z1)*(y3-y1);
	B = -( (x2-x1)*(z3-z1)-(z2-z1)*(x3-x1) );
	C = (x2-x1)*(y3-y1)-(y2-y1)*(x3-x1);
	
	double DummyP2[3], RealP2[3];
	DummyP2[0]=x1+1.0;
	DummyP2[1]=(B/A)*(DummyP2[0]-x1)+y1;
	DummyP2[2]=(C/A)*(DummyP2[0]-x1)+z1;
	RealP2[0]=x1, RealP2[1]=y1, RealP2[2]=z1;	
	
	//translate RealP0-DummyP0 to RealP2-DummyP2 by overlap RealP0 to RealP2
	double translation[3];
	translation[0]=RealP0[0]-RealP2[0];
	translation[1]=RealP0[1]-RealP2[1];
	translation[2]=RealP0[2]-RealP2[2];
	
	double aftertrans[3];
	aftertrans[0]=DummyP0[0]-translation[0];
	aftertrans[1]=DummyP0[1]-translation[1];
	aftertrans[2]=DummyP0[2]-translation[2];
	
	dih = Angle(aftertrans, RealP2, DummyP2);
	angle[2]=dih;
	

	//==calculate the angle between b-A(i) and b-A(i+2), as the 4th angle==//
	double ang0 = Angle(xyz[0], xyz[1], xyz[2]);
	double ang2 = Angle(xyz[8], xyz[9], xyz[10]);
	
	
	vector<double> point, axis1, axis2;
	point.push_back(xyz[0][0]);
	point.push_back(xyz[0][1]);
	point.push_back(xyz[0][2]);
	
	axis1.push_back(RealP0[0]);
	axis1.push_back(RealP0[1]);
	axis1.push_back(RealP0[2]);

	axis2.push_back(DummyP0[0]);
	axis2.push_back(DummyP0[1]);
	axis2.push_back(DummyP0[2]);

	vector<double> DummyB0;
	CoordinateRotation(point, axis1, axis2, ang0/2, DummyB0);


	point.clear();
	axis1.clear();
	axis2.clear();
	point.push_back(xyz[8][0]);
	point.push_back(xyz[8][1]);
	point.push_back(xyz[8][2]);
	
	axis1.push_back(RealP2[0]);
	axis1.push_back(RealP2[1]);
	axis1.push_back(RealP2[2]);

	axis2.push_back(DummyP2[0]);
	axis2.push_back(DummyP2[1]);
	axis2.push_back(DummyP2[2]);
	
	vector<double> DummyB2;
	CoordinateRotation(point, axis1, axis2, ang2/2, DummyB2);
	
	aftertrans[0]=DummyB0[0]-translation[0];
	aftertrans[1]=DummyB0[1]-translation[1];
	aftertrans[2]=DummyB0[2]-translation[2];
	
	double B2[3];
	B2[0]=DummyB2[0];
	B2[1]=DummyB2[1];
	B2[2]=DummyB2[2];
	
	dih = Angle(aftertrans, RealP2, B2);
	angle[3]=dih;	
	
	point.clear();
	axis1.clear();
	axis2.clear();
	DummyB0.clear();
	DummyB2.clear();
	
	return 0;
}
*/



//=============================================================================
// According to the threading objective function, pre-calculate
// 	all Eij scores of each site-pairs of two aligned sequence;
// The whole score-matrix is the N*M set of all single Eij;
// If there are several sets in the objective function, such as
// 	E0,E1,E2,E3,et al., so use which set and the weight of 
// 	set can be controlled by the factor[], the
// 	factor also include the gap open and extension penalties.
// Objective Function:
//	Eij = f0*Eseq,seq + f1*E2nd + f2*Etrp + f3*Eseq,str + f4
//=============================================================================
int GetSelfSM(vector<vector<double > >& selfsm, SeqStr& queueseq, SeqStr& targstr, double factor[8])
{
	int i=0, j=0, k=0;
	
	double f0, f1, f2, f3;
	double f4, f5, f6, f7;
	f0 = factor[0];
	f1 = factor[1];
	f2 = factor[2];
	f3 = factor[3];
	f4 = factor[4];

	f5 = factor[5];
	f6 = factor[6];
	f7 = factor[7];

	

	double Extra=1.0e-20;
	vector<vector<double > > E0;
	vector<vector<double > > E1;
	vector<vector<double > > E2;
	vector<vector<double > > E3;

	vector<vector<double > > E5;
	vector<vector<double > > E6;
	vector<vector<double > > E7;
	
	double score=0.0;
	vector<double> scoreline;
	double N=queueseq.resinum;
	double M=targstr.resinum;
	
	//get E0, the first energy set Eseq,seq
	if(f0>Extra || f0<-Extra)
	{
		for(i=0; i<N; i++)
		{
			scoreline.clear();
			for(j=0; j<M; j++)
			{
				score = Eseqseq(queueseq.freq[i], targstr.prof[j]);
				
				scoreline.push_back(score);
			}
			
			E0.push_back(scoreline);
		}
	}

	
	//get E1, the second structure comparative set E2nd
	if(f1>Extra || f1<-Extra)
	{
		for(i=0; i<N; i++)
		{
			scoreline.clear();
			for(j=0; j<M; j++)
			{
				score = E2nd(queueseq.second[i], targstr.second[j]);
				
				scoreline.push_back(score);
			}
			
			E1.push_back(scoreline);
		}
	}
	
	
	//get E2, the third set Etrp
	if(f2>Extra || f2<-Extra)
	{
		for(i=0; i<N; i++)
		{
			scoreline.clear();
			for(j=0; j<M; j++)
			{
				score = Etrpbin(i, j, queueseq.resiindex, targstr.angleindex);
				
				scoreline.push_back(score);
			}
			
			E2.push_back(scoreline);
		}
	}
	
	
	//get E3, the four set Eseq,str
	if(f3>Extra || f3<-Extra)
	{
		for(i=0; i<N; i++)
		{
			scoreline.clear();
			for(j=0; j<M; j++)
			{
				score = Eseqstr(queueseq.prof[i], targstr.strfreq[j]);
				
				scoreline.push_back(score);
			}
			
			E3.push_back(scoreline);
		}
	}

	//get E5, the five set HCA
	if(f5>Extra || f5<-Extra)
	{
		for(i=0; i<N; i++)
		{
			scoreline.clear();
			for(j=0; j<M; j++)
			{
				score = HCA(queueseq.sequence[i], targstr.sequence[j]);
				
				scoreline.push_back(score);
			}
			
			E5.push_back(scoreline);
		}
	}

	
	
	//get the final score matrix --- selfsm according to E0, E1, E2, E3 and all factors
	bool flag0=false, flag1=false, flag2=false, flag3=false;
	bool flag5=false, flag6=false, flag7=false;

	if( !E0.empty() ) flag0=true;
	if( !E1.empty() ) flag1=true;
	if( !E2.empty() ) flag2=true;
	if( !E3.empty() ) flag3=true;

	if( !E5.empty() ) flag5=true;
	if( !E6.empty() ) flag6=true;
	if( !E7.empty() ) flag7=true;
	
	if( !flag0 && !flag1 && !flag2 && !flag3 && !flag5 && !flag6 && !flag7)
	{
		cout << "\nGetSelfSM()--Warning: all energy sets have no elements!\n";
		return -1;
	}
	
	for(i=0; i<N; i++)
	{
		scoreline.clear();
		for(j=0; j<M; j++)
		{
			score = f4;
			if( flag0 ) score += f0*E0[i][j]; 
			if( flag1 ) score += f1*E1[i][j];
			if( flag2 ) score += f2*E2[i][j];
			if( flag3 ) score += f3*E3[i][j];
			if( flag5 ) score += f5*E5[i][j];
			if( flag6 ) score += f6*E6[i][j];
			if( flag7 ) score += f7*E7[i][j];
			//cout << score << ' ';
			scoreline.push_back(score);
		}
		//cout << endl;
		selfsm.push_back(scoreline);
	}
	
	
	scoreline.clear();
	
	E0.clear();
	E1.clear();
	E2.clear();
	E3.clear();
	E5.clear();
	E6.clear();
	E7.clear();
	
	return 0;
}





//================================================================================
// Calculate the similarity between two positions with their
//	Profile-information;
// freq vector is the position 1's frequency-distribution among
//	20 kinds of amino acides;
// prof vector is the position 2's profile-distribution among
//	20 kinds of amino acides;
//================================================================================
double Eseqseq(vector<double>& freq, vector<double>& prof)
{
	if(freq.size()!=20)
	{
		cout << "Eseqseq()--Wrong: the size of input vector<double>& freq != 20" << endl;
		exit(-1);
	}

	if(prof.size()!=20)
	{
		cout << "Eseqseq()--Wrong: the size of input vector<double>& prof != 20" << endl;
		exit(-2);
	}
	
	int i=0;	
	double score=0.0;
	
	for(i=0; i<20; i++)
	{
		score += freq[i]*prof[i];
	}
	
	
	return score;
}



//================================================================================
// Calculate the second structure aligned score of two positions
// 
//================================================================================
double E2nd(char origsecond, char targsecond)
{
	double score=0.0;
	
	if( origsecond == targsecond )
	{
		score = 1.0;
	}
	else
	{
		score = -1.0;
	}
	
	return score;
}

//================================================================================
// Calculate the Hydrophic scoring matrix of two positions
/* Referenfence:Author: Pedro J. Silva
Reference: Silva, P.J.
(2007) "Assessing the reliability of sequence similarities detected through hydrophobic sequence analysis"
 Proteins: Structure, Function and Bioinformatics, 70, 1588-1594. 
*/
//================================================================================
double HCA(char origseq, char targseq)
{
	double score=0.0;
 	string type1,type2;
//VILFYWM
// string type:"Hydrophobic" or "No"
	if(origseq =='V'||origseq =='I'||origseq =='L'||origseq =='F'||origseq =='Y'||origseq =='W'||origseq =='M')
	{
		type1 ="Hydrophobic";
	}
	else
	{
		type1 ="No";
	}

	if(targseq =='V'||targseq =='I'||targseq =='L'||targseq =='F'||targseq =='Y'||targseq =='W'||targseq =='M')
	{
		type2 ="Hydrophobic";
	}
	else
	{
		type2 ="No";
	}

	if(type1 =="Hydrophobic" && type2 =="Hydrophobic" )
	{
		score = 1.0;
		return score;
	}
	
	if(origseq ==targseq)
	{
		if( origseq =='P' )
		{
			score = 1.0;
			return score;
		}
		else
		{
			score = 0.7;
			return score;
		}
	}
	return score;
}



//================================================================================
// Get the sequence and conformation dependant local triplet-residues 
// 	energy;
// ATTENTION: at the N-terminal and C-terminal of target structure,
// 	there are two 2-residues local fragments, have to deal them
//	additionly. 
//================================================================================
double Etrpbin(int qpos, int tpos, vector<int>& resiindex, vector<vector<int > >& angleindex)
{
	int a1=-1, a2=-1, a3=-1, a4=-1;
	double score=0.0;
	
	//deal with the N-terminal and C-terminal of queue sequence and target structure
	if(qpos==0 || qpos==resiindex.size()-1 || tpos==0 || tpos==angleindex.size()-1)
	{
		score = 0.0;
		return score;
	}
	else	//general calculation
	{
		if(angleindex[tpos][0]==-1)	//if there are atom deletions in residue, set the score as ZERO
		{
			score=0.0;
			return score;
		}
		
		a1=resiindex[qpos-1];
		a2=resiindex[qpos];
		a3=resiindex[qpos+1];
		
		a4=6*6*6*angleindex[tpos][0]+6*6*angleindex[tpos][1]+6*angleindex[tpos][2]+angleindex[tpos][3];
		
		score = TRPTABLE[a1][a2][a3][a4];	//TRPTABLE[][][][] is a global vector stored all Etrp parameters
		
		if(score>1.0-1.0e-10)	//ATTENTATION: if set the 1.0 bins in TRPTABLE to reference state (value=0.0)
		{
			score=0.0;
		}
	}
	
	
	return score;
}



//================================================================================
// Calculate the similarity between two positions with their
//	sequence and structure profile-information;
// prof vector is the position 1's sequence profile-distribution
//	among 20 kinds of amino acides;
// strfreq vector if the frequence-distribution based local
//	9-resudues' structure information;
//================================================================================
double Eseqstr(vector<double>& prof, vector<double>& strfreq)
{
	if(prof.size()!=20)
	{
		cout << "Eseqseq()--Wrong: the size of input vector<double>& prof != 20" << endl;
		exit(-1);
	}

	if(strfreq.size()!=20)
	{
		cout << "Eseqseq()--Wrong: the size of input vector<double>& strfreq != 20" << endl;
		exit(-2);
	}
	
	int i=0;	
	double score=0.0;
	
	for(i=0; i<20; i++)
	{
		score += prof[i]*strfreq[i];
	}
	
	
	return score;
	
}




//================================================================================
// Read the standard aligned file in ProSup dataset
// 
//================================================================================
int ReadProSupFile(vector<PosPair>& stdalign, const char* prosupfile)
{
	int i=0, j=0;
	
	ifstream fin(prosupfile);
	if( !fin.is_open() )
	{
		cout << "\nReadProSupFile()--Error: can not open " << prosupfile << endl;
		exit(-1);
	}
	
	PosPair pospair;
	char buf[100];
	string record[5];
	bool flag=false;
	while( !fin.eof() )
	{
		if( !flag )
		{
			fin.getline(buf, 100);
			if(strncmp(buf, "@begin equivalences for 1", 25)==0)
			{
				flag=true;
				continue;
			}
		}
		else
		{
			fin >> record[0] >> record[1] >> record[2] >> record[3] >> record[4];
			
			if(record[0]=="@end")
			{
				break;
			}
			else
			{
				pospair.id1 = atoi( record[1].c_str() );
				pospair.id2 = atoi( record[3].c_str() );
				
				const char* chtemp=record[4].c_str();
				
				pospair.rname1 = chtemp[2];
				pospair.rname2 = chtemp[4];
				
				stdalign.push_back(pospair);
			}
		}
	}
	fin.close();
	
	
	return 0;
}



//================================================================================
// Read the standard aligned file in ProSup dataset
// read the first index alignment in the largest cluster!
//	 (maybe is not the alignment 1)
//================================================================================
int ReadProSupFilePro(vector<PosPair>& stdalign, const char* prosupfile)
{
	int i=0, j=0;
	
	ifstream fin(prosupfile);
	if( !fin.is_open() )
	{
		cout << "\nReadProSupFilePro()--Error: can not open " << prosupfile << endl;
		exit(-1);
	}


	char ch='1';
	char buf[100];	
	string record[5];
	while( !fin.eof() )
	{
		fin.getline(buf, 100);
		if( strncmp(buf, "Cluster: Alignment", 18)==0 )
		{
			fin >> record[0] >> record[1];
			if(record[0]=="1:")
			{
				const char* tempch=record[1].c_str();
				ch=tempch[0];	
			}
		}
		
	}
	fin.close();
	
	char posch[30]="@begin equivalences for 1";
	posch[24]=ch;
/*原始程序读下面的比对的时候会出错误
,因为前面几个分为5部分，后面的都是分为四部分，
我修改了比对，人为地加入了空格.
42:32:{  36 }:{1009 }:M:M改为
42:32:{  36 }:{ 1009 }:M:M


26:16:{  20 }:{ 993 }:G:R
27:17:{  21 }:{ 994 }:P:E
28:18:{  22 }:{ 995 }:R:K
29:19:{  23 }:{ 996 }:Y:I
30:20:{  24 }:{ 997 }:T:T
42:32:{  36 }:{1009 }:M:M
43:33:{  37 }:{1010 }:V:V
44:34:{  38 }:{1011 }:C:Y
45:35:{  39 }:{1012 }:S:E
46:36:{  40 }:{1013 }:A:G
47:37:{  41 }:{1014 }:Y:N
48:38:{  42 }:{1015 }:D:A
49:39:{  43 }:{1016 }:N:R
50:40:{  44 }:{1017 }:L:D
*/	
	PosPair pospair;
	bool flag=false;
	ifstream fin2;
	fin2.open(prosupfile);
	while( !fin2.eof() )
	{
		if( !flag )
		{
			fin2.getline(buf, 100);
			if(strncmp(buf, posch, 25)==0)
			{
				flag=true;
				continue;
			}
		}
		else
		{
			fin2 >> record[0] >> record[1] >> record[2] >> record[3] >> record[4];
			
			if(record[0]=="@end")
			{
				break;
			}
			else
			{
				pospair.id1 = atoi( record[1].c_str() );
				pospair.id2 = atoi( record[3].c_str() );
				
				const char* chtemp=record[4].c_str();
				
				pospair.rname1 = chtemp[2];
				pospair.rname2 = chtemp[4];
				
				stdalign.push_back(pospair);
			}
		}
	}
	fin2.close();
	
	
	return 0;
}




//================================================================================
// Compare the temporary aligned result to the standard aligned 
//	in ProSup;
// Output the percentage of correct aligned resudue-pairs as the
//	the measure of the input aligned result;
//================================================================================
int CompAlign_Local_details(AlignResult& alignresult, vector<PosPair>& stdalign, int initpos1, int initpos2,int max_i,int max_j,
		vector<int>& rserial1, vector<int>& rserial2)
{
	int i=0, j=0, k=0;
	int stdnum=stdalign.size();
	int alignednum=0, alignednum1=0, alignednum2=0;
	//cout << "stdnum" << stdnum <<"**" <<endl;
      cout << "alignment result.gapqseq.size = " << alignresult.gapqseq.size() << endl;
      cout << "alignment result.gaptseq.size = " << alignresult.gaptseq.size() << endl;
	const char* chstr1=alignresult.gapqseq.c_str();
	const char* chstr2=alignresult.gaptseq.c_str();
	int pos=0;
     	for(int i= 0;i< rserial2.size();i++)
      	{
        	//cout << "i = "<< i << "," <<rserial1[i] <<" ";//rserial in the checked PDB file
      	 }
        int L_gapqseq = alignresult.gapqseq.size();
        int L_gaptseq = alignresult.gaptseq.size();
        //cout << "L_gapqseq  " << L_gapqseq <<";L_gaptseq  " << L_gaptseq << endl;
    
        //caculate the gap numbers of the gapqseq and gaptseq
        int G_gapqseq = 0,G_gaptseq = 0;
        for(j=0; j<L_gapqseq; j++)
	{
			if(chstr1[j]=='-') G_gapqseq++;
			if(chstr2[j]=='-') G_gaptseq++;    
        }
        cout << "G_gapqseq  " <<  G_gapqseq <<";G_gaptseq  " << G_gaptseq << endl;

       //caculate the start rserial[] of the local aligment
       int gapqseq_start ,gaptseq_start ;
       //gapqseq_start = max_i-(L_gapqseq-G_gapqseq)+1;
       //gaptseq_start = max_j-(L_gaptseq-G_gaptseq)+1;
       gapqseq_start = max_i-(L_gapqseq-G_gapqseq);
       gaptseq_start = max_j-(L_gaptseq-G_gaptseq);
    	cout << "max_i  " << max_i<<"  ,max_j  " << max_j <<endl; 
    	cout << "gapqseq_start  " << gapqseq_start <<"  ,gaptseq_start  " << gaptseq_start <<endl; 
         int pos2 =0;
       cout << endl;
	for(i=0; i<stdnum; i++)
	{
		pos=stdalign[i].id1-initpos1+1;
            pos2=stdalign[i].id2-initpos2+1;//定位到每一个标准比对的残基号码,
		//cout<< "pos = " << pos << "*****"<<endl;
                //cout<< "pos2 = " << pos2 << "*****"<<endl;
		int i1= -1, i2= -1, alli=0;
		for(j=0; j<alignresult.gapqseq.size(); j++)
		{
			if(chstr1[j]!='-') i1++;
			if(chstr2[j]!='-') i2++;
			if( (rserial1[i1]+gapqseq_start)==pos )
			{
                                 //cout << "in the If loop i1 = " << i1 ;
                                 //cout << "rserial1[i1]:" << rserial1[i1]<<"--pos=" << pos << endl;
                                //cout << "in the If loop i2 = " << i2; 
                                //cout << "  rserial2[i2]:" << rserial1[i2]<<"--pos2=" << pos2 <<endl;
                                //cout << "chstr2[alli]= " <<chstr2[alli]<<","<<" std.rname2 "<<stdalign[i].rname2 ;
                                 //cout << "i2 = " << i2 ;
                               //cout << ",rserial2[ii2]= " << rserial2[i2]+gaptseq_start <<", std.POS2::" <<stdalign[i].id2-initpos2+1<<endl;
				if(chstr2[alli]==stdalign[i].rname2 && (rserial2[i2]+gaptseq_start)==pos2) 
				{
                              cout << "ch2[alli]= " <<chstr2[alli]<<","<<" std.r2 "<<stdalign[i].rname2 ;
                              cout << ",rs2[ii2]= " << rserial2[i2]+gaptseq_start <<", std.POS2 Norm::" <<stdalign[i].id2-initpos2+1<<", std.POS2 Org::" <<stdalign[i].id2<<endl;
					alignednum++;
					break;
				}
				
		/*		for(k=-4; k<=4; k++)
				{
					if(i+k>=0 && i+k<stdnum && chstr2[alli]==stdalign[i+k].rname2 && rserial2[i2]==(stdalign[i+k].id2-initpos2+1))
					{
						alignednum1++;
						break;
					}
				}
				
				//break;
			}
			
			if( rserial2[i2]==stdalign[i].id2-initpos2+1 )
			{
				for(k=-4; k<=4; k++)
				{
					if(i+k>=0 && i+k<stdnum && chstr1[alli]==stdalign[i+k].rname1 && rserial1[i1]==pos)
					{
						alignednum2++;
						
						break;
					}
				}
				
				//break;	
		*/
			}
			
			
			alli++;
		}
		
	}

	cout << alignednum << '\t' << stdnum << '\t' << alignednum*1.0/stdnum << endl;	
	
	//alignednum = alignednum1+alignednum2-alignednum;
	//cout << alignednum << '\t' << stdnum << '\t' << alignednum*1.0/stdnum << endl;
	
	return alignednum;
}
//================================================================================
// Compare the temporary aligned result to the standard aligned 
//	in ProSup;
// Output the percentage of correct aligned resudue-pairs as the
//	the measure of the input aligned result;
//================================================================================
int CompAlign_Local(AlignResult& alignresult, vector<PosPair>& stdalign, int initpos1, int initpos2,int max_i,int max_j,
		vector<int>& rserial1, vector<int>& rserial2)
{
	int i=0, j=0, k=0;
	int stdnum=stdalign.size();
	int alignednum=0, alignednum1=0, alignednum2=0;
	//cout << "stdnum" << stdnum <<"**" <<endl;
        //cout << "alignment result.gapqseq.size = " << alignresult.gapqseq.size() << endl;
        //cout << "alignment result.gaptseq.size = " << alignresult.gaptseq.size() << endl;
	const char* chstr1=alignresult.gapqseq.c_str();
	const char* chstr2=alignresult.gaptseq.c_str();
	int pos=0;
     // for(int i= 0;i< rserial2.size();i++)
       // {
        //cout << "i = "<< i << "," <<rserial1[i] <<" ";//rserial in the checked PDB file
        //}
        int L_gapqseq = alignresult.gapqseq.size();
        int L_gaptseq = alignresult.gaptseq.size();
        //cout << "L_gapqseq  " << L_gapqseq <<";L_gaptseq  " << L_gaptseq << endl;
    
        //caculate the gap numbers of the gapqseq and gaptseq
        int G_gapqseq = 0,G_gaptseq = 0;
        for(j=0; j<L_gapqseq; j++)
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
       //cout << "gapqseq_start  " << gapqseq_start <<"  ,gaptseq_start  " << gaptseq_start <<endl; need to change
         int pos2 =0;
       cout << endl;
	for(i=0; i<stdnum; i++)
	{
		pos=stdalign[i].id1-initpos1+1;
                pos2=stdalign[i].id2-initpos2+1;//定位到每一个标准比对的残基号码,
		//cout<< "pos = " << pos << "*****"<<endl;
                //cout<< "pos2 = " << pos2 << "*****"<<endl;
		int i1= -1, i2= -1, alli=0;
		for(j=0; j<alignresult.gapqseq.size(); j++)
		{
			if(chstr1[j]!='-') i1++;
			if(chstr2[j]!='-') i2++;
			if( (rserial1[i1]+gapqseq_start)==pos )
			{
                                 //cout << "in the If loop i1 = " << i1 ;
                                 //cout << "rserial1[i1]:" << rserial1[i1]<<"--pos=" << pos << endl;
                                //cout << "in the If loop i2 = " << i2; 
                                //cout << "  rserial2[i2]:" << rserial1[i2]<<"--pos2=" << pos2 <<endl;
                                //cout << "chstr2[alli]= " <<chstr2[alli]<<","<<" std.rname2 "<<stdalign[i].rname2 ;
                                 //cout << "i2 = " << i2 ;
                               //cout << ",rserial2[ii2]= " << rserial2[i2]+gaptseq_start <<", std.POS2::" <<stdalign[i].id2-initpos2+1<<endl;
				if(chstr2[alli]==stdalign[i].rname2 && (rserial2[i2]+gaptseq_start)==pos2) alignednum++;
				break;
				
		/*		for(k=-4; k<=4; k++)
				{
					if(i+k>=0 && i+k<stdnum && chstr2[alli]==stdalign[i+k].rname2 && rserial2[i2]==(stdalign[i+k].id2-initpos2+1))
					{
						alignednum1++;
						break;
					}
				}
				
				//break;
			}
			
			if( rserial2[i2]==stdalign[i].id2-initpos2+1 )
			{
				for(k=-4; k<=4; k++)
				{
					if(i+k>=0 && i+k<stdnum && chstr1[alli]==stdalign[i+k].rname1 && rserial1[i1]==pos)
					{
						alignednum2++;
						
						break;
					}
				}
				
				//break;	
		*/
			}
			
			
			alli++;
		}
		
	}

	cout << alignednum << '\t' << stdnum << '\t' << alignednum*1.0/stdnum << endl;	
	
	//alignednum = alignednum1+alignednum2-alignednum;
	//cout << alignednum << '\t' << stdnum << '\t' << alignednum*1.0/stdnum << endl;
	
	return alignednum;
}

//================================================================================
// Compare the temporary aligned result to the standard aligned 
//	in ProSup;
// Output the percentage of correct aligned resudue-pairs as the
//	the measure of the input aligned result;
//================================================================================
int CompAlign_Global(AlignResult& alignresult, vector<PosPair>& stdalign, int initpos1, int initpos2,
		vector<int>& rserial1, vector<int>& rserial2)
{
	int i=0, j=0, k=0;
	int stdnum=stdalign.size();
	int alignednum=0;
	
	const char* chstr1=alignresult.gapqseq.c_str();
	const char* chstr2=alignresult.gaptseq.c_str();
	int pos=0;
	for(i=0; i<stdnum; i++)
	{
		pos=stdalign[i].id1-initpos1+1;
		
		int i1=-1, i2=-1, alli=0;
		for(j=0; j<alignresult.gapqseq.size(); j++)
		{
			if(chstr1[j]!='-') i1++;
			if(chstr2[j]!='-') i2++;

			if( rserial1[i1]==pos )
			{
				if(chstr2[alli]==stdalign[i].rname2 && rserial2[i2]==(stdalign[i].id2-initpos2+1)) alignednum++;
				
				break;
			}
			
			alli++;
		}
		
	}

	//cout << alignednum << '\t' << stdnum << '\t' << alignednum*1.0/stdnum << endl;	
	
	return alignednum;
}




