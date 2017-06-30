				/*********************************************************
				Copyright	    	 : Edited by Wuaiping, 2006
				Program name	 	 : PTADEE
				Program version	 : 0.2
				Begin time		 : Dec. 28 2005
				Refined time 	 : Jun. 14 2007
				Email    		 : wuaiping@moon.ibp.ac.cn
				*********************************************************/
				/**********************************************************
				Copyright(c)2012, IBP, CHINESE ACADEMY OF SCIENCE
				All rights reserved.
				NAME:		IOprotein.cpp
				ABSTRACT:	THREADING IOPROTEIN
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
//Define class PDB

#include "IOprotein.hpp"
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <iomanip>

using namespace std;


/******************************************************************************
//Get the number of atoms in the input pdb file
******************************************************************************/
int PDB::AtomNum(const char* PDBfile)
{
	atomnum=0;
	char buf[100];
	
	ifstream infile;
	infile.open(PDBfile);
	if( !infile.is_open() )
	{
		cout << "\nCan not open " << PDBfile << endl;
		exit(-1);
	}
	
	while( !infile.eof() )
	{
		infile.getline(buf, 99);
		if( strncmp(buf, "ATOM", 4)==0 )
		{
			atomnum++;
		}
	}
	
	infile.close();
	
	return atomnum;
}


/******************************************************************************
//Get the number of residues in the input pdb file
******************************************************************************/
int PDB::ResiNum(const char* PDBfile)
{
	resinum=0;
	char buf[100];
	int serial1=-1000, serial2;
	int i=22;
	
	ifstream infile;
	infile.open(PDBfile);
	if( !infile.is_open() )
	{
		cout << "\nCan not open " << PDBfile << endl;
		exit(-1);
	}
	
	while( !infile.eof() )
	{
		i=22;
		infile.getline(buf, 99);
		if( strncmp(buf, "ATOM", 4)==0 )
		{
			while(buf[i]==' ')
			{ i++; }
			serial2=atoi(&buf[i]);
			
			if(serial1!=serial2)
			{
				resinum++;
				serial1=serial2;
			}			
			
		}
	}
	
	infile.close();
	
	return resinum;
}
/******************************************************************************
//READ A FASTA FILE
******************************************************************************/
std::string ReadFasta(char* Fasfile)
{
	std::string ss;
	std::string buf;
	std::ifstream infile;
	infile.open(Fasfile);
	if(!infile){std::cerr<<"\nCan not open: "<<Fasfile<<'\n';exit(0);}
	while(!infile.eof())
	{
		buf.clear();
		getline(infile,buf);
		if(buf.empty())break;
		if(buf[0]=='>')continue;
		ss.append(buf);
	};infile.close();
	return ss;
}

/******************************************************************************
//Get the atom-names or residue-names of the input pdb file
******************************************************************************/
void PDB::ReadName(const char* PDBfile, string name[], string datalist)
{
	char buf[100];
	int i,j=0;
	if(datalist=="ATOMNAME")
		{ i=13; }
	else if(datalist=="RESINAME")
		{ i=17; }
	else
		{ cout << "Wrong: the final parameter is ATOMNAME or RESINAME." << endl; }
	
	ifstream infile;
	infile.open(PDBfile);
	if( !infile.is_open() )
	{
		cout << "\nCan not open " << PDBfile << endl;
		exit(-1);
	}
	
	while( !infile.eof() )
	{
		infile.getline(buf, 99);
		if( strncmp(buf, "ATOM", 4)==0 )
		{
			for(int k=i; k<i+4; k++)
			{
				if(buf[k]!=' ')
					name[j]= name[j]+buf[k];
			}
			j++;
		}
	}
	
	infile.close();

}


/******************************************************************************
//Get the XYZ-coordinates of the input pdb file
******************************************************************************/
void PDB::ReadCoord(const char* PDBfile, double coordinate[], string datalist)
{
	char buf[100];
	int i,j=0;
	int ti=0;
	if(datalist=="XCOORD")
		{ i=30; }
	else if(datalist=="YCOORD")
		{ i=38; }
	else if(datalist=="ZCOORD")
		{ i=46; }
	else
		{ cout << "Wrong: the final parameter is XCOORD, YCOORD or ZCOORD." << endl; }
	
	ifstream infile;
	infile.open(PDBfile);
	if( !infile.is_open() )
	{
		cout << "\nCan not open " << PDBfile << endl;
		exit(-1);
	}
	
	while( !infile.eof() )
	{
		infile.getline(buf, 99);
		if( strncmp(buf, "ATOM", 4)==0 )
		{
			ti=i;
			while(buf[ti]==' ')
				ti++;
			coordinate[j]=atof( &buf[ti] );
			
			j++;
		}
	}
	
	infile.close();

}
/******************************************************************************
//MEASURE THE HP TYPE FOR A AMINO ACID
******************************************************************************/
bool hca(char aa)
{
	switch(aa)
	{
		case 'V':return true;
		case 'I':return true;
		case 'L':return true;
		case 'F':return true;
		case 'Y':return true;
		case 'W':return true;
		case 'M':return true;
		default:return false;
	}
	return false;
}

/******************************************************************************
//Get the atom-serials or residue-serials of the input pdb file
******************************************************************************/
void PDB::ReadSerial(const char* PDBfile, int serial[], string datalist)
{
	char buf[100];
	int i,j=0;
	int ti=0;
	if(datalist=="ATOMSERIAL")
		{ i=7; }
	else if(datalist=="RESISERIAL")
		{ i=23; }
	else
		{ cout << "Wrong: the final parameter is ATOMSERIAL or RESISERIAL." << endl; }

	
	ifstream infile;
	infile.open(PDBfile);
	if( !infile.is_open() )
	{
		cout << "\nCan not open " << PDBfile << endl;
		exit(-1);
	}
	
	while( !infile.eof() )
	{
		infile.getline(buf, 99);
		if( strncmp(buf, "ATOM", 4)==0 )
		{
			ti=i;
			while(buf[ti]==' ')
				ti++;
			serial[j]=atoi( &buf[ti] );
			
			j++;
		}
	}
	
	infile.close();

}


/******************************************************
//Output the input data into a PDB file
//
******************************************************/
int OutProtein(const char* outname, int atomnum, char chainname, int resiserial[], double x[], double y[], double z[], string aname[], string rname[])
{
	//open a output file
	ofstream outfile(outname, ios::out | ios::trunc);
	if(!outfile)
	{
		cerr << "\nError: unable to open output file!\n";
		return -1;
	}
	
	char atomtype[atomnum];
	for(int i=0; i<atomnum; i++)
	{
		const char* tempatomname=aname[i].c_str();
		atomtype[i]=tempatomname[0];
	}
	
	//output all the data into file by the pdb file format
	for(int i=0; i<atomnum; i++)
	{
		outfile.seekp(81*i+0,ios::beg);
		outfile << "ATOM";

		outfile.seekp(81*i+4,ios::beg);
		char index[10];
		sprintf(index, "%d", i+1);
		outfile << setw(7) << index;
		outfile.setf(ios::right);
		
		outfile.seekp(81*i+11,ios::beg);
		outfile << "  ";
		
		outfile.seekp(81*i+13,ios::beg);
		if(aname[i].length()==1)
		{
			outfile << aname[i]<<"   ";
		}
		else if(aname[i].length()==2)
		{
			outfile << aname[i]<<"  ";
		}
		else if(aname[i].length()==3)
		{
			outfile << aname[i]<<" ";
		}
		else if(aname[i].length()==4)
		{
			outfile << aname[i];
		}
		
		outfile.seekp(81*i+17,ios::beg);
		const char* temprname=rname[i].c_str();
		if(strncmp(temprname,"HSD",3)==0)
		{	
			outfile << setw(3)<< "HIS";
		}
		else
			outfile << setw(3)<< rname[i];
		outfile.setf(ios::right);

		outfile.seekp(81*i+20,ios::beg);
		outfile << setw(2) << chainname;
		outfile.setf(ios::right);

		outfile.seekp(81*i+22,ios::beg);
		char resindex[10];
		sprintf(resindex, "%d", resiserial[i]);
		outfile << setw(4) << resindex;
		outfile.setf(ios::right);

		char xcoor[10];
		sprintf(xcoor, "%.3f", x[i]);
		outfile.seekp(81*i+26,ios::beg);
		outfile << setw(12) << xcoor;
		outfile.setf(ios::right);
		
		char ycoor[10];
		sprintf(ycoor, "%.3f", y[i]);
		outfile.seekp(81*i+38,ios::beg);
		outfile << setw(8)<< ycoor;
		outfile.setf(ios::right);

		char zcoor[10];
		sprintf(zcoor, "%.3f", z[i]);
		outfile.seekp(81*i+46,ios::beg);
		outfile << setw(8) << zcoor;
		outfile.setf(ios::right);

		outfile.seekp(81*i+54,ios::beg);
		outfile << setw(6)<<"1.00";
		outfile.setf(ios::right);

		outfile.seekp(81*i+60,ios::beg);
		outfile << setw(18)<<atomtype[i];
		outfile.setf(ios::right);

		outfile.seekp(81*i+78,ios::beg);
		outfile << setw(3) <<'\n';
		outfile.setf(ios::right);
	}

	outfile.close();
	
	return 0;
}

/**/
void AAname1to3(char name1, string& name3)
{

	switch( name1 )  //switch语句   switch后面括号只能接受整型，包括char，short，int，long等，不能是浮点数。double不行
	{
		case 'A': name3="ALA"; break;
		case 'R': name3="ARG"; break;
		case 'N': name3="ASN"; break;
		case 'D': name3="ASP"; break;
		case 'C': name3="CYS"; break;
		case 'Q': name3="GLN"; break;
		case 'E': name3="GLU"; break;
		case 'G': name3="GLY"; break;
		case 'H': name3="HIS"; break;
		case 'I': name3="ILE"; break;
		case 'L': name3="LEU"; break;
		case 'K': name3="LYS"; break;
		case 'M': name3="MET"; break;
		case 'F': name3="PHE"; break;
		case 'P': name3="PRO"; break;
		case 'S': name3="SER"; break;
		case 'T': name3="THR"; break;
		case 'W': name3="TRP"; break;
		case 'Y': name3="TYR"; break;
		case 'V': name3="VAL"; break;
		default: break;
	}

}

bool AAname3to1( string name3,char& name1)
{

		if( name3 == "ALA"){ name1='A';}
		else if( name3 == "ARG"){ name1='R';}
		else if( name3 == "ASN"){ name1='N';}
		else if( name3 == "ASP"){ name1='D';}
		else if( name3 == "CYS"){ name1='C';}
		else if( name3 == "GLN"){ name1='Q';}
		else if( name3 == "GLU"){ name1='E';}
		else if( name3 == "GLY"){ name1='G';}
		else if( name3 == "HIS"){ name1='H';}
		else if( name3 == "ILE"){ name1='I';}
		else if( name3 == "LEU"){ name1='L';}
		else if( name3 == "LYS"){ name1='K';}
		else if( name3 == "MET"){ name1='M';}
		else if( name3 == "PHE"){ name1='F';}
		else if( name3 == "PRO"){ name1='P';}
		else if( name3 == "SER"){ name1='S';}
		else if( name3 == "THR"){ name1='T';}
		else if( name3 == "TRP"){ name1='W';}
		else if( name3 == "TYR"){ name1='Y';}
		else if( name3 == "VAL"){ name1='V';}
            else {cout << "In IOprotein.cpp() AA name3 error! " << endl;return false;}
		return true;

}

