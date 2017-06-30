				/*********************************************************
				Copyright	    	 : Edited by Wuaiping, 2006
				Program name	 	 : PTADEE
				Program version	 : 0.2
				Begin time		 : Dec. 28 2005
				Refined time 	 : Jun. 14 2007
				Email    		 : wuaiping@moon.ibp.ac.cn
				*********************************************************/

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
	string bufline;
	
	ifstream infile;
	infile.open(PDBfile);
	if( !infile.is_open() )
	{
		cout << "\nCan not open " << PDBfile << endl;
		exit(-1);
	}
	
	while( getline(infile,bufline) )
	{
		if( bufline.compare(0,4,"ATOM")==0 )
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
	string bufline;
	int serial1=-1000, serial2;
	int i=22;
	
	ifstream infile;
	infile.open(PDBfile);
	if( !infile.is_open() )
	{
		cout << "\nCan not open " << PDBfile << endl;
		exit(-1);
	}
	
	while( getline(infile,bufline) )
	{
		i=22;		
		if( bufline.compare(0,4,"ATOM")==0 )
		{
			while(bufline[i]==' ')
			{ i++; }
			serial2=atoi(&bufline[i]);
			
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
//Get the atom-names or residue-names of the input pdb file
******************************************************************************/
void PDB::ReadName(const char* PDBfile, string name[], string datalist)
{	
	string bufline;
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
	
	while( getline(infile,bufline) )
	{		
		if( bufline.compare(0,4,"ATOM")==0 )
		{
			//for(int k=i; k<i+4; k++)
			//{
				//if(buf[k]!=' ')
					//name[j]= name[j]+buf[k];
			//}
			//j++;
			if(i==13)
			{
				if(bufline[12]!=' ')
				{
					for(int k=i-1; k<i+3; k++)
					{
						if(bufline[k]!=' ')
							name[j]= name[j]+bufline[k];
					}
					j++;
				}
				else
				{
					for(int k=i; k<i+3; k++)
					{
						if(bufline[k]!=' ')
							name[j]= name[j]+bufline[k];
					}
					j++;
				}
			}
			else
			{
				for(int k=i; k<i+3; k++)
				{
					if(bufline[k]!=' ')
						name[j]= name[j]+bufline[k];
				}
				j++;
			}
		}
	}
	
	infile.close();

}



/******************************************************************************
//Get the XYZ-coordinates of the input pdb file
******************************************************************************/
void PDB::ReadCoord(const char* PDBfile, double coordinate[], string datalist)
{
	string bufline;
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
	
	while( getline(infile,bufline) )
	{		
		if( bufline.compare(0,4,"ATOM")==0 )
		{
			ti=i;
			while(bufline[ti]==' ')
				ti++;
			coordinate[j]=atof( &bufline[ti] );
			
			j++;
		}
	}
	
	infile.close();

}


/******************************************************************************
//Get the atom-serials or residue-serials of the input pdb file
******************************************************************************/
void PDB::ReadSerial(const char* PDBfile, int serial[], string datalist)
{
	string bufline;
	int i,j=0;
	int ti=0;
	if(datalist=="ATOMSERIAL")
		{ i=7; }
	else if(datalist=="RESISERIAL")
		{ i=22; }
	else
		{ cout << "Wrong: the final parameter is ATOMSERIAL or RESISERIAL." << endl; }

	
	ifstream infile;
	infile.open(PDBfile);
	if( !infile.is_open() )
	{
		cout << "\nCan not open " << PDBfile << endl;
		exit(-1);
	}
	
	while( getline(infile,bufline) )
	{
		if( bufline.compare(0,4,"ATOM")==0 )
		{
			ti=i;
			while(bufline[ti]==' ')
				ti++;
			serial[j]=atoi( &bufline[ti] );
			
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
	
	char* atomtype=new char[atomnum];
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
	delete[] atomtype;
	
	return 0;
}

