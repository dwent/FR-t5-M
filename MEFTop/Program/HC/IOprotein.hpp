				/*********************************************************
				Copyright	    	 : Edited by Wuaiping, 2006
				Program name	 	 : PTADEE
				Program version	 : 0.2
				Begin time		 : Dec. 28 2005
				Refined time 	 : Jun. 14 2007
				Email    		 : wuaiping@moon.ibp.ac.cn
				*********************************************************/

#ifndef IOPROTEIN_H
#define IOPROTEIN_H

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>


using namespace std;

class PDB
{
public:
	
	int AtomNum(char* PDBfile);
	int ResiNum(char* PDBfile);
	
	void ReadName(char* PDBfile, string name[], string datalist);
	void ReadCoord(char* PDBfile, double coordinate[], string datalist);
	void ReadSerial(char* PDBfile, int serial[], string datalist);
	

private:
	int atomnum;
	int resinum;

};

int OutProtein(char* outname, int atomnum, char chainname, int resiserial[], double x[], double y[], double z[], string aname[], string rname[]);


#endif
