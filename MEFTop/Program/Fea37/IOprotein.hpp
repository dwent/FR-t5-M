				/*********************************************************
				Copyright	    	 : Edited by Wuaiping, 2006
				Program name	 	 : PTADEE
				Program version	 : 0.2
				Begin time		 : Dec. 28 2005
				Refined time 	 : Jun. 14 2007
				Email    		 : wuaiping@moon.ibp.ac.cn

				Edited by WangXuan on 2011.5.5 see readme for detail
				mainly related to VLA of C language and the forgotten <string>
				*********************************************************/

#ifndef IOPROTEIN_H
#define IOPROTEIN_H

#include <iostream>
#include <fstream>
#include <cstring>
#include <string>

using namespace std;

class PDB
{
public:
	
	int AtomNum(const char* PDBfile);
	int ResiNum(const char* PDBfile);
	
	void ReadName(const char* PDBfile, string name[], string datalist);
	void ReadCoord(const char* PDBfile, double coordinate[], string datalist);
	void ReadSerial(const char* PDBfile, int serial[], string datalist);
	

private:
	int atomnum;
	int resinum;

};
int OutProtein(const char* outname, int atomnum, char chainname, int resiserial[], double x[], double y[], double z[], string aname[], string rname[]);


#endif
