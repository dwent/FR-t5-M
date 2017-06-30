//===========================================================================================================
//Edited by Wu Aiping (Email: wuaiping@moon.ibp.ac.cn)
//Energy function based residue-five-beads model (N,CA,C,O, and centroid of side-chain)
//Energy Function: V = Econ + a*Etrp + b*Ehp
//	1)Econ : atom-pairs contact MIU-potential;
//	2)Etrp : local 3-residues' sequence-dependent conformation MIU-potential;
//	3)Ehp  : hydrogen-bonding potential, realized by MIU-poteintial
//	4)a, b and c in Ehp: three weigth factors to balance different energy sets; 
//===========================================================================================================

#include "FiveBeadEnergy.hpp"


using namespace std;



//=============================================================================================================
//Read all parameter tables used in the program, each table means as following:
//1.vdwradtable, the minimum distance of each atom-type-pairwise gotten from PDB database;
//2.contable, parameter table of Econ (atom contact set) energy set;
//3.trptable, parameter table of Etrp (local 3-residues conformation energy) energy set;
//4.hbtable, parameter table of Ehb (hydrogen bonding energy) energy set;
//5.centable, parameter table of PHI-PSI-dependent centroid atom positions;
//=============================================================================================================
void ReadAllPar(DV2& vdwradtable, const char* vdwfile, DV2& contable, const char* confile, DV4& trptable, const char* trpfile, DV6& hbtable, const char* hbfile, DV4& centable, const char* centroidINFOfile)
{
	int i=0, j=0, k=0, l=0, m=0, n=0, p=0;

	//Read INTERRESI and vdwradtable[][]
	int ATYPE=99;
	int INTERRESI=-1;
	vdwradtable.resize(ATYPE);
	for(i=0; i<ATYPE; i++)
	{
		vdwradtable[i].resize(ATYPE);
	}
		
	ReadVdwRadius(INTERRESI, vdwradtable, vdwfile);
	
	
	//Read Econ parameter table
	contable.resize(ATYPE);
	for(i=0; i<ATYPE; i++)
	{
		contable[i].resize(ATYPE);
	}
	
	ReadConPar(contable, confile);
	
	
	//Read Etrp parameter table
	int RNUM=20;
	int LEVEL4=1296;
	trptable.resize(RNUM);
	for(i=0; i<RNUM; i++)
	{
		trptable[i].resize(RNUM);
		for(j=0; j<RNUM; j++)
		{
			trptable[i][j].resize(RNUM);
			for(k=0; k<RNUM; k++)
			{
				trptable[i][j][k].resize(LEVEL4);
			}
		}
	}
	
	ReadTrpPar(trptable, trpfile);
	
	
	//Read Ehp parameter table
	int BIN=9;
	hbtable.resize(BIN);
	for(i=0; i<BIN; i++)
	{
		hbtable[i].resize(BIN);
		for(j=0; j<BIN; j++)
		{
			hbtable[i][j].resize(BIN);
			for(k=0; k<BIN; k++)
			{
				hbtable[i][j][k].resize(BIN);
				for(l=0; l<BIN; l++)
				{
					hbtable[i][j][k][l].resize(BIN);
					for(m=0; m<BIN; m++)
					{
						hbtable[i][j][k][l][m].resize(BIN);
					}
				}
			}
		}
	}
	
	ReadHbPar(hbtable, hbfile);
	
	
	
	//Read backbone-dependent centroid of side-chain internal parameter table (19*36*36*6)
	centable.resize(19);
	for(i=0; i<19; i++)
	{
		centable[i].resize(36);
		for(j=0; j<36; j++)
		{
			centable[i][j].resize(36);
			for(k=0; k<36; k++)
			{
				centable[i][j][k].resize(6);
			}
		}
	}
	
	ReadCentroid(centable, centroidINFOfile);
	
}

//=====================================================================================================
//READ THE PSIPRED FILE INTO A STRING;

//=====================================================================================================
std::string ReadHoriz(char* ssfile)
{
	std::string ss,line;
	std::ifstream infile(ssfile);
	if(!infile){std::cerr<<"\nCan not open: "<<ssfile<<'\n';exit(0);}
	while(!infile.eof())
	{
		line.clear();
		getline(infile,line);
		if(line.empty())continue;
		if(line.substr(0,4)=="Pred")ss.append(line.substr(6,line.length()-6));
	};infile.close();
	return ss;
}


//=============================================================================================================
//Read INTERRESI value and VDWRADIUS[99][99] matrix that calculated from
//	"../bin/getvdwmatrix99_99 ......"
//INTERRESI means unless how many residues between two atoms calculated as atom-atom pair
//VDWRADIUS[99][99] is the minimum distance matrix statified with INTERRESI condition
//=============================================================================================================
void ReadVdwRadius(int& INTERRESI, DV2& VDWRADIUS, const char* vdwfile)
{
	int i=0, j=0;
	int NUM=VDWRADIUS.size();
	
	string s1,s2;
	string record[NUM];
	
	ifstream fin(vdwfile);
	if( !fin.is_open() )
	{
		cout << "\nReadVdwRadius()--Error: can not open " << vdwfile << endl;
		exit(-1);
	}
	
	bool flag=false;
	int index=0;
	while( !fin.eof() )
	{
		if( !flag )
		{
			fin >> s1 >> s2;
			
			if(s1=="INTERRESI")
			{
				INTERRESI = atoi(s2.c_str());
				flag=true;
				continue;
			}
		}	
		else if(flag && index<NUM)
		{
			for(i=0; i<NUM; i++)
			{
				fin >> record[i];
				
				VDWRADIUS[index][i] = atof(record[i].c_str());
			}
			
			index++;
		}
		else
		{
			break;
		}
		
	}
	
	fin.close();
}




//=====================================================================================================
//Read parameter file of 23*23 kinds of atom-pairs' contact table 
//23 atom types included N, CA, C, O and 19 kinds of centroid of residue side-chain
//	except GLY;
//=====================================================================================================
int ReadConPar(DV2& contable, const char* parfile)
{
	int i=0;
	int NUM=contable.size();
	
	ifstream fin(parfile);
	if( !fin.is_open() )
	{
		cout << "\nReadConPar()--Error: can not open " << parfile << endl;
		exit(-1);
	}
	
	string data[NUM];
	for(int index=0; index<NUM; index++)
	{
		for(i=0; i<NUM; i++)
		{
			fin >> data[i];
				
			contable[index][i] = atof(data[i].c_str());
		}
	}
	fin.close();
	
	return 0;
}




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




//=====================================================================================================
//Read parameter file of 9*9*9*9*9*9 bins of Ehb table;
//For two 3-residues fragment, each fragment have 6 vectors (namely, each residue
//	have 2 vectors --b(i) and P(i)), and 6 angles can be substracted from these
//	6-6 vector-pairs;
//For each angle, cut it into 9 bins by 20-degrees, so there are 9*9*9*9*9*9 bins
//=====================================================================================================
int ReadHbPar(DV6& hbtable, const char* parfile)
{
	int a=0, b=0, c=0, d=0, e=0;
	int i=0, j=0;
	
	ifstream fin(parfile);
	if( !fin.is_open() )
	{
		cout << "\nReadHbPar()--Error: can not open " << parfile << endl;
		exit(-1);
	}
	
	
	string data[9];
	int num5=9*9*9*9*9;
	int num4=9*9*9*9;
	int num3=9*9*9;
	int num2=9*9;
	int num1=9;
	for(i=0; i<num5; i++)
	{
		a=i/num4;
		b=(i-a*num4)/num3;
		c=(i-a*num4-b*num3)/num2;
		d=(i-a*num4-b*num3-c*num2)/num1;
		e=i-a*num4-b*num3-c*num2-d*num1;
		
		
		for(j=0; j<9; j++)
		{
			fin >> data[j];
			
			hbtable[a][b][c][d][e][j] = atof(data[j].c_str());
		}
	}
	fin.close();
	
	return 0;
}




//=====================================================================================================
//Read parameter file of 19*36*36*6 bins of centroid table;
//
//=====================================================================================================
int ReadCentroid(DV4& centable, const char* cenfile)
{
	int a=0, b=0, c=0;
	int i=0, j=0;
	
	ifstream fin(cenfile);
	if( !fin.is_open() )
	{
		cout << "\nReadCentroid()--Error: can not open " << cenfile << endl;
		exit(-1);
	}
	
	
	int num1=19*36*36;
	int num2=36*36;
	int num3=36;
	string data[10];
	for(i=0; i<num1; i++)
	{	
		a=i/num2;
		b=(i-a*num2)/num3;
		c=i-a*num2-b*num3;
		
		for(j=0; j<10; j++)
		{
			fin >> data[j];
		}
		
		for(j=0; j<6; j++)
		{
			centable[a][b][c][j] = atof(data[j+4].c_str());
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




//=====================================================================================================
//Calculate the atom-pairs' contact energy of residue's five-Beads model (99*99 types)
//
//=====================================================================================================
double EconOfFiveBead99_99(ChainCoord& chain, DV2& VDWRADIUS, DV2& contable)
{
	int i=0, j=0, k=0, l=0, m=0, n=0;
	double outE=0.0;
	
	
	//=========================================//
	string RNAME[20]={"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "HIS", "ILE",
			     "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "GLY"};
	
	
	double LOWLIMIT=1.0, HIGHLIMIT=1.9;
	double LowD2[99][99], HighD2[99][99];
	for(i=0; i<99; i++)
	{
		for(j=0; j<99; j++)
		{
			LowD2[i][j]=LOWLIMIT*VDWRADIUS[i][j];
			LowD2[i][j]=LowD2[i][j]*LowD2[i][j];
			
			HighD2[i][j]=HIGHLIMIT*VDWRADIUS[i][j];
			HighD2[i][j]=HighD2[i][j]*HighD2[i][j];
		}
	}
	//=========================================//
	
	vector<int> RI;
	for(i=0; i<chain.resinum; i++)
	{
		for(j=0; j<20; j++)
		{
			if(chain.resiname[i]==RNAME[j])
			{
				RI.push_back(5*j);
				break;
			}
		}
	}
	
	
	double x1, y1, z1, x2, y2, z2;
	double d=0.0;
	int i1=-1, i2=-1;
	int len=chain.xyz.size();
	
	//get the distance distribution, when distance*distance of two atoms large than thresh, not calculated their interaction
	double thresh=16*16;
	vector<vector<double > > DISTANCE;
	vector<double> tempDIS;
	for(j=0; j<len-1; j++)
	{
		x1=chain.xyz[j][0][0];
		y1=chain.xyz[j][0][1];
		z1=chain.xyz[j][0][2];
			
		tempDIS.clear();
		for(l=j+1; l<len; l++)
		{
			x2=chain.xyz[l][0][0];
			y2=chain.xyz[l][0][1];
			z2=chain.xyz[l][0][2];
						
			d = (x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1);
			
			tempDIS.push_back(d);
		}
			
		DISTANCE.push_back(tempDIS);
	}
	
	
	//calculate the Econ score
	for(j=0; j<len-1; j++)
	{
		for(l=j+1; l<len; l++)
		{	
			if(DISTANCE[j][l-j-1]<=thresh)
			{
				for(k=0; k<chain.xyz[j].size(); k++)
				{
					x1=chain.xyz[j][k][0];
					y1=chain.xyz[j][k][1];
					z1=chain.xyz[j][k][2];
				
					i1=RI[j]+k;
				
					for(m=0; m<chain.xyz[l].size(); m++)
					{
						x2=chain.xyz[l][m][0];
						y2=chain.xyz[l][m][1];
						z2=chain.xyz[l][m][2];
						
						d = (x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1);
							
						i2=RI[l]+m;
	
						if(d<HighD2[i1][i2] && d>=LowD2[i1][i2])
						{
							outE += contable[i1][i2];
						}
					}
				}
			}
		}
	}
	
	RI.clear();
	
	return outE;
}





//=============================================================================================
//Read protein structure pdb file
//Substract the infomation of five-beads (N, CA, C, O and all heavy atoms 
//	of side-chain) of each residue
//=============================================================================================
int ReadPdb(ChainCoord& chain, const char* pdbfile)
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

	
	vector<string> tempaname;
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
				chain.atomname.push_back(tempaname);
				chain.xyz.push_back(resi);
				resi.clear();
				tempaname.clear();
			}
			
			tempaname.push_back(aname[i]);
			
			atom.clear();
			atom.push_back(x[i]);
			atom.push_back(y[i]);
			atom.push_back(z[i]);
			
			resi.push_back(atom);
						
			chain.resiname.push_back(rname[i]);
			chain.resiserial.push_back(rserial[i]);
		}
		else
		{
			tempaname.push_back(aname[i]);
			
			atom.clear();
			atom.push_back(x[i]);
			atom.push_back(y[i]);
			atom.push_back(z[i]);
			
			resi.push_back(atom);

		}

		
		if(i==atomnum-1)
		{
			chain.atomname.push_back(tempaname);
			chain.xyz.push_back(resi);
		}
		
	}


	return 0;	
	
}



//===================================================================================
//Translate atom name of a residue into an index in the 84 atom types
//	to construct 84*84 atom-pairs matrix
//===================================================================================
int Aname2Index(string rname, string aname)
{
	int index=-1;
	
	if(aname=="N") index=79;
	else if(aname=="CA" && rname!="GLY") index=80;
	else if(aname=="C") index=81;
	else if(aname=="O") index=82;
	else if(aname=="OXT" || aname=="OCT") index=83;
	else if(rname=="ALA")
	{
		if(aname=="CB") index=0;
	}
	else if(rname=="ARG")
	{
		if(aname=="CB") index=1;
		else if(aname=="CG") index=2;
		else if(aname=="CD") index=3;
		else if(aname=="NE") index=4;	
		else if(aname=="CZ") index=5;	
		else if(aname=="NH1" || aname=="NH2") index=6;	
		
	}
	else if(rname=="ASN")
	{
		if(aname=="CB") index=7;
		else if(aname=="CG") index=8;
		else if(aname=="OD1") index=9;
		else if(aname=="ND2") index=10;
	
	}
	else if(rname=="ASP")
	{
		if(aname=="CB") index=11;
		else if(aname=="CG") index=12;
		else if(aname=="OD1" || aname=="OD2") index=13;
	
	}
	else if(rname=="CYS")
	{
		if(aname=="CB") index=14;
		else if(aname=="SG") index=15;
	
	}
	else if(rname=="GLN")
	{
		if(aname=="CB") index=16;
		else if(aname=="CG") index=17;
		else if(aname=="CD") index=18;
		else if(aname=="OE1") index=19;
		else if(aname=="NE2") index=20;
	
	}
	else if(rname=="GLU")
	{
		if(aname=="CB") index=21;
		else if(aname=="CG") index=22;
		else if(aname=="CD") index=23;
		else if(aname=="OE1" || aname=="OE2") index=24;
	
	}
	else if(rname=="HIS")
	{
		if(aname=="CB") index=25;
		else if(aname=="CG") index=26;
		else if(aname=="ND1") index=27;
		else if(aname=="CD2") index=28;
		else if(aname=="CE1") index=29;
		else if(aname=="NE2") index=30;
	
	}
	else if(rname=="ILE")
	{
		if(aname=="CB") index=31;
		else if(aname=="CG1") index=32;
		else if(aname=="CG2") index=33;
		else if(aname=="CD1") index=34;
	
	}
	else if(rname=="LEU")
	{
		if(aname=="CB") index=35;
		else if(aname=="CG") index=36;
		else if(aname=="CD1" || aname=="CD2") index=37;
	
	}
	else if(rname=="LYS")
	{
		if(aname=="CB") index=38;
		else if(aname=="CG") index=39;
		else if(aname=="CD") index=40;
		else if(aname=="CE") index=41;
		else if(aname=="NZ") index=42;
	
	}
	else if(rname=="MET")
	{
		if(aname=="CB") index=43;
		else if(aname=="CG") index=44;
		else if(aname=="SD") index=45;
		else if(aname=="CE") index=46;
	
	}
	else if(rname=="PHE")
	{
		if(aname=="CB") index=47;
		else if(aname=="CG") index=48;
		else if(aname=="CD1" || aname=="CD2") index=49;
		else if(aname=="CE1" || aname=="CE2") index=50;
		else if(aname=="CZ") index=51;
		
	}
	else if(rname=="PRO")
	{
		if(aname=="CB") index=52;
		else if(aname=="CG") index=53;
		else if(aname=="CD") index=54;
	
	}
	else if(rname=="SER")
	{
		if(aname=="CB") index=55;
		else if(aname=="OG") index=56;
	
	}
	else if(rname=="THR")
	{
		if(aname=="CB") index=57;
		else if(aname=="OG1") index=58;
		else if(aname=="CG2") index=59;
	
	}	
	else if(rname=="TRP")
	{
		if(aname=="CB") index=60;
		else if(aname=="CG") index=61;
		else if(aname=="CD1") index=62;
		else if(aname=="CD2") index=63;
		else if(aname=="NE1") index=64;
		else if(aname=="CE2") index=65;
		else if(aname=="CE3") index=66;
		else if(aname=="CZ2") index=67;
		else if(aname=="CZ3") index=68;
		else if(aname=="CH2") index=69;
	
	}
	else if(rname=="TYR")
	{
		if(aname=="CB") index=70;
		else if(aname=="CG") index=71;
		else if(aname=="CD1" || aname=="CD2") index=72;
		else if(aname=="CE2" || aname=="CE2") index=73;
		else if(aname=="CZ") index=74;
		else if(aname=="OH") index=75;
	
	}
	else if(rname=="VAL")
	{
		if(aname=="CB") index=76;
		else if(aname=="CG1") index=77;
		else if(aname=="CG2") index=78;
	
	}
	else if(rname=="GLY")
	{
		if(aname=="CA") index=79;
	}

	
	
	return index;	
}





//=====================================================================================================
//Calculate the atom-pairs' contact energy of residue's all-atoms model (84*84 types)
//
//=====================================================================================================
double EconOfAllAtom84_84(ChainCoord& chain, int INTERRESI, DV2& VDWRADIUS, DV2& contable)
{
	int i=0, j=0, k=0, l=0, m=0, n=0;
	double outE=0.0;
	
	
	//=========================================//
	double LOWLIMIT=1.0, HIGHLIMIT=1.9;
	double LowD[84][84], HighD[84][84];
	for(i=0; i<84; i++)
	{
		for(j=0; j<84; j++)
		{
			LowD[i][j]=LOWLIMIT*VDWRADIUS[i][j];
			HighD[i][j]=HIGHLIMIT*VDWRADIUS[i][j];
		}
	}
	//=========================================//
	
	
	
	double x1, y1, z1, x2, y2, z2;
	double d=0.0;

	for(j=0; j<chain.xyz.size()-1; j++)
	{
		int rs1=chain.resiserial[j];
		
		for(k=0; k<chain.xyz[j].size(); k++)
		{
			x1=chain.xyz[j][k][0];
			y1=chain.xyz[j][k][1];
			z1=chain.xyz[j][k][2];
				
			int i1=Aname2Index(chain.resiname[j], chain.atomname[j][k]);			
				
			if(i1<0 || i1>83)
			{
				//cout << "\nError: data wrong! i1= " << i1 << endl;
				continue;
			}				
									
				
			for(l=j+1; l<chain.xyz.size(); l++)
			{
				int rs2=chain.resiserial[l];
				if(rs2-rs1<=INTERRESI)
				{
					continue;
				}
					
				for(m=0; m<chain.xyz[l].size(); m++)
				{
					x2=chain.xyz[l][m][0];
					y2=chain.xyz[l][m][1];
					z2=chain.xyz[l][m][2];
						
					d = (x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1);
					d = sqrt(d);
						
					int i2=Aname2Index(chain.resiname[l], chain.atomname[l][m]);					
					
					if(i2<0 || i2>83)
					{
						//cout << "\nError: data wrong! i2= " << i2 << endl;
						continue;
					}
					

					if(d<HighD[i1][i2] && d>=LowD[i1][i2])
					{
						outE += contable[i1][i2];
						//cout << i1 << ' ' << i2 << '\t' << contable[i1][i2] << endl;
					}
				}
			}
		}
	}
	
	
	return outE;
}




//=========================================================================================
//Substract 2 dihedrals and 2 angles from a local 3-residues' fragment
//For fragment A(i)A(i+1)A(i+2), 2 dihedrals are Phi(i+1) and Psi(i+1),
//	2 angles are angle between P(i) and P(i+2), angle between b(i) and
//	b(i+2)
//  A SEVERE BUG IN THIS FUNCTION UN-CORRECTABLE
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
	
	
	return 0;
}




//=========================================================================================
//For an input *.3resi file, get all the 4-angles information of all the 
//	3-residues fragment included in this file
//=========================================================================================
int Chain2LocalFrag(vector<LocalFrag>& localfrag, ChainCoord& chain)
{
	int i=0, j=0, k=0, m=0;
	

	//substract and record the information of all local fragments
	LocalFrag frag;
	double xyz[12][3];

	for(i=0; i<chain.resinum-2; i++)
	{
		frag.rname[0]=chain.resiname[i];
		frag.rname[1]=chain.resiname[i+1];
		frag.rname[2]=chain.resiname[i+2];
	
	
		xyz[0][0]=chain.xyz[i][0][0], xyz[0][1]=chain.xyz[i][0][1], xyz[0][2]=chain.xyz[i][0][2];
		xyz[1][0]=chain.xyz[i][1][0], xyz[1][1]=chain.xyz[i][1][1], xyz[1][2]=chain.xyz[i][1][2];
		xyz[2][0]=chain.xyz[i][2][0], xyz[2][1]=chain.xyz[i][2][1], xyz[2][2]=chain.xyz[i][2][2];
		xyz[3][0]=chain.xyz[i][3][0], xyz[3][1]=chain.xyz[i][3][1], xyz[3][2]=chain.xyz[i][3][2];

		xyz[4][0]=chain.xyz[i+1][0][0], xyz[4][1]=chain.xyz[i+1][0][1], xyz[4][2]=chain.xyz[i+1][0][2];
		xyz[5][0]=chain.xyz[i+1][1][0], xyz[5][1]=chain.xyz[i+1][1][1], xyz[5][2]=chain.xyz[i+1][1][2];
		xyz[6][0]=chain.xyz[i+1][2][0], xyz[6][1]=chain.xyz[i+1][2][1], xyz[6][2]=chain.xyz[i+1][2][2];
		xyz[7][0]=chain.xyz[i+1][3][0], xyz[7][1]=chain.xyz[i+1][3][1], xyz[7][2]=chain.xyz[i+1][3][2];
		
		xyz[8][0]=chain.xyz[i+2][0][0], xyz[8][1]=chain.xyz[i+2][0][1], xyz[8][2]=chain.xyz[i+2][0][2];
		xyz[9][0]=chain.xyz[i+2][1][0], xyz[9][1]=chain.xyz[i+2][1][1], xyz[9][2]=chain.xyz[i+2][1][2];
		xyz[10][0]=chain.xyz[i+2][2][0], xyz[10][1]=chain.xyz[i+2][2][1], xyz[10][2]=chain.xyz[i+2][2][2];
		xyz[11][0]=chain.xyz[i+2][3][0], xyz[11][1]=chain.xyz[i+2][3][1], xyz[11][2]=chain.xyz[i+2][3][2];

		
		double tempangle[4];
		Get4Angle(tempangle, xyz);
		
		frag.angle[0]=tempangle[0];
		frag.angle[1]=tempangle[1];
		frag.angle[2]=tempangle[2];
		frag.angle[3]=tempangle[3];
	
		localfrag.push_back(frag);
	}

	
	return 0;
}




//=============================================================================================
//Calculate the set of Etrp energy---local triplet sequence-dependant enregy
//MIU=0.991, statisfied from 6000-single chains database
//=============================================================================================
double Etrp(ChainCoord& chain, DV4& trptable)
{
	int i=0, j=0, k=0, l=0, m=0, n=0, p=0;
	double outE=0.0;
	
	
	//read all local 3-residues' fragments
	vector<LocalFrag> allfrag;
	Chain2LocalFrag(allfrag, chain);
	
	
	//statify all local 3-residues fragments into bins
	int fragnum=allfrag.size();

	
	string RNAME[20]={"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
			     "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL" };//serials from AAindex

	int ANGLE=6;
	int ANGLE2=ANGLE*ANGLE, ANGLE3=ANGLE2*ANGLE, ANGLE4=ANGLE3*ANGLE;

	
	int angpos=0;	
	int a1=-1, a2=-1, a3=-1;
	int d1=-1, d2=-1, d3=-1, d4=-1;
	for(i=0; i<fragnum; i++)
	{
		a1=-1, a2=-1, a3=-1;
		d1=-1, d2=-1, d3=-1, d4=-1;
		
		for(j=0; j<20; j++)
		{
			if(allfrag[i].rname[0]==RNAME[j])
			{
				a1=j;
			}
			
			if(allfrag[i].rname[1]==RNAME[j])
			{
				a2=j;
			}
			
			if(allfrag[i].rname[2]==RNAME[j])
			{
				a3=j;
			}
		}
		
		if(a1==-1 || a2==-1 || a3==-1)
		{
			continue;
		}
		

		d1 = (int)floor(allfrag[i].angle[0]+180)/60;
		d2 = (int)floor(allfrag[i].angle[1]+180)/60;
		d3 = (int)floor(allfrag[i].angle[2])/30;
		d4 = (int)floor(allfrag[i].angle[3])/30;
		
		if(d1<0 || d1>5 || d2<0 || d2>5 || d3<0 || d3>5 || d4<0 || d4>5)
		{
			continue;
		}
		else
		{
			angpos=d1*ANGLE3+d2*ANGLE2+d3*ANGLE+d4;
		}
		
		outE += trptable[a1][a2][a3][angpos];
	}
	

	
	return outE;
}





//=========================================================================================
//Calculate the 6 angles between two 3-residues' fragments
//These 6 angles are between b(i-1)-b(j+1), b(i)-b(j), b(i+1)-b(j-1) and 
//	P(i-1)-P(j+1), P(i)-P(j), P(i+1)-P(j-1)
//=========================================================================================
int Get6Angle(Frag6Vector& frag1, Frag6Vector& frag2, double angle[6])
{
	if(frag2.rindex-frag1.rindex<2)
	{
		return -1;
	}
	
	
	angle[0] = VectorAngle(frag1.b1, frag2.b3);
	angle[1] = VectorAngle(frag1.b2, frag2.b2);
	angle[2] = VectorAngle(frag1.b3, frag2.b1);
	angle[3] = VectorAngle(frag1.P1, frag2.P3);
	angle[4] = VectorAngle(frag1.P2, frag2.P2);
	angle[5] = VectorAngle(frag1.P3, frag2.P1);
	
	if(angle[0]>180 || angle[1]>180 || angle[2]>180 || angle[3]>180 || angle[4]>180 || angle[5]>180)
	{
		return -2;
	}
	
	return 0;
} 




//=========================================================================================
//Get the angle between two space vector
//=========================================================================================
double VectorAngle(double v1[2][3], double v2[2][3])
{
	int i=0, j=0;
	double angle=0.0;
	
	double translation[3];
	translation[0] = v1[0][0]-v2[0][0];
	translation[1] = v1[0][1]-v2[0][1];
	translation[2] = v1[0][2]-v2[0][2];
	
	double point[3];
	point[0] = v1[1][0]-translation[0];
	point[1] = v1[1][1]-translation[1];
	point[2] = v1[1][2]-translation[2];
	
	angle = Angle(point, v2[0], v2[1]);
	
	return angle;
}





//=========================================================================================
//Substract 3 b-vectors and 3 P-vectors from a local 2-residues' fragment
//For fragment A(i)A(i+1), there are 6-vectors to represent their structures
//=========================================================================================
void Get6Vector_2resi(Frag6Vector& frag6vector, double xyz[12][3], int rindex)
{

	frag6vector.rindex=rindex;
	
	frag6vector.Oatom[0]=xyz[3][0];
	frag6vector.Oatom[1]=xyz[3][1];
	frag6vector.Oatom[2]=xyz[3][2];
	
	double translation[3];
	translation[0]=xyz[10][0]-xyz[11][0];
	translation[1]=xyz[10][1]-xyz[11][1];
	translation[2]=xyz[10][2]-xyz[11][2];
	
	double dco=translation[0]*translation[0]+translation[1]*translation[1]+translation[2]*translation[2];
	dco = sqrt(dco);
	
	frag6vector.Hatom[0] = xyz[8][0]+translation[0]/dco;
	frag6vector.Hatom[1] = xyz[8][1]+translation[1]/dco;
	frag6vector.Hatom[2] = xyz[8][2]+translation[2]/dco;
	
	
	double resi[4][3];
		
	resi[0][0]=xyz[0][0], resi[0][1]=xyz[0][1], resi[0][2]=xyz[0][2];
	resi[1][0]=xyz[1][0], resi[1][1]=xyz[1][1], resi[1][2]=xyz[1][2];
	resi[2][0]=xyz[2][0], resi[2][1]=xyz[2][1], resi[2][2]=xyz[2][2];
	resi[3][0]=xyz[3][0], resi[3][1]=xyz[3][1], resi[3][2]=xyz[3][2];
	
	Resi2Vector(frag6vector.b1, frag6vector.P1, resi);
	
	resi[0][0]=xyz[4][0], resi[0][1]=xyz[4][1], resi[0][2]=xyz[4][2];
	resi[1][0]=xyz[5][0], resi[1][1]=xyz[5][1], resi[1][2]=xyz[5][2];
	resi[2][0]=xyz[6][0], resi[2][1]=xyz[6][1], resi[2][2]=xyz[6][2];
	resi[3][0]=xyz[7][0], resi[3][1]=xyz[7][1], resi[3][2]=xyz[7][2];
	
	Resi2Vector(frag6vector.b2, frag6vector.P2, resi);
		
	resi[0][0]=xyz[8][0], resi[0][1]=xyz[8][1], resi[0][2]=xyz[8][2];
	resi[1][0]=xyz[9][0], resi[1][1]=xyz[9][1], resi[1][2]=xyz[9][2];
	resi[2][0]=xyz[10][0], resi[2][1]=xyz[10][1], resi[2][2]=xyz[10][2];
	resi[3][0]=xyz[11][0], resi[3][1]=xyz[11][1], resi[3][2]=xyz[11][2];
	
	Resi2Vector(frag6vector.b3, frag6vector.P3, resi);
	
}






//=================================================================================================
//Substract a b-vector and a P-vector from a residue include 4 atoms (N,CA,C,O)
//For residue A, there has two vectors to simplify, they are b-vector to bisect
//	angle N-CA-C and P-vector (the normal vector of the plane N-CA-C)
//================================================================================================
void Resi2Vector(double bv[2][3], double Pv[2][3], double resi[4][3])
{	
	
	int i=0, j=0, k=0, m=0;
	
	//get the P-A(i) vector---Pv[2][3]
	Pv[0][0]=resi[1][0];
	Pv[0][1]=resi[1][1];
	Pv[0][2]=resi[1][2];
	
	double x1, y1, z1, x2, y2, z2, x3, y3, z3;
	x1=resi[1][0], y1=resi[1][1], z1=resi[1][2];
	x2=resi[0][0], y2=resi[0][1], z2=resi[0][2];
	x3=resi[2][0], y3=resi[2][1], z3=resi[2][2];
	
	double A, B, C;
	A = (y2-y1)*(z3-z1)-(z2-z1)*(y3-y1);
	B = -( (x2-x1)*(z3-z1)-(z2-z1)*(x3-x1) );
	C = (x2-x1)*(y3-y1)-(y2-y1)*(x3-x1);
	
	Pv[1][0]=x1+1.0;
	Pv[1][1]=(B/A)*(Pv[1][0]-x1)+y1;
	Pv[1][2]=(C/A)*(Pv[1][0]-x1)+z1;
	
	
	//get the b-A(i) vector---bv[2][3]
	double ang = Angle(resi[0], resi[1], resi[2]);
	
	vector<double> point, axis1, axis2, outpoint;
	point.push_back(resi[0][0]);
	point.push_back(resi[0][1]);
	point.push_back(resi[0][2]);
	
	axis1.push_back(Pv[0][0]);
	axis1.push_back(Pv[0][1]);
	axis1.push_back(Pv[0][2]);

	axis2.push_back(Pv[1][0]);
	axis2.push_back(Pv[1][1]);
	axis2.push_back(Pv[1][2]);

	CoordinateRotation(point, axis1, axis2, ang/2, outpoint);
	
	point.clear();
	axis1.clear();
	axis2.clear();
	
	bv[0][0]=resi[1][0];
	bv[0][1]=resi[1][1];
	bv[0][2]=resi[1][2];
	
	bv[1][0]=outpoint[0];
	bv[1][1]=outpoint[1];
	bv[1][2]=outpoint[2];
	
}



//=========================================================================================
//For the input peptide coordinates information, get all the 6-vectors 
//	information of all the 2-residues fragment included in this file
//=========================================================================================
int Chain2FragVector_2resi(vector<Frag6Vector>& chain6v, ChainCoord& chain)
{
	int i=0, j=0, k=0, m=0;
	
	
	//substract and record the information of all local fragments
	Frag6Vector frag6v;
	double xyz[12][3];
	for(i=0; i<chain.resinum-2; i++)
	{
		xyz[0][0]=chain.xyz[i][0][0], xyz[0][1]=chain.xyz[i][0][1], xyz[0][2]=chain.xyz[i][0][2];
		xyz[1][0]=chain.xyz[i][1][0], xyz[1][1]=chain.xyz[i][1][1], xyz[1][2]=chain.xyz[i][1][2];
		xyz[2][0]=chain.xyz[i][2][0], xyz[2][1]=chain.xyz[i][2][1], xyz[2][2]=chain.xyz[i][2][2];
		xyz[3][0]=chain.xyz[i][3][0], xyz[3][1]=chain.xyz[i][3][1], xyz[3][2]=chain.xyz[i][3][2];

		xyz[4][0]=chain.xyz[i][1][0], xyz[4][1]=chain.xyz[i][1][1], xyz[4][2]=chain.xyz[i][1][2];
		xyz[5][0]=chain.xyz[i][2][0], xyz[5][1]=chain.xyz[i][2][1], xyz[5][2]=chain.xyz[i][2][2];
		xyz[6][0]=chain.xyz[i+1][0][0], xyz[6][1]=chain.xyz[i+1][0][1], xyz[6][2]=chain.xyz[i+1][0][2];
		xyz[7][0]=chain.xyz[i+1][1][0], xyz[7][1]=chain.xyz[i+1][1][1], xyz[7][2]=chain.xyz[i+1][1][2];
		
		xyz[8][0]=chain.xyz[i+1][0][0], xyz[8][1]=chain.xyz[i+1][0][1], xyz[8][2]=chain.xyz[i+1][0][2];
		xyz[9][0]=chain.xyz[i+1][1][0], xyz[9][1]=chain.xyz[i+1][1][1], xyz[9][2]=chain.xyz[i+1][1][2];
		xyz[10][0]=chain.xyz[i+1][2][0], xyz[10][1]=chain.xyz[i+1][2][1], xyz[10][2]=chain.xyz[i+1][2][2];
		xyz[11][0]=chain.xyz[i+1][3][0], xyz[11][1]=chain.xyz[i+1][3][1], xyz[11][2]=chain.xyz[i+1][3][2];
		
				
		int rindex=chain.resiserial[i+1];
		
		Get6Vector_2resi(frag6v, xyz, rindex);
		
		chain6v.push_back(frag6v);
	}

	
	return 0;
}


//=========================================================================================
//Calculate the set of Ehp energy---hydrogen-bondding enregy
//	,and decide the value of MIU, make the net sum nearly 0
//=========================================================================================
void OxyHyDis(double d[2], Frag6Vector& frag1, Frag6Vector& frag2)
{
	double doh=0.0, dho=0.0;
	
	doh = (frag1.Oatom[0]-frag2.Hatom[0])*(frag1.Oatom[0]-frag2.Hatom[0])+(frag1.Oatom[1]-frag2.Hatom[1])*(frag1.Oatom[1]-frag2.Hatom[1])+(frag1.Oatom[2]-frag2.Hatom[2])*(frag1.Oatom[2]-frag2.Hatom[2]);
	doh = sqrt(doh);

	dho = (frag1.Hatom[0]-frag2.Oatom[0])*(frag1.Hatom[0]-frag2.Oatom[0])+(frag1.Hatom[1]-frag2.Oatom[1])*(frag1.Hatom[1]-frag2.Oatom[1])+(frag1.Hatom[2]-frag2.Oatom[2])*(frag1.Hatom[2]-frag2.Oatom[2]);
	dho = sqrt(doh);
	
	d[0]=doh;
	d[1]=dho;
}





//=========================================================================================
//Calculate the set of Ehp_2resi energy---hydrogen-bondding enregy
//MIU gotten from 6000 single chains database
//=========================================================================================
double Ehb_2resi(ChainCoord& chain, DV6& hbtable)
{
	int i=0, j=0, k=0, l=0, m=0, n=0, p=0;
	double outE=0.0;
	
	
	//read all local 3-residues' fragments
	vector<Frag6Vector> allfrag;
	Chain2FragVector_2resi(allfrag, chain);

	
	int WIDTH=20;
	double HBDIS=2.5;
	double Doh=0.0;
	double tempang[6];
	int a0=-1, a1=-1, a2=-1, a3=-1, a4=-1, a5=-1;
	for(j=0; j<allfrag.size()-1; j++)
	{
		for(k=j+1; k<allfrag.size(); k++)
		{
			
			if(allfrag[k].rindex-allfrag[j].rindex<=2)
			{
				continue;
			}
			
			double D[2];
			OxyHyDis(D, allfrag[j], allfrag[k]);
			if(D[0]>HBDIS && D[1]>HBDIS)
			{
				continue;
			}
				
			int angflag=Get6Angle(allfrag[j], allfrag[k], tempang);
			if(angflag<0)
			{
				continue;
			}
			
	
			a0 = (int)floor(tempang[0])/WIDTH;
			a1 = (int)floor(tempang[1])/WIDTH;
			a2 = (int)floor(tempang[2])/WIDTH;
			a3 = (int)floor(tempang[3])/WIDTH;
			a4 =(int) floor(tempang[4])/WIDTH;
			a5 =(int) floor(tempang[5])/WIDTH;
			
			outE += hbtable[a0][a1][a2][a3][a4][a5];
		}
		
	}

	
	return outE;
}




//=========================================================================================
//According to the 4-atoms (N,CA,C,O) information and the centroid table,
//	add each residues' centroid atoms
//=========================================================================================
void AddCentroid(ChainCoord& chain, DV2& tolerance, DV4& centable)
{
	int i=0, j=0, k=0;
	
	string RNAME[19]={"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "HIS", "ILE",
			     "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL" };//serials from AAindex, except GLY
	
	//get value from centable (19*36*36*6)
	int a, b, c;
	double p1[3], p2[3], p3[3], p4[3];
	double PHI=0.0, PSI=0.0;
	vector<double> vtemp;
	for(i=0; i<chain.resinum; i++)
	{
		a=-1, b=-1, c=-1;
		
		//get a (residue type index)
		if(chain.resiname[i]=="GLY")
		{
			vtemp.clear();
			vtemp.push_back(0);
			vtemp.push_back(0);
			vtemp.push_back(0);
			tolerance.push_back(vtemp);
			
			continue;
		}
		else
		{
			for(j=0; j<19; j++)
			{
				if(chain.resiname[i]==RNAME[j])
				{
					a=j;
					break;
				}
			}
		}
		
		//get b and c (bins' index of PHI and PSI)
		if(i==0)
		{
			p1[0]=chain.xyz[i][0][0];p1[1]=chain.xyz[i][0][1];p1[2]=chain.xyz[i][0][2];
			p2[0]=chain.xyz[i][1][0];p2[1]=chain.xyz[i][1][1];p2[2]=chain.xyz[i][1][2];
			p3[0]=chain.xyz[i][2][0];p3[1]=chain.xyz[i][2][1];p3[2]=chain.xyz[i][2][2];
			p4[0]=chain.xyz[i+1][0][0];p4[1]=chain.xyz[i+1][0][1];p4[2]=chain.xyz[i+1][0][2];
			
			PSI = Dihedral(p1, p2, p3, p4);
			
			c = (int)(floor(PSI)+180)/10;
			
			for(j=0; j<36; j++)
			{
				if(centable[a][j][c][0]>1.0)
				{
					b=j;
					break;
				}
			}
			
			if(b==-1)
			{
				cout << "\nError: can not add centroid for residue " << i << endl;
				continue;
			}
			
		}
		else if(i==chain.resinum-1)
		{
			p1[0]=chain.xyz[i-1][2][0];p1[1]=chain.xyz[i-1][2][1];p1[2]=chain.xyz[i-1][2][2];		
			p2[0]=chain.xyz[i][0][0];p2[1]=chain.xyz[i][0][1];p2[2]=chain.xyz[i][0][2];
			p3[0]=chain.xyz[i][1][0];p3[1]=chain.xyz[i][1][1];p3[2]=chain.xyz[i][1][2];
			p4[0]=chain.xyz[i][2][0];p4[1]=chain.xyz[i][2][1];p4[2]=chain.xyz[i][2][2];
			
			
			PHI = Dihedral(p1, p2, p3, p4);
			
			b = (int)(floor(PHI)+180)/10;
			
			for(j=0; j<36; j++)
			{
				if(centable[a][b][j][0]>1.0)
				{
					c=j;
					break;
				}
			}
			
			if(c==-1)
			{
				cout << "\nError: can not add centroid for residue " << i << endl;
				continue;
			}
			
		}
		else
		{
			p1[0]=chain.xyz[i-1][2][0];p1[1]=chain.xyz[i-1][2][1];p1[2]=chain.xyz[i-1][2][2];		
			p2[0]=chain.xyz[i][0][0];p2[1]=chain.xyz[i][0][1];p2[2]=chain.xyz[i][0][2];
			p3[0]=chain.xyz[i][1][0];p3[1]=chain.xyz[i][1][1];p3[2]=chain.xyz[i][1][2];
			p4[0]=chain.xyz[i][2][0];p4[1]=chain.xyz[i][2][1];p4[2]=chain.xyz[i][2][2];
			
			
			PHI = Dihedral(p1, p2, p3, p4);
			
			b = (int)(floor(PHI)+180)/10;


			p1[0]=chain.xyz[i+1][0][0];p1[1]=chain.xyz[i+1][0][1];p1[2]=chain.xyz[i+1][0][2];
			
			PSI = Dihedral(p2, p3, p4, p1);
			
			c = (int)(floor(PSI)+180)/10;
			
		}
		
		double inter[3];
		inter[0]=centable[a][b][c][0];
		inter[1]=centable[a][b][c][2];
		inter[2]=centable[a][b][c][4];
		
		vtemp.clear();
		vtemp.push_back(centable[a][b][c][1]);
		vtemp.push_back(centable[a][b][c][3]);
		vtemp.push_back(centable[a][b][c][5]);
		
		tolerance.push_back(vtemp);
		
		
		p1[0]=chain.xyz[i][0][0];p1[1]=chain.xyz[i][0][1];p1[2]=chain.xyz[i][0][2];
		p2[0]=chain.xyz[i][2][0];p2[1]=chain.xyz[i][2][1];p2[2]=chain.xyz[i][2][2];
		p3[0]=chain.xyz[i][1][0];p3[1]=chain.xyz[i][1][1];p3[2]=chain.xyz[i][1][2];
		
		internal2cartesian(p1, p2, p3, inter, p4);
		
		vtemp.clear();
		vtemp.push_back(p4[0]);
		vtemp.push_back(p4[1]);
		vtemp.push_back(p4[2]);
		
		chain.xyz[i].push_back(vtemp);
	}
	
	
}


/*
//=========================================================================================
//Translate Internal data into Cartesian data for 
//	the last of the four points
//=========================================================================================
bool internal2cartesian (double *c1, double *c2, double *c3, double *p, double *c4)
{
   double G_EXTRA=1.0e-20;
   double d1[3];
   double d2[3];
   double xp[3];

   subtract(c1, c2, d1);
   subtract(c3, c2, d2);
   crossproduct(d1, d2, xp);

   if( (xp[0]<G_EXTRA && xp[0]<-G_EXTRA) && (xp[1]<G_EXTRA && xp[1]>-G_EXTRA) && (xp[2]<G_EXTRA && xp[2]>-G_EXTRA) )
    {
        cout<<"Error! Points 1, 2, 3 are in one line & no plane!\n";

        return false;        
    }

   double d3[3];
   double yp[3];
   double r[3] ;
   double ypp[3];
   double tmp1[3];

   crossproduct(d2, xp, yp);
   double ang1 = deg2rad(p[2]);
   norm(xp);
   norm(yp);
   multi(cos(ang1), xp);
   multi(sin(ang1), yp);
   vectorsum( xp, yp, r);

   crossproduct(d2, r, ypp);
   double ang2 = deg2rad(p[1]);
   norm(d2 );
   norm(ypp);
   multi(-cos(ang2), d2);
   multi(sin(ang2), ypp);
   vectorsum( d2, ypp, d3 );
   multi(p[0], d3);
   vectorsum(c3, d3, tmp1);

   c4[0] = tmp1[0];
   c4[1] = tmp1[1];
   c4[2] = tmp1[2];


   return true;

}
*/















