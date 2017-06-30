//===================================================================================
//
// Calculate the NACACO score as feature for MS
//
// AUTHORS: daiwentao
// E-mail: daiwentao@moon.ibp.ac.cn
// Date: 2012.06.18
//              
//
//===================================================================================

//===========================================================================================================
//Edited by Wentao Dai (Email: daiwentao@moon.ibp.ac.cn)
//Energy function based residue-five-beads model (N,CA,C,O, and centroid of side-chain)
//	1)Econ : atom-pairs contact MIU-potential;
//	2)Etrp : local 3-residues' sequence-dependent conformation MIU-potential;
//	3)Esolve  : solvation potential, realized by MIU-poteintial
//	4)a, b and c in Ehp: three weigth factors to balance different energy sets; 
//===========================================================================================================

#include "MathTool.hpp"
#include "IOprotein.hpp"
#include "NCACO.hpp"

using namespace std;

void ReadAllPar(DV2& vdwradtable, const string vdwfile, DV2& contable, const string confile, DV4& trptable, const string trpfile, DV2& solvetable, const string solvefile, DV4& centable, const string centroidINFOfile, DV4&struenergyset, const string stru_energy_file)
{
	int i,j,k;
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Begin to Read Energy parameters~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
	//Read INTERRESI and VDWRADIUS[][]
	int ATYPE=99;
	int INTERRESI;
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
	int LEVEL4=12*12*12*12;
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
	
	
	//Read Esol parameter table
	solvetable.resize(20);
	for(i=0; i<20; i++)
	{
		solvetable[i].resize(100);
	}
	
	ReadSolvPar(solvetable, solvefile);
	
	
	
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
	
	//read frag contact energy table; angle1 angle2 dihedral rotate_angle
	ReadStruEnergyTable(stru_energy_file,struenergyset);
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~End of Reading Energy parameters~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
}

//=============================================================================================================
//Read INTERRESI value and VDWRADIUS[23][23] matrix that calculated from
//	"../bin/getvdwmatrix ......"
//INTERRESI means unless how many residues between two atoms calculated as atom-atom pair
//VDWRADIUS[23][23] is the minimum distance matrix statified with INTERRESI condition
//=============================================================================================================
void ReadVdwRadius(int& INTERRESI, DV2& VDWRADIUS, const string vdwfile)
{
	int i=0, j=0;
	int NUM=VDWRADIUS.size();
	
	string s1,s2;
	string record[NUM];
	
	ifstream fin(vdwfile.c_str());
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
int ReadConPar(DV2& contable, const string parfile)
{
	int i=0;
	int NUM=contable.size();
	
	ifstream fin(parfile.c_str());
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
int ReadTrpPar(DV4& trptable, const string parfile)
{
	int a=0, b=0, c=0;
	int i=0, j=0;
	
	ifstream fin(parfile.c_str());
	if( !fin.is_open() )
	{
		cout << "\nReadTrpPar()--Error: can not open " << parfile << endl;
		exit(-1);
	}
	
	
	string data[12*12*12*12];
	for(i=0; i<8000; i++)
	{	
		a=i/400;
		b=(i-400*a)/20;
		c=i-400*a-20*b;
		
		for(j=0; j<12*12*12*12; j++)
		{
			fin >> data[j];
			if(data[j]=="inf") trptable[a][b][c][j] = 2.0;
			else trptable[a][b][c][j] = atof(data[j].c_str());
		}
	}
	fin.close();
	
	
	return 0;
}



//=====================================================================================================
//Read parameter Esol table;
//=====================================================================================================
int ReadSolvPar(DV2 &solvetable, const string solvefile)
{
	double dbuf;
	int ibuf;
	string sbuf;
	int i=0, j=0,k=0;
	
	ifstream fin(solvefile.c_str());
	if( !fin.is_open() )
	{
		cout << "\nReadSolvPar()--Error: can not open " << solvefile << endl;
		exit(-1);
	}
	
	
	for(j=0; j<20; j++)
	{
		fin>>sbuf;
		
		for(k=0;k<100;k++)
		{
			fin >>sbuf>>sbuf>>sbuf>>ibuf>>sbuf>>sbuf>>sbuf;
			fin >>sbuf;
			if(sbuf=="inf" || sbuf=="nan") solvetable[j][k]=4.0;
			else if(ibuf<10) solvetable[j][k]=0;
			else solvetable[j][k] = atof(sbuf.c_str());
		}
	}
	fin.close();
	
	return 0;
}



//=====================================================================================================
//Read parameter file of 19*36*36*6 bins of centroid table;
//
//=====================================================================================================
int ReadCentroid(DV4& centable, const string cenfile)
{
	int a=0, b=0, c=0;
	int i=0, j=0;
	
	ifstream fin(cenfile.c_str());
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

void ReadStruEnergyTable(const string energy_file,DV4 &StruEnergySet)
{
	ifstream infile(energy_file.c_str());
	if(!infile) cout<<"can not open stru energy file"<<endl;
		
	string sbuf;
	double dbuf;
	int i0,i1,i2,i3,i4,i5,i6;

	StruEnergySet.resize(4);
	for(i1=0;i1<StruEnergySet.size();i1++)
	{
		StruEnergySet[i1].resize(4);
		for(i2=0;i2<StruEnergySet[i1].size();i2++)
		{
			StruEnergySet[i1][i2].resize(6);
			for(i3=0;i3<StruEnergySet[i1][i2].size();i3++)
			{
				StruEnergySet[i1][i2][i3].resize(6);
			}
		}
	}

	getline(infile,sbuf);
	do
	{
		sbuf.clear();
		infile>>sbuf;
		if(sbuf.empty()) break;
			
		//angle1
		if(sbuf=="60") i0=0;
		else if(sbuf=="90") i0=1;
		else if(sbuf=="120") i0=2;
		else if(sbuf=="150") i0=3;
		else cout<<"error_stru_index1"<<endl;
		infile>>sbuf;

		//angle2
		infile>>sbuf;
		if(sbuf=="60") i1=0;
		else if(sbuf=="90") i1=1;
		else if(sbuf=="120") i1=2;
		else if(sbuf=="150") i1=3;
		else cout<<"error_stru_index2"<<endl;
		infile>>sbuf;
		
		//dihedral
		infile>>sbuf;
		if(sbuf=="0") i2=0;
		else if(sbuf=="30") i2=1;
		else if(sbuf=="60") i2=2;
		else if(sbuf=="90") i2=3;
		else if(sbuf=="120") i2=4;
		else if(sbuf=="150") i2=5;
		else cout<<"error_stru_index3"<<endl;
		infile>>sbuf;
		
		//rotate_angle
		infile>>sbuf;
		if(sbuf=="0") i3=0;
		else if(sbuf=="30") i3=1;
		else if(sbuf=="60") i3=2;
		else if(sbuf=="90") i3=3;
		else if(sbuf=="120") i3=4;
		else if(sbuf=="150") i3=5;
		else cout<<"error_stru_index4"<<endl;
		infile>>sbuf;
				
		//energy
		infile>>sbuf>>sbuf>>sbuf;
		
		if(sbuf=="inf") StruEnergySet[i0][i1][i2][i3]=10.0;	
		else StruEnergySet[i0][i1][i2][i3]=atof(sbuf.c_str());
		
	}while(!infile.eof());
	infile.close();
}

//=============================================================================================
//Read protein structure pdb file
//Substract the infomation of four-beads (N, CA, C, O) of each residue
//=============================================================================================
int Read4AtomPdb(ChainCoord& chain, char* pdbfile)
{
	int i=0, j=0, k=0, l=0;

	ifstream fin(pdbfile);
	if( !fin.is_open() )
	{
		cout << "\nError: can not open " << pdbfile << endl;
		return -1;
	}

/*
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
*/	
	int atomindex=0;
	string tmpstr;
	while( getline(fin,tmpstr) )
	{	
		
		if(tmpstr.length()>10)
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
	//cout<<resinum<<" "<<atomnum<<endl;
	
	
	double x[atomnum], y[atomnum], z[atomnum];
	string aname[atomnum], rname[atomnum];
	int rserial[atomnum];
	
	pdbobj.ReadCoord(pdbfile, x, "XCOORD");
	pdbobj.ReadCoord(pdbfile, y, "YCOORD");
	pdbobj.ReadCoord(pdbfile, z, "ZCOORD");
	pdbobj.ReadName(pdbfile, aname, "ATOMNAME");
	pdbobj.ReadName(pdbfile, rname, "RESINAME");
	pdbobj.ReadSerial(pdbfile, rserial, "RESISERIAL");

	//cout<<"good!"<<endl;
	vector<vector<double > > resi;
	vector<double> atom;
	int weight=1;
	double radius=0.0;
	double cc[3];

	//bool hasOyet=false;//for each residue, has read O atom already? 
	const int ONAMENUM=2;
	//define a priority list
	string Onamels[ONAMENUM]={"O","O1"};
	//because O atom may be replaced more than once, so...
	vector<double> Oatom;
	Oatom.resize(3);
	//CurrentOAtomRank should be larger than ONAMENUM
	int CurrentOAtomRank=ONAMENUM+5;
	for(i=0; i<atomnum; i++)
	{
		atom.clear();
		if(aname[i]=="N")
		{
			if(i>0)
			{
				resi.push_back(Oatom);
				chain.xyz.push_back(resi);
				resi.clear();

				//because every residue contains only one O atom, so all O atom assistant variable should be restore
				//hasOyet=false;
				Oatom.clear();
				Oatom.resize(3);
				CurrentOAtomRank=ONAMENUM+5;
				
			}
			
			atom.clear();
			atom.push_back(x[i]);
			atom.push_back(y[i]);
			atom.push_back(z[i]);
			
			resi.push_back(atom);
						
			chain.resiname.push_back(rname[i]);
			chain.resiserial.push_back(rserial[i]);
			
		}
		else if(aname[i]=="CA" || aname[i]=="C")
		{
			atom.push_back(x[i]);
			atom.push_back(y[i]);
			atom.push_back(z[i]);
			
			resi.push_back(atom);
		}
		else if (aname[i].compare(0,1,"O")==0)
		{
			int ORank=ONAMENUM+5;
			//find aname location in priority list
			for(int lci=0;lci<ONAMENUM;lci++)
			{
				if(Onamels[lci]==aname[i])
				{
					ORank=lci;
					break;
				}				
			}
			//if an O atom has a higher priority, then replace it
			if(ORank<CurrentOAtomRank)
			{
				Oatom[0]=x[i];
				Oatom[1]=y[i];
				Oatom[2]=z[i];
				CurrentOAtomRank=ORank;
			}
		}
		
		if(i==atomnum-1)
		{
			resi.push_back(Oatom);
			chain.xyz.push_back(resi);
			resi.clear();
		}
		
	}

	return 0;	
	
}

//===============================================================================
//According to the 4-atoms (N,CA,C,O) information and the centroid table,
//	add each residues' centroid atoms
//===============================================================================
void AddCentroid(ChainCoord& chain, DV4& centable)
{
	int i=0, j=0, k=0;
	
	string RNAME[19]={"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "HIS", "ILE",
				 "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL" };//serials from AAindex, except GLY
	
	
	//get value from centable (19*36*36*6)
	int a=-1, b=-1, c=-1;
	double p1[3], p2[3], p3[3], p4[3];
	double PHI=0.0, PSI=0.0;
	vector<double> vtemp;
	for(i=0; i<chain.resiname.size(); i++)
	{
		a=-1;
		b=-1;
		c=-1;

		
		//get a (residue type index)
		for(j=0; j<19; j++)
		{
			if(chain.resiname[i]==RNAME[j])
			{
				a=j;
				break;
			}
		}
		
		if(a==-1) continue;
		
		
		
		//get b and c (bins' index of PHI and PSI)
		if(i==0)
		{
			p1[0]=chain.xyz[i][0][0];p1[1]=chain.xyz[i][0][1];p1[2]=chain.xyz[i][0][2];
			p2[0]=chain.xyz[i][1][0];p2[1]=chain.xyz[i][1][1];p2[2]=chain.xyz[i][1][2];
			p3[0]=chain.xyz[i][2][0];p3[1]=chain.xyz[i][2][1];p3[2]=chain.xyz[i][2][2];
			p4[0]=chain.xyz[i+1][0][0];p4[1]=chain.xyz[i+1][0][1];p4[2]=chain.xyz[i+1][0][2];
			
			PSI = Dihedral(p1, p2, p3, p4);
			
			c = (floor(PSI)+180)/10;
			
			if(c==36) c=35;
			
			if(c<0 || c>35) continue;
			
			
			for(j=0; j<36; j++)
			{
				if(centable[a][j][c][0]>1.0)
				{
					b=j;
					break;
				}
			}
			
		}
		else if(i==chain.resinum-1)
		{
			p1[0]=chain.xyz[i-1][2][0];p1[1]=chain.xyz[i-1][2][1];p1[2]=chain.xyz[i-1][2][2];		
			p2[0]=chain.xyz[i][0][0];p2[1]=chain.xyz[i][0][1];p2[2]=chain.xyz[i][0][2];
			p3[0]=chain.xyz[i][1][0];p3[1]=chain.xyz[i][1][1];p3[2]=chain.xyz[i][1][2];
			p4[0]=chain.xyz[i][2][0];p4[1]=chain.xyz[i][2][1];p4[2]=chain.xyz[i][2][2];
			
			
			PHI = Dihedral(p1, p2, p3, p4);
			
			b = (floor(PHI)+180)/10;
			
			if(b==36) b=35;
			
			if(b<0 || b>35) continue;

			
			for(j=0; j<36; j++)
			{
				if(centable[a][b][j][0]>1.0)
				{
					c=j;
					break;
				}
			}
		
		}
		else
		{	//cout<<"Error isn't here!"<<endl;
			//cout<<chain.xyz[i][2][0]<<" "<<endl;
			p1[0]=chain.xyz[i-1][2][0];p1[1]=chain.xyz[i-1][2][1];p1[2]=chain.xyz[i-1][2][2];		
			p2[0]=chain.xyz[i][0][0];p2[1]=chain.xyz[i][0][1];p2[2]=chain.xyz[i][0][2];
			p3[0]=chain.xyz[i][1][0];p3[1]=chain.xyz[i][1][1];p3[2]=chain.xyz[i][1][2];
			p4[0]=chain.xyz[i][2][0];p4[1]=chain.xyz[i][2][1];p4[2]=chain.xyz[i][2][2];
			
			
			
			PHI = Dihedral(p1, p2, p3, p4);
			
			b = (floor(PHI)+180)/10;


			p1[0]=chain.xyz[i+1][0][0];p1[1]=chain.xyz[i+1][0][1];p1[2]=chain.xyz[i+1][0][2];
			
			PSI = Dihedral(p2, p3, p4, p1);
			
			c = (floor(PSI)+180)/10;
			
			if(b==36) b=35;
			if(c==36) c=35;
			
		}
		//cout<<"Residue "<<i<<" "<<RNAME[a]<<endl;
		if(b<0 || b>35 || c<0 || c>35)
		{
			continue;
		}
		else
		{
			double inter[3];
			inter[0]=centable[a][b][c][0];
			inter[1]=centable[a][b][c][2];
			inter[2]=centable[a][b][c][4];
				
		
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
	
}

//=====================================================================================================
//Calculate the atom-pairs' contact energy of residue's five-Beads model (99*99 types)
//
//=====================================================================================================
double EconOfFiveBead99_99(ChainCoord& chain, DV2& VDWRADIUS, DV2& contable, int &collidetnum)
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
		x1=chain.xyz[j][1][0];
		y1=chain.xyz[j][1][1];
		z1=chain.xyz[j][1][2];
			
		tempDIS.clear();
		for(l=j+1; l<len; l++)
		{
			x2=chain.xyz[l][1][0];
			y2=chain.xyz[l][1][1];
			z2=chain.xyz[l][1][2];
						
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
						else if(d<LowD2[i1][i2]) collidetnum++;
					}
				}
			}
		}
	}

	RI.clear();
	
	return outE;
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

	int ANGLE=12;
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
		

		d1 = floor(allfrag[i].angle[0]+180)/30;
		d2 = floor(allfrag[i].angle[1]+180)/30;
		d3 = floor(allfrag[i].angle[2]+180)/30;
		d4 = floor(allfrag[i].angle[3]+180)/30;
		
		if(d1<0 || d1>11 || d2<0 || d2>11 || d3<0 || d3>11 || d4<0 || d4>11)
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

//==========================================
//calculate Esolve
//==========================================
double Esolve(ChainCoord& chain,const vector<vector<double> >& solvetable)
{
	const double cutoff=12.5;
	
	int i,j;
	
	//transfer FragCartesian to vector<residueInfor>
	vector<residueInfor> residue;
	for(i=0;i<chain.resinum;i++)
	{
		residueInfor singleResi;
		singleResi.resid=chain.resiname[i];
		singleResi.index=chain.resiserial[i];
		singleResi.CAcoordinate=chain.xyz[i][1];
		
		residue.push_back(singleResi);
	}
	
	string RNAME[20]={"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
				 "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL" };//serials from AAindex
	
	double Esol_energy=0;
	for(i=0;i<residue.size();i++)
	{
		int itemp=0;
		for(j=0;j<residue.size();j++)
		{
			if(i==j) continue;
			if(MyDistance(residue[i].CAcoordinate,residue[j].CAcoordinate)<cutoff) itemp++;
		}
		
		if(itemp>99) {cout<<"ERROR too much contact"<<endl;exit(0);}
		
		int a1=-1;
		for(j=0; j<20; j++)
		{
			if(residue[i].resid==RNAME[j])
			{
				a1=j;
				break;
			}
		}
		
		Esol_energy+=solvetable[a1][itemp];
	}
	
	return Esol_energy;
}

//===========================================
//Calculate the beta contact score 
//===========================================
double Estru(ChainCoord& chain, const vector<vector<vector<vector<double > > > > &StruEnergySet)
{
	double StruEnergy=0;
	
	int i,j;
	
	map<string,char> residueName;
	residueName[string("ALA")]='A';
	residueName[string("ARG")]='R';
	residueName[string("ASN")]='N';
	residueName[string("ASP")]='D';
	residueName[string("CYS")]='C';
	residueName[string("GLN")]='Q';
	residueName[string("GLU")]='E';
	residueName[string("GLY")]='G';
	residueName[string("HIS")]='H';
	residueName[string("ILE")]='I';
	residueName[string("LEU")]='L';
	residueName[string("LYS")]='K';
	residueName[string("MET")]='M';
	residueName[string("PHE")]='F';
	residueName[string("PRO")]='P';
	residueName[string("SER")]='S';
	residueName[string("THR")]='T';
	residueName[string("TRP")]='W';
	residueName[string("TYR")]='Y';
	residueName[string("VAL")]='V';

	//transfer FragCartesian to vector<residueInfor>
	vector<residueInfor> residue;
	for(i=0;i<chain.resinum;i++)
	{
		residueInfor singleResi;
		singleResi.resid=chain.resiname[i];
		singleResi.index=chain.resiserial[i];
		singleResi.CAcoordinate=chain.xyz[i][1];
		
		residue.push_back(singleResi);
	}

	//whether continue frag
	vector<int> continueIndex;
	map<int,vector<double> > index2coor;
	map<int,char> index2aa;
	
	for(i=0;i<residue.size()-2;i++)
	{
		if((residue[i+1].index-residue[i].index==1)&&(residue[i+2].index-residue[i+1].index==1))
		{
			if((MyDistance(residue[i+1].CAcoordinate,residue[i].CAcoordinate)<4)&&
			   (MyDistance(residue[i+2].CAcoordinate,residue[i+1].CAcoordinate)<4))
			continueIndex.push_back(residue[i].index);
		}
	}

	for(i=0;i<residue.size();i++)
	{
		index2coor[residue[i].index]=residue[i].CAcoordinate;
		index2aa[residue[i].index]=residueName[residue[i].resid];
	}

	//stat pair info
	int frag1A,frag1B,frag1C,frag2A,frag2B,frag2C;//residue index for 2 frag 
	
	double dis1,dis2;
	
	//derive contact cutoff
	const double ContactDisCutoff=5.6;//if the distance between two frag is <ContactDisCutoff, define the two frag contact

	//item
	double angle1,angle2,dihedral,rotate_angle;
	double d1,d2,d3;
	string frag1seq,frag2seq;
	frag1seq.resize(3);frag2seq.resize(3);
	
	for(i=0;i<continueIndex.size();i++)
	{
		for(j=i+1;j<continueIndex.size();j++)
		{		
			if(continueIndex[j]-continueIndex[i]<6) continue;//not stat nearby frag pair, must be divided by >=3 residue
					
			frag1A=continueIndex[i],frag1B=continueIndex[i]+1,frag1C=continueIndex[i]+2;
			frag2A=continueIndex[j],frag2B=continueIndex[j]+1,frag2C=continueIndex[j]+2;
			
			//determine whether contact and contact pattern
			bool whetherContact1=false;//pattern1
			bool whetherContact2=false;//pattern2
			
			vector<vector<double> > vdtemp1,vdtemp2;
			vdtemp1.push_back(index2coor[frag1A]);vdtemp1.push_back(index2coor[frag1B]);vdtemp1.push_back(index2coor[frag1C]);
			vdtemp2.push_back(index2coor[frag2A]);vdtemp2.push_back(index2coor[frag2B]);vdtemp2.push_back(index2coor[frag2C]);
			whetherContact1=CalContact(ContactDisCutoff,vdtemp1,vdtemp2,dis1);
			
			vdtemp2.clear();
			vdtemp2.push_back(index2coor[frag2C]);vdtemp2.push_back(index2coor[frag2B]);vdtemp2.push_back(index2coor[frag2A]);
			whetherContact2=CalContact(ContactDisCutoff,vdtemp1,vdtemp2,dis2);
			
			if(whetherContact1||whetherContact2)//contact pair
			{
				if(whetherContact1&&whetherContact2)
				{
					if(dis1>dis2)//change pattern
					{
						frag2A=frag2A+2;
						frag2C=frag2C-2;
					}
				}
				else if((!whetherContact1)&&whetherContact2)//change pattern
				{
					frag2A=frag2A+2;
					frag2C=frag2C-2;
				}
				
				//define frag angle
				angle1=MyAngle(index2coor[frag1A],index2coor[frag1B],index2coor[frag1C]);
				angle2=MyAngle(index2coor[frag2A],index2coor[frag2B],index2coor[frag2C]);
				
				//define frag pair dihedral
				dihedral=MyDihedral(index2coor[frag1A],index2coor[frag1B],index2coor[frag1C],
									index2coor[frag2A],index2coor[frag2B],index2coor[frag2C]);
				
				//define frag pair rotate angle
				rotate_angle=MyRotateAngle(index2coor[frag1A],index2coor[frag1B],index2coor[frag2A],index2coor[frag2B]);
			}
			
			//add energy for contact energy
			if(whetherContact1||whetherContact2)
			{
				//for stru energy
				int i0,i1,i2,i3;
				
				if(angle1<90) i0=0;
				else if(angle1<120) i0=1;
				else if(angle1<150) i0=2;
				else i0=3;
				
				if(angle2<90) i1=0;
				else if(angle2<120) i1=1;
				else if(angle2<150) i1=2;
				else i1=3;
					
				if(dihedral<30) i2=0;
				else if(dihedral<60) i2=1;
				else if(dihedral<90) i2=2;
				else if(dihedral<120) i2=3;
				else if(dihedral<150) i2=4;
				else i2=5;
					
				if(rotate_angle<30) i3=0;
				else if(rotate_angle<60) i3=1;
				else if(rotate_angle<90) i3=2;
				else if(rotate_angle<120) i3=3;
				else if(rotate_angle<150) i3=4;
				else i3=5;
				
				StruEnergy+=StruEnergySet[i0][i1][i2][i3];
			}
		}
	}
	
	return StruEnergy;
}

bool CalContact(const double ContactDisCutoff,const vector<vector<double> > &CAcoordinate1,
				const vector<vector<double> > & CAcoordinate2,double &dis)
{
	int i,j;
	double dtemp=0;
	vector<double> vdtemp;
	dis=0;
	bool btemp=true;
	
	for(i=0;i<3;i++)
	{	
		dtemp=0;
		for(j=0;j<3;j++)
		{
			dtemp+=(CAcoordinate1[i][j]-CAcoordinate2[i][j])*(CAcoordinate1[i][j]-CAcoordinate2[i][j]);
		}
		dtemp=sqrt(dtemp);
			
		if(dtemp>ContactDisCutoff)
		{
			btemp=false;
		}
		
		dis+=dtemp;
	}
	
	dis=dis/3;
	return btemp;
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

//=========================================================================================
//Substract 2 dihedrals and 2 angles from a local 3-residues' fragment
//For fragment A(i)A(i+1)A(i+2), 4 dihedrals are Psi(i), Phi(i+1), Psi(i+1), Phi(i+2)
//=========================================================================================
int Get4Angle(double angle[4], double xyz[12][3])
{
	int i=0, j=0, k=0, m=0;
	
	double dih=0.0;
	
	dih = Dihedral(xyz[0], xyz[1], xyz[2], xyz[4]);
	angle[0]=dih;
	
	dih = Dihedral(xyz[2], xyz[4], xyz[5], xyz[6]);
	angle[1]=dih;
	
	dih = Dihedral(xyz[4], xyz[5], xyz[6], xyz[8]);
	angle[2]=dih;
	
	dih = Dihedral(xyz[6], xyz[8], xyz[9], xyz[10]);
	angle[3]=dih;
		
	return 0;
}

//=====================================================================
//Translate Internal data into Cartesian data for 
//	the last of the four points
//====================================================================
bool internal2cartesian (double *c1, double *c2, double *c3, double *p, double *c4)
{
	double	G_EXTRA=1.0e-10;
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

//=============================================================
// Global shuffling for the sequence
//=============================================================
void randomShuffle(vector<string>& array)  
{  
	int len=array.size();
	for(int i=0;i<len;i++)  
	{  
		srand ( time(NULL) );
		int j=rand()%(len-i)+i; 
		array[i].swap(array[j]); 
	}  
} 

//=============================================================================================
//Calculate the set of Etrp energy---local triplet sequence-dependant enregy
//MIU=0.991, statisfied from 6000-single chains database
//=============================================================================================
double Etrp_n(ChainCoord& chain, DV4& trptable)
{
	int i=0, j=0, k=0, l=0, m=0, n=0, p=0;
	double outE=0.0,norE=0.0;
	
	
	//read all local 3-residues' fragments
	vector<LocalFrag> allfrag;
	Chain2LocalFrag(allfrag, chain);
	
	
	//statify all local 3-residues fragments into bins
	int fragnum=allfrag.size();

	
	string RNAME[20]={"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
				 "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL" };//serials from AAindex

	int ANGLE=12;
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
		

		d1 = floor(allfrag[i].angle[0]+180)/30;
		d2 = floor(allfrag[i].angle[1]+180)/30;
		d3 = floor(allfrag[i].angle[2]+180)/30;
		d4 = floor(allfrag[i].angle[3]+180)/30;
		
		norE += trptable[a1][a2][a3][0];
		
		if(d1<0 || d1>11 || d2<0 || d2>11 || d3<0 || d3>11 || d4<0 || d4>11)
		{
			continue;
		}
		else
		{
			angpos=d1*ANGLE3+d2*ANGLE2+d3*ANGLE+d4;
		}
		
		outE += trptable[a1][a2][a3][angpos];
	}
	outE = exp(log(2)*(-(2*(outE-norE)/(norE+outE))*(2*(outE-norE)/(norE))));
	
	return outE;
	//return norE;
}

//===============================================================================
//Generate the distance matrix from residues' centroid atoms
//===============================================================================
void GenCENmatrix(ChainCoord& chain, vector<vector<double> >& Cmatrix)
{
		for(int i=0; i<chain.xyz.size(); i++)
		{
			cout<<i<<" "<<chain.xyz[i].size()<<endl;
			vector<double> ldist; 
			for(int j=0; j<chain.xyz.size(); j++)
			{
				//cout<<i<<" "<<chain.xyz[i].size()<<" "<<j<<" "<<chain.xyz[j].size()<<endl;
				double dist=0.0;
				if(chain.xyz[i].size()==5 && chain.xyz[j].size()==5)
				{
					dist= sqrt( (chain.xyz[i][4][0]-chain.xyz[j][4][0])*(chain.xyz[i][4][0]-chain.xyz[j][4][0])
									+(chain.xyz[i][4][1]-chain.xyz[j][4][1])*(chain.xyz[i][4][1]-chain.xyz[j][4][1])
									+(chain.xyz[i][4][2]-chain.xyz[j][4][2])*(chain.xyz[i][4][2]-chain.xyz[j][4][2]));
					ldist.push_back(dist);
				}
				else if(chain.xyz[i].size()==5 && chain.xyz[j].size()!=5)
				{
					dist= sqrt( (chain.xyz[i][4][0]-chain.xyz[j][1][0])*(chain.xyz[i][4][0]-chain.xyz[j][1][0])
									+(chain.xyz[i][4][1]-chain.xyz[j][1][1])*(chain.xyz[i][4][1]-chain.xyz[j][1][1])
									+(chain.xyz[i][4][2]-chain.xyz[j][1][2])*(chain.xyz[i][4][2]-chain.xyz[j][1][2]));
					ldist.push_back(dist);
				}
				else if(chain.xyz[i].size()!=5 && chain.xyz[j].size()==5)
				{
					dist= sqrt( (chain.xyz[i][1][0]-chain.xyz[j][4][0])*(chain.xyz[i][1][0]-chain.xyz[j][4][0])
									+(chain.xyz[i][1][1]-chain.xyz[j][4][1])*(chain.xyz[i][1][1]-chain.xyz[j][4][1])
									+(chain.xyz[i][1][2]-chain.xyz[j][4][2])*(chain.xyz[i][1][2]-chain.xyz[j][4][2]));
					ldist.push_back(dist);
				}
				else
				{
					dist= sqrt( (chain.xyz[i][1][0]-chain.xyz[j][1][0])*(chain.xyz[i][1][0]-chain.xyz[j][1][0])
									+(chain.xyz[i][1][1]-chain.xyz[j][1][1])*(chain.xyz[i][1][1]-chain.xyz[j][1][1])
									+(chain.xyz[i][1][2]-chain.xyz[j][1][2])*(chain.xyz[i][1][2]-chain.xyz[j][1][2]));
					ldist.push_back(dist);
				}
			}
			Cmatrix.push_back(ldist);
		}
}
