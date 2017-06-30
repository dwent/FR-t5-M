//===================================================================================
//
// Description:
//	1. cluster all hydrophobic residues and decide the hydrophobic-cores;
//	2. match the hydro-faces into hydro-cores;
//	3. define the geometry relation between hydro-faces in the same core.
//
// Contact    : wuaiping@moon.ibp.ac.cn, 2010_04_06
//
//===================================================================================


#include "HydroCluster.hpp"
#include "IOprotein.hpp"
#include "MathTool.hpp"
#include <algorithm>


using namespace std;



//=================================================================
//Read parameter file of 19*36*36*6 bins of centroid table;
//=================================================================
int ReadCentroid(DV4& centable, char* cenfile)
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



//==========================================================================
// Read residue, second structure and ACC infomation from DSSP file;
//==========================================================================
int ReadDSSP(vector<char>& resi, vector<char>& second, vector<double>& acc, const char* secondfile)
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
		if(stemp[10]!=' ') continue;
		
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
		
		resi.push_back(stemp[13]);
		
		string str="";
		for(int i=35; i<38; i++)
		{
			if(stemp[i]!=' ') str.push_back(stemp[i]);
		}
		acc.push_back(atof(str.c_str()));
		
		
		index++;
		
	}while(!infileDssp.eof());
	infileDssp.close();

	
	return 0;
}




//=============================================================================
//Read protein structure pdb file
//Substract the infomation of four-beads (N, CA, C, O) of each residue
//=============================================================================
int Read4AtomPdb(ChainCoord& chain, char* pdbfile)
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
		{
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


//=============================================================
// Translate residue names from 3-lette to 1-letter form
//=============================================================
void AAname3to1(char& name1, string name3)
{

	if(name3=="ALA") name1='A';
	else if(name3=="ARG") name1='R';
	else if(name3=="ASN") name1='N';
	else if(name3=="ASP") name1='D';
	else if(name3=="CYS") name1='C';
	else if(name3=="GLN") name1='Q';
	else if(name3=="GLU") name1='E';
	else if(name3=="GLY") name1='G';
	else if(name3=="HIS") name1='H';
	else if(name3=="ILE") name1='I';
	else if(name3=="LEU") name1='L';
	else if(name3=="LYS") name1='K';
	else if(name3=="MET") name1='M';
	else if(name3=="PHE") name1='F';
	else if(name3=="PRO") name1='P';
	else if(name3=="SER") name1='S';
	else if(name3=="THR") name1='T';
	else if(name3=="TRP") name1='W';
	else if(name3=="TYR") name1='Y';
	else if(name3=="VAL") name1='V';
	else name1='?';

}


//=============================================================
// Calculate the distance matrix and normalized vectors
//	matrix among all hydrophobic residues' centroids.
//=============================================================
void CentroidMatrix(DV2& Dcent, DV3& NVcent, ChainCoord chain, vector<int> hi)
{
	int i=0, j=0;
	double dist=0.0;
	int num=hi.size();
	vector<double> vtemp;
	
	Dcent.resize(num);
	NVcent.resize(num);
	for(i=0; i<num; i++)
	{
		Dcent[i].resize(num);
		NVcent[i].resize(num);
		for(j=0; j<num; j++) NVcent[i][j].resize(3);
	}
	
	for(i=0; i<num-1; i++)
	{
		int ti=hi[i];
		vector<double> vt1=chain.xyz[ti][chain.xyz[ti].size()-1];
		
		for(j=i+1; j<num; j++)
		{
			int tj=hi[j];
			vector<double> vt2=chain.xyz[tj][chain.xyz[tj].size()-1];
			
			dist=VDistance(vt1, vt2);
			Dcent[i][j]=dist;
			Dcent[j][i]=dist;
			
			vtemp.clear();
			vtemp.push_back( (vt2[0]-vt1[0])/dist );
			vtemp.push_back( (vt2[1]-vt1[1])/dist );
			vtemp.push_back( (vt2[2]-vt1[2])/dist );
			NVcent[i][j]=vtemp;
			vtemp[0] *= -1.0;
			vtemp[1] *= -1.0;
			vtemp[2] *= -1.0;
			NVcent[j][i]=vtemp;
		}
	}
	
}



//=============================================================
// Calculate the distance matrix and normalized vectors
//	matrix among all hydrophobic residues' centroids.
//=============================================================
void ClusterMap(vector<vector<int > >& map, 
		  const DV2& Dcent, 
		  const DV3& NVcent, 
		  const ChainCoord& chain, 
		  const vector<int>& hi)
{
	int i=0, j=0, k=0, m=0;
	double dist=0.0;
	double Dcut=2.5; // to be tested!
	int num=Dcent.size();
	map.resize(num);
	for(i=0; i<num; i++) map[i].resize(num);
	
	//map[i][j]=1 means residue i and j maybe belong to a same core
	
	// according to the distance, generate some points between Cc-Cc
	vector<double> point;
	vector<vector<double > > allp;
	for(i=0; i<num-1; i++)
	{
		int ti=hi[i];
		vector<double> vt1=chain.xyz[ti][chain.xyz[ti].size()-1];
		
		for(j=i+1; j<num; j++)
		{
			int tj=hi[j];
			vector<double> nv=NVcent[i][j];
			
		
			//==== new added constraint, 2010_05_14 ====//
			double dca=VDistance(chain.xyz[hi[i]][1], chain.xyz[hi[j]][1]);
			if(Dcent[i][j]>dca) continue;
			//==============================================//
		
			
			int num=0;
			if (Dcent[i][j]>15.0) continue;	//to be tested!
			else if(Dcent[i][j]<5.0) num=2;
			else if(Dcent[i][j]<10.0) num=5;
			else num=10;
			
			allp.clear();
			double nd=Dcent[i][j]/num;
			for(k=1; k<num; k++)
			{
				point.clear();
				point.push_back(vt1[0]+nv[0]*k*nd);
				point.push_back(vt1[1]+nv[1]*k*nd);
				point.push_back(vt1[2]+nv[2]*k*nd);
				
				allp.push_back(point);
			}
			
			//== 
			bool flag=true;
			for(k=0; k<allp.size(); k++)
			{
				for(m=0; m<chain.resinum; m++)
				{
					double d1=VDistance(allp[k], chain.xyz[m][0]);
					double d2=VDistance(allp[k], chain.xyz[m][1]);
					double d3=VDistance(allp[k], chain.xyz[m][2]);
					if(d1<Dcut || d2<Dcut || d3<Dcut) {flag=false; break;}
				}
				
				if(!flag) break;
			}
			
			if(flag)
			{
				map[i][j]=1;
				map[j][i]=1;
			}
		}
	}
}


//===============================================================
// Cluster all hydrophobic residues into hydrophobic cores
//===============================================================
void HydroCluster(vector<vector<int > >& cluster, 
		  vector<vector<int > > map, 
		  const DV2& Dcent,
		  const vector<int>& hi)
{	
	//== 疏水核心识别方法：
	//== 1.生成起始种子：选取有最大连接度的残基作为第一个类的起始；
	//== 2.从和起始残基相连的所有残基中选取发生contact具有最大连接度的第二个残基进行扩展；
	//== 3.记录剩余残基中和前两个残基中一个有contact的残基，加入类；
	//== 4.重复上述1-3步，直到收敛。	
	
	
	int i=0, j=0, k=0, m=0, n=0;
	vector<int> vtemp;
	double Dcut=7.0; //to be tested!
	
	vector<int> left=hi;
	vector<int> LI;
	for(i=0; i<left.size(); i++) LI.push_back(i);
	
	bool bflag=true;
	while(bflag)
	{
		//== generate the initiative residue
		vector<int> lnum;
		for(i=0; i<left.size(); i++)
		{
			n=0;
			for(j=0; j<left.size(); j++)
			{
				if(i==j) continue;
				else if(map[LI[i]][LI[j]]==1) n++;
			}
			
			lnum.push_back(n);
		}
		int ti=(int) (max_element(lnum.begin(), lnum.end())-lnum.begin());
		vtemp.push_back(left[ti]);
		
		
		//== get the second hydro-residue
		lnum.clear();
		for(i=0; i<left.size(); i++)
		{
			n=0;
			if(map[LI[ti]][LI[i]]==0 || i==ti) {lnum.push_back(n); continue;}
			else
			{
				for(j=0; j<left.size(); j++)
				{
					if(i==j) continue;
					else if(map[LI[i]][LI[j]]==1 && Dcent[LI[i]][LI[j]]<Dcut) n++;
				}
			}
			
			lnum.push_back(n);
		}
		int tj=(int) (max_element(lnum.begin(), lnum.end())-lnum.begin());
		if(lnum[tj]>0) vtemp.push_back(left[tj]);
		else tj=ti;
		
		
		//== find the contact hydrophobic residues
		for(i=0; i<left.size(); i++)
		{
			if(i==ti || i==tj) continue;
			bool flag=false;
			for(j=0; j<vtemp.size(); j++)
			{
				for(k=0; k<hi.size(); k++)
				{
					if(vtemp[j]==hi[k]) break;
				}
				
				if(Dcent[LI[i]][k]<Dcut) {flag=true; break;}
			}
			
			if(flag) vtemp.push_back(left[i]);
		}
		
		//== the next round or break 
		if(vtemp.size()<=3) break;
		else
		{
			cluster.push_back(vtemp);
			
			for(i=left.size()-1; i>=0; i--)
			{
				bool flag=false;
				for(j=0; j<vtemp.size(); j++)
				{
					if(left[i]==vtemp[j]) {flag=true; break;}
				}
				
				if(flag)
				{
					left.erase(left.begin()+i);
					LI.erase(LI.begin()+i);
				}
			}
			
			vtemp.clear();
		}
	}
	
}


//===================================================================
// Recognize the secondary elements of the input sequence
//===================================================================
void InitialSSE(vector<ELE>& allele, vector<char> Second, ChainCoord& chain)
{
	int i,j,k,len,num;
	len =Second.size();

	vector<int> Flag;
	for(i=0; i<len; i++) Flag.push_back(0);
	
	vector<char> type;
	vector<int> initpos;
	vector<int> length;

	for(i=0; i<len; i++)
	{
		if(Second[i]!='C' && Flag[i] == 0)
		{
			num=0;
			char Type =Second[i];
			type.push_back(Type);
			initpos.push_back(i);
			for(j=0; j<len-i;j++)
			{
				if(Second[i+j]==Type )
			       {	
					num ++;
					Flag[i+j] =1;
					if((i+j) == len-1)
                                	{
						length.push_back(num);
						break;
					}
				}
				else
                                {
					length.push_back(num);
					break;
				}
			}		
		}
	}
	
	for(i=0; i<type.size(); i++)
	{
		ELE ele;
		ele.type=type[i];
		ele.pos=initpos[i];
		ele.len=length[i];
		
		if(ele.type=='H' && ele.len>=7) allele.push_back(ele);
		else if(ele.type=='E' && ele.len>=4) allele.push_back(ele);
	}
	
	
	//== push CA data into allele
	for(i=0; i<allele.size(); i++)
	{
		vector<vector<double > > CA;
		for(j=0; j<allele[i].len; j++)
		{
			CA.push_back(chain.xyz[allele[i].pos+j][1]);
		}
		allele[i].CA=CA;
	}
}



//===================================================================
// Cluster all hydrophobic residues into hydrophobic cores based
//	on the hydrophobic-faces information on elements
//===================================================================
void FaceCluster(vector<vector<int > >& cluster,
		  vector<vector<int > > map, 
		  const DV2& Dcent,
		  const vector<int>& hi,
		  vector<ELE>& allele)
{
	int i=0, j=0, k=0, m=0, n=0;
	vector<int> vtemp;
	double Dcut=7.0; //to be tested!
	

	//== 疏水核心识别方法(有别于HydroCluster()函数，加入了hydro-face的约束)：
	//== 1.同一个beta-strand上的两个faces不能属于同一个hydro-cluster；(还未使用)
	//== 2.同一个hydro-face上的所有残基应包含在同一个hydro-cluster中。
	
//cout << "Pos 1" << endl;	
	//== construct the face-index vector for each hydrophobic residue
	vector<vector<int > > VF;
	for(i=0; i<hi.size(); i++)
	{
		vtemp.clear();
		bool flag=false;
		for(j=0; j<allele.size(); j++)
		{
			if(allele[j].hface.empty()) continue;
			
			for(k=0; k<allele[j].hface.size(); k++)
			{
				for(m=0; m<allele[j].hface[k].ri.size(); m++)
				{
					if(hi[i]==allele[j].hface[k].ri[m]+allele[j].pos)
					{
						vtemp=allele[j].hface[k].ri;
						vtemp.erase(vtemp.begin()+m);
						for(n=0; n<vtemp.size(); n++) vtemp[n]+=allele[j].pos;
						
						flag=true;
						break;
					}
				}
				if(flag) break;
			}
			if(flag) break;
		}
		VF.push_back(vtemp);
	}
	
//cout << "Pos 2" << endl;	
	//==================================================
	vector<int> left=hi;
	vector<int> LI;
	for(i=0; i<left.size(); i++) LI.push_back(i);
//cout << "Pos 3" << endl;	
	bool bflag=true;
	while(bflag)
	{
		//== generate the initiative residue
		vector<int> lnum;
		for(i=0; i<left.size(); i++)
		{
			n=0;
			for(j=0; j<left.size(); j++)
			{
				if(i==j) continue;
				else if(map[LI[i]][LI[j]]==1) n++;
			}
			
			lnum.push_back(n);
		}
		
		if(lnum.empty()) break;
		
		int ti=(int) (max_element(lnum.begin(), lnum.end())-lnum.begin());
		vtemp.push_back(left[ti]);
		if(!VF[LI[ti]].empty())
		{
			for(i=0; i<VF[LI[ti]].size(); i++)
			{
				bool rflag=Redundant(VF[LI[ti]][i], vtemp);
				if(!rflag) vtemp.push_back(VF[LI[ti]][i]);
			}
		}
		
	//cout << "Pos 4: " << vtemp.size() << endl;	
		//== get the second hydro-residue
		int tj=ti;
		if(vtemp.size()==1)
		{
			lnum.clear();
			for(i=0; i<left.size(); i++)
			{
				n=0;
				if(map[LI[ti]][LI[i]]==0 || i==ti) {lnum.push_back(n); continue;}
				else
				{
					for(j=0; j<left.size(); j++)
					{
						if(i==j) continue;
						else if(map[LI[i]][LI[j]]==1 && Dcent[LI[i]][LI[j]]<Dcut) n++;
					}
				}
		
				lnum.push_back(n);
			}
			
			if(!lnum.empty())
			{
				tj=(int) (max_element(lnum.begin(), lnum.end())-lnum.begin());
				if(lnum[tj]>0)
				{
					bool rflag=Redundant(left[tj], vtemp);
					if(!rflag) vtemp.push_back(left[tj]);
				}
			}
		}
		
	//cout << "Pos 5: " << vtemp.size() << endl;	
		//== find the contact hydrophobic residues
		for(i=0; i<left.size(); i++)
		{
			if(i==ti || i==tj) continue;
			bool flag=false;
			for(j=0; j<vtemp.size(); j++)
			{
				for(k=0; k<hi.size(); k++)
				{
					if(vtemp[j]==hi[k]) break;
				}
				
				if(Dcent[LI[i]][k]<Dcut) {flag=true; break;}
			}
			
			if(flag)
			{
				///////// restraints, added 2010_05_14 //////////
				if(vtemp.size()>=4)
				{
					int Nmark=0;
					for(j=0; j<vtemp.size(); j++)
					{
						for(k=0; k<hi.size(); k++)
						{
							if(vtemp[j]==hi[k]) break;
						}
					
						if(map[LI[i]][k]==1) Nmark++;
					}
					if(Nmark<2) continue;
				}
				else if(vtemp.size()>=8)
				{
					int Nmark=0;
					for(j=0; j<vtemp.size(); j++)
					{
						for(k=0; k<hi.size(); k++)
						{
							if(vtemp[j]==hi[k]) break;
						}
					
						if(map[LI[i]][k]==1) Nmark++;
					}
					if(Nmark<3) continue;
				}
				//////////////////////////////////////////////////////
				
			
				bool rflag=Redundant(left[i], vtemp);
				if(!rflag) vtemp.push_back(left[i]);
				if(!VF[LI[i]].empty())
				{
					for(j=0; j<VF[LI[i]].size(); j++)
					{
						rflag=Redundant(VF[LI[i]][j], vtemp);
						if(!rflag) vtemp.push_back(VF[LI[i]][j]);
					}
				}
			}
		}
		
	//cout << "Pos 6: " << vtemp.size() << endl;	
		//== the next round or break 
		if(vtemp.size()<=3) break;
		else
		{
			cluster.push_back(vtemp);
			
			for(i=left.size()-1; i>=0; i--)
			{
				bool flag=false;
				for(j=0; j<vtemp.size(); j++)
				{
					if(left[i]==vtemp[j]) {flag=true; break;}
				}
				
				if(flag)
				{
					left.erase(left.begin()+i);
					LI.erase(LI.begin()+i);
				}
			}
			
			vtemp.clear();
		}
		
		if(left.size()<=3) break;
	
	//cout << "Pos 7: " << vtemp.size() << endl;	
	}
	
	//cout << "Pos 8" << endl;	
}


//================================================
// Check if the id is included in the set
//================================================
bool Redundant(int id, vector<int> set)
{
	bool flag=false;
	
	for(int i=0; i<set.size(); i++)
	{
		if(id==set[i]) {flag=true; break;}
	}
	
	return flag;
}


//===============================================================================
// Calculate the number of hydro-faces that included in the hydro-cluster
//===============================================================================
int FaceInCluster(const vector<ELE>& allele, const vector<int>& cluster)
{
	int num=0;
	int i=0, j=0, k=0;
	
	for(i=0; i<allele.size(); i++)
	{
		for(j=0; j<allele[i].hface.size(); j++)
		{
			bool flag=false;
			for(k=0; k<allele[i].hface[j].ri.size(); k++)
			{
				for(int l=0; l<cluster.size(); l++)
				{
				if(cluster[l]==allele[i].pos+allele[i].hface[j].ri[k]) 
				{
					flag=true;
					break;
				}
				}
				if(flag) 
				{
					num++;
					break;
				}
			}
			
			//if(flag) num++;
		}
	}
	
	
	return num;
}


//===============================================================================
// Calculate the centroid and max radius of the sphere that include all 
//	residues in the hydrophobic cluster.
//===============================================================================
double ClusterCentroid(vector<double>& cent, 
			const vector<int>& cluster,
			const ChainCoord& chain)
{
	int i=0, ti=0;
	
	double x=0.0, y=0.0, z=0.0;
	for(i=0; i<cluster.size(); i++)
	{
		ti=cluster[i];
		
		if(chain.resiname[ti]=="GLY")
		{
			x += chain.xyz[ti][1][0];
			y += chain.xyz[ti][1][1];
			z += chain.xyz[ti][1][2];
		}
		else
		{
			x += chain.xyz[ti][4][0];
			y += chain.xyz[ti][4][1];
			z += chain.xyz[ti][4][2];
		}
	}
	
	cent.push_back(x/cluster.size());
	cent.push_back(y/cluster.size());
	cent.push_back(z/cluster.size());
	
	
	double Rmax=0.0;
	for(i=0; i<cluster.size(); i++)
	{
		ti=cluster[i];
		double tempR=0.0;
		if(chain.resiname[ti]=="GLY") tempR=VDistance(cent, chain.xyz[ti][1]);
		else tempR=VDistance(cent, chain.xyz[ti][4]);
		
		if(tempR>Rmax) Rmax=tempR;
	}
	
	
	return Rmax;
}



//===============================================================================
// Calculate the direction vector between all hydrophobic residues and
//	the centroid of the hydrophobic cluster.
//===============================================================================
void SphereVector(vector<vector<double > >& svect,
		const vector<double>& cent, 
		const vector<int>& cluster,
		const ChainCoord& chain)
{
	int i=0, j=0, k=0;
	vector<double> vtemp;
	
	for(i=0; i<cluster.size(); i++)
	{
		j=cluster[i];
		
		vtemp.clear();
		if(chain.resiname[j]=="GLY") k=1;
		else k=4;
		
		vtemp.push_back(cent[0]-chain.xyz[j][k][0]);
		vtemp.push_back(cent[1]-chain.xyz[j][k][1]);
		vtemp.push_back(cent[2]-chain.xyz[j][k][2]);
		
		svect.push_back(vtemp);
	}
}


//======================================================================================
// Calculate the number of residues of each quadrant around the pointed residue
//======================================================================================
void ResiQuadrant(vector<int>& Nquad, int ri, const ChainCoord& chain)
{
	int i=0, j=0;
	double Dcut=15.0;
	double D=0.0;
	
	vector<double> vtemp;	//Centroid atom
	if(chain.resiname[ri]=="GLY") vtemp=chain.xyz[ri][1];
	else vtemp=chain.xyz[ri][4];
	
	vector<int> num;
	num.resize(8);
	for(i=0; i<chain.resinum; i++)
	{
		if(i==ri) continue;
		
		vector<double> v2;
		if(chain.resiname[i]=="GLY") v2=chain.xyz[i][1];	//CA atom
		else v2=chain.xyz[i][1];
		D=VDistance(vtemp, v2);
		if(D<Dcut)
		{
			vector<double> vd;
			vd.push_back(v2[0]-vtemp[0]);
			vd.push_back(v2[1]-vtemp[1]);
			vd.push_back(v2[2]-vtemp[2]);
			
			j=0;
			if(vd[0]>0.0) j+=1;
			if(vd[1]>0.0) j+=2;
			if(vd[2]>0.0) j+=4;
			
			num[j]++;
		}
	}
	
	Nquad=num;
	
}



//======================================================================================
// Decide if the residue is buried or exposed according to its contact residues
//	in its eight quadrants.
//======================================================================================
bool IfResiBury(const vector<int>& Nquad)
{
	bool flag=true;
	
	for(int i=0; i<Nquad.size(); i++)
	{
		if(Nquad[i]<=1)
		{
			flag=false;
			break;
		}
	}
	
	return flag;
}


//====================================================
// Decide if the hydrophobic-face is buried or 
//	exposed according to its residues
//====================================================
bool IfFaceBury(int initpos, HFace& hface, const ChainCoord& chain)
{
	bool flag=true;
	int num=0;
	
	vector<int> Nquad;
	for(int i=0; i<hface.ri.size(); i++)
	{
		Nquad.clear();
		ResiQuadrant(Nquad, initpos+hface.ri[i], chain);
		bool tf=IfResiBury(Nquad);
		
		if(tf) num++;
	}
	
	if(num==0) flag=false;	//be carefule! "num==0" means the hydrophobic-face is buried even if only a residue is buried.
		
	
	return flag;
}


//========================================================
// Decide if the hydrophobic-core is buried or 
//	exposed according to its residues infomation
//========================================================
bool IfCoreBury(vector<int>& cluster, const ChainCoord& chain)
{
	bool flag=true;
	int num=0;
	
	vector<int> Nquad;
	for(int i=0; i<cluster.size(); i++)
	{
		int ti=cluster[i];
		Nquad.clear();
		
		ResiQuadrant(Nquad, ti, chain);
		bool tf=IfResiBury(Nquad);
		
		if(tf) num++;
	}
	
	if(num<cluster.size()*0.5 && num<10) flag=false;
		
	
	return flag;	
}



