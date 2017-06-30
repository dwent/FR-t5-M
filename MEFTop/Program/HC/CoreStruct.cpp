//=======================================================
//
// Description:
//	Seperate hydrophobic-core from protein and
//		output its structure.
//
// Contact: wuaiping@moon.ibp.ac.cn, 2010_05_19
//
//=======================================================


#include "CoreStruct.hpp"
#include "IOprotein.hpp"


using namespace std;



//===============================================================
// Collect relate information about hydro-core and seperate
//	it from the protein structure.
//===============================================================
void GetHcoreStr(HcoreStr& hc, 
		  const ChainCoord& chain, 
		  const vector<int>& hcid, 
		  const vector<char>& second, 
		  const vector<ELE>& allele,
		  int Fnum)
{
	int i=0, j=0, k=0, m=0;
	
	
	//== based on hcid[], collect all residues around hydrophobic-core (include many polar residues);
	vector<int> ti;
	for(i=0; i<hcid.size(); i++)
	{
		// check if hcid[i] belongs to a hydro-face
		bool flag=false;
		for(j=0; j<allele.size(); j++)
		{
			if(hcid[i]>=allele[j].pos && hcid[i]<allele[j].pos+allele[j].len)
			{
				for(k=allele[j].pos; k<allele[j].pos+allele[j].len; k++)
				{
					ti.push_back(k);
				}
				
				flag=true;
			}
			
			if(flag) break;
		}
		
		
		int cutoff=3;
		if(!flag)
		{
			if(hcid[i]>=cutoff)
			{
				for(j=hcid[i]-cutoff; j<=hcid[i]; j++) ti.push_back(j);
			}
			else for(j=0; j<=hcid[i]; j++) ti.push_back(j);
			
			if(hcid[i]<chain.resinum-cutoff)
			{
				for(j=hcid[i]+1; j<=hcid[i]+cutoff; j++) ti.push_back(j);
			}
			else for(j=hcid[i]+1; j<=chain.resinum-1; j++) ti.push_back(j);	
		}
	}
	
	//== delete redundant IDs and sort IDs
	for(i=ti.size()-1; i>0; i--)
	{
		for(j=i-1; j>=0; j--)
		{
			if(ti[i]==ti[j]) {ti.erase(ti.begin()+i); break;}
		}
	}
	sort(ti.begin(), ti.end());
	
	
	//== record related data for residues ti[] in struct HcoreStr.
	vector<int> hcmark;
	string gseq="";
	string g2nd="";
	vector<char> rname;
	vector<int> rserial;
	vector<vector<vector<double > > > xyz;
	
		// get hcmark[]
	for(i=0; i<ti.size(); i++)
	{
		bool flag=false;
		for(j=0; j<hcid.size(); j++)
		{
			if(ti[i]==hcid[j]) {flag=true; break;}
		}
		
		if(flag) hcmark.push_back(1);
		else hcmark.push_back(0);
	}
	
		// get gseq, g2nd, rname[], rserial[], xyz[][][];
	for(i=0; i<ti.size(); i++)
	{
		char ch=' ';
		AAname3to1(ch, chain.resiname[ti[i]]); //Modify the order
		gseq.push_back(ch);
		g2nd.push_back(second[ti[i]]);
		if(i!=ti.size()-1 && ti[i+1]-ti[i]!=1)
		{
			gseq.push_back('-'); 
			g2nd.push_back('-');
		}
		

		rname.push_back(ch);
		rserial.push_back(chain.resiserial[ti[i]]);
		xyz.push_back(chain.xyz[ti[i]]);
	}
	
	
	
	//typedef struct _HcoreStr
	//{
	//	vector<int> hcmark;		//hcmark[i]=1 or 0, "1" meas hydrophobic-core positon, "0" means not.
	//	
	//	string gapseq;			//sequence with gap mark "-", means chain with break.
	//	string gap2nd;			//second structure with gap mark "-", means chain with break.
	//	
	//	vector<char> rname;		//residue names, 1-letter.
	//	vector<int> rserial;		//residue serials, not begin from zero.
	//	
	//	vector<vector<vector<double > > > xyz;	//coordinates of each residue: N, CA, C, O, four atoms.
	//	
	//} HcoreStr;
	
	hc.hcmark=hcmark;
	hc.gapseq=gseq;
	hc.gap2nd=g2nd;
	hc.rname=rname;
	hc.rserial=rserial;
	hc.xyz=xyz;
	hc.Facenum=Fnum;
}


//===============================================================
// Call OutProtein() function, output hydro-core structures.
//===============================================================
void OutputHCore(string outname, const HcoreStr& hc)
{
	int i=0, j=0;
	const char* chname=outname.c_str();
	int resinum=hc.hcmark.size();
	int atomnum=resinum*4;
	char chainname=' ';
	int rserial[atomnum];
	double x[atomnum], y[atomnum], z[atomnum];
	string aname[atomnum], rname[atomnum];
	for(i=0; i<resinum; i++)
	{
		rserial[4*i+0]=hc.rserial[i];
		rserial[4*i+1]=hc.rserial[i];
		rserial[4*i+2]=hc.rserial[i];
		rserial[4*i+3]=hc.rserial[i];
		
		for(j=0; j<4; j++)
		{
			x[4*i+j]=hc.xyz[i][j][0];
			y[4*i+j]=hc.xyz[i][j][1];
			z[4*i+j]=hc.xyz[i][j][2];
		}
		
		aname[4*i+0]="N";
		aname[4*i+1]="CA";
		aname[4*i+2]="C";
		aname[4*i+3]="O";

		string str="";
		AAname1to3(hc.rname[i], str);
		rname[4*i+0]=str;
		rname[4*i+1]=str;
		rname[4*i+2]=str;
		rname[4*i+3]=str;
	}
	
	char fname[500];
	for(i=0; i<strlen(chname); i++) fname[i]=chname[i];
	fname[strlen(chname)]='\0';
	OutProtein(fname, atomnum, chainname, rserial, x, y, z, aname, rname);
	
	
	//== output record information (except coordinates)
	ofstream fout;
	char infname[500];
	for(i=0; i<strlen(chname)-3; i++) infname[i]=chname[i];
	infname[strlen(chname)-3]='i';
	infname[strlen(chname)-2]='n';
	infname[strlen(chname)-1]='f';			
	infname[strlen(chname)]='\0';	
	
	fout.open(infname);
	fout << "SEQ    " << hc.gapseq << endl;
	fout << "2ND    " << hc.gap2nd << endl;
	fout << "MARK   ";
	j=0;
	for(i=0; i<hc.gapseq.size(); i++)
	{
		if(hc.gapseq.substr(i,1)=="-") fout << '-';
		else {fout << hc.hcmark[j]; j++;}
	}
	fout << endl;
	fout << "Face   "<<hc.Facenum<<endl;
	fout << "Radius "<<hc.radius<<endl;
	
	
	fout.close();
}



//=============================================================
// Translate residue names from 1-lette to 3-letter form
//=============================================================
void AAname1to3(char name1, string& name3)
{

	switch( name1 )
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


