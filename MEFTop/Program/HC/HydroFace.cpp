//======================================================================
//
// Description:
//	Generate hydrophobic-face on alpha-helix based on
//	 center-axis and on beta-strand based on residue-position.
//
// Contact: wuaiping@moon.ibp.ac.cn, 2010_05_11
//
//======================================================================


#include "HydroFace.hpp"


using namespace std;



//===============================================================
// Translate 20 amino acids into Hydrophobic or Polar (HP) 
//===============================================================
void Seq2HP(vector<char>& HPseq, vector<char> seq)
{
	HPseq.clear();
	for(int i=0; i<seq.size(); i++)
	{
		if(seq[i]=='A' || seq[i]=='I' || seq[i]=='L' || seq[i]=='V' || seq[i]=='G' || seq[i]=='F' || seq[i]=='M' || seq[i]=='C') HPseq.push_back('H');
		else HPseq.push_back('P');
	}
	
}



//===============================================================
// Generate the hydrophobic faces on the secondary element
//===============================================================
bool HydroFace(ELE& ele)
{
	int i=0;
	bool flag=false;
	const char* chseq=ele.HPseq.c_str();
	int len=ele.len;
	
	
	if(ele.type=='H')	//== recognize hydro-face on alpha-helix
	{
		vector<HFace> allhf;
		
		HelixHydroFace(allhf, ele.CA, ele.HPseq);
		
		if(!allhf.empty()) ele.hface=allhf;
	}
	else if(ele.type=='E')	//== recognize hydro-face on beta-strand
	{
		vector<int> hi1, hi2;
		for(i=0; i<len; i++)
		{
			if(i%2==0 && chseq[i]=='H') hi1.push_back(i);
			else if(i%2==1 && chseq[i]=='H') hi2.push_back(i);
		}
		
		if(hi1.size()>=2)
		{
			HFace hf;
			hf.ri=hi1;
			ele.hface.push_back(hf);
		}
		
		if(hi2.size()>=2)
		{
			HFace hf;
			hf.ri=hi2;
			ele.hface.push_back(hf);
		}
	}
	
	
	if(!ele.hface.empty()) flag=true;
	return flag;
}



//===============================================================
// Generate the hydrophobic faces on the alpha-helix
//===============================================================
void HelixHydroFace(vector<HFace>& hface,
		     const vector<vector<double > >& CA,
		     const string HPseq)
{
	int i=0, j=0, k=0, m=0;


	vector<int> hi;	//the index of each hydrophobic residue
	for(i=0; i<HPseq.size(); i++)
	{
		if(HPseq.substr(i,1)=="H") hi.push_back(i);
	}
	if(hi.size()<2) return;
		
		
	//== get center-axis
	vector<vector<double > > axis;
	HelixAxis(axis, CA);
		
	
	//== find the first hydro-face with maximum hydro-residues
	vector<vector<int > > faceindex;
	for(i=0; i<hi.size(); i++)
	{
		int ti=hi[i];
		
		vector<double> pa;		//point on axis;
		if(ti==0) pa=axis[0];
		else if(ti==HPseq.size()-1) pa=axis[axis.size()-1];
		else pa=axis[ti-1];
		
		vector<double> ps=CA[ti];	//point on surface relate to pa[];
		
		vector<double> vd;
		vd.push_back(ps[0]-pa[0]);
		vd.push_back(ps[1]-pa[1]);
		vd.push_back(ps[2]-pa[2]);
		
		// calculate a surface-axis that pass point CA[ti]
		vector<vector<double > > faceaxis;
		vector<double> vtemp;
		for(j=0; j<axis.size(); j++)
		{
			vtemp.clear();
			vtemp.push_back(axis[j][0]+vd[0]);
			vtemp.push_back(axis[j][1]+vd[1]);
			vtemp.push_back(axis[j][2]+vd[2]);
			
			faceaxis.push_back(vtemp);
		}
		
		// calculate the number of hydrophobic residues around this surface-axis
		double scope=1.6;
		vector<int> index;
		for(j=0; j<CA.size(); j++)
		{
			if(j==0)
			{
				if(CalDistance(CA[j], faceaxis[0])<scope)
				{
					if(HPseq.substr(j,1)=="H") index.push_back(j);
				}
			}
			else if(j==CA.size()-1)
			{
				if(CalDistance(CA[j], faceaxis[faceaxis.size()-1])<scope)
				{
					if(HPseq.substr(j,1)=="H") index.push_back(j);
				}
			}
			else
			{
				if(CalDistance(CA[j], faceaxis[j-1])<scope)
				{
					if(HPseq.substr(j,1)=="H") index.push_back(j);
				}			
			}
		}
		if(index.size()>=2) faceindex.push_back(index);
	}
	
	
	//== delete redundant faceindex[]
	vector<int> mark(faceindex.size());
	for(i=faceindex.size()-1; i>=1; i--)
	{
		bool flag=false;
		for(j=i-1; j>=0; j--)
		{
			vector<int> DI, MI;
			for(k=0; k<faceindex[i].size(); k++)
			{
				bool tf=false;
				for(m=0; m<faceindex[j].size(); m++)
				{
					if(faceindex[i][k]==faceindex[j][m]) {DI.push_back(faceindex[i][k]); tf=true; break;}
				}
				
				if(!tf) MI.push_back(faceindex[i][k]);
			}
			
			if(DI.size()==faceindex[i].size()) {flag=true; break;}
		/*	else if(DI.size()>=2 && mark[j]==0)
			{
				for(k=0; k<MI.size(); k++) faceindex[j].push_back(MI[k]);
				flag=true;
				mark[j]=1;
				break;
			}
		*/	else if(DI.size()!=0)
			{
				if(faceindex[i].size()>=faceindex[j].size())
				{
					for(k=0; k<DI.size(); k++)
					{
						for(m=faceindex[j].size()-1; m>=0; m--)
						{
							if(DI[k]==faceindex[j][m]) faceindex[j].erase(faceindex[j].begin()+m);
						}
					}
					
				}
				else
				{
					for(k=0; k<DI.size(); k++)
					{
						for(m=faceindex[i].size()-1; m>=0; m--)
						{
							if(DI[k]==faceindex[i][m]) faceindex[i].erase(faceindex[i].begin()+m);
						}
					}
					
				}
			}
		}
		
		if(flag) faceindex.erase(faceindex.begin()+i);
	}
	
	for(i=faceindex.size()-1; i>=0; i--)
	{
		if(faceindex[i].size()<2) faceindex.erase(faceindex.begin()+i);
	}
	

	//for(i=0; i<faceindex.size(); i++)
	//{
	//	cout << "Helix " << i+1 << " : ";
	//	for(j=0; j<faceindex[i].size(); j++)
	//		cout << faceindex[i][j] << ' ';
	//	cout << endl;
	//}
	
	
	for(i=0; i<faceindex.size(); i++)
	{
		HFace hf;
		hf.ri=faceindex[i];
		hface.push_back(hf);
	}
		
}



