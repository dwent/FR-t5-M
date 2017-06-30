//===================================================================================
//
// Read information from files
//
// AUTHORS: daiwentao
// E-mail: daiwentao@moon.ibp.ac.cn
// Date: 2012.01.04
//              
//
//===================================================================================


#include "ReadInfor.hpp"


//===================================================================
// Read all IDs from the ID-list file
//===================================================================
int ReadID(vector<string>& allid, string IDfile)
{
	ifstream fcid;
	fcid.open(IDfile.c_str());
	if( !fcid.is_open() )
	{
		cout << "Can not open " << IDfile << endl;
		exit(-1); 
	}


	string stemp="";
	while( !fcid.eof() )
	{
		stemp="";
		fcid >> stemp;
		if(stemp!="") allid.push_back(stemp);
	}
	fcid.close();
	
	return 0;	
}

//==========================================================================
// Read sequence from fastafile;
//==========================================================================
void ReadSequence(const char* fastafile,string& sequence)
{
	ifstream fin(fastafile);
	if( !fin.is_open() )
	{
		cout << "\nCan not open fasta file in ReadFasta()!!\n";
		exit(-1);	
	}
	
	char buf[5000]; 
	char seq[5000]; 
	strcpy(seq, ""); 
	int seqnum=0;
	while( !fin.eof() )
	{
		strcpy(buf, "");
		fin.getline(buf, 5000);
		
		if(buf[0]!='>' && strlen(buf)>0)
		{
			strcat(seq, buf);
		}
	}
	
	seqnum=strlen(seq); 
	for(int i=0; i<seqnum; i++)
	{
		sequence.push_back(seq[i]);
	}
	
	fin.close(); 
}

//==========================================================================
// Read sequence from fastafile;
//==========================================================================
int ReadCM(string& cmfile,string& ss,string& rs, vector< vector<double> >& CMA, string& seq)
{
	ifstream fin;
	fin.open(cmfile.c_str());
	if( !fin.is_open() )
	{
		cout << "Can not open " << cmfile<< endl;
		return -1; 
	}
	
	int lnu=0;
	vector<double> lcm;
	while( !fin.eof() )
	{
		string stemp="";
		fin >> stemp;
		if(lnu<2) lnu++;
		else if(lnu==2)
		{
			//cout <<stemp<< endl;
			lnu++;
			ss=stemp;
		}
		else if(lnu==3)
		{
			lnu++;
			rs=stemp;
		}
		else lcm.push_back(atof(stemp.c_str()));
	}
	fin.close();
	//cout<< ss<<'\n'<<rs<<'\n'<<lcm.size()<<endl;
	//cout<<lcm[0]<<" "<<lcm[17688]<<endl;
	if(ss.size()==seq.size() && rs.size()==seq.size())
	{
		int lenP=ss.size();
		int rn=0;
		for(int i=0; i<lenP; i++)
		{
			vector<double> rcm;
			for(int j=rn*lenP; j<(rn+1)*lenP; j++) rcm.push_back(lcm[j]);
			CMA.push_back(rcm);
			rn++;		
		}
	}
	else cout<<"CM file doesn't match sequence"<<endl;
	
	return 0;
}

//==========================================================================
// Read sse and surface from *.sssa file;
//==========================================================================
int Readsssa(string& cmfile,string& ss,string& rs)
{
	ifstream fin;
	fin.open(cmfile.c_str());
	if( !fin.is_open() )
	{
		cout << "Can not open " << cmfile<< endl;
		return -1; 
	}
	
	int lnu=0;
	string lss="";
	string lrs="";
	while( !fin.eof() )
	{
		string stemp="";
		fin >> stemp;
		if(lnu<2) lnu++;
		else if(lnu==2)
		{
			//cout <<stemp<< endl;
			lnu++;
			lss=stemp;
		}
		else if(lnu==3)
		{
			lnu++;
			lrs=stemp;
		}		
	}
	
	if(lss.size() != ss.size() || lrs.size()!= rs.size())
	{
		cout<<"*.sssa file can't match!"<<endl;
	}
	else
	{
		for(int i=0; i<lss.size(); i++)
		{
			if(lrs[i]=='b') lrs[i]='-';
			if( (lss[i] == ss[i]) && (lrs[i] == rs[i]) ) continue;
			else
			{
				if(lss[i]!=ss[i])
				{
					if(lss[i] == 'H') ss[i]='h';
					else if(lss[i] == 'E') ss[i]='e';
					else ss[i]='c';
				}
				if(lrs[i]!=rs[i])
				{
					if(lrs[i]=='e') rs[i]='E';
					else rs[i]='b';
				}
			}
		}
	}
	
}

//===================================================================
// Read all XYZ-coordinates of atoms in *.pdb file
// Residue coordinates of N,CA,C,O,CB and et.
//===================================================================
int Readpdbcoord(vector<vector< vector<double> > >& coord, 
		    string name)
{
	int i=0, j=0, k=0;
	
	ifstream fin(name.c_str());
	if( !fin.is_open() )
	{
		cout << "Error: can not open " << name << endl;
		return -2;
	}
	
	vector<vector<double> > Rescoord;
	char buf[200];
	while( !fin.eof() )
	{
		fin.getline(buf, 200);
		string head="";
		head.push_back(buf[0]);
		head.push_back(buf[1]);
		head.push_back(buf[2]);
		head.push_back(buf[3]);
		
		if(head=="ATOM")
		{
			string xyz="";
			string smark="";
			smark.push_back(buf[13]);
			smark.push_back(buf[14]);
			vector<double> Atocoord;
			
			if(smark=="N ")
			{
				if(Rescoord.size() >= 4)
				{
					coord.push_back(Rescoord);
					Rescoord.clear();
				}
				xyz="";
				for(i=30; i<39; i++)
				{
					if(buf[i]!=' ') xyz.push_back(buf[i]);
				}
				Atocoord.push_back( atof(xyz.c_str()) );
				
				xyz="";
				for(i=39; i<47; i++)
				{
					if(buf[i]!=' ') xyz.push_back(buf[i]);
				}
				Atocoord.push_back( atof(xyz.c_str()) );				

				xyz="";
				for(i=47; i<55; i++)
				{
					if(buf[i]!=' ') xyz.push_back(buf[i]);
				}
				Atocoord.push_back( atof(xyz.c_str()) );	
				
				Rescoord.push_back(Atocoord);			
			}
			else
			{
				xyz="";
				for(i=30; i<39; i++)
				{
					if(buf[i]!=' ') xyz.push_back(buf[i]);
				}
				Atocoord.push_back( atof(xyz.c_str()) );
				
				xyz="";
				for(i=39; i<47; i++)
				{
					if(buf[i]!=' ') xyz.push_back(buf[i]);
				}
				Atocoord.push_back( atof(xyz.c_str()) );				

				xyz="";
				for(i=47; i<55; i++)
				{
					if(buf[i]!=' ') xyz.push_back(buf[i]);
				}
				Atocoord.push_back( atof(xyz.c_str()) );	
				
				Rescoord.push_back(Atocoord);	
			}
		}
	}
	fin.close();
	coord.push_back(Rescoord);

	
	return 0;
}

//===================================================================
// Read all XYZ-coordinates of atoms in *.pdb file
// Residue coordinates of N,CA,C,O.
//===================================================================
int ReadMCcoord(vector<vector< vector<double> > >& coord, 
		    string name)
{
	int i=0, j=0, k=0;
	
	ifstream fin(name.c_str());
	if( !fin.is_open() )
	{
		cout << "Error: can not open " << name << endl;
		return -2;
	}
	
	vector<vector<double> > Rescoord;
	char buf[200];
	while( !fin.eof() )
	{
		fin.getline(buf, 200);
		string head="";
		head.push_back(buf[0]);
		head.push_back(buf[1]);
		head.push_back(buf[2]);
		head.push_back(buf[3]);
		
		if(head=="ATOM")
		{
			string xyz="";
			string smark="";
			smark.push_back(buf[13]);
			smark.push_back(buf[14]);
			vector<double> Atocoord;
			
			if(buf[12] == ' ')
			{
			if(smark=="N ")
			{
				if(Rescoord.size() >= 4)
				{
					coord.push_back(Rescoord);
					Rescoord.clear();
				}
				xyz="";
				for(i=30; i<39; i++)
				{
					if(buf[i]!=' ') xyz.push_back(buf[i]);
				}
				Atocoord.push_back( atof(xyz.c_str()) );
				
				xyz="";
				for(i=39; i<47; i++)
				{
					if(buf[i]!=' ') xyz.push_back(buf[i]);
				}
				Atocoord.push_back( atof(xyz.c_str()) );				

				xyz="";
				for(i=47; i<55; i++)
				{
					if(buf[i]!=' ') xyz.push_back(buf[i]);
				}
				Atocoord.push_back( atof(xyz.c_str()) );	
				
				Rescoord.push_back(Atocoord);			
			}
			else if(smark=="CA" || smark=="C " || smark=="O " )
			{
				xyz="";
				for(i=30; i<39; i++)
				{
					if(buf[i]!=' ') xyz.push_back(buf[i]);
				}
				Atocoord.push_back( atof(xyz.c_str()) );
				
				xyz="";
				for(i=39; i<47; i++)
				{
					if(buf[i]!=' ') xyz.push_back(buf[i]);
				}
				Atocoord.push_back( atof(xyz.c_str()) );				

				xyz="";
				for(i=47; i<55; i++)
				{
					if(buf[i]!=' ') xyz.push_back(buf[i]);
				}
				Atocoord.push_back( atof(xyz.c_str()) );	
				
				Rescoord.push_back(Atocoord);	
			}
		
			}
		}
	}
	fin.close();
	coord.push_back(Rescoord);

	
	return 0;
}

//==========================================================================
// Read residue, second structure and ACC infomation from DSSP file;
//==========================================================================
int ReadDSSP(string& second, vector<int>& acc, string secondfile)
{
	string stemp;
	
	ifstream infileDssp(secondfile.c_str());
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
		
		//resi.push_back(stemp[13]);
		
		string str="";
		for(int i=35; i<38; i++)
		{
			if(stemp[i]!=' ') str.push_back(stemp[i]);
		}
		acc.push_back(atoi(str.c_str()));
		
		
		index++;
		
	}while(!infileDssp.eof());
	infileDssp.close();

	
	return 0;
}

//====================================================
//Get the surface residues
//====================================================
int SurfRes(string res, vector<int>& acc, string& SIs )
{
	int maxArea[20];
	maxArea[0] = 106; maxArea[1] = 142; maxArea[2] = 135; maxArea[3] = 163;
	maxArea[4] = 194; maxArea[5] = 197; maxArea[6] =  84; maxArea[7] = 184;
	maxArea[8] = 169; maxArea[9] = 227; maxArea[10]= 205; maxArea[11]= 164;
	maxArea[12]= 188; maxArea[13]= 157; maxArea[14]= 222; maxArea[15]= 136;
	maxArea[16]= 198; maxArea[17]= 248; maxArea[18]= 130; maxArea[19]= 142;
	
	for(int i=0; i<res.size(); i++)
	{
		//cout<<res[i];
		int j=0;
		if(res[i]=='V' || res[i]=='W' || res[i]=='Y' || res[i]=='X' )
		{
			
			if(res[i]=='V')  j='B'-65;
			if(res[i]=='W')  j='J'-65;
			if(res[i]=='Y')  j='O'-65;
			if(res[i]=='X')  j='H'-65;
				//cout<<j<<endl;
		}
		else
		j=res[i]-65;
		
		if(acc[i]>0.1*maxArea[j])
		{
			SIs.push_back('e');
		}
		else SIs.push_back('-');
	}
}

//====================================================
//Get the CA distance matrix
//====================================================
int GenCAmatrix(vector<vector<double> >& matrix, vector<vector< vector<double> > > xyz)
{
	for(int i=0; i<xyz.size(); i++)
	{
		//cout<<xyz[i][1][0]<<" "<<xyz[i][1][1]<<" "<<xyz[i][1][2]<<endl;
		vector<double> ldist;
		for(int j=0; j<xyz.size(); j++)
		{
			double dist = sqrt( (xyz[i][1][0]-xyz[j][1][0])*(xyz[i][1][0]-xyz[j][1][0])
									+(xyz[i][1][1]-xyz[j][1][1])*(xyz[i][1][1]-xyz[j][1][1])
									+(xyz[i][1][2]-xyz[j][1][2])*(xyz[i][1][2]-xyz[j][1][2]));
			ldist.push_back(dist);
		}
		matrix.push_back(ldist);
	}
}

//======================================================
//Clear struct QUE
//======================================================
void CleQUE(QUE& querry )
{
	querry.seq.clear();
	querry.ss.clear();
	querry.rs.clear();
	querry.CM8a.clear();
	querry.CM12a.clear();
	querry.sse.clear();
	querry.radius=0.0;
}

//=================================================
//Caculate the radius of model
//=================================================
double RGvalue(vector<vector< vector<double> > > xyz)
{
	int Resinum = xyz.size();
	double Center[3];
	double CAcoor[Resinum][3];
	double SumX =0.0;
	double SumY =0.0;
	double SumZ =0.0;
	
	for(int i=0; i<Resinum; i++)
	{
		CAcoor[i][0]=xyz[i][1][0];
		CAcoor[i][1]=xyz[i][1][1];
		CAcoor[i][2]=xyz[i][1][2];

		SumX +=CAcoor[i][0];
		SumY +=CAcoor[i][1];
		SumZ +=CAcoor[i][2];
	} 
	Center[0] =SumX/Resinum;
	Center[1] =SumY/Resinum;
	Center[2] =SumZ/Resinum;

	SumX =0.0;SumY =0.0;SumZ =0.0;
	
	for(int i=0; i<Resinum; i++)
	{
		SumX +=(CAcoor[i][0]-Center[0])*(CAcoor[i][0]-Center[0]);
		SumY +=(CAcoor[i][1]-Center[1])*(CAcoor[i][1]-Center[1]);
		SumZ +=(CAcoor[i][2]-Center[2])*(CAcoor[i][2]-Center[2]);
	} 

	double Sum = (SumX + SumY +SumZ)/Resinum;
	
	return pow(Sum,0.5);
}

//=================================================
//Caculate the radius of model
//=================================================
double PredRG(string sequence)
{
	return 2.2 * pow(sequence.size(),0.38);
}

//===============================
//Read TMscore
//===============================
double ReadTM(string name)
{
	ifstream fin(name.c_str());
	if( !fin.is_open() )
	{
		cout << "Error: can not open " << name << endl;
		return -2;
	}
	
	char buf[200];
	while( !fin.eof() )
	{
		fin.getline(buf, 200);
		string head="";
		head.push_back(buf[0]);
		head.push_back(buf[1]);
		head.push_back(buf[2]);
		head.push_back(buf[3]);
		head.push_back(buf[4]);
		head.push_back(buf[5]);
		head.push_back(buf[6]);
		head.push_back(buf[7]);
		
		if(head=="TM-score")
		{
			string score="";
			for(int i=14; i<20; i++) score.push_back(buf[i]);
			return atof(score.c_str());
		}
		
	}
	fin.close();
}
