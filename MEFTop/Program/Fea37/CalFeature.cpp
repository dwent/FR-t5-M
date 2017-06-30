//===================================================================================
//
// Caculate the feature score based on the querry and model
//
// AUTHORS: daiwentao
// E-mail: daiwentao@moon.ibp.ac.cn
// Date: 2012.03.12
//              
//
//===================================================================================


#include "CalFeature.hpp"

//================================================
//Secondary structure state of residue
//================================================
double resSS(string qss, string mss)
{
	//int hm=0,em=0,cm=0;
	int mn=0;
	for(int i=0; i<qss.size(); i++)
	{
		if(qss[i]=='h') qss[i]='H';	
		else if(qss[i]=='e') qss[i]='E';
		else if(qss[i]=='c') qss[i]='C';
		
		if(qss[i]==mss[i]) mn++;
		
	}
	return double(mn)/double(mss.size());

}

//============================================
//Buried states of residue
//============================================
double resBE(string qrs, string mrs)
{
	int mn=0;
	for(int i=0; i<qrs.size(); i++)
	{
		if(qrs[i]=='b') qrs[i]='-';	
		else if(qrs[i]=='E') qrs[i]='e';
		
		if(qrs[i]==mrs[i]) mn++;		
	}
	return double(mn)/double(mrs.size());
}

//================================================
//Secondary structure state of residue
//================================================
double resSS1(string qss, string mss)
{
	//int hm=0,em=0,cm=0;
	int mn=0,hn=0;
	for(int i=0; i<qss.size(); i++)
	{
		if(qss[i]==mss[i]) mn++;
		else if(qss[i]=='h') 
		{
			if(mss[i]=='H')
			{
				//mn++;
				hn++;
			}	
		}
		else if(qss[i]=='e') 
		{
			if(mss[i]=='E')
			{
				//mn++;
				hn++;
			}
		}
		else if(qss[i]=='c')
		{
			 if(mss[i]=='C')
			{
				//mn++;
				hn++;
			}
		}
		
		
		
	}
	return (double(mn)+0.5*hn)/double(mss.size());

}

//============================================
//Buried states of residue
//============================================
double resBE1(string qrs, string mrs)
{
	int mn=0,hn=0;
	for(int i=0; i<qrs.size(); i++)
	{
		if(qrs[i]==mrs[i]) mn++;	
		else if(qrs[i]=='b') 
		{
			if(mrs[i]=='-') hn++;
		}	
		else if(qrs[i]=='E')
		{
			if(mrs[i]=='e') hn++;
		}
			
	}
	return (double(mn)+0.5*hn)/double(mrs.size());
}

//============================================
//1D feature pair score
//============================================
void pairSSA(vector<double>& score, QUE querry, MOD model)
{
	//double hm=0.0,em=0.0,cm=0.0,hq=0.0,eq=0.0,cq=0.0,
	//sm=0.0,bm=0.0,sq=0.0,bq=0.0;
	vector<double> num;
	num.resize(10, 0.0);
	
	for(int i=0; i<querry.ss.size(); i++)
	{
		switch(querry.ss[i])
		{
			case 'H':
				num[0] = num[0]+1.0;
			break;
			case 'h':
				num[0] = num[0]+1.0;
			break;
			case 'E':
				num[1] = num[1]+1.0;
			break;
			case 'e':
				num[1] = num[1]+1.0;
			break;
			case 'C':
				num[2] = num[2]+1.0;
			break;
			case 'c':
				num[2] = num[2]+1.0;
			break;
		}
		switch(model.ss[i])
		{
			case 'H':
				num[5] = num[5]+1.0;
			break;
			case 'E':
				num[6] = num[6]+1.0;
			break;
			case 'C':
				num[7] = num[7]+1.0;
			break;
		}
		switch(querry.rs[i])
		{
			case '-':
				num[3] = num[3]+1.0;
			break;
			case 'b':
				num[3] = num[3]+1.0;
			break;
			case 'E':
				num[4] = num[4]+1.0;
			break;
			case 'e':
				num[4] = num[4]+1.0;
			break;
		}
		switch(model.rs[i])
		{
			case '-':
				num[8] = num[8]+1.0;
			break;
			case 'b':
				num[8] = num[8]+1.0;
			break;
			case 'E':
				num[9] = num[9]+1.0;
			break;
			case 'e':
				num[9] = num[9]+1.0;
			break;
		}
		
	}
	
	for(int i=0; i<num.size(); i++)
	{
		num[i] /= querry.ss.size();
		score[i]=num[i];
	}
	score[num.size()]=cosine(num);
	score[num.size()+1]=correl(num);
	score[num.size()+2]=expdist(num);
	score[num.size()+3]=dotprod(num);
}

//============================================
//cosine score
//============================================
double cosine(vector<double> ratio)
{
	double q_ave=0.0,m_ave=0.0,qm_sum=0.0,
	q_sqr=0.0,m_sqr=0.0, score=0.0;
	for(int i=0; i<ratio.size()/2; i++)
	{
		qm_sum += ratio[i]*ratio[i+ratio.size()/2];
		q_sqr += ratio[i]*ratio[i];
		m_sqr += ratio[i+ratio.size()/2]*ratio[i+ratio.size()/2];
	}
	q_sqr = sqrt(q_sqr);
	m_sqr = sqrt(m_sqr);
	if(q_sqr*m_sqr == 0) score=0.0;
	else
	{
		score = qm_sum/(q_sqr*m_sqr);
		if(score<0) score=0.0;
		else if(score>1) score=1.0;
	}
	return score;
}

//============================================
//Correlation score
//============================================
double correl(vector<double> ratio)
{
	double q_ave=0.0,m_ave=0.0,qm_sum=0.0,
	q_sqr=0.0,m_sqr=0.0, score=0.0;
	for(int i=0; i<ratio.size()/2; i++)
	{
		q_ave += ratio[i];
		m_ave += ratio[i+ratio.size()/2];
	}
	q_ave /= (ratio.size()/2);
	m_ave /= (ratio.size()/2);
	for(int i=0; i<ratio.size()/2; i++)
	{
		qm_sum += (ratio[i]-q_ave)*(ratio[i+ratio.size()/2]-m_ave);
		q_sqr += (ratio[i]-q_ave)*(ratio[i]-q_ave);
		m_sqr += (ratio[i+ratio.size()/2]-m_ave)*(ratio[i+ratio.size()/2]-m_ave);
	}
	q_sqr = sqrt(q_sqr);
	m_sqr = sqrt(m_sqr);
	if(q_sqr*m_sqr == 0) score=0.0;
	else
	{
		score = qm_sum/(q_sqr*m_sqr);
		if(score<-1) score=-1.0;
		else if(score>1) score=1.0;
	}
	return score;
}

//============================================
//Expdist score
//============================================
double expdist(vector<double> ratio)
{
	double score=0.0;
	for(int i=0; i<ratio.size()/2; i++)
	{
		score += (ratio[i]-ratio[i+ratio.size()/2])*(ratio[i]-ratio[i+ratio.size()/2]);
	}
	return exp(-sqrt(score));
}

//============================================
//Dotproduct score
//============================================
double dotprod(vector<double> ratio)
{
	double score=0.0;
	for(int i=0; i<ratio.size()/2; i++)
	{
		score += ratio[i]*ratio[i+ratio.size()/2];
	}
	return score;
}

//============================================
//Contact map score
//============================================
void cmFeature(vector<double>& score, vector< vector<double> > cm, 
vector< vector<double> > ca, double cutoff, int nu)
{
	int len=cm.size(),cnu=0;
	vector<double> connum;
	connum.resize(len*2, 0.0);
	vector<double> conord;
	conord.resize(len*2, 0.0);
	double normcon=0.0;
	for(int i=0; i<len; i++)
	{
		for(int j=0; j<len; j++)
		{
			if((j-i)>=6)
			{
				connum[i] += cm[i][j];
				conord[i] += cm[i][j]*(j-i);
				if(ca[i][j]<cutoff)
				{
					normcon += cm[i][j];
					cnu++;
					connum[i+len] += 1.0;
					conord[i+len] += j-i;
				}
			}
			
		}
	}
	if(cnu==0) normcon=0;
	else normcon /= cnu;
	score[nu]=normcon;
	score[nu+1]=cosine(connum);
	score[nu+2]=correl(connum);
	score[nu+3]=cosine(conord);
	score[nu+4]=correl(conord);
	
}

//============================================
//Contact map score of index
//============================================
void cmFeatureIn(vector<int> si, vector<double>& score, vector< vector<double> > cm, 
vector< vector<double> > ca, double cutoff, int nu)
{
	int len=cm.size(),cnu=0;
	vector<double> connum;
	connum.resize(len*2, 0.0);
	vector<double> conord;
	conord.resize(len*2, 0.0);
	double normcon=0.0;
	for(int i=0; i<len; i++)
	{
		for(int j=0; j<len; j++)
		{
			if((j-i)>=6)
			{
				connum[i] += cm[i][j];
				conord[i] += cm[i][j]*(j-i);
				if(si[i]!=-1 && si[j]!=-1) //This is index judgement
				{
				if(ca[si[i]][si[j]]<cutoff)
				{
					normcon += cm[i][j];
					cnu++;
					connum[i+len] += 1.0;
					conord[i+len] += j-i;
				}
				}
				else continue;
			}
			
		}
	}
	if(cnu==0) normcon=0;
	else normcon /= cnu;
	score[nu]=normcon;
	score[nu+1]=cosine(connum);
	score[nu+2]=correl(connum);
	score[nu+3]=cosine(conord);
	score[nu+4]=correl(conord);
	
}

//======================================================
//SSE length feature score and contact score
//======================================================
void slenFeature(vector<double>& score, QUE querry, MOD model, int nu)
{
	vector<double> Slen;
	Slen.resize(model.sse.size(), 0.0);
	int mindex=0;
	double wei=0.0, sl=0.0, slt=0.0;
	vector<int> SMindex;
	SMindex.resize(model.sse.size(), -1);
	for(int i=0; i<model.sse.size(); i++)
	{
		int delta=100;
		if(i==0)
		{
			for(int j=mindex; j<querry.sse.size(); j++)
			{
			if((i+1)<model.sse.size()) //Some structure has only 1 sse
			{
			if(abs(model.sse[i].pos-querry.sse[j].pos) < delta && abs(model.sse[i].pos-querry.sse[j].pos) < abs(model.sse[i+1].pos-querry.sse[j].pos))
			{
				delta=abs(model.sse[i].pos-querry.sse[j].pos);
				mindex=j;
				SMindex[i]=mindex;
			}
			}
			else
			{
			if(abs(model.sse[i].pos-querry.sse[j].pos) < delta )
			{
				delta=abs(model.sse[i].pos-querry.sse[j].pos);
				mindex=j;
				SMindex[i]=mindex;
			}
			}
			}
		}
		else
		{
			if(i<model.sse.size()-1)
			{
			for(int j=mindex+1; j<querry.sse.size(); j++)
			{
			if(abs(model.sse[i].pos-querry.sse[j].pos) < delta && abs(model.sse[i].pos-querry.sse[j].pos) < abs(model.sse[i+1].pos-querry.sse[j].pos) )
			{
				delta=abs(model.sse[i].pos-querry.sse[j].pos);
				mindex=j;
				SMindex[i]=mindex; //Align index
				break;
			}
			}
			}
			else
			{
			for(int j=mindex+1; j<querry.sse.size(); j++)
			{
			if(abs(model.sse[i].pos-querry.sse[j].pos) < delta  )
			{
				delta=abs(model.sse[i].pos-querry.sse[j].pos);
				mindex=j;
				SMindex[i]=mindex; //Align index
				break;
			}
			}
			}	
		}
		if(delta < 90)
		{
		Slen[i]=exp(log(2)*(-((model.sse[i].len-querry.sse[mindex].len)/(0.25*querry.sse[mindex].len))*((model.sse[i].len-querry.sse[mindex].len)/(0.25*querry.sse[mindex].len)))) ;
		if(model.sse[i].type != querry.sse[mindex].type) slt += 0.5*Slen[i]*model.sse[i].len; //SSE type penalty
		else slt += Slen[i]*model.sse[i].len;
		}
		else
		Slen[i]=0.0; //SSE string gap penalty
		wei += model.sse[i].len;
		sl += Slen[i]*model.sse[i].len;
	}
	if(wei==0.0)
	{
		score[nu]=0.0;
		score[nu+1] = 0.0;
	}
	else
	{
		score[nu]=sl/wei;
		score[nu+1] = slt/wei;
	}
	
	//for(int i=0; i<SMindex.size(); i++) /*if(SMindex[i]<0)*/ cout<<i<<" : "<<SMindex[i]<<endl;
	score[nu+2]=ssCon(querry, model, SMindex);
	//cout<<score[nu]<<" "<<score[nu+1]<<" "<<score[nu+2]<<endl;
}

//===================================================================
//SSE length feature score and contact score of index
//===================================================================
void slenFeatureIn(vector<int> ri, vector<double>& score, QUE querry, MOD model, int nu)
{
	vector<double> Slen;
	Slen.resize(model.sse.size(), 0.0);
	int mindex=0;
	double wei=0.0, sl=0.0, slt=0.0;
	vector<int> SMindex;
	SMindex.resize(model.sse.size(), -1);

	//if(model.sse.size()==0 || querry.sse.size()==0) cout<<"Here"<<endl;
	//cout<<"slen "<<model.sse.size()<<" "<<querry.sse.size()<<endl;
	for(int i=0; i<model.sse.size(); i++)
	{
		int delta=100; 
		if(i==0)
		{
			for(int j=mindex; j<querry.sse.size(); j++)
			{
			if((i+1)<model.sse.size()) //Some structure has only 1 sse
			{
			if(abs(ri[model.sse[i].pos]-querry.sse[j].pos-1) < delta && abs(ri[model.sse[i].pos]-querry.sse[j].pos-1) < abs(ri[model.sse[i+1].pos]-querry.sse[j].pos-1))
			{
				delta=abs(ri[model.sse[i].pos]-1-querry.sse[j].pos);
				mindex=j;
				SMindex[i]=mindex;				
			}
			}
			else
			{
			if(abs(model.sse[i].pos-querry.sse[j].pos) < delta )
			{
				delta=abs(model.sse[i].pos-querry.sse[j].pos);
				mindex=j;
				SMindex[i]=mindex;
			}
			}			
			}
		}
		else
		{
			if(i<model.sse.size()-1)
			{
			for(int j=mindex+1; j<querry.sse.size(); j++)
			{
			if(abs(ri[model.sse[i].pos]-querry.sse[j].pos-1) < delta && abs(ri[model.sse[i].pos]-querry.sse[j].pos-1) < abs(ri[model.sse[i+1].pos]-querry.sse[j].pos-1))
			{
				delta=abs(ri[model.sse[i].pos]-1-querry.sse[j].pos);
				mindex=j;
				SMindex[i]=mindex; //Align index
				break;
			}
			}
			}
			else
			{
				for(int j=mindex+1; j<querry.sse.size(); j++)
				{
					if(abs(ri[model.sse[i].pos]-querry.sse[j].pos-1) < delta )
					{
						delta=abs(ri[model.sse[i].pos]-1-querry.sse[j].pos);
						mindex=j;
						SMindex[i]=mindex; //Align index
						break;
					}
				}
			}		
		}
		if(delta < 90)
		{
		Slen[i]=exp(log(2)*(-((model.sse[i].len-querry.sse[mindex].len)/(0.25*querry.sse[mindex].len))*((model.sse[i].len-querry.sse[mindex].len)/(0.25*querry.sse[mindex].len)))) ;
		if(model.sse[i].type != querry.sse[mindex].type) slt += 0.5*Slen[i]*model.sse[i].len; //SSE type penalty
		else slt += Slen[i]*model.sse[i].len;
		}
		else
		Slen[i]=0.0; //SSE string gap penalty
		wei += model.sse[i].len;
		sl += Slen[i]*model.sse[i].len;
	}
	if(wei==0.0)
	{
		score[nu]=0.0;
		score[nu+1] = 0.0;
	}
	else
	{
		score[nu]=sl/wei;
		score[nu+1] = slt/wei;
	}
	
	//for(int i=0; i<SMindex.size(); i++) /*if(SMindex[i]<0)*/ cout<<i<<" : "<<SMindex[i]<<endl;
	score[nu+2]=ssCon(querry, model, SMindex);
	//cout<<score[nu]<<" "<<score[nu+1]<<" "<<score[nu+2]<<endl;
}

//============================================
//SSE contact feature score
//============================================
double ssCon(QUE querry, MOD model, vector<int> smindex)
{
	vector<double> qsi,msi;
	vector<int> qw,mw;
	double s1=0.0,s2=0.0,s3=0.0;
	for(int i=0; i<smindex.size(); i++)
	{
		if(smindex[i] == -1) continue;
		for(int j=0; j<smindex.size(); j++)
		{
			if((smindex[j] == -1) || i==j) continue;
			else
			{
				double mindist=20.0;
				int mcn=0,qcn=0;
				for(int k=model.sse[i].pos; k<(model.sse[i].pos+model.sse[i].len); k++)
				{
					for(int l=model.sse[j].pos; l<(model.sse[j].pos+model.sse[j].len); l++)
					{
						if(model.CAmatrix[k][l]<8.5)
						{
							mcn++;
							if(mindist>model.CAmatrix[k][l]) mindist=model.CAmatrix[k][l];
						}
					}
				}
				if(mindist<8.5)
				{
					double maxpro=0.0;
					for(int k=querry.sse[smindex[i]].pos; k<querry.sse[smindex[i]].pos+querry.sse[smindex[i]].len; k++)
					{
						for(int l=querry.sse[smindex[j]].pos; l<querry.sse[smindex[j]].pos+querry.sse[smindex[j]].len; l++)
						{
							if(querry.CM8a[k][l]>0.1)
							{
								qcn++;
								if(maxpro<querry.CM8a[k][l]) maxpro=querry.CM8a[k][l];
							}
						}
					}
					if(maxpro<0.1) maxpro=16.0; //Double declining balance method
					else
					{
						maxpro = 4.67612-5.51673*log(maxpro+0.083631);
						if(maxpro<3.8) maxpro=3.8;
					}
					qsi.push_back(maxpro);
					msi.push_back(mindist);
					qw.push_back(qcn); // weight could try len?
					mw.push_back(mcn);
				}
			}
			
		}
	}
	for(int i=0; i<qsi.size(); i++)
	{
		s1 += exp(log(2)*(-((msi[i]-qsi[i])/(qsi[i]))*((msi[i]-qsi[i])/(qsi[i]))))*qw[i];
		s2 += exp(log(2)*(-((msi[i]-qsi[i])/(msi[i]))*((msi[i]-qsi[i])/(msi[i]))))*mw[i];
		s3 += (qw[i]+mw[i]);
	}
	//cout<<s1<<" "<<s2<<" "<<s3<<endl;
	if(s3 == 0.0) return 0.0;
	else
	return ((s1+s2)/s3);
}

//============================================
//SSE contact feature score
//============================================
void sconFeature(vector<double>& score, QUE querry, MOD model, int nu)
{
	int len=model.sse.size(),cnu=0;
	vector<double> connum;
	connum.resize(len*2, 0.0);
	vector<double> conord;
	conord.resize(len*2, 0.0);
	double normcon=0.0;
	for(int i=0; i<len; i++)
	{
		for(int j=0; j<len; j++)
		{
			int rcn=0;
			vector<double> ssp; 
			ssp.resize(model.sse[i].len, 0.0);
			double mindist=20.0;
			vector<int> ssm;
			ssm.resize(model.sse[i].len, 0);
			if(i==j) continue;
			else
			{
				
				for(int k=model.sse[i].pos; k<(model.sse[i].pos+model.sse[i].len); k++)
				{
					int ri=k-model.sse[i].pos;
					for(int l=model.sse[j].pos; l<(model.sse[j].pos+model.sse[j].len); l++)
					{
						if(querry.CM8a[k][l]>0.1) 
						{
							rcn++;
							if(ssp[ri]<querry.CM8a[k][l]) ssp[ri]=querry.CM8a[k][l];
						}
						if(model.CAmatrix[k][l]<8.5)
						{
							ssm[ri]=1;
							if(mindist>model.CAmatrix[k][l]) mindist=model.CAmatrix[k][l];
						}
					}
				}
				if(rcn>=4)
				{
				double cn=0.0;
				for(int k=0; k<ssp.size(); k++)
				{
					cn += ssp[k];
				}
				connum[i] = cn/model.sse[i].len;
				conord[i] = (cn/model.sse[i].len)*abs(model.sse[i].ei-model.sse[j].ei);	
				}
				else 
				{
					connum[i]=0;
					conord[i]=0;
				}
				if(mindist<8.5)
				{
					double cn=0.0;
					for(int k=0; k<ssm.size(); k++) cn += ssm[k];
					if((cn/4.0)>1.0) cn=1.0;
					else cn /= 4.0;
					connum[i+len] = cn;
					conord[i+len] = cn*abs(model.sse[i].ei-model.sse[j].ei);
					normcon += connum[i];
					cnu++;
				}
				else
				{
					connum[i+len]=0;
					conord[i+len]=0;
				}
			}			
		}
	}
	if(cnu==0) normcon=0;
	else normcon /= cnu;
	score[nu]=normcon;
	score[nu+1]=cosine(connum);
	score[nu+2]=correl(connum);
	/*Order is equal for ever*/
	//score[nu+3]=cosine(conord);
	//score[nu+4]=correl(conord);
	
}

//============================================
//SSE contact feature score of index
//============================================
void sconFeatureIn(vector<int> resi, vector<double>& score, QUE querry, MOD model, int nu)
{
	int len=model.sse.size(),cnu=0;
	vector<double> connum;
	connum.resize(len*2, 0.0);
	vector<double> conord;
	conord.resize(len*2, 0.0);
	double normcon=0.0;
	for(int i=0; i<len; i++)
	{
		for(int j=0; j<len; j++)
		{
			int rcn=0;
			vector<double> ssp; 
			ssp.resize(model.sse[i].len, 0.0);
			double mindist=20.0;
			vector<int> ssm;
			ssm.resize(model.sse[i].len, 0);
			if(i==j) continue;
			else
			{
				
				for(int k=model.sse[i].pos; k<(model.sse[i].pos+model.sse[i].len); k++)
				{
					int ri=k-model.sse[i].pos;
					for(int l=model.sse[j].pos; l<(model.sse[j].pos+model.sse[j].len); l++)
					{
						if(querry.CM8a[resi[k]][resi[l]]>0.1) 
						{
							rcn++;
							if(ssp[ri]<querry.CM8a[resi[k]][resi[l]]) ssp[ri]=querry.CM8a[resi[k]][resi[l]];
						}
						
						if(model.CAmatrix[k][l]<8.5)
						{
							ssm[ri]=1;
							if(mindist>model.CAmatrix[k][l]) mindist=model.CAmatrix[k][l];
						}
						
					}
				}
				if(rcn>=4)
				{
				double cn=0.0;
				for(int k=0; k<ssp.size(); k++)
				{
					cn += ssp[k];
				}
				connum[i] = cn/model.sse[i].len;
				conord[i] = (cn/model.sse[i].len)*abs(model.sse[i].ei-model.sse[j].ei);	
				}
				else 
				{
					connum[i]=0;
					conord[i]=0;
				}
				if(mindist<8.5)
				{
					double cn=0.0;
					for(int k=0; k<ssm.size(); k++) cn += ssm[k];
					if((cn/4.0)>1.0) cn=1.0;
					else cn /= 4.0;
					connum[i+len] = cn;
					conord[i+len] = cn*abs(model.sse[i].ei-model.sse[j].ei);
					normcon += connum[i];
					cnu++;
				}
				else
				{
					connum[i+len]=0;
					conord[i+len]=0;
				}
			}			
		}
	}
	if(cnu==0) normcon=0;
	else normcon /= cnu;
	score[nu]=normcon;
	score[nu+1]=cosine(connum);
	score[nu+2]=correl(connum);
	/*Order is equal for ever*/
	//score[nu+3]=cosine(conord);
	//score[nu+4]=correl(conord);
	
}

//================================================
//Radius feature
//================================================
void radFeature(vector<double>& score, double rad1, double rad2)
{
	score[score.size()-2]=exp(log(2)*(-(2*(rad2-rad1)/(rad1+rad2))*(2*(rad2-rad1)/(rad1+rad2))));
	score[score.size()-1]=exp(log(2)*(-(2*(rad2-rad1)/(rad1))*(2*(rad2-rad1)/(rad1))));
}

//===================================================================
// Read infromation from the inf file
//===================================================================
int ReadInf(string& cseq, string& chf, string& cmark, int& num, double& rad, const char* inffile)
{
	string stemp;
	
		ifstream infileinf(inffile);
		if(!infileinf)
		{
			cerr<<"ReadInf()---error: can not open inf file " << inffile << endl;
			exit(-1);
		}
		
		do
		{
			infileinf>>stemp;
			
			if(stemp=="SEQ") 
			{
				getline(infileinf, stemp);
				//cout<<stemp[4]<<endl;
				cseq.append(stemp.begin()+4, stemp.end());
			}
			else if(stemp=="2ND")
			{
				getline(infileinf, stemp);
				chf.append(stemp.begin()+4, stemp.end());
			} 
			else if(stemp=="MARK")
			{
				getline(infileinf, stemp);
				cmark.append(stemp.begin()+3, stemp.end());
			}
			else if(stemp=="Face")
			{
				getline(infileinf, stemp);
				num=atoi(stemp.c_str());
				//cout<<num<<endl;
			}
			else if(stemp=="Radius")
			{
				getline(infileinf, stemp);
				rad=atof(stemp.c_str());
				//cout<<num<<endl;
			}
		}while(!infileinf.eof());
		
	
}

//================================================
//HC feature
//================================================
void HCFeature(vector<double>& score, MOD model, vector<vector<double> > hc)
{
	double gRad=1.0,hfn=0.0,hcs=0.0; //cout<<"OK!"<<endl;
	for(int i=0; i<hc.size(); i++)
	{
		gRad = gRad*hc[i][0];
		hfn += hc[i][1];
		hcs += hc[i][2];
	} 
	gRad=pow(gRad, 1.0/double(hc.size()));
	double pface=0.0;
	for(int i=0; i<model.sse.size(); i++)
	{
		if(model.sse[i].type =='H') pface += 4.0;
		else if(model.sse[i].type =='E') pface += 2.0;
		else continue;
	}
	if(gRad <=1.0) score[score.size()-3]=0.0;
	else
	score[score.size()-3]=gRad/model.radius;
	if(pface ==0.0 || hfn == 0.0) score[score.size()-2]=0;
	else
	{
	score[score.size()-2]=hfn/pface; 
	}
	if(model.ss.size() == 0) score[score.size()-1]=0;
	else score[score.size()-1]=hcs/model.ss.size();
	
}
