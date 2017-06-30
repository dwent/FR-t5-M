//===================================================================================
//
// Generate SSE information and index
//
// AUTHORS: daiwentao
// E-mail: daiwentao@moon.ibp.ac.cn
// Date: 2012.01.05
//              
//
//===================================================================================

#include "GenSSE.hpp"

//===================================================================
// Recognize the secondary elements of the input sequence
//===================================================================
void InitialSSE(vector<ELE>& allele, string Second)
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
		if(Second[i]=='h') Second[i]='H';
		else if(Second[i]=='e') Second[i]='E';
		else if(Second[i]=='c') Second[i]='C';
		else continue;
	}

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
		
		if(ele.type=='H' && ele.len>=4) allele.push_back(ele);
		else if(ele.type=='E' && ele.len>=4) allele.push_back(ele);
	}
}
