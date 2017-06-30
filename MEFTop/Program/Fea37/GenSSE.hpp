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


#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <list>
#include <iostream>
#include <cstring>
#include <vector>
#include <fstream>

using namespace std;

typedef struct _ElE
{
	char type;		//'H' or 'E'
	int ei;					//index of element;
	int pos;		//the initpos of this element in sequence
	int len;		//the length of this element
}ELE;

void InitialSSE(vector<ELE>& allele, string Second);
