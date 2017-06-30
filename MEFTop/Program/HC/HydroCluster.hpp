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


#include "HydroFace.hpp"


using namespace std;


//=====================================================
typedef struct _HydroCluster
{
	vector<int> ri;
	vector<HFace> allhf;
} HClust;
//=====================================================



//=====================================================
typedef struct _ChainCoord	//to store chain's five-beads model information
{
	int resinum;
	int atomnum;

	double beginxyz[3];
	double endxyz[6];
	
	vector<vector<string > > atomname;
	vector<string> resiname;
	vector<int> resiserial;
	vector<vector<vector<double > > > xyz;

} ChainCoord;

typedef vector<vector<double > > DV2;
typedef vector<vector<vector<double > > > DV3;
typedef vector<vector<vector<vector<double > > > > DV4;


int ReadCentroid(DV4& centable, char* cenfile);
int ReadDSSP(vector<char>& resi, 
		vector<char>& second, 
		vector<double>& acc, 
		const char* secondfile);
int Read4AtomPdb(ChainCoord& chain, char* pdbfile);

//add centroid atom to 4-atom's residue (N,CA,C,O)
void AddCentroid(ChainCoord& chain, DV4& centable);
//=====================================================

void AAname3to1(char& name1, string name3);
void InitialSSE(vector<ELE>& allele, vector<char> Second, ChainCoord& chain);


void CentroidMatrix(DV2& Dcent, DV3& NVcent, ChainCoord chain, vector<int> hi);
void ClusterMap(vector<vector<int > >& map, 
		  const DV2& Dcent, 
		  const DV3& NVcent, 
		  const ChainCoord& chain, 
		  const vector<int>& hi);
void HydroCluster(vector<vector<int > >& cluster, 
		  vector<vector<int > > map, 
		  const DV2& Dcent,
		  const vector<int>& hi);

void FaceCluster(vector<vector<int > >& cluster,
		  vector<vector<int > > map, 
		  const DV2& Dcent,
		  const vector<int>& hi,
		  vector<ELE>& allele);
bool Redundant(int id, vector<int> set);
int FaceInCluster(const vector<ELE>& allele, const vector<int>& cluster);



//===========================================================
// 两个问题：
// 	一是从得到的疏水内核中识别出内核的封闭程度；
// 	二是基于开放的疏水内核，定出待组装位置；
// 对第一个问题的解决方案：
//	1. 以每个残基为侧链质心为中心，以一个距离阈值（如15A）为
//		半径建立一个球形，然后分析球形空间内八个象限内其
//		他残基侧链质心的数量分布；（识别暴露残基）
//	2. 根据内核中的暴露残基个数判定内核封闭与否（识别暴露内核）
//
// 对第二个问题的解决方案：
//	1. 以每个暴露残基为基准，在空白象限中定义出几个潜在坐标；
//	2. 将所有暴露残基的潜在坐标去冗余并聚类或计算平均轴，作为
//		候选的元件组装目标轴。
//
//===========================================================

double ClusterCentroid(vector<double>& cent, 
			const vector<int>& cluster,
			const ChainCoord& chain);
void SphereVector(vector<vector<double > >& svect,
		const vector<double>& cent, 
		const vector<int>& cluster,
		const ChainCoord& chain);
void ResiQuadrant(vector<int>& Nquad, int ri, const ChainCoord& chain);

bool IfResiBury(const vector<int>& Nquad);
bool IfFaceBury(int initpos, HFace& hface, const ChainCoord& chain);
bool IfCoreBury(vector<int>& cluster, const ChainCoord& chain);




