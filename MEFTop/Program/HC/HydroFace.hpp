//======================================================================
//
// Description:
//	Generate hydrophobic-face on alpha-helix based on
//	 center-axis and on beta-strand based on residue-position.
//
// Contact: wuaiping@moon.ibp.ac.cn, 2010_05_11
//
//======================================================================


#ifndef HYDROFACE_H
#define HYDROFACE_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include "HelixAxis.hpp"

using namespace std;



//==================================================================================
// Hydrophobic-face recognition method of alpha-helix:
// 1. calculate the center-axis of alpha-helix; (call HelixAxis())
// 2. translate axis to each hydrophobic CA-atom, form surface-axises;
// 3. calculate the hydro-residues around every surface-axis;
// 4. select the surface-axis with maximum hydro-resi number;
// 5. rotate the surface-axis with 120 and 240 degree, get two new faces;
// 6. deredundant of these three faces, and decide if they are hydro-faces. 
//==================================================================================




typedef struct _HFace
{
	vector<int> ri;	//from 0 to ri.size()-1
} HFace;



typedef struct _ELE
{
	char type;		//'H' or 'E'
	int pos;		//the initpos of this element in sequence
	int len;		//the length of this element
	string seq;		//the sequence of this element
	string HPseq;		//the Hydro/Polar property, such as "HPHPHH..."
	vector<vector<double > > CA;		//all CA coordinates
	
	vector<HFace> hface;	//all hydrophobic faces on this element;
} ELE;



void Seq2HP(vector<char>& HPseq, vector<char> seq);

bool HydroFace(ELE& ele);

void HelixHydroFace(vector<HFace>& hface,
		     const vector<vector<double > >& CA,
		     const string HPseq);




#endif



