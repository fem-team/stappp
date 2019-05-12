/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Material.h"

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

//	Read material data from stream Input
bool CBarMaterial::Read(ifstream& Input, unsigned int mset)
{
	Input >> nset;	// Number of property set

	if (nset != mset + 1)
	{
		cerr << "*** Error *** Material sets must be inputted in order !" << endl 
			 << "    Expected set : " << mset + 1 << endl
			 << "    Provided set : " << nset << endl;

		return false;
	}

	Input >> E >> Area >> rho_3D;	// Young's modulus and section area

	return true;
}

bool CPlateMaterial::Read(ifstream& Input, unsigned int mset)
{
	Input >> nset;	// Number of property set

	if (nset != mset + 1)
	{
		cerr << "*** Error *** Material sets must be inputted in order !" << endl
			<< "    Expected set : " << mset + 1 << endl
			<< "    Provided set : " << nset << endl;

		return false;
	}

	Input >> E >> nu >> h;	// Young's modulus and section area

	return true;
}

bool CIEMMaterial::Read(ifstream& Input, unsigned int mset)
{
	Input >> nset;	// Number of property set

	Input >> nset;	// Number of property set

	if (nset != mset + 1)
	{
		cerr << "*** Error *** Material sets must be inputted in order !" << endl
			<< "    Expected set : " << mset + 1 << endl
			<< "    Provided set : " << nset << endl;

		return false;
	}

	Input >> E >> nu;	// Young's modulus and section area

	return true;
}

bool CQ4Material::Read(ifstream& Input, unsigned int mset)
{
	Input >> nset;	// Number of property set

	if (nset != mset + 1)
	{
		cerr << "*** Error *** Material sets must be inputted in order !" << endl
			<< "    Expected set : " << mset + 1 << endl
			<< "    Provided set : " << nset << endl;

		return false;
	}

	Input >> E >> Nu;	// Young's modulus and section area

	return true;
}

bool CQ5Material::Read(ifstream& Input, unsigned int mset)
{
	Input >> nset;	// Number of property set

	if (nset != mset + 1)
	{
		cerr << "*** Error *** Material sets must be inputted in order !" << endl
			<< "    Expected set : " << mset + 1 << endl
			<< "    Provided set : " << nset << endl;

		return false;
	}

	Input >> E >> nu;	// Young's modulus and section area

	return true;
}

bool CBeamMaterial::Read(ifstream& Input, unsigned int mset)
{
	Input >> nset;	// Number of property set

	if (nset != mset + 1)
	{
		cerr << "*** Error *** Material sets must be inputted in order !" << endl
			<< "    Expected set : " << mset + 1 << endl
			<< "    Provided set : " << nset << endl;

		return false;
	}

	Input >> E >> Nu >> a >> b >> x1 >> x2 >> y1 >> y2 >> n1 >> n2 >> n3;	// Young's modulus and I

	return true;
}

//	Write material data to Stream
void CBarMaterial::Write(COutputter& output, unsigned int mset)
{
	output << setw(5) << mset+1 << setw(16) << E << setw(16) << Area << endl;
}




void CQ4Material::Write(COutputter& output, unsigned int mset)
{
	output << setw(5) << mset+1 << setw(16) << E << setw(16) << Nu << endl;
}

void CQ5Material::Write(COutputter& output, unsigned int mset)
{
	output << setw(5) << mset + 1 << setw(16) << E << setw(16) << nu << endl;
}


void CBeamMaterial::Write(COutputter& output, unsigned int mset)
{
	output << setw(5) << mset+1 << setw(16) << E << setw(16) << Nu << setw(16) <<a << setw(16) <<b<< endl;
}

void CPlateMaterial::Write(COutputter& output, unsigned int mset)
{
	output << setw(5) << mset + 1 << setw(16) << E << setw(16) << nu << setw(16) << h << endl;
}

void CIEMMaterial::Write(COutputter& output, unsigned int mset)
{
	output << setw(5) << mset + 1 << setw(16) << E << setw(16) << nu << endl;
}
