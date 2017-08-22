/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Domain.h"
#include "Bar.h"
#include "Outputter.h"

#include <string>
#include <iostream>

using namespace std;

int main(int argc, char *argv[])
{
	if (argc != 2) //  Print help message
	{
	    cout << "Usage: stap++ InputFileName\n";
		exit(1);
	}

	string filename(argv[1]);
	string InFile = filename + ".dat";
	string OutFile = filename + ".out";

	CDomain* FEMData = CDomain::Instance();

	if (!FEMData->ReadData(InFile, OutFile))
	{
		cout << "*** Error *** Data input failed!" << endl;
		exit(1);
	}

	FEMData->AllocateMatrices();
	FEMData->AssembleStiffnessMatrix();
	
	CLDLTSolver* S = new CLDLTSolver(FEMData);
	S->Solve();
	
	COutputter* Output = COutputter::Instance();
	Output->OutputElementStress();

	return 0;
}