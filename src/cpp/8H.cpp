/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "8H.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

//	Constructor
C8H::C8H()
{
	NEN = 8;
	nodes = new CNode*[NEN];
    
    ND = 24;
    LocationMatrix = new unsigned int[ND];

	ElementMaterial = nullptr;
}

//	Desconstructor
C8H::~C8H()
{
	delete [] nodes;
    delete [] LocationMatrix;
}
//	Read element data from stream Input
bool C8H::Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList)
{
	unsigned int N;
	Input >> N;	// element number
	if (N != Ele + 1)
	{
		cerr << "*** Error *** Elements must be inputted in order !" << endl 
			 << "    Expected element : " << Ele + 1 << endl
			 << "    Provided element : " << N << endl;
		return false;
	}
	unsigned int MSet;
	unsigned int N1, N2, N3 ,N4 ,N5 ,N6 ,N7 ,N8;
	Input >> N1 >> N2 >> N3 >> N4 >> N5 >> N6 >> N7 >> N8 >> MSet;
    ElementMaterial = dynamic_cast<C8HMaterial*>(MaterialSets) + MSet - 1;
	nodes[0] = &NodeList[N1 - 1];
	nodes[1] = &NodeList[N2 - 1];
	nodes[2] = &NodeList[N3 - 1];
	nodes[3] = &NodeList[N4 - 1];
	nodes[4] = &NodeList[N5 - 1];
	nodes[5] = &NodeList[N6 - 1];
	nodes[6] = &NodeList[N7 - 1];
	nodes[7] = &NodeList[N8 - 1];
	return true;
}
//	Write element data to stream
void C8H::Write(COutputter& output, unsigned int Ele)
{
  output << setw(5) << Ele+1 << setw(11) << nodes[0]->NodeNumber 
		 << setw(9) << nodes[1]->NodeNumber 
		 << setw(9) << nodes[2]->NodeNumber 
		 << setw(9) << nodes[3]->NodeNumber
		 << setw(9) << nodes[4]->NodeNumber
		 << setw(9) << nodes[5]->NodeNumber
		 << setw(9) << nodes[6]->NodeNumber
		 << setw(9) << nodes[7]->NodeNumber
		 << setw(12) << ElementMaterial->nset << endl;
}

//  Generate location matrix: the global equation number that corresponding to each DOF of the element
//	Caution:  Equation number is numbered from 1 !
void C8H::GenerateLocationMatrix()
{
    unsigned int i = 0;
    for (unsigned int N = 0; N < NEN; N++)
        for (unsigned int D = 0; D < 3; D++)
            LocationMatrix[i++] = nodes[N]->bcode[D];
}

//	Return the size of the element stiffness matrix (stored as an array column by column)
//	For Plate element, element stiffness is a 12x12 matrix, whose upper triangular part
//	has 78 elements
unsigned int C8H::SizeOfStiffnessMatrix() { return 300; }

//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void C8H::ElementStiffness(double* Matrix)
{
	clear(Matrix, SizeOfStiffnessMatrix());

	C8HMaterial* material = dynamic_cast<C8HMaterial*>(ElementMaterial);
	double E = material->E;
	double nv = material->nv;
	double mu = E/(1+nv);
	double lambda = nv*E/((1+nv)*(1-2*nv));
	double shape[2];
	shape[0] = -1/sqrt(3);
	shape[1] = 1/sqrt(3);
	for (unsigned int m = 0; m < 2; m++)//shape function.
	{
		for (unsigned int n = 0; n < 2; n++)
		{
			for (unsigned int o = 0; o < 2; o++)
			{
				double xi = shape[m];
				double eta = shape[n];
				double zeta = shape[o];
				double GN[12];
				GN[0] = 0.125*(1 - eta)*(1 - zeta);
				GN[1] = 0.125*(1 + eta)*(1 - zeta);
				GN[2] = 0.125*(1 - eta)*(1 + zeta);
				GN[3] = 0.125*(1 + eta)*(1 + zeta);
				GN[4] = 0.125*(1 - xi )*(1 - zeta);
				GN[5] = 0.125*(1 + xi )*(1 - zeta);
				GN[6] = 0.125*(1 - xi )*(1 + zeta);
				GN[7] = 0.125*(1 + xi )*(1 + zeta);
				GN[8] = 0.125*(1 - xi )*(1 - eta);
				GN[9] = 0.125*(1 + xi )*(1 - eta);
				GN[10] = 0.125*(1 + xi)*(1 + eta);
				GN[11] = 0.125*(1 - xi)*(1 + eta);

				double J[9];
				J[0] = GN[0]*(nodes[1]->XYZ[0] - nodes[0]->XYZ[0]) + GN[1]*(nodes[2]->XYZ[0] - nodes[3]->XYZ[0]) + GN[2]*(nodes[5]->XYZ[0] - nodes[4]->XYZ[0]) + GN[3]*(nodes[6]->XYZ[0] - nodes[7]->XYZ[0]);
				J[1] = GN[0]*(nodes[1]->XYZ[1] - nodes[0]->XYZ[1]) + GN[1]*(nodes[2]->XYZ[1] - nodes[3]->XYZ[1]) + GN[2]*(nodes[5]->XYZ[1] - nodes[4]->XYZ[1]) + GN[3]*(nodes[6]->XYZ[1] - nodes[7]->XYZ[1]);
				J[2] = GN[0]*(nodes[1]->XYZ[2] - nodes[0]->XYZ[2]) + GN[1]*(nodes[2]->XYZ[2] - nodes[3]->XYZ[2]) + GN[2]*(nodes[5]->XYZ[2] - nodes[4]->XYZ[2]) + GN[3]*(nodes[6]->XYZ[2] - nodes[7]->XYZ[2]);
				J[3] = GN[4]*(nodes[3]->XYZ[0] - nodes[0]->XYZ[0]) + GN[5]*(nodes[2]->XYZ[0] - nodes[1]->XYZ[0]) + GN[6]*(nodes[7]->XYZ[0] - nodes[4]->XYZ[0]) + GN[7]*(nodes[6]->XYZ[0] - nodes[5]->XYZ[0]);
				J[4] = GN[4]*(nodes[3]->XYZ[1] - nodes[0]->XYZ[1]) + GN[5]*(nodes[2]->XYZ[1] - nodes[1]->XYZ[1]) + GN[6]*(nodes[7]->XYZ[1] - nodes[4]->XYZ[1]) + GN[7]*(nodes[6]->XYZ[1] - nodes[5]->XYZ[1]);
				J[5] = GN[4]*(nodes[2]->XYZ[1] - nodes[0]->XYZ[2]) + GN[5]*(nodes[2]->XYZ[2] - nodes[1]->XYZ[2]) + GN[6]*(nodes[7]->XYZ[2] - nodes[4]->XYZ[2]) + GN[7]*(nodes[6]->XYZ[2] - nodes[5]->XYZ[2]);
				J[6] = GN[8]*(nodes[2]->XYZ[0] - nodes[0]->XYZ[0]) + GN[9]*(nodes[5]->XYZ[0] - nodes[1]->XYZ[0]) + GN[10]*(nodes[6]->XYZ[0] - nodes[2]->XYZ[0]) + GN[11]*(nodes[7]->XYZ[0] - nodes[3]->XYZ[0]);
				J[7] = GN[8]*(nodes[4]->XYZ[1] - nodes[0]->XYZ[1]) + GN[9]*(nodes[5]->XYZ[1] - nodes[1]->XYZ[1]) + GN[10]*(nodes[6]->XYZ[1] - nodes[2]->XYZ[1]) + GN[11]*(nodes[7]->XYZ[1] - nodes[3]->XYZ[1]);
				J[8] = GN[8]*(nodes[4]->XYZ[2] - nodes[0]->XYZ[2]) + GN[9]*(nodes[5]->XYZ[2] - nodes[1]->XYZ[2]) + GN[10]*(nodes[6]->XYZ[2] - nodes[2]->XYZ[2]) + GN[11]*(nodes[7]->XYZ[2] - nodes[3]->XYZ[2]);
				double DetJ1 = J[0]*J[4]*J[8] - J[0]*J[5]*J[7] - J[1]*J[3]*J[8] + J[1]*J[5]*J[6] + J[2]*J[3]*J[7] - J[2]*J[4]*J[6];
				double DetJ = abs(DetJ1);

				double InvJ[9];
				InvJ[0] = (J[4]*J[8]-J[5]*J[7])/DetJ1;
			    InvJ[1] = -(J[1]*J[8]-J[2]*J[7])/DetJ1;
			    InvJ[2] = (J[1]*J[5]-J[2]*J[4])/DetJ1;
			    InvJ[3] = -(J[3]*J[8]-J[5]*J[6])/DetJ1;
			    InvJ[4] = (J[0]*J[8]-J[2]*J[6])/DetJ1;
			    InvJ[5] = -(J[0]*J[5]-J[2]*J[3])/DetJ1;
			    InvJ[6] = (J[3]*J[7]-J[4]*J[6])/DetJ1;
			    InvJ[7] = -(J[0]*J[7]-J[1]*J[6])/DetJ1;
			    InvJ[8] = (J[0]*J[4]-J[1]*J[3])/DetJ1;

				double B1[24];
				B1[0] = -GN[0]*InvJ[0]-GN[4]*InvJ[1]-GN[8]*InvJ[2];
			    B1[1] = GN[0]*InvJ[0]-GN[5]*InvJ[1]-GN[9]*InvJ[2];
			    B1[2] = GN[1]*InvJ[0]+GN[5]*InvJ[1]-GN[10]*InvJ[2];
			    B1[3] = -GN[1]*InvJ[0]+GN[4]*InvJ[1]-GN[11]*InvJ[2];
				B1[4] = -GN[2]*InvJ[0]-GN[6]*InvJ[1]+GN[8]*InvJ[2];
				B1[5] = GN[2]*InvJ[0]-GN[7]*InvJ[1]+GN[9]*InvJ[2];
				B1[6] = GN[3]*InvJ[0]+GN[7]*InvJ[1]+GN[10]*InvJ[2];
				B1[7] = -GN[3]*InvJ[0]+GN[6]*InvJ[1]+GN[11]*InvJ[2];
			    B1[8] = -GN[0]*InvJ[3]-GN[4]*InvJ[4]-GN[8]*InvJ[5];
			    B1[9] = GN[0]*InvJ[3]-GN[5]*InvJ[4]-GN[9]*InvJ[5];
			    B1[10] = GN[1]*InvJ[3]+GN[5]*InvJ[4]-GN[10]*InvJ[5];
			    B1[11] = -GN[1]*InvJ[3]+GN[4]*InvJ[4]-GN[11]*InvJ[5];
				B1[12] = -GN[2]*InvJ[3]-GN[6]*InvJ[4]+GN[8]*InvJ[5];
				B1[13] = GN[2]*InvJ[3]-GN[7]*InvJ[4]+GN[9]*InvJ[5];
				B1[14] = GN[3]*InvJ[3]+GN[7]*InvJ[4]+GN[10]*InvJ[5];
				B1[15] = -GN[3]*InvJ[3]+GN[6]*InvJ[4]+GN[11]*InvJ[5];
				B1[16] = -GN[0]*InvJ[6]-GN[4]*InvJ[7]-GN[8]*InvJ[8];
			    B1[17] = GN[0]*InvJ[6]-GN[5]*InvJ[7]-GN[9]*InvJ[8];
			    B1[18] = GN[1]*InvJ[6]+GN[5]*InvJ[7]-GN[10]*InvJ[8];
			    B1[19] = -GN[1]*InvJ[6]+GN[4]*InvJ[7]-GN[11]*InvJ[8];
				B1[20] = -GN[2]*InvJ[6]-GN[6]*InvJ[7]+GN[8]*InvJ[8];
				B1[21] = GN[2]*InvJ[6]-GN[7]*InvJ[7]+GN[9]*InvJ[8];
				B1[22] = GN[3]*InvJ[6]+GN[7]*InvJ[7]+GN[10]*InvJ[8];
				B1[23] = -GN[3]*InvJ[6]+GN[6]*InvJ[7]+GN[11]*InvJ[8];


			       Matrix[0] += (B1[0]*B1[0]*(mu+lambda)+B1[8]*mu*B1[8]+B1[16]*mu*B1[16])*DetJ;
				   Matrix[1] += (B1[8]*B1[8]*(mu+lambda)+B1[0]*mu*B1[0]+B1[16]*mu*B1[16])*DetJ;
				   Matrix[2] += (B1[0]*mu*B1[8]+B1[8]*lambda*B1[0])*DetJ;
				   Matrix[3] += (B1[16]*B1[16]*(mu+lambda)+B1[0]*mu*B1[0]+B1[8]*mu*B1[8])*DetJ;
				   Matrix[4] += (B1[8]*mu*B1[16]+B1[16]*lambda*B1[8])*DetJ;
				   Matrix[5] += (B1[0]*mu*B1[16]+B1[16]*lambda*B1[0])*DetJ;
				   Matrix[6] += (B1[1]*B1[1]*(mu+lambda)+B1[9]*mu*B1[9]+B1[17]*mu*B1[17])*DetJ;
				   Matrix[7] += (B1[17]*mu*B1[0]+B1[1]*lambda*B1[16])*DetJ;
				   Matrix[8] += (B1[9]*mu*B1[0]+B1[1]*lambda*B1[8])*DetJ;
				   Matrix[9] += (B1[1]*B1[0]*(mu+lambda)+B1[9]*mu*B1[8]+B1[17]*mu*B1[16])*DetJ;
				   Matrix[10] += (B1[9]*B1[9]*(mu+lambda)+B1[1]*mu*B1[1]+B1[17]*mu*B1[17])*DetJ;
				   Matrix[11] += (B1[1]*mu*B1[9]+B1[9]*lambda*B1[1])*DetJ;
				   Matrix[12] += (B1[17]*mu*B1[8]+B1[9]*lambda*B1[16])*DetJ;
				   Matrix[13] += (B1[9]*B1[8]*(mu+lambda)+B1[1]*mu*B1[0]+B1[17]*mu*B1[16])*DetJ;
				   Matrix[14] += (B1[1]*mu*B1[8]+B1[9]*lambda*B1[0])*DetJ;
				   Matrix[15] += (B1[17]*B1[17]*(mu+lambda)+B1[1]*mu*B1[1]+B1[9]*mu*B1[9])*DetJ;
				   Matrix[16] += (B1[9]*mu*B1[17]+B1[17]*lambda*B1[9])*DetJ;
				   Matrix[17] += (B1[1]*mu*B1[17]+B1[17]*lambda*B1[1])*DetJ;
				   Matrix[18] += (B1[17]*B1[16]*(mu+lambda)+B1[1]*mu*B1[0]+B1[9]*mu*B1[8])*DetJ;
				   Matrix[19] += (B1[9]*mu*B1[16]+B1[17]*lambda*B1[8])*DetJ;
				   Matrix[20] += (B1[1]*mu*B1[16]+B1[17]*lambda*B1[0])*DetJ;
				   Matrix[21] += (B1[2]*B1[2]*(mu+lambda)+B1[10]*mu*B1[10]+B1[18]*mu*B1[18])*DetJ;
				   Matrix[22] += (B1[18]*mu*B1[1]+B1[2]*lambda*B1[17])*DetJ;
				   Matrix[23] += (B1[10]*mu*B1[1]+B1[2]*lambda*B1[9])*DetJ;
				   Matrix[24] += (B1[2]*B1[1]*(mu+lambda)+B1[10]*mu*B1[9]+B1[18]*mu*B1[17])*DetJ;
				   Matrix[25] += (B1[18]*mu*B1[0]+B1[2]*lambda*B1[16])*DetJ;
				   Matrix[26] += (B1[10]*mu*B1[0]+B1[2]*lambda*B1[8])*DetJ;
				   Matrix[27] += (B1[2]*B1[0]*(mu+lambda)+B1[10]*mu*B1[8]+B1[18]*mu*B1[16])*DetJ;
				   Matrix[28] += (B1[10]*B1[10]*(mu+lambda)+B1[2]*mu*B1[2]+B1[18]*mu*B1[18])*DetJ;
				   Matrix[29] += (B1[2]*mu*B1[10]+B1[10]*lambda*B1[2])*DetJ;
				   Matrix[30] += (B1[18]*mu*B1[9]+B1[10]*lambda*B1[17])*DetJ;
				   Matrix[31] += (B1[10]*B1[9]*(mu+lambda)+B1[2]*mu*B1[1]+B1[18]*mu*B1[17])*DetJ;
				   Matrix[32] += (B1[2]*mu*B1[9]+B1[10]*lambda*B1[1])*DetJ;
				   Matrix[33] += (B1[18]*mu*B1[8]+B1[10]*lambda*B1[16])*DetJ;
				   Matrix[34] += (B1[10]*B1[8]*(mu+lambda)+B1[2]*mu*B1[0]+B1[18]*mu*B1[16])*DetJ;
				   Matrix[35] += (B1[2]*mu*B1[8]+B1[10]*lambda*B1[0])*DetJ;
				   Matrix[36] += (B1[18]*B1[18]*(mu+lambda)+B1[2]*mu*B1[2]+B1[10]*mu*B1[10])*DetJ;
				   Matrix[37] += (B1[10]*mu*B1[18]+B1[18]*lambda*B1[10])*DetJ;
				   Matrix[38] += (B1[2]*mu*B1[18]+B1[18]*lambda*B1[2])*DetJ;
				   Matrix[39] += (B1[18]*B1[17]*(mu+lambda)+B1[2]*mu*B1[1]+B1[10]*mu*B1[9])*DetJ;
				   Matrix[40] += (B1[10]*mu*B1[17]+B1[18]*lambda*B1[9])*DetJ;
				   Matrix[41] += (B1[2]*mu*B1[17]+B1[18]*lambda*B1[1])*DetJ;
				   Matrix[42] += (B1[18]*B1[16]*(mu+lambda)+B1[2]*mu*B1[0]+B1[10]*mu*B1[8])*DetJ;
				   Matrix[43] += (B1[10]*mu*B1[16]+B1[18]*lambda*B1[8])*DetJ;
				   Matrix[44] += (B1[2]*mu*B1[16]+B1[18]*lambda*B1[0])*DetJ;
				   Matrix[45] += (B1[3]*B1[3]*(mu+lambda)+B1[11]*mu*B1[11]+B1[19]*mu*B1[19])*DetJ;
				   Matrix[46] += (B1[19]*mu*B1[2]+B1[3]*lambda*B1[18])*DetJ;
				   Matrix[47] += (B1[11]*mu*B1[2]+B1[3]*lambda*B1[10])*DetJ;
				   Matrix[48] += (B1[3]*B1[2]*(mu+lambda)+B1[11]*mu*B1[10]+B1[19]*mu*B1[18])*DetJ;
				   Matrix[49] += (B1[19]*mu*B1[1]+B1[3]*lambda*B1[17])*DetJ;
				   Matrix[50] += (B1[11]*mu*B1[1]+B1[3]*lambda*B1[9])*DetJ;
				   Matrix[51] += (B1[3]*B1[1]*(mu+lambda)+B1[11]*mu*B1[9]+B1[19]*mu*B1[17])*DetJ;
				   Matrix[52] += (B1[19]*mu*B1[0]+B1[3]*lambda*B1[16])*DetJ;
				   Matrix[53] += (B1[11]*mu*B1[0]+B1[3]*lambda*B1[8])*DetJ;
				   Matrix[54] += (B1[3]*B1[0]*(mu+lambda)+B1[11]*mu*B1[8]+B1[19]*mu*B1[16])*DetJ;
				   Matrix[55] += (B1[11]*B1[11]*(mu+lambda)+B1[3]*mu*B1[3]+B1[19]*mu*B1[19])*DetJ;
				   Matrix[56] += (B1[3]*mu*B1[11]+B1[11]*lambda*B1[3])*DetJ;
				   Matrix[57] += (B1[19]*mu*B1[10]+B1[11]*lambda*B1[18])*DetJ;
				   Matrix[58] += (B1[11]*B1[10]*(mu+lambda)+B1[3]*mu*B1[2]+B1[19]*mu*B1[18])*DetJ;
				   Matrix[59] += (B1[3]*mu*B1[10]+B1[11]*lambda*B1[2])*DetJ;
				   Matrix[60] += (B1[19]*mu*B1[9]+B1[11]*lambda*B1[17])*DetJ;
				   Matrix[61] += (B1[11]*B1[9]*(mu+lambda)+B1[3]*mu*B1[1]+B1[19]*mu*B1[17])*DetJ;
				   Matrix[62] += (B1[3]*mu*B1[9]+B1[11]*lambda*B1[1])*DetJ;
				   Matrix[63] += (B1[19]*mu*B1[8]+B1[11]*lambda*B1[16])*DetJ;
				   Matrix[64] += (B1[11]*B1[8]*(mu+lambda)+B1[3]*mu*B1[0]+B1[19]*mu*B1[16])*DetJ;
				   Matrix[65] += (B1[3]*mu*B1[8]+B1[11]*lambda*B1[0])*DetJ;
				   Matrix[66] += (B1[19]*B1[19]*(mu+lambda)+B1[3]*mu*B1[3]+B1[11]*mu*B1[11])*DetJ;
				   Matrix[67] += (B1[11]*mu*B1[19]+B1[19]*lambda*B1[11])*DetJ;
				   Matrix[68] += (B1[3]*mu*B1[19]+B1[19]*lambda*B1[3])*DetJ;
				   Matrix[69] += (B1[19]*B1[18]*(mu+lambda)+B1[3]*mu*B1[2]+B1[11]*mu*B1[10])*DetJ;
				   Matrix[70] += (B1[11]*mu*B1[18]+B1[19]*lambda*B1[10])*DetJ;
				   Matrix[71] += (B1[3]*mu*B1[18]+B1[19]*lambda*B1[2])*DetJ;
				   Matrix[72] += (B1[19]*B1[17]*(mu+lambda)+B1[3]*mu*B1[1]+B1[11]*mu*B1[9])*DetJ;
				   Matrix[73] += (B1[11]*mu*B1[17]+B1[19]*lambda*B1[9])*DetJ;
				   Matrix[74] += (B1[3]*mu*B1[17]+B1[19]*lambda*B1[1])*DetJ;
				   Matrix[75] += (B1[19]*B1[16]*(mu+lambda)+B1[3]*mu*B1[0]+B1[11]*mu*B1[8])*DetJ;
				   Matrix[76] += (B1[11]*mu*B1[16]+B1[19]*lambda*B1[8])*DetJ;
				   Matrix[77] += (B1[3]*mu*B1[16]+B1[19]*lambda*B1[0])*DetJ;
				   Matrix[78] += (B1[4]*B1[4]*(mu+lambda)+B1[12]*mu*B1[12]+B1[20]*mu*B1[20])*DetJ;
				   Matrix[79] += (B1[20]*mu*B1[3]+B1[4]*lambda*B1[19])*DetJ;
				   Matrix[80] += (B1[12]*mu*B1[3]+B1[4]*lambda*B1[11])*DetJ;
				   Matrix[81] += (B1[4]*B1[3]*(mu+lambda)+B1[12]*mu*B1[11]+B1[20]*mu*B1[19])*DetJ;
				   Matrix[82] += (B1[20]*mu*B1[2]+B1[4]*lambda*B1[18])*DetJ;
				   Matrix[83] += (B1[12]*mu*B1[2]+B1[4]*lambda*B1[10])*DetJ;
				   Matrix[84] += (B1[4]*B1[2]*(mu+lambda)+B1[12]*mu*B1[10]+B1[20]*mu*B1[18])*DetJ;
				   Matrix[85] += (B1[20]*mu*B1[1]+B1[4]*lambda*B1[17])*DetJ;
				   Matrix[86] += (B1[12]*mu*B1[1]+B1[4]*lambda*B1[9])*DetJ;
				   Matrix[87] += (B1[4]*B1[1]*(mu+lambda)+B1[12]*mu*B1[9]+B1[20]*mu*B1[17])*DetJ;
				   Matrix[88] += (B1[20]*mu*B1[0]+B1[4]*lambda*B1[16])*DetJ;
				   Matrix[89] += (B1[12]*mu*B1[0]+B1[4]*lambda*B1[8])*DetJ;
				   Matrix[90] += (B1[4]*B1[0]*(mu+lambda)+B1[12]*mu*B1[8]+B1[20]*mu*B1[16])*DetJ;
				   Matrix[91] += (B1[12]*B1[12]*(mu+lambda)+B1[4]*mu*B1[4]+B1[20]*mu*B1[20])*DetJ;
				   Matrix[92] += (B1[4]*mu*B1[12]+B1[12]*lambda*B1[4])*DetJ;
				   Matrix[93] += (B1[20]*mu*B1[11]+B1[12]*lambda*B1[19])*DetJ;
				   Matrix[94] += (B1[12]*B1[11]*(mu+lambda)+B1[4]*mu*B1[3]+B1[20]*mu*B1[19])*DetJ;
				   Matrix[95] += (B1[4]*mu*B1[11]+B1[12]*lambda*B1[3])*DetJ;
				   Matrix[96] += (B1[20]*mu*B1[10]+B1[12]*lambda*B1[18])*DetJ;
				   Matrix[97] += (B1[12]*B1[10]*(mu+lambda)+B1[4]*mu*B1[2]+B1[20]*mu*B1[18])*DetJ;
				   Matrix[98] += (B1[4]*mu*B1[10]+B1[12]*lambda*B1[2])*DetJ;
				   Matrix[99] += (B1[20]*mu*B1[9]+B1[12]*lambda*B1[17])*DetJ;
				   Matrix[100] += (B1[12]*B1[9]*(mu+lambda)+B1[4]*mu*B1[1]+B1[20]*mu*B1[17])*DetJ;
				   Matrix[101] += (B1[4]*mu*B1[9]+B1[12]*lambda*B1[1])*DetJ;
				   Matrix[102] += (B1[20]*mu*B1[8]+B1[12]*lambda*B1[16])*DetJ;
				   Matrix[103] += (B1[12]*B1[8]*(mu+lambda)+B1[4]*mu*B1[0]+B1[20]*mu*B1[16])*DetJ;
				   Matrix[104] += (B1[4]*mu*B1[8]+B1[12]*lambda*B1[0])*DetJ;
				   Matrix[105] += (B1[20]*B1[20]*(mu+lambda)+B1[4]*mu*B1[4]+B1[12]*mu*B1[12])*DetJ;
				   Matrix[106] += (B1[12]*mu*B1[20]+B1[20]*lambda*B1[12])*DetJ;
				   Matrix[107] += (B1[4]*mu*B1[20]+B1[20]*lambda*B1[4])*DetJ;
				   Matrix[108] += (B1[20]*B1[19]*(mu+lambda)+B1[4]*mu*B1[3]+B1[12]*mu*B1[11])*DetJ;
				   Matrix[109] += (B1[12]*mu*B1[19]+B1[20]*lambda*B1[11])*DetJ;
				   Matrix[110] += (B1[4]*mu*B1[19]+B1[20]*lambda*B1[3])*DetJ;
				   Matrix[111] += (B1[20]*B1[18]*(mu+lambda)+B1[4]*mu*B1[2]+B1[12]*mu*B1[10])*DetJ;
				   Matrix[112] += (B1[12]*mu*B1[18]+B1[20]*lambda*B1[10])*DetJ;
				   Matrix[113] += (B1[4]*mu*B1[18]+B1[20]*lambda*B1[2])*DetJ;
				   Matrix[114] += (B1[20]*B1[17]*(mu+lambda)+B1[4]*mu*B1[1]+B1[12]*mu*B1[9])*DetJ;
				   Matrix[115] += (B1[12]*mu*B1[17]+B1[20]*lambda*B1[9])*DetJ;
				   Matrix[116] += (B1[4]*mu*B1[17]+B1[20]*lambda*B1[1])*DetJ;
				   Matrix[117] += (B1[20]*B1[16]*(mu+lambda)+B1[4]*mu*B1[0]+B1[12]*mu*B1[8])*DetJ;
				   Matrix[118] += (B1[12]*mu*B1[16]+B1[20]*lambda*B1[8])*DetJ;
				   Matrix[119] += (B1[4]*mu*B1[16]+B1[20]*lambda*B1[0])*DetJ;
				   Matrix[120] += (B1[5]*B1[5]*(mu+lambda)+B1[13]*mu*B1[13]+B1[21]*mu*B1[21])*DetJ;
				   Matrix[121] += (B1[21]*mu*B1[4]+B1[5]*lambda*B1[20])*DetJ;
				   Matrix[122] += (B1[13]*mu*B1[4]+B1[5]*lambda*B1[12])*DetJ;
				   Matrix[123] += (B1[5]*B1[4]*(mu+lambda)+B1[13]*mu*B1[12]+B1[21]*mu*B1[20])*DetJ;
				   Matrix[124] += (B1[21]*mu*B1[3]+B1[5]*lambda*B1[19])*DetJ;
				   Matrix[125] += (B1[13]*mu*B1[3]+B1[5]*lambda*B1[11])*DetJ;
				   Matrix[126] += (B1[5]*B1[3]*(mu+lambda)+B1[13]*mu*B1[11]+B1[21]*mu*B1[19])*DetJ;
				   Matrix[127] += (B1[21]*mu*B1[2]+B1[5]*lambda*B1[18])*DetJ;
				   Matrix[128] += (B1[13]*mu*B1[2]+B1[5]*lambda*B1[10])*DetJ;
				   Matrix[129] += (B1[5]*B1[2]*(mu+lambda)+B1[13]*mu*B1[10]+B1[21]*mu*B1[18])*DetJ;
				   Matrix[130] += (B1[21]*mu*B1[1]+B1[5]*lambda*B1[17])*DetJ;
				   Matrix[131] += (B1[13]*mu*B1[1]+B1[5]*lambda*B1[9])*DetJ;
				   Matrix[132] += (B1[5]*B1[1]*(mu+lambda)+B1[13]*mu*B1[9]+B1[21]*mu*B1[17])*DetJ;
				   Matrix[133] += (B1[21]*mu*B1[0]+B1[5]*lambda*B1[16])*DetJ;
				   Matrix[134] += (B1[13]*mu*B1[0]+B1[5]*lambda*B1[8])*DetJ;
				   Matrix[135] += (B1[5]*B1[0]*(mu+lambda)+B1[13]*mu*B1[8]+B1[21]*mu*B1[16])*DetJ;
				   Matrix[136] += (B1[13]*B1[13]*(mu+lambda)+B1[5]*mu*B1[5]+B1[21]*mu*B1[21])*DetJ;
				   Matrix[137] += (B1[5]*mu*B1[13]+B1[13]*lambda*B1[5])*DetJ;
				   Matrix[138] += (B1[21]*mu*B1[12]+B1[13]*lambda*B1[20])*DetJ;
				   Matrix[139] += (B1[13]*B1[12]*(mu+lambda)+B1[5]*mu*B1[4]+B1[21]*mu*B1[20])*DetJ;
				   Matrix[140] += (B1[5]*mu*B1[12]+B1[13]*lambda*B1[4])*DetJ;
				   Matrix[141] += (B1[21]*mu*B1[11]+B1[13]*lambda*B1[19])*DetJ;
				   Matrix[142] += (B1[13]*B1[11]*(mu+lambda)+B1[5]*mu*B1[3]+B1[21]*mu*B1[19])*DetJ;
				   Matrix[143] += (B1[5]*mu*B1[11]+B1[13]*lambda*B1[3])*DetJ;
				   Matrix[144] += (B1[21]*mu*B1[10]+B1[13]*lambda*B1[18])*DetJ;
				   Matrix[145] += (B1[13]*B1[10]*(mu+lambda)+B1[5]*mu*B1[2]+B1[21]*mu*B1[18])*DetJ;
				   Matrix[146] += (B1[5]*mu*B1[10]+B1[13]*lambda*B1[2])*DetJ;
				   Matrix[147] += (B1[21]*mu*B1[9]+B1[13]*lambda*B1[17])*DetJ;
				   Matrix[148] += (B1[13]*B1[9]*(mu+lambda)+B1[5]*mu*B1[1]+B1[21]*mu*B1[17])*DetJ;
				   Matrix[149] += (B1[5]*mu*B1[9]+B1[13]*lambda*B1[1])*DetJ;
				   Matrix[150] += (B1[21]*mu*B1[8]+B1[13]*lambda*B1[16])*DetJ;
				   Matrix[151] += (B1[13]*B1[8]*(mu+lambda)+B1[5]*mu*B1[0]+B1[21]*mu*B1[16])*DetJ;
				   Matrix[152] += (B1[5]*mu*B1[8]+B1[13]*lambda*B1[0])*DetJ;
				   Matrix[153] += (B1[21]*B1[21]*(mu+lambda)+B1[5]*mu*B1[5]+B1[13]*mu*B1[13])*DetJ;
				   Matrix[154] += (B1[13]*mu*B1[21]+B1[21]*lambda*B1[13])*DetJ;
				   Matrix[155] += (B1[5]*mu*B1[21]+B1[21]*lambda*B1[5])*DetJ;
				   Matrix[156] += (B1[21]*B1[20]*(mu+lambda)+B1[5]*mu*B1[4]+B1[13]*mu*B1[12])*DetJ;
				   Matrix[157] += (B1[13]*mu*B1[20]+B1[21]*lambda*B1[12])*DetJ;
				   Matrix[158] += (B1[5]*mu*B1[20]+B1[21]*lambda*B1[4])*DetJ;
				   Matrix[159] += (B1[21]*B1[19]*(mu+lambda)+B1[5]*mu*B1[3]+B1[13]*mu*B1[11])*DetJ;
				   Matrix[160] += (B1[13]*mu*B1[19]+B1[21]*lambda*B1[11])*DetJ;
				   Matrix[161] += (B1[5]*mu*B1[19]+B1[21]*lambda*B1[3])*DetJ;
				   Matrix[162] += (B1[21]*B1[18]*(mu+lambda)+B1[5]*mu*B1[2]+B1[13]*mu*B1[10])*DetJ;
				   Matrix[163] += (B1[13]*mu*B1[18]+B1[21]*lambda*B1[10])*DetJ;
				   Matrix[164] += (B1[5]*mu*B1[18]+B1[21]*lambda*B1[2])*DetJ;
				   Matrix[165] += (B1[21]*B1[17]*(mu+lambda)+B1[5]*mu*B1[1]+B1[13]*mu*B1[9])*DetJ;
				   Matrix[166] += (B1[13]*mu*B1[17]+B1[21]*lambda*B1[9])*DetJ;
				   Matrix[167] += (B1[5]*mu*B1[17]+B1[21]*lambda*B1[1])*DetJ;
				   Matrix[168] += (B1[21]*B1[16]*(mu+lambda)+B1[5]*mu*B1[0]+B1[13]*mu*B1[8])*DetJ;
				   Matrix[169] += (B1[13]*mu*B1[16]+B1[21]*lambda*B1[8])*DetJ;
				   Matrix[170] += (B1[5]*mu*B1[16]+B1[21]*lambda*B1[0])*DetJ;
				   Matrix[171] += (B1[6]*B1[6]*(mu+lambda)+B1[14]*mu*B1[14]+B1[22]*mu*B1[22])*DetJ;
				   Matrix[172] += (B1[22]*mu*B1[5]+B1[6]*lambda*B1[21])*DetJ;
				   Matrix[173] += (B1[14]*mu*B1[5]+B1[6]*lambda*B1[13])*DetJ;
				   Matrix[174] += (B1[6]*B1[5]*(mu+lambda)+B1[14]*mu*B1[13]+B1[22]*mu*B1[21])*DetJ;
				   Matrix[175] += (B1[22]*mu*B1[4]+B1[6]*lambda*B1[20])*DetJ;
				   Matrix[176] += (B1[14]*mu*B1[4]+B1[6]*lambda*B1[12])*DetJ;
				   Matrix[177] += (B1[6]*B1[4]*(mu+lambda)+B1[14]*mu*B1[12]+B1[22]*mu*B1[20])*DetJ;
				   Matrix[178] += (B1[22]*mu*B1[3]+B1[6]*lambda*B1[19])*DetJ;
				   Matrix[179] += (B1[14]*mu*B1[3]+B1[6]*lambda*B1[11])*DetJ;
				   Matrix[180] += (B1[6]*B1[3]*(mu+lambda)+B1[14]*mu*B1[11]+B1[22]*mu*B1[19])*DetJ;
				   Matrix[181] += (B1[22]*mu*B1[2]+B1[6]*lambda*B1[18])*DetJ;
				   Matrix[182] += (B1[14]*mu*B1[2]+B1[6]*lambda*B1[10])*DetJ;
				   Matrix[183] += (B1[6]*B1[2]*(mu+lambda)+B1[14]*mu*B1[10]+B1[22]*mu*B1[18])*DetJ;
				   Matrix[184] += (B1[22]*mu*B1[1]+B1[6]*lambda*B1[17])*DetJ;
				   Matrix[185] += (B1[14]*mu*B1[1]+B1[6]*lambda*B1[9])*DetJ;
				   Matrix[186] += (B1[6]*B1[1]*(mu+lambda)+B1[14]*mu*B1[9]+B1[22]*mu*B1[17])*DetJ;
				   Matrix[187] += (B1[22]*mu*B1[0]+B1[6]*lambda*B1[16])*DetJ;
				   Matrix[188] += (B1[14]*mu*B1[0]+B1[6]*lambda*B1[8])*DetJ;
				   Matrix[189] += (B1[6]*B1[0]*(mu+lambda)+B1[14]*mu*B1[8]+B1[22]*mu*B1[16])*DetJ;
				   Matrix[190] += (B1[14]*B1[14]*(mu+lambda)+B1[6]*mu*B1[6]+B1[22]*mu*B1[22])*DetJ;
				   Matrix[191] += (B1[6]*mu*B1[14]+B1[14]*lambda*B1[6])*DetJ;
				   Matrix[192] += (B1[22]*mu*B1[13]+B1[14]*lambda*B1[21])*DetJ;
				   Matrix[193] += (B1[14]*B1[13]*(mu+lambda)+B1[6]*mu*B1[5]+B1[22]*mu*B1[21])*DetJ;
				   Matrix[194] += (B1[6]*mu*B1[13]+B1[14]*lambda*B1[5])*DetJ;
				   Matrix[195] += (B1[22]*mu*B1[12]+B1[14]*lambda*B1[20])*DetJ;
				   Matrix[196] += (B1[14]*B1[12]*(mu+lambda)+B1[6]*mu*B1[4]+B1[22]*mu*B1[20])*DetJ;
				   Matrix[197] += (B1[6]*mu*B1[12]+B1[14]*lambda*B1[4])*DetJ;
				   Matrix[198] += (B1[22]*mu*B1[11]+B1[14]*lambda*B1[19])*DetJ;
				   Matrix[199] += (B1[14]*B1[11]*(mu+lambda)+B1[6]*mu*B1[3]+B1[22]*mu*B1[19])*DetJ;
				   Matrix[200] += (B1[6]*mu*B1[11]+B1[14]*lambda*B1[3])*DetJ;
				   Matrix[201] += (B1[22]*mu*B1[10]+B1[14]*lambda*B1[18])*DetJ;
				   Matrix[202] += (B1[14]*B1[10]*(mu+lambda)+B1[6]*mu*B1[2]+B1[22]*mu*B1[18])*DetJ;
				   Matrix[203] += (B1[6]*mu*B1[10]+B1[14]*lambda*B1[2])*DetJ;
				   Matrix[204] += (B1[22]*mu*B1[9]+B1[14]*lambda*B1[17])*DetJ;
				   Matrix[205] += (B1[14]*B1[9]*(mu+lambda)+B1[6]*mu*B1[1]+B1[22]*mu*B1[17])*DetJ;
				   Matrix[206] += (B1[6]*mu*B1[9]+B1[14]*lambda*B1[1])*DetJ;
				   Matrix[207] += (B1[22]*mu*B1[8]+B1[14]*lambda*B1[16])*DetJ;
				   Matrix[208] += (B1[14]*B1[8]*(mu+lambda)+B1[6]*mu*B1[0]+B1[22]*mu*B1[16])*DetJ;
				   Matrix[209] += (B1[6]*mu*B1[8]+B1[14]*lambda*B1[0])*DetJ;
				   Matrix[210] += (B1[22]*B1[22]*(mu+lambda)+B1[6]*mu*B1[6]+B1[14]*mu*B1[14])*DetJ;
				   Matrix[211] += (B1[14]*mu*B1[22]+B1[22]*lambda*B1[14])*DetJ;
				   Matrix[212] += (B1[6]*mu*B1[22]+B1[22]*lambda*B1[6])*DetJ;
				   Matrix[213] += (B1[22]*B1[21]*(mu+lambda)+B1[6]*mu*B1[5]+B1[14]*mu*B1[13])*DetJ;
				   Matrix[214] += (B1[14]*mu*B1[21]+B1[22]*lambda*B1[13])*DetJ;
				   Matrix[215] += (B1[6]*mu*B1[21]+B1[22]*lambda*B1[5])*DetJ;
				   Matrix[216] += (B1[22]*B1[20]*(mu+lambda)+B1[6]*mu*B1[4]+B1[14]*mu*B1[12])*DetJ;
				   Matrix[217] += (B1[14]*mu*B1[20]+B1[22]*lambda*B1[12])*DetJ;
				   Matrix[218] += (B1[6]*mu*B1[20]+B1[22]*lambda*B1[4])*DetJ;
				   Matrix[219] += (B1[22]*B1[19]*(mu+lambda)+B1[6]*mu*B1[3]+B1[14]*mu*B1[11])*DetJ;
				   Matrix[220] += (B1[14]*mu*B1[19]+B1[22]*lambda*B1[11])*DetJ;
				   Matrix[221] += (B1[6]*mu*B1[19]+B1[22]*lambda*B1[3])*DetJ;
				   Matrix[222] += (B1[22]*B1[18]*(mu+lambda)+B1[6]*mu*B1[2]+B1[14]*mu*B1[10])*DetJ;
				   Matrix[223] += (B1[14]*mu*B1[18]+B1[22]*lambda*B1[10])*DetJ;
				   Matrix[224] += (B1[6]*mu*B1[18]+B1[22]*lambda*B1[2])*DetJ;
				   Matrix[225] += (B1[22]*B1[17]*(mu+lambda)+B1[6]*mu*B1[1]+B1[14]*mu*B1[9])*DetJ;
				   Matrix[226] += (B1[14]*mu*B1[17]+B1[22]*lambda*B1[9])*DetJ;
				   Matrix[227] += (B1[6]*mu*B1[17]+B1[22]*lambda*B1[1])*DetJ;
				   Matrix[228] += (B1[22]*B1[16]*(mu+lambda)+B1[6]*mu*B1[0]+B1[14]*mu*B1[8])*DetJ;
				   Matrix[229] += (B1[14]*mu*B1[16]+B1[22]*lambda*B1[8])*DetJ;
				   Matrix[230] += (B1[6]*mu*B1[16]+B1[22]*lambda*B1[0])*DetJ;
				   Matrix[231] += (B1[7]*B1[7]*(mu+lambda)+B1[15]*mu*B1[15]+B1[23]*mu*B1[23])*DetJ;
				   Matrix[232] += (B1[23]*mu*B1[6]+B1[7]*lambda*B1[22])*DetJ;
				   Matrix[233] += (B1[15]*mu*B1[6]+B1[7]*lambda*B1[14])*DetJ;
				   Matrix[234] += (B1[7]*B1[6]*(mu+lambda)+B1[15]*mu*B1[14]+B1[23]*mu*B1[22])*DetJ;
				   Matrix[235] += (B1[23]*mu*B1[5]+B1[7]*lambda*B1[21])*DetJ;
				   Matrix[236] += (B1[15]*mu*B1[5]+B1[7]*lambda*B1[13])*DetJ;
				   Matrix[237] += (B1[7]*B1[5]*(mu+lambda)+B1[15]*mu*B1[13]+B1[23]*mu*B1[21])*DetJ;
				   Matrix[238] += (B1[23]*mu*B1[4]+B1[7]*lambda*B1[20])*DetJ;
				   Matrix[239] += (B1[15]*mu*B1[4]+B1[7]*lambda*B1[12])*DetJ;
				   Matrix[240] += (B1[7]*B1[4]*(mu+lambda)+B1[15]*mu*B1[12]+B1[23]*mu*B1[20])*DetJ;
				   Matrix[241] += (B1[23]*mu*B1[3]+B1[7]*lambda*B1[19])*DetJ;
				   Matrix[242] += (B1[15]*mu*B1[3]+B1[7]*lambda*B1[11])*DetJ;
				   Matrix[243] += (B1[7]*B1[3]*(mu+lambda)+B1[15]*mu*B1[11]+B1[23]*mu*B1[19])*DetJ;
				   Matrix[244] += (B1[23]*mu*B1[2]+B1[7]*lambda*B1[18])*DetJ;
				   Matrix[245] += (B1[15]*mu*B1[2]+B1[7]*lambda*B1[10])*DetJ;
				   Matrix[246] += (B1[7]*B1[2]*(mu+lambda)+B1[15]*mu*B1[10]+B1[23]*mu*B1[18])*DetJ;
				   Matrix[247] += (B1[23]*mu*B1[1]+B1[7]*lambda*B1[17])*DetJ;
				   Matrix[248] += (B1[15]*mu*B1[1]+B1[7]*lambda*B1[9])*DetJ;
				   Matrix[249] += (B1[7]*B1[1]*(mu+lambda)+B1[15]*mu*B1[9]+B1[23]*mu*B1[17])*DetJ;
				   Matrix[250] += (B1[23]*mu*B1[0]+B1[7]*lambda*B1[16])*DetJ;
				   Matrix[251] += (B1[15]*mu*B1[0]+B1[7]*lambda*B1[8])*DetJ;
				   Matrix[252] += (B1[7]*B1[0]*(mu+lambda)+B1[15]*mu*B1[8]+B1[23]*mu*B1[16])*DetJ;
				   Matrix[253] += (B1[15]*B1[15]*(mu+lambda)+B1[7]*mu*B1[7]+B1[23]*mu*B1[23])*DetJ;
				   Matrix[254] += (B1[7]*mu*B1[15]+B1[15]*lambda*B1[7])*DetJ;
				   Matrix[255] += (B1[23]*mu*B1[14]+B1[15]*lambda*B1[22])*DetJ;
				   Matrix[256] += (B1[15]*B1[14]*(mu+lambda)+B1[7]*mu*B1[6]+B1[23]*mu*B1[22])*DetJ;
				   Matrix[257] += (B1[7]*mu*B1[14]+B1[15]*lambda*B1[6])*DetJ;
				   Matrix[258] += (B1[23]*mu*B1[13]+B1[15]*lambda*B1[21])*DetJ;
				   Matrix[259] += (B1[15]*B1[13]*(mu+lambda)+B1[7]*mu*B1[5]+B1[23]*mu*B1[21])*DetJ;
				   Matrix[260] += (B1[7]*mu*B1[13]+B1[15]*lambda*B1[5])*DetJ;
				   Matrix[261] += (B1[23]*mu*B1[12]+B1[15]*lambda*B1[20])*DetJ;
				   Matrix[262] += (B1[15]*B1[12]*(mu+lambda)+B1[7]*mu*B1[4]+B1[23]*mu*B1[20])*DetJ;
				   Matrix[263] += (B1[7]*mu*B1[12]+B1[15]*lambda*B1[4])*DetJ;
				   Matrix[264] += (B1[23]*mu*B1[11]+B1[15]*lambda*B1[19])*DetJ;
				   Matrix[265] += (B1[15]*B1[11]*(mu+lambda)+B1[7]*mu*B1[3]+B1[23]*mu*B1[19])*DetJ;
				   Matrix[266] += (B1[7]*mu*B1[11]+B1[15]*lambda*B1[3])*DetJ;
				   Matrix[267] += (B1[23]*mu*B1[10]+B1[15]*lambda*B1[18])*DetJ;
				   Matrix[268] += (B1[15]*B1[10]*(mu+lambda)+B1[7]*mu*B1[2]+B1[23]*mu*B1[18])*DetJ;
				   Matrix[269] += (B1[7]*mu*B1[10]+B1[15]*lambda*B1[2])*DetJ;
				   Matrix[270] += (B1[23]*mu*B1[9]+B1[15]*lambda*B1[17])*DetJ;
				   Matrix[271] += (B1[15]*B1[9]*(mu+lambda)+B1[7]*mu*B1[1]+B1[23]*mu*B1[17])*DetJ;
				
			 }
		}
	}
}

void C8H::ElementGravity(double* bodyforce, double Gravity)
{
	clear(bodyforce,24);
}

//	Calculate element stress 
void C8H::ElementStress(double* stress4, double* Displacement)
{	
	C8HMaterial* material = dynamic_cast<C8HMaterial*>(ElementMaterial);
	double E = material->E;
	double nv = material->nv;
	double mu = E/(1+nv);
	double lambda = nv*E/((1+nv)*(1-2*nv));

	double disp[24];
	for (unsigned int i = 0; i < 24; i++)
	{
		if (LocationMatrix[i])
		{
			disp[i] = Displacement[LocationMatrix[i] - 1];
		}
		else
		{
			disp[i] = 0;
		}
	}

	double shape[2];
	shape[0] = -1/sqrt(3);
	shape[1] = 1/sqrt(3);
	for (unsigned int m = 0; m < 2; m++)
	{
		for (unsigned int n = 0; n < 2; n++)
		{
			for (unsigned int o = 0; o < 2; o++)
			{
				double xi = shape[m];
				double eta = shape[n];
				double zeta = shape[o];
				double GN[12];
				GN[0] = 0.125*(1 - eta)*(1 - zeta);
				GN[1] = 0.125*(1 + eta)*(1 - zeta);
				GN[2] = 0.125*(1 - eta)*(1 + zeta);
				GN[3] = 0.125*(1 + eta)*(1 + zeta);
				GN[4] = 0.125*(1 - xi)*(1 - zeta);
				GN[7] = 0.125*(1 + xi)*(1 + zeta);
				GN[8] = 0.125*(1 - xi)*(1 - eta);
				GN[9] = 0.125*(1 + xi)*(1 - eta);
				GN[10] = 0.125*(1 + xi)*(1 + eta);
				GN[11] = 0.125*(1 - xi)*(1 + eta);

				double J[9];
				J[0] = GN[0]*(nodes[1]->XYZ[0] - nodes[0]->XYZ[0]) + GN[1]*(nodes[2]->XYZ[0] - nodes[3]->XYZ[0]) + GN[2]*(nodes[5]->XYZ[0] - nodes[4]->XYZ[0]) + GN[3]*(nodes[6]->XYZ[0] - nodes[7]->XYZ[0]);
				J[1] = GN[0]*(nodes[1]->XYZ[1] - nodes[0]->XYZ[1]) + GN[1]*(nodes[2]->XYZ[1] - nodes[3]->XYZ[1]) + GN[2]*(nodes[5]->XYZ[1] - nodes[4]->XYZ[1]) + GN[3]*(nodes[6]->XYZ[1] - nodes[7]->XYZ[1]);
				J[2] = GN[0]*(nodes[1]->XYZ[2] - nodes[0]->XYZ[2]) + GN[1]*(nodes[2]->XYZ[2] - nodes[3]->XYZ[2]) + GN[2]*(nodes[5]->XYZ[2] - nodes[4]->XYZ[2]) + GN[3]*(nodes[6]->XYZ[2] - nodes[7]->XYZ[2]);
				J[5] = GN[4]*(nodes[3]->XYZ[2] - nodes[0]->XYZ[2]) + GN[5]*(nodes[2]->XYZ[2] - nodes[1]->XYZ[2]) + GN[6]*(nodes[7]->XYZ[2] - nodes[4]->XYZ[2]) + GN[7]*(nodes[6]->XYZ[2] - nodes[5]->XYZ[2]);
				J[6] = GN[8]*(nodes[4]->XYZ[0] - nodes[0]->XYZ[0]) + GN[9]*(nodes[5]->XYZ[0] - nodes[1]->XYZ[0]) + GN[10]*(nodes[6]->XYZ[0] - nodes[2]->XYZ[0]) + GN[11]*(nodes[7]->XYZ[0] - nodes[3]->XYZ[0]);
				J[7] = GN[8]*(nodes[4]->XYZ[1] - nodes[0]->XYZ[1]) + GN[9]*(nodes[5]->XYZ[1] - nodes[1]->XYZ[1]) + GN[10]*(nodes[6]->XYZ[1] - nodes[2]->XYZ[1]) + GN[11]*(nodes[7]->XYZ[1] - nodes[3]->XYZ[1]);
				J[8] = GN[8]*(nodes[4]->XYZ[2] - nodes[0]->XYZ[2]) + GN[9]*(nodes[5]->XYZ[2] - nodes[1]->XYZ[2]) + GN[10]*(nodes[6]->XYZ[2] - nodes[2]->XYZ[2]) + GN[11]*(nodes[7]->XYZ[2] - nodes[3]->XYZ[2]);
				double DetJ = J[0]*J[4]*J[8] - J[0]*J[5]*J[7] - J[1]*J[3]*J[8] + J[1]*J[5]*J[6] + J[2]*J[3]*J[7] - J[2]*J[4]*J[6];

				double InvJ[9];
				InvJ[0] = (J[4]*J[8]-J[5]*J[7])/DetJ;
			    InvJ[1] = -(J[1]*J[8]-J[2]*J[7])/DetJ;
			    InvJ[2] = (J[1]*J[5]-J[2]*J[4])/DetJ;
	
				double B1[24];
				B1[0] = -GN[0]*InvJ[0]-GN[4]*InvJ[1]-GN[8]*InvJ[2];
			    B1[1] = GN[0]*InvJ[0]-GN[5]*InvJ[1]-GN[9]*InvJ[2];
			    B1[2] = GN[1]*InvJ[0]+GN[5]*InvJ[1]-GN[10]*InvJ[2];
			    B1[3] = -GN[1]*InvJ[0]+GN[4]*InvJ[1]-GN[11]*InvJ[2];
				B1[4] = -GN[2]*InvJ[0]-GN[6]*InvJ[1]+GN[8]*InvJ[2];
				B1[5] = GN[2]*InvJ[0]-GN[7]*InvJ[1]+GN[9]*InvJ[2];
				B1[6] = GN[3]*InvJ[0]+GN[7]*InvJ[1]+GN[10]*InvJ[2];
				B1[7] = -GN[3]*InvJ[0]+GN[6]*InvJ[1]+GN[11]*InvJ[2];
			    B1[8] = -GN[0]*InvJ[3]-GN[4]*InvJ[4]-GN[8]*InvJ[5];
			    B1[9] = GN[0]*InvJ[3]-GN[5]*InvJ[4]-GN[9]*InvJ[5];
			    B1[10] = GN[1]*InvJ[3]+GN[5]*InvJ[4]-GN[10]*InvJ[5];
			    B1[11] = -GN[1]*InvJ[3]+GN[4]*InvJ[4]-GN[11]*InvJ[5];
				B1[12] = -GN[2]*InvJ[3]-GN[6]*InvJ[4]+GN[8]*InvJ[5];
	
				double EPS[6];
				clear(EPS,6);
				EPS[0] = B1[0]*disp[0]+B1[1]*disp[3]+B1[2]*disp[6]+B1[3]*disp[9]+B1[4]*disp[12]+B1[5]*disp[15]+B1[6]*disp[18]+B1[7]*disp[21];
		   	    EPS[1] = B1[8]*disp[1]+B1[9]*disp[4]+B1[10]*disp[7]+B1[11]*disp[10]+B1[12]*disp[13]+B1[13]*disp[16]+B1[14]*disp[19]+B1[15]*disp[22];
			    EPS[2] = B1[16]*disp[2]+B1[17]*disp[5]+B1[18]*disp[8]+B1[19]*disp[11]+B1[20]*disp[14]+B1[21]*disp[17]+B1[22]*disp[20]+B1[23]*disp[23];
	
				stress4[0] = EPS[1]*lambda+EPS[2]*lambda+EPS[0]*(mu+lambda);
			    stress4[1] = EPS[0]*lambda+EPS[2]*lambda+EPS[1]*(mu+lambda);
			    stress4[2] = EPS[0]*lambda+EPS[1]*lambda+EPS[2]*(mu+lambda);
			    stress4[3] = EPS[3]*mu;
			    stress4[4] = EPS[4]*mu;
			    stress4[5] = EPS[5]*mu;

			}
		}
	}
}
