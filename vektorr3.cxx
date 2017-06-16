
#include "vektorr3.h"

// Elementfunktionen von VektorR3
/*
const VektorR3 VektorR3::CrossProduct(const VektorR3 &a, const VektorR3 &b)
{
	VektorR3 r;

	r[0] = a[1]*b[2] - a[2]*b[1];
	r[1] = a[2]*b[0] - a[0]*b[2];
	r[2] = a[0]*b[1] - a[1]*b[0];

	return r;
} // VektorR3::CrossProduct
*/
const VektorR3 VektorR3::A(VektorR3 *A) const
{
	VektorR3 b;

	b[0] = A[0][0]*c[0] + A[1][0]*c[1] + A[2][0]*c[2];
	b[1] = A[0][1]*c[0] + A[1][1]*c[1] + A[2][1]*c[2];
	b[2] = A[0][2]*c[0] + A[1][2]*c[1] + A[2][2]*c[2];

	return(b);
} // VektorR3::A

void VektorR3::SetRotMatrix(VektorR3 *Matrix, double Angle_x, double Angle_y, double Angle_z)
{
	double a_x = (double)Angle_x*(M_PI/180.0);
	double a_y = (double)Angle_y*(M_PI/180.0);
	double a_z = (double)Angle_z*(M_PI/180.0);

	Matrix[0][0] = (double)(cos(a_y)*cos(a_z));
	Matrix[0][1] = (double)(sin(a_x)*sin(a_y)*cos(a_z)+cos(a_x)*sin(a_z));
	Matrix[0][2] = (double)(cos(a_x)*sin(a_y)*cos(a_z)-sin(a_x)*sin(a_z));

	Matrix[1][0] = (double)(-cos(a_y)*sin(a_z));
	Matrix[1][1] = (double)(-sin(a_x)*sin(a_y)*sin(a_z)+cos(a_x)*cos(a_z));
	Matrix[1][2] = (double)(-cos(a_x)*sin(a_y)*sin(a_z)-sin(a_x)*cos(a_z));

	Matrix[2][0] = (double)-sin(a_y);
	Matrix[2][1] = (double)(sin(a_x)*cos(a_y));
	Matrix[2][2] = (double)(cos(a_x)*cos(a_y));
} // VektorR3::SetRotMatrix

void VektorR3::SetRotMatrix(VektorR3 *Matrix, VektorR3 v1, VektorR3 v2)
{
	VektorR3 SpinCS[3], SpinCSInv[3], RotMat[3], TmpMat[3], test;
	double    Angle, sign1, sign2, tf, MinDist = 100000.0;
	int		iteration = 0;

	v1 /= v1.Norm2(); v2 /= v2.Norm2();

	if(v1 == v2)
	{
		Matrix[0] = VektorR3(1.0, 0.0, 0.0);
		Matrix[1] = VektorR3(0.0, 1.0, 0.0);
		Matrix[2] = VektorR3(0.0, 0.0, 1.0);

      return;
	} // if

	do
	{
		switch(iteration)
		{
			case 0:
				sign1 = 1.0; sign2 = 1.0;
				break;

			case 1:
				sign1 = -1.0; sign2 = 1.0;
				break;

			case 2:
				sign1 = 1.0; sign2 = -1.0;
				break;

			case 3:
				sign1 = -1.0; sign2 = -1.0;
				break;
		} // switch
		iteration++;

		// Lokales (Dreh-) KS berechnen, in dem um die Z-Achse gedreht wird

		SpinCS[0] = v1;
		SpinCS[2] = CrossProduct(v1, v2); SpinCS[2] /= sign1*SpinCS[2].Norm2();
		SpinCS[1] = CrossProduct(v1, SpinCS[2]); SpinCS[1] /= sign2*SpinCS[1].Norm2();

		// Inverses DrehKS berechnen

		CalcInverse(SpinCS, SpinCSInv);
//		MatrixMult(TmpMat, SpinCS, SpinCSInv);
		tf = v1*v2;
		Angle = ((double)acos((double)tf))*((double)(180.0/M_PI));

		SetRotMatrix(RotMat, 0.0, 0.0, Angle);

		MatrixMult(TmpMat, RotMat, SpinCSInv);
		MatrixMult(Matrix, SpinCS, TmpMat);

		test = v1.A(Matrix);
		tf = (v2-test).Norm2();
		if(tf < MinDist) MinDist = tf;
	} while((tf > 0.001) && (iteration < 4));

	if(iteration == 4)
	  cerr  << "(VektorR3::SetRotMatrix): Die Rotationsmatrix konnte nicht berechnet werden.";
} // VektorR3::SetRotMatrix

void VektorR3::MatrixMult(VektorR3 *pParaE, VektorR3 *A, VektorR3 *B)
{
	VektorR3 E[3];

	// Matrixmultiplikation: E = A*B

	for(int zeile = 0; zeile < 3; zeile++)
	{
		for(int spalte = 0; spalte < 3; spalte++)
		{
			E[spalte][zeile] = 0.0;

			for(int i = 0; i < 3; i++)
				E[spalte][zeile] += A[i][zeile] * B[spalte][i];
		} // for
	} // for

	for(int i = 0; i < 3; i++)
		pParaE[i] = E[i];
} // VektorR3::MatrixMult

void VektorR3::CalcInverse(const VektorR3 const *SourceMat, VektorR3 *DestMat)
{
	Gauss LGS(9);

	// Besetzen des 9x9 LGS mit den Elementen der Quellmatrix

	LGS.A(1, 1, SourceMat[0][0]); LGS.A(1, 2, SourceMat[0][1]); LGS.A(1, 3, SourceMat[0][2]);
	LGS.A(2, 1, SourceMat[1][0]); LGS.A(2, 2, SourceMat[1][1]); LGS.A(2, 3, SourceMat[1][2]);
	LGS.A(3, 1, SourceMat[2][0]); LGS.A(3, 2, SourceMat[2][1]); LGS.A(3, 3, SourceMat[2][2]);

	LGS.A(4, 4, SourceMat[0][0]); LGS.A(4, 5, SourceMat[0][1]); LGS.A(4, 6, SourceMat[0][2]);
	LGS.A(5, 4, SourceMat[1][0]); LGS.A(5, 5, SourceMat[1][1]); LGS.A(5, 6, SourceMat[1][2]);
	LGS.A(6, 4, SourceMat[2][0]); LGS.A(6, 5, SourceMat[2][1]); LGS.A(6, 6, SourceMat[2][2]);

	LGS.A(7, 7, SourceMat[0][0]); LGS.A(7, 8, SourceMat[0][1]); LGS.A(7, 9, SourceMat[0][2]);
	LGS.A(8, 7, SourceMat[1][0]); LGS.A(8, 8, SourceMat[1][1]); LGS.A(8, 9, SourceMat[1][2]);
	LGS.A(9, 7, SourceMat[2][0]); LGS.A(9, 8, SourceMat[2][1]); LGS.A(9, 9, SourceMat[2][2]);

	// Ergebnis-Vektor definieren

	LGS.b(1, 1.0); LGS.b(5, 1.0); LGS.b(9, 1.0);
	LGS.Calc(NULL); // Keine Speicherung im Stream

	// Zielmatrix besetzen

	DestMat[0][0] = LGS.x(1); DestMat[1][0] = LGS.x(2); DestMat[2][0] = LGS.x(3);
	DestMat[0][1] = LGS.x(4); DestMat[1][1] = LGS.x(5); DestMat[2][1] = LGS.x(6);
	DestMat[0][2] = LGS.x(7); DestMat[1][2] = LGS.x(8); DestMat[2][2] = LGS.x(9);
} // VektorR3::CalcInverse

void VektorR3::SaveVektorR3(ostream &strm)
{
	strm.write((char*)c, sizeof(VektorR3));
} // VektorR3::Save

void VektorR3::LoadVektorR3(istream &strm)
{
	strm.read((char*)c, sizeof(VektorR3));
} // VektorR3::Save
