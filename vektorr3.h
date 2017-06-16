/*
*/

#ifndef _VektorR3_h
#define _VektorR3_h
#include <iostream>
#include <math.h>

using namespace std;

class VektorR3
{
 private:
  double c[3]; // Drei Komponenten des dreidim. Vektors
  
 public:
  VektorR3() {}
  VektorR3(double x, double y, double z) {c[0] = x; c[1] = y; c[2] = z;}
  VektorR3(VektorR3& Vektor) {c[0] = Vektor.c[0]; c[1] = Vektor.c[1]; c[2] = Vektor.c[2];}
  VektorR3(VektorR3 const &Vektor) {c[0] = Vektor.c[0]; c[1] = Vektor.c[1]; c[2] = Vektor.c[2];}
  
  double&		  			 operator [] (const int);
  double		  			 	 operator [] (const int) const;
  const bool				 operator == (const VektorR3&) const;
  const VektorR3 		 operator +  (const VektorR3&) const;
  const VektorR3& 		 operator += (const VektorR3&);
  const VektorR3        operator -  (const VektorR3&) const;
  const VektorR3& 		 operator -= (const VektorR3&);
  const double    		 operator *  (const VektorR3&) const;
  friend const VektorR3 operator *  (const double, const VektorR3&);
  const VektorR3        operator *  (const double);
  const VektorR3& 	operator *= (const double);
  const VektorR3	operator /  (const double) const;
  const VektorR3&	operator /= (const double);
  static const VektorR3 CrossProduct(const VektorR3&, const VektorR3&);
  
  void print() const {cout << c[0] << "," << c[1] << "," << c[2];};
  double 		Norm2() const;
  const double          Normalize();
  const VektorR3	A(VektorR3*) const;
  static void		SetRotMatrix(VektorR3*, double, double, double);
  static void		SetRotMatrix(VektorR3*, VektorR3, VektorR3);
  static void		MatrixMult(VektorR3*, VektorR3*, VektorR3*);
  // static void	CalcInverse(const VektorR3 const*, VektorR3*);
  void			LoadVektorR3(istream&);
  void			SaveVektorR3(ostream&);
}; // class VektorR3

// Inline-Elementfunktionen von VektorR3

inline double& VektorR3::operator [] (const int c_index)
{
	return(c[c_index]);
} // VektorR3::operator []

inline double VektorR3::operator [] (const int c_index) const
{
	return(c[c_index]);
} // VektorR3::operator []

inline const bool VektorR3::operator == (const VektorR3& Vektor) const
{
	return((c[0] == Vektor.c[0]) && (c[1] == Vektor.c[1]) && (c[2] == Vektor.c[2]));
} // VektorR3::operator ==

inline const VektorR3 VektorR3::operator + (const VektorR3& Vektor) const
{
	return VektorR3(c[0] + Vektor.c[0], c[1] + Vektor.c[1], c[2] + Vektor.c[2]);
} // VektorR3::operator+

inline const VektorR3& VektorR3::operator += (const VektorR3& Vektor)
{
	c[0] += Vektor.c[0];
	c[1] += Vektor.c[1];
	c[2] += Vektor.c[2];

	return *this;
} // VektorR3::operator +=

inline const VektorR3 VektorR3::operator - (const VektorR3& Vektor) const
{
	return VektorR3(c[0] - Vektor.c[0], c[1] - Vektor.c[1], c[2] - Vektor.c[2]);
} // VektorR3::operator -

inline const VektorR3& VektorR3::operator -= (const VektorR3& Vektor)
{
	c[0] -= Vektor.c[0];
	c[1] -= Vektor.c[1];
	c[2] -= Vektor.c[2];

	return *this;
} // VektorR3::operator -=

inline const double VektorR3::operator * (const VektorR3& Vektor) const
{
	return (c[0]*Vektor.c[0] + c[1]*Vektor.c[1] + c[2]*Vektor.c[2]);
} // VektorR3::operator *

inline const VektorR3 operator * (const double s, const VektorR3& Vektor)
{
	return VektorR3(s*Vektor.c[0], s*Vektor.c[1], s*Vektor.c[2]);
} // friend operator * (const double, const VektorR3&);

inline const VektorR3 VektorR3::operator * (const double s)
{
	return VektorR3(s*c[0], s*c[1], s*c[2]);
} // VektorR3::operator * (const double, const VektorR3&);

inline const VektorR3& VektorR3::operator *= (const double s)
{
	c[0] *= s;
	c[1] *= s;
	c[2] *= s;

	return *this;
} // VektorR3::opertor *= (const double&)

inline const VektorR3 VektorR3::operator / (const double s) const
{
	return VektorR3(c[0]/s, c[1]/s, c[2]/s);
} // VektorR3::operator /

inline const VektorR3& VektorR3::operator /= (const double s)
{
	c[0] /= s;
	c[1] /= s;
	c[2] /= s;

	return *this;
} // VektorR3::operator /=

inline double VektorR3::Norm2() const
{
	return ((double)(sqrt((double)(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]))));
} // VektorR3::Norm2

inline const double VektorR3::Normalize() 
{
   double length;

   if((length = Norm2()) != 0.0)
   {
      *this /= length;
      return length;
   } // if

   return -1.0;
} // VektorR3::Normalize

// Globale Operatoren
/*
inline ostream& operator << (ostream &OutStrm, VektorR3 &v)
{
	v.Save(OutStrm);
	return(OutStrm);
} // operator <<

inline istream& operator >> (istream &InStrm, VektorR3 &v)
{
	v.Load(InStrm);
	return(InStrm);
} // operator >>
*/

inline const VektorR3 VektorR3::CrossProduct(const VektorR3 &a, const VektorR3 &b)
{
	VektorR3 r;

	r[0] = a[1]*b[2] - a[2]*b[1];
	r[1] = a[2]*b[0] - a[0]*b[2];
	r[2] = a[0]*b[1] - a[1]*b[0];

	return r;
} // VektorR3::CrossProduct

inline const VektorR3 VektorR3::A(VektorR3 *A) const
{
	VektorR3 b;

	b[0] = A[0][0]*c[0] + A[1][0]*c[1] + A[2][0]*c[2];
	b[1] = A[0][1]*c[0] + A[1][1]*c[1] + A[2][1]*c[2];
	b[2] = A[0][2]*c[0] + A[1][2]*c[1] + A[2][2]*c[2];

	return(b);
} // VektorR3::A

inline void VektorR3::SetRotMatrix(VektorR3 *Matrix, double Angle_x, double Angle_y, double Angle_z)
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

inline void VektorR3::SetRotMatrix(VektorR3 *Matrix, VektorR3 v1, VektorR3 v2)
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

		// CalcInverse(SpinCS, SpinCSInv);
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
/*
	if(iteration == 4)
		MessageBox(NULL, "Die Rotationsmatrix konnte nicht berechnet werden.", "VektorR3::SetRotMatrix", MB_OK|MB_ICONQUESTION);
*/
} // VektorR3::SetRotMatrix

inline void VektorR3::MatrixMult(VektorR3 *pParaE, VektorR3 *A, VektorR3 *B)
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
/*
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
*/
inline void VektorR3::SaveVektorR3(ostream &strm)
{
	strm.write((char*)c, sizeof(VektorR3));
} // VektorR3::Save

inline void VektorR3::LoadVektorR3(istream &strm)
{
	strm.read((char*)c, sizeof(VektorR3));
} // VektorR3::Save

#endif
