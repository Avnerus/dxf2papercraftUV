#ifndef _VektorR2_h
#define _VektorR2_h
#include <iostream>
#include <math.h> 

// #include <stdlib.h>
// #include <fstream>
// #include <iomanip>

using namespace std;

class VektorR2
{
private:
		double c[2]; // Komponenten des zweidim. Vektors
		
public:
		VektorR2() {}
		VektorR2(double x, double y) {c[0] = x; c[1] = y;}
		VektorR2(const VektorR2& Vektor) {c[0] = Vektor.c[0]; c[1] = Vektor.c[1];}

		double&		operator [] (const int); 
		double   	operator [] (const int) const;
		const bool	operator == (const VektorR2&) const;
		VektorR2        operator +  (const VektorR2&) const;
		VektorR2&       operator += (const VektorR2&);
		VektorR2        operator -  (const VektorR2&) const;
		VektorR2&       operator -= (const VektorR2&);
		double           operator *  (const VektorR2&) const;
		friend VektorR2 operator *  (const double, const VektorR2&);
		VektorR2        operator *  (const double);
		VektorR2& 	operator *= (const double);
		VektorR2	operator /  (const double) const;
		VektorR2&	operator /= (const double);
		static const VektorR2 CrossProduct(VektorR2&, VektorR2&);
		
		void            print() const {cout << c[0] << "," << c[1];};
		double 		Norm2() const;
		double		Normalize();

		VektorR2	Rotate(const double) const;
}; // class VektorR2

// Inline-Elementfunktionen von VektorR2

inline double& VektorR2::operator [] (const int c_index)
{
  return(c[c_index]);
} // VektorR2::operator []

inline double VektorR2::operator [] (const int c_index) const
{
	return(c[c_index]);
} // VektorR2::operator []

inline const bool VektorR2::operator == (const VektorR2& Vektor) const
{
	return((c[0] == Vektor.c[0]) && (c[1] == Vektor.c[1]));
} // VektorR2::operator ==

inline VektorR2 VektorR2::operator + (const VektorR2& Vektor) const
{
	return VektorR2(c[0] + Vektor.c[0], c[1] + Vektor.c[1]);
} // VektorR2::operator+

inline VektorR2& VektorR2::operator += (const VektorR2& Vektor)
{
	c[0] += Vektor.c[0];
	c[1] += Vektor.c[1];
	
	return *this;
} // VektorR2::operator +=

inline VektorR2 VektorR2::operator - (const VektorR2& Vektor) const
{
	return VektorR2(c[0] - Vektor.c[0], c[1] - Vektor.c[1]);
} // VektorR2::operator -

inline VektorR2& VektorR2::operator -= (const VektorR2& Vektor)
{
	c[0] -= Vektor.c[0];
	c[1] -= Vektor.c[1];
	
	return *this;
} // VektorR2::operator -=

inline double VektorR2::operator * (const VektorR2& Vektor) const
{
	return (c[0]*Vektor.c[0] + c[1]*Vektor.c[1]);
} // VektorR2::operator *

inline VektorR2 operator * (const double s, const VektorR2& Vektor)
{
	return VektorR2(s*Vektor.c[0], s*Vektor.c[1]);
} // friend operator * (const double, const VektorR2&);

inline VektorR2 VektorR2::operator * (const double s)
{
	return VektorR2(s*c[0], s*c[1]);
} // VektorR2::operator * (const double, const VektorR2&);

inline VektorR2& VektorR2::operator *= (const double s)
{
	c[0] *= s;
	c[1] *= s;
	
	return *this;
} // VektorR2::opertor *= (const double&)

inline VektorR2 VektorR2::operator / (const double s) const
{
	return VektorR2(c[0]/s, c[1]/s);
} // VektorR2::operator /

inline VektorR2& VektorR2::operator /= (const double s)
{
	c[0] /= s;
	c[1] /= s;
	
	return *this;
} // VektorR2::operator /=

inline double VektorR2::Norm2() const
{
	return ((double)(sqrt((double)(c[0]*c[0] + c[1]*c[1]))));
} // VektorR2::Norm2

inline double VektorR2::Normalize() {
  double length = Norm2();
  
  if(length > 0.0) *this /= length;
  return length;
} // VektorR2::Normalize

inline VektorR2 VektorR2::Rotate(const double angle) const {
	
	double   current_angle;
	double    length;
	VektorR2 v = *this;

	length = v.Norm2();
	v /= length;

	current_angle = acos((double)v[0]);
	if(v[1] < 0.0) current_angle = 2.0*M_PI - current_angle;
	current_angle += (double)angle;

	return length*VektorR2((double)cos(current_angle), (double)sin(current_angle));
} // VektorR2::Rotate
#endif
