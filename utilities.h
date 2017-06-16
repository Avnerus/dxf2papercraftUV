
#include "vektorr2.h"
#include "vektorr3.h"
#include "turtle.h"

bool PointBetweenLineEnds(VektorR2, VektorR2, VektorR2);
void writePoint(ofstream*, VektorR2);
bool IntersectLineLineParam(VektorR2, VektorR2, VektorR2, VektorR2, double&, double&);
bool IntersectLineLine(VektorR2, VektorR2, VektorR2, VektorR2);
bool PointBetweenLineEnds(VektorR2, VektorR2, VektorR2);
bool PointOnLine(VektorR2, VektorR2, VektorR2, VektorR2);
bool PointOnLine(VektorR3, VektorR3, VektorR3, VektorR3);
void writePoint(ofstream*, VektorR2);
void writeDigit(ofstream*, VektorR2, VektorR2, VektorR2, double, int);
void writeNumber(ofstream*, VektorR2, VektorR2, VektorR2, double, long);
void writeLine(ofstream*, VektorR2, VektorR2);
int getNextFace(int, int, VektorR2&, VektorR2&, int&, int&, VektorR2&);
int getPointIndex(VektorR3, int, double);
/*
inline int index_save(int index, int max)
{
  if(index < 0) return max+index;
  return index % max;
} // index_underflow
*/
