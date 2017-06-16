
#ifndef _Turtle_h
#define _Turtle_h

#include <fstream>
#include "vektorr2.h"

using namespace std;

class Turtle {
  ofstream* file;
  VektorR2  loc, dir;
  VektorR2  x_axis, y_axis;
  double     size;

public:
  Turtle(ofstream* DXFFile, VektorR2 x, VektorR2 y, double s) : file(DXFFile), dir(x), x_axis(x), y_axis(y), size(s) 
  {loc = VektorR2(0.0, 0.0);};

  void set_global(VektorR2);
  void set_local(VektorR2);
  void move(double);
  void rot(double);
}; // class Turtle

#endif
