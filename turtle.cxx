
#include "turtle.h"
#include "utilities.h"

void Turtle::set_global(VektorR2 g_offset)
{
  loc = g_offset;
  //  cout << "set: " << loc[0] << "," << loc[1] << endl;
} // void Turtle::move

void Turtle::set_local(VektorR2 l_offset)
{
  loc += size*l_offset[0]*x_axis + size*l_offset[1]*y_axis;
  // cout << "set: " << loc[0] << "," << loc[1] << endl;
  writePoint(file, loc);
} // void Turtle::move

void Turtle::move(double width)
{
  loc += size*width*dir;
  writePoint(file, loc);
} // void Turtle::move

void Turtle::rot(double angle)
{
  VektorR2 dir_rot;

  dir_rot[0] = cos(angle*2*M_PI/360.0)*dir[0] - sin(angle*2*M_PI/360.0)*dir[1];
  dir_rot[1] = sin(angle*2*M_PI/360.0)*dir[0] + cos(angle*2*M_PI/360.0)*dir[1];

  dir = dir_rot;
  dir.Normalize();
} // void Turtle::rot
