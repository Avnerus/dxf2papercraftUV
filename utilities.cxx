#include <fstream>
#include <stdlib.h>
#include "utilities.h"
#include "paperface.h"

extern double     epsilon;
extern ofstream*  DXFFile;
extern VektorR3*  pointArray;
extern paperFace* faceArray;

bool IntersectLineLineParam(VektorR2 a1, VektorR2 a2, VektorR2 b1, VektorR2 b2, double& p1, double& p2)
{
  VektorR2 r, s, n;

  s = b2-b1;
  r = a2-a1;
  n = VektorR2(-s[1], s[0]);

  double tmp = r*n;

  if(fabs(tmp) < 0.0001) return false;

  p1 = ((b1-a1)*n)/tmp;

  n = VektorR2(-r[1], r[0]);
  tmp = s*n;

  if(fabs(tmp) < 0.0001) return false;

  p2 = ((a1-b1)*n)/tmp;

  return true;
}

bool IntersectLineLine(VektorR2 a1, VektorR2 a2, VektorR2 b1, VektorR2 b2)
{
  double p1, p2;

  if(IntersectLineLineParam(a1, a2, b1, b2, p1, p2) == true) {
    if((p1 < 0.0001) || (p1 > 0.9999)) return false;
    if((p2 < 0.0001) || (p2 > 0.9999)) return false;

    //    cout << "    p1 = " << p1 << "   p2 = " << p2 << endl;

    return true;
  } // if

  return false;
} // IntersectLineLine

bool PointBetweenLineEnds(VektorR2 p, VektorR2 x1, VektorR2 x2)
{
  // We assume, the meeting the vertex does
  // not count as lieing on the line

  if((p-x1).Norm2() < epsilon) return false;
  if((p-x2).Norm2() < epsilon) return false;

  VektorR2 dir = x2-x1;

  for(int i = 0; i < 2; i++) {
    if(dir[i] > 0.0001) {
      double r = (p[i]-x1[i])/dir[i];
      
      if((r > epsilon) && (r < (1.0-epsilon))) return true;
    } // if
  } // for

  return false;
} // PointBetweenLineEnds

bool PointOnLine(VektorR2 p, VektorR2 direction, VektorR2 x1, VektorR2 x2)
{
  // We assume, the meeting the vertex does
  // not count as lieing on the line

  if((p-x1).Norm2() < epsilon) return false;
  if((p-x2).Norm2() < epsilon) return false;

  VektorR2 dir = x2-x1;
  VektorR2 n = VektorR2(-dir[1], dir[0]);
  VektorR2 diff = (p-x1);

  diff.Normalize();
  n.Normalize();

  if(fabs(diff*n) < 0.001) {

    // p lies between x1 and x2. But the 
    // lines only overlap if their point
    // into the same direction

    dir.Normalize();
    direction.Normalize();

    if((dir*direction) > 0.999) return true;
  } // if

  return false;
} // PointOnLine

bool PointOnLine(VektorR3 p, VektorR3 direction, VektorR3 x1, VektorR3 x2)
{
  // We assume, the meeting the vertex does
  // not count as lieing on the line

  if((p-x1).Norm2() < epsilon) return false;
  if((p-x2).Norm2() < epsilon) return false;

  VektorR3 dir = x2-x1;
  double    param[3];
  bool     match[3];

  for(int i = 0; i < 3; i++) {

    match[i] = false;

    if(fabs(dir[i]) > epsilon) {
      param[i] = (p[i] - x1[i])/dir[i];
      if((param[i] < epsilon) || (param[i] > (1.0-epsilon))) return false;
      else {
	match[i] == true;
      } // else
    } // if
  } // for

  double sum = 0.0;

  for(int i = 0; i < 3; i++) {
    if(match[i] == true) {
      sum += param[i];
    } // if
  } // for

  sum /= 3.0;

  for(int i = 0; i < 3; i++) {
    if(match[i] == true) {
      if(fabs(sum-param[i]) > 0.0001) return false;
    } // if
  } // for

  return true;
} // PointOnLine

void writePoint(ofstream* DXFFile, VektorR2 v)
{
  (*DXFFile) << "  0\nVERTEX" << endl;
  (*DXFFile) << " 10" << endl;
  (*DXFFile) << v[0] << endl;
  
  (*DXFFile) << " 20" << endl;
  (*DXFFile) << v[1] << endl;
} // writePoint

void writeLine(ofstream* DXFFile, VektorR2 p1, VektorR2 p2)
{
  (*DXFFile) << "  0\nPOLYLINE\n 70\n     0" << endl;

  writePoint(DXFFile, p1);
  writePoint(DXFFile, p2);

  (*DXFFile) << "  0\nSEQEND" << endl;  
} // writeLine

void writeDigit(ofstream* DXFFile, VektorR2 offset, VektorR2 x_axis, VektorR2 y_axis, double size, int digit)
{
  Turtle turtle(DXFFile, x_axis, y_axis, size);

  (*DXFFile) << "  0\nPOLYLINE\n 70\n     0" << endl;

  turtle.set_global(offset);

  switch(digit) {
  case 0:
    turtle.set_local(VektorR2(0.0, 1.0));
    turtle.move(1); turtle.rot(90);
    turtle.move(1); turtle.rot(90);
    turtle.move(1); turtle.rot(90);
    turtle.move(1); turtle.rot(90);
    break;

  case 1:
    turtle.set_local(VektorR2(0.0, 0.0));
    turtle.move(0.5);
    turtle.rot(-90); turtle.move(1);
    turtle.rot(-90); turtle.move(0.5);
    turtle.rot(180); turtle.move(1.0);
    turtle.rot(90); turtle.move(1.0/3.0);    
    break;

  case 2:
    turtle.set_local(VektorR2(0.0, 0.0));
    turtle.move(1); turtle.rot(-90);
    turtle.move(0.5); turtle.rot(-90);
    turtle.move(1); turtle.rot(90);
    turtle.move(0.5); turtle.rot(90);
    turtle.move(1);
    break;

  case 3:
    turtle.set_local(VektorR2(0.0, 0.0));
    turtle.move(1); turtle.rot(-90);
    turtle.move(0.5); turtle.rot(-90);
    turtle.move(0.5); turtle.rot(180);
    turtle.move(0.5); turtle.rot(-90);
    turtle.move(0.5); turtle.rot(-90);
    turtle.move(1); 
    break;

  case 4:
    turtle.set_local(VektorR2(0.0, 0.0));
    turtle.rot(-90); turtle.move(0.5); 
    turtle.rot(90); turtle.move(1.0); 
    turtle.rot(90); turtle.move(0.25); 
    turtle.rot(180); turtle.move(0.75); 
    break;

  case 5:
    turtle.set_local(VektorR2(1.0, 0.0));
    turtle.rot(180); turtle.move(1.0); 
    turtle.rot(90); turtle.move(0.5); 
    turtle.rot(90); turtle.move(1.0); 
    turtle.rot(-90); turtle.move(0.5); 
    turtle.rot(-90); turtle.move(1.0); 
    break;

  case 6:
    turtle.set_local(VektorR2(1.0/3.0, 0.0));
    turtle.rot(180); turtle.move(1.0/3.0); 
    turtle.rot(90); turtle.move(2.0/3.0); 
    turtle.rot(90); turtle.move(1.0); 
    turtle.rot(-90); turtle.move(1.0/3.0); 
    turtle.rot(-90); turtle.move(1.0); 
    turtle.rot(-90); turtle.move(1.0/3.0); 
    break;

  case 7:
    turtle.set_local(VektorR2(0.0, 1.0/3.0));
    turtle.rot(90); turtle.move(1.0/3.0);
    turtle.rot(-90);     
    turtle.move(1.0); turtle.rot(-90); 
    turtle.move(1.0/3.0); turtle.rot(-(90.0-asin((1.0/3.0)/0.5)*360.0/(2.0*M_PI))); 
    turtle.move(sqrt(1.0/9.0+0.25)); turtle.rot((90.0-asin((1.0/3.0)/0.5)*360.0/(2.0*M_PI)));
    turtle.move(1.0/3.0);
    break;

  case 8:
    turtle.set_local(VektorR2(0.75, 0.5));
    turtle.rot(90); turtle.move(0.5); 
    turtle.rot(90); turtle.move(0.5); 
    turtle.rot(90); turtle.move(0.5); 
    turtle.rot(-90); turtle.move(0.25); 
    turtle.rot(90); turtle.move(0.5); 
    turtle.rot(90); turtle.move(1.0); 
    turtle.rot(90); turtle.move(0.5); 
    turtle.rot(90); turtle.move(0.75); 
    break;

  case 9:
    turtle.set_local(VektorR2(1.0, 1.0/3.0));
    turtle.rot(90); turtle.move(1.0/3.0); 
    turtle.rot(90); turtle.move(1.0); 
    turtle.rot(90); turtle.move(1.0/3.0); 
    turtle.rot(90); turtle.move(1.0); 
    turtle.rot(-90); turtle.move(2.0/3.0); 
    turtle.rot(-90); turtle.move(1.0/3.0); 
    break;

  } // switch

  (*DXFFile) << "  0\nSEQEND" << endl;
} // writeDigit

void writeNumber(ofstream* DXFFile, VektorR2 offset, VektorR2 x_axis, VektorR2 y_axis, double size, long number)
{
  int digit[21], counter = 0;
  
  x_axis.Normalize();
  y_axis.Normalize();

  while((number > 0) && (counter < 20)) {
    digit[counter++] = number % 10;
    number -= (number % 10);
    number /= 10;
  } // while

  for(; counter >= 0; counter--) {
    writeDigit(DXFFile, offset, x_axis, y_axis, size, digit[counter]);
    offset += 1.5*size*x_axis;
  } // for
} // writeNumber

int getPointIndex(VektorR3 point, int no_points, double threshold)
{
  for(int i = 0; i < no_points; i++) {
    if((pointArray[i] - point).Norm2() < threshold) return i;
    // if((pointArray[i] - point).Norm2() == 0.0) return i;
  } // for
} // getPointIndex

int getNextFace(int index_last_face, int face_count, VektorR2& offset, VektorR2& dir, int& start_with_index_returned_face,
                int& start_with_index_old_face, VektorR2& normal)
{
  for(int i = 0; i < face_count; i++) {
    if((faceArray[i].processed == false) && (i != index_last_face)) {

      int cp = faceArray[i].no_points;

      if(faceArray[index_last_face].sharesSomeEdge(faceArray[i], start_with_index_old_face, start_with_index_returned_face) == true) {
        offset = faceArray[index_last_face].projection[start_with_index_old_face];
        dir = faceArray[index_last_face].projection[faceArray[index_last_face].index_save(start_with_index_old_face+1)] - faceArray[index_last_face].projection[start_with_index_old_face];

        normal = faceArray[index_last_face].calcNormal(start_with_index_old_face);

        faceArray[index_last_face].adherent_connected[start_with_index_old_face] = true;
        faceArray[i].adherent_connected[start_with_index_returned_face] = true;

        if(faceArray[i].adherent == true) {
          cout << "error: getNextFace() returned adherent face" << endl;
          exit(-1);
        } // if

        return i;
      } // if

    } // if
  } // for

  return -1;
} // getNextFace
