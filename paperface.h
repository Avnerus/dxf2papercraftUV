
#define MAX_POINTS_PER_FACE        300
#define ADHERENT_SIZE              0.15

// enclosing angle in degrees between two neighboring
// faces which is considered to be zero

#define NEGLIGIBLE_ENCLOSING_ANGLE 0.025

#include "vektorr2.h"
#include "vektorr3.h"
#include "json/json.h"

class paperFace {
private:
  double cost; // cost value needed in generateCutoutSheet

public:
  long       point[MAX_POINTS_PER_FACE];
  VektorR2   projection[MAX_POINTS_PER_FACE];
  bool       adherent_connected[MAX_POINTS_PER_FACE];
  int        adherent_ID[MAX_POINTS_PER_FACE];
  paperFace* neighbor[MAX_POINTS_PER_FACE];
  paperFace* my_father; // used for generating the cutout sheet

  int  no_points;
  int  ID;
  bool processed, drawn, adherent, visited;

  bool     writeToDXF(ofstream*, Json::Value &, VektorR2, VektorR2, int, VektorR2, int, bool, bool);  
  void     writeAdherent(ofstream*, double, int&);
  void     printProjection(int);
  void     print3DVertices(int);
  VektorR2 calcNormal(int);
  VektorR3 calc3DNormal();

  bool     mergeFace(paperFace&);
  void     copy(paperFace&);
  bool     sharesSomeEdge(paperFace&, int&, int&);
  int      getSharedEdgeID(paperFace&, int);
  bool     SharesProjectedVertex(int, VektorR2, VektorR2);
  void     clearConnectedFlags();
  bool     EdgeOnEdge(int, VektorR2, VektorR2);
  bool     EdgeOnEdge(int, VektorR3, VektorR3);
  void     printProjection();
  void     plotProjection();
  void     plot3DVertices();
  bool     CheckColinear();
  int      getConnectedFace(int, int);
  void     removeColinearPoints();
  bool     IsInnerEdge(int, int);
  bool     PointInside(VektorR2, VektorR2);
  void     initialize();
  void     generateConnectedFaceGraph(int);
  bool     generateCutoutSheet(ofstream*, Json::Value &, VektorR2, VektorR2, int, VektorR2, int, double, bool, bool);
  bool     sharesProjectedEdge(int, int);
  VektorR2 getMidPoint();
  int      index_save(int);

  static double sumMinimalEnclosingAngles(int);

 private:
  void     swapVertexData(int, int);
  void     swapEdgeData(int, int);
  bool     sharesThisEdge(paperFace&, int, int&, bool);
  bool     intersectsOtherFace(int);
  double   estimateExtend();
  double   getCumulatedDistance(int);
  bool     getAdjacentEdge(int, int, VektorR2&, VektorR2&, VektorR2&, bool);
  double   getEnclosedAngle(int, int, VektorR2&, bool);
  bool     calcCostForCutoutSheet(int, VektorR2, int);
  double   getCost() {return cost;};
  double   getMinimalEnclosingAngle(int);

}; // class paperFace

inline int paperFace::index_save(int index)
{
  if(index < 0) return no_points+index;
  return index % no_points;
}
