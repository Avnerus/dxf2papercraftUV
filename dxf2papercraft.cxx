#include <dime/entities/Entity.h>
#include <dime/Input.h>
#include <dime/Model.h>
#include <dime/State.h>
#include <dime/util/BSPTree.h>
#include <stdlib.h>
#include <fstream>
#include <iomanip>
#include <getopt.h>

#include "utilities.h"
#include "turtle.h"
#include "paperface.h"

#include "json/json.h"

const int MAX_NO_SPLIT_FACES = 100;

VektorR3*  pointArray = NULL;
paperFace* faceArray = NULL;
double     epsilon;
int        global_adherent_counter = 0;
int        global_counter = 0;
ofstream*  debugFile = NULL;
bool       SHOW_DEBUG;
int        DEBUG_COUNTER = 0;
int        SUPPRESSED_ADHERENTS = 0;
int        STRATEGY;
bool       MERGE_FACES;
bool       ALLOW_INTERSECTION;
bool       GENERATE_GLUE_TABS;
bool       CONNECT_FACES;
bool       NUMBER_FACES;
bool       FORCE_ADHERENTS;
int        indexSplitFaces;
int        arraySplitFaces[MAX_NO_SPLIT_FACES];

class faceset_data {
public:
  faceset_data() {
  }
  ~faceset_data() {
  }
  void clean() {
    bspTree.clear();
    indices.makeEmpty();
  }
public:
  dimeArray <int> indices;
  dimeArray <int> lineIndices;
  dimeBSPTree bspTree;
}; 

// some static arrays 
static dimeArray <int> indexarray;
static dimeArray <int> tmparray;
static dimeArray <dimeVec3f> vertexarray;
static faceset_data *facesets;

static bool my_dime_callback(const dimeState * const state, 
			     dimeEntity *entity, void *);
static void convertModel(faceset_data *data, const int colnum, const char* const);

ofstream* openDXFFile(const char* const p_filename)
{
  ofstream* DXFfile = new ofstream(p_filename);

  (*DXFfile) << "  0\nSECTION\n  2\nENTITIES" << endl;
  
  return DXFfile;
} // openDXFFile

void closeDXFFile(ofstream* DXFFile)
{
  (*DXFFile) << "   0\nENDSEQ\n  0\nEOF" << endl;
  DXFFile->close();
  delete DXFFile;
} // closeDXFFile

//
//  Call this to initialize static members.
//
static void
init()
{
  if (facesets) delete [] facesets;
  facesets = new faceset_data[255];
}

//
//  Call this to free static data.
//
static void
end()
{
  delete [] facesets;
  facesets = NULL;
  indexarray.makeEmpty();
  vertexarray.makeEmpty();

  if(pointArray) delete [] pointArray;
  if(faceArray) delete [] faceArray;
}

static bool
convertDxf2Wrl(const char * const dxffilename,
	       const char * const wrlfilename)
{
  dimeInput in;
  if (!in.setFile(dxffilename)) {
    return false;
  }
  
  dimeModel model;
  if (!model.read(&in)) {
    fprintf(stderr,"Read error in line: %d\n", in.getFilePosition());
    return false;
  }

  FILE *out = fopen("delete_me_later", "w");
  if (!out) return false;
  
  fprintf(out,
	  "#VRML V1.0 ascii\n\n"
	  "Separator {\n");
  
  model.traverseEntities(my_dime_callback, &model, 
			 false, true, false);
  
  for (int i = 0; i < 255; i++) {
    convertModel(&facesets[i], i+1, wrlfilename);
  }

  fprintf(out,"}\n");


  return true;
}

void calcExtend(int facecnt, double& extend_x, double& extend_y)
{
  double min_x, min_y, max_x, max_y;
  double extend;
  bool  inside = false;

  min_x = min_y = 1000000;
  max_x = max_y = -1000000;

  for(int face = 0; face < facecnt; face++) {

    if(faceArray[face].drawn == true) {

      inside = true;

      for(int edge = 0; edge < faceArray[face].no_points; edge++) {
	
	if(min_x > faceArray[face].projection[edge][0]) min_x = faceArray[face].projection[edge][0];
	if(max_x < faceArray[face].projection[edge][0]) max_x = faceArray[face].projection[edge][0];
	
	if(min_y > faceArray[face].projection[edge][1]) min_y = faceArray[face].projection[edge][1];
	if(max_y < faceArray[face].projection[edge][1]) max_y = faceArray[face].projection[edge][1];
	
      } // for
    } // if
  } // for

  if(inside) {
    extend_x = max_x - min_x;
    extend_y = max_y - min_y;
  } // if
  else {
    extend_x = 0.0;
    extend_y = 0.0;  
  } // else
} // calcExtend

static void 
normalizeUvs(Json::Value &faces) {
    int min_x = 1000;
    int max_x = -1000;
    int min_y = 1000;
    int max_y = -1000;

    // First get min max
    for (auto itr : faces) {
        Json::Value uvs = itr["uvs"];
        for (auto uv_itr : uvs) {
            int u = uv_itr["u"].asInt();
            int v = uv_itr["v"].asInt();

            if (u < min_x) {
                min_x = u;
            }
            if (v < min_y) {
                min_y = v;
            }
            if (u > max_x) {
                max_x = u;
            }
            if (v > max_y) {
                max_y = v;
            }
        }
    }

    std::cout << "Min X: " << min_x << " Min Y: " << min_y << " Max X: " << max_x << " Max Y: " << max_y << std::endl;
    int width = max_x - min_x;
    int height = max_y - min_y;
    std::cout << "Width: " << width << ",Height" << height << std::endl;

    // Make it square
    if (width < height) {
        width = height;
    } else {
        height = width;
    }
    // Offset to start at 0
    int offset_x = width / 2;
    int offset_y = height / 2;

    std::cout << "Fixed Width: " << width << ",Height" << height << std::endl;
    std::cout << "Offset Width: " << offset_x << ",Offset Height" << offset_y << std::endl;

    for (auto & itr : faces) {
        Json::Value &uvs = itr["uvs"];
        for (auto & uv_itr : uvs) {
            float u = uv_itr["u"].asFloat();
            float v = uv_itr["v"].asFloat();
            uv_itr["u"] = Json::Value((u + offset_x) / width);
            uv_itr["v"] = Json::Value((v + offset_y) / height);
        }
    }
}

static void convertModel(faceset_data *data, const int colnum, const char* const filename)
{
  int i;
  int facecnt = 1;

  int icnt = data->indices.count(); // dxf-stuff I don't understand...
  int lcnt = data->lineIndices.count();
  int vcnt = data->bspTree.numPoints();


  Json::Value root;
  Json::Value faces(Json::arrayValue);

  double min_x, max_x, min_y, max_y, min_z, max_z;

  if (vcnt == 0 || (icnt == 0 && lcnt == 0)) return;

  // If not present, append a -1 at the end of the index array
  if (icnt && data->indices[icnt-1] >= 0) data->indices.append(-1); // dxf-stuff I don't understand
  if (lcnt && data->lineIndices[lcnt-1] >= 0) data->lineIndices.append(-1);
  
  dimeBSPTree &bsptree = data->bspTree; // dxf-stuff I don't understand

  float r, g, b;
  dimeLayer::colorToRGB(colnum, r, g, b); // dxf-stuff I don't understand, probably not needed

  int first = 1;

  // points are stored in an array. the faces' vertices
  // are only indexes into this point array. we need to
  // know the number of points as 'vcnt' (VertexCouNT).

  pointArray = new VektorR3 [vcnt];

  // also obtain the min- and -max coordinates
  // for later layout optimization

  min_x = min_y = min_z = 1000000;
  max_x = max_y = max_z = -1000000;

  for (i = 0; i < vcnt; i++) {
    dimeVec3f tmp;
    bsptree.getPoint(i, tmp);

    if(min_x > tmp[0]) min_x = tmp[0];
    if(max_x < tmp[0]) max_x = tmp[0];

    if(min_y > tmp[1]) min_y = tmp[1];
    if(max_y < tmp[1]) max_y = tmp[1];

    if(min_z > tmp[2]) min_z = tmp[2];
    if(max_z < tmp[2]) max_z = tmp[2];

    pointArray[i] = VektorR3(tmp[0], tmp[1], tmp[2]);
  }

  double extend3D = -1;

  if((max_x - min_x) > extend3D) extend3D = max_x - min_x;
  if((max_y - min_y) > extend3D) extend3D = max_y - min_y;
  if((max_z - min_z) > extend3D) extend3D = max_z - min_z;
  
  bsptree.clear();

  // now we know all points which will be loaded. it comes
  // in handy that VRML also stores the vertices used in
  // the beginning of the file.

  // count the number of faces as 'facecnt'
  
  if (icnt) {

    // this is only intended to count the number of faces

    int *ptr = data->indices.arrayPointer(); // lib stuff I don't understand
    facecnt = 1; // 1 and not zero because the last face has no index == -1
    
    for (i = 0; i < icnt-1; i++) {
      int index = *ptr++;

      if (index < 0) // index == -1 codes: "now comes the next face"
	facecnt++;

    }

    // allocate the necessary number of faces 'facecnt'.
    // each side of the face might have a glue tab (adherent).
    // these could possibly be stored later. that's why
    // we allocate facecnt*MAX_POINTS_PER_FACE faces to
    // be on the save side.

    faceArray = new paperFace [facecnt*MAX_POINTS_PER_FACE];

    ptr = data->indices.arrayPointer();
    facecnt = 0;
    int pointIndex = 0;

    // initilize first face

    faceArray[facecnt].initialize();
    
    // now we load all the vertices of a face

    for (i = 0; i < icnt-1; i++) {
      
      // the 'index' is a vertex of the current face faceArray[facecnt].
      // 'index' points into the array pointArray[]

      int index = *ptr++;      

      if(index >= 0) {

	// if (index >= 0) we have to store it as a new vertex of the current face

	if(pointIndex >= MAX_POINTS_PER_FACE) {
	  cout << "error: Number of point exceeds static limit. Please adjust MAX_POINTS_PER_FACE and recompile." << endl;
	  exit(-1);
	} // if

	// why don't we simply store 'faceArray[facecnt].point[pointIndex] = index'? Because blender
	// might create points which are very close to one another. logically they should be the same
	// but numerically they may differ. So lets try to find a previously stored point be could use.
	// in most cases the routing getPointIndex() will only find a single match anyway.

	faceArray[facecnt].point[pointIndex] = getPointIndex(pointArray[index], vcnt, extend3D/1000.0);
	pointIndex++;
	faceArray[facecnt].no_points = pointIndex; // store number of vertices known so far

      } // if
      else {

	// if (index < 0) meants that the current face is finished
	// and we store the vertices of the next face

	pointIndex = 0;

	if(faceArray[facecnt].CheckColinear() == false) {
	  facecnt++;                                       // next vertices go into next face in faceArray
	  faceArray[facecnt].initialize();
	} // if
	else
	  cout << "WARNING: skiped colinear face" << endl;
      } // else
    }

    facecnt++; // this is necessary to store the last face in the list
               // because for the last face we get no index = -1 as
               // for the previous faces

    data->indices.makeEmpty(); // free memory. lib stuff I don't understand
  }

  epsilon = extend3D/10000.0;

  // loading faces is finished here

  cout << "faces loaded: " << facecnt << endl;
    
  // merge faces within the same polygon

  if(MERGE_FACES == true) {

    bool action;
    int  pass_counter = 0;

    //    cout << ">> merging faces pass ";

    do {
      action = false;
      
      //      cout << "#" << pass_counter++ << " " << flush;
      
      for(int i = 0; i < facecnt-1; i++) {	

	for(int j = i+1; j < facecnt; j++) {
	  if(faceArray[i].mergeFace(faceArray[j]) == true) {

	    //	    cout << "that was a merging success with faces #" << i << " and #" << j << endl;

	    action = true;
	    
	    if(j < facecnt-1) {

	      // fill the space of face j with the
	      // last face of the array and reduce
	      // array-size by 1

	      //	      cout << "copy face " << facecnt-1 << " into " << j << endl;
	      faceArray[j].copy(faceArray[facecnt-1]);
	    } // if

	    facecnt--;

	    //	    cout << "-----------------------------------------------" << endl;
	    
	  } // if
	  else {
	  } // else
	} // for

      } // for
              
      for(int r = 0; r < facecnt; r++) {
	if(faceArray[r].CheckColinear() == true) {
	  //	  cout << "delete colinear face" << endl;
	  //	  cout << "copy face " << facecnt-1 << " into " << r << endl;
	  faceArray[r].copy(faceArray[facecnt-1]);
	  facecnt--;
	  r--;
	} // if
      } // for
           	
    } while(action == true);

    cout << endl;
  } // if    

  if (lcnt) {
    cout << "error: line sets are not yet processed" << endl;

    /*
    fprintf(out,"IndexedLineSet {\n"
	    "coordIndex [\n");
    */
    int *ptr = data->lineIndices.arrayPointer();
    int linecnt = 1;
    
    for (i = 0; i < lcnt-1; i++) {
      int index = *ptr++;
      // fprintf(out, "%d,", index);
      if (index < 0) {
	linecnt++;
	//	if ((linecnt & 5) == 0) fprintf(out,"\n"); 
      }
    }
    // fprintf(out, "-1\n]\n}\n"); 
    data->lineIndices.makeEmpty(); // free memory
  }

  // remove colinear points on polygon

  for(int i = 0; i < facecnt; i++) {
    faceArray[i].removeColinearPoints();
  } // for
  
  for(int i = 0; i < facecnt; i++) {
    if(faceArray[i].CheckColinear() == true) {
      faceArray[i].copy(faceArray[facecnt-1]);
      facecnt--;
      i--;
    } // if
  } // for

  // number faces

  for(int i = 0; i < facecnt; i++) {
    faceArray[i].ID = i;
  } // for

  if(facecnt > 0) {
   faceArray[0].generateConnectedFaceGraph(facecnt);
  } // if
  // generate papercraft model

  ofstream* DXFFile = openDXFFile(filename);
  VektorR2 last_offset(0.0, 0.0), last_dir(1.0, 0.0);
  int      last_face = 0;
  double   extend, extend_x, extend_y; 

  int*     gen_face_array = new int [facecnt];
  int      gen_face_index = 0;

  VektorR2 normal_vec(0.55, 0.77);
  normal_vec.Normalize();
  
  int start_with_index = 0, current_index, counter = 0;
  int initial_facecnt = facecnt;

  cout << "faces to generate: " << initial_facecnt << endl;

  bool action;

  for(int offset = 3; offset < 4; offset++) {
    
    int start_off = 0;
    SUPPRESSED_ADHERENTS = 0;
    
    for(int i = 0; i < facecnt; i++) {
      for(int edge = 0; edge < faceArray[i].no_points; edge++) {
	faceArray[i].adherent_connected[edge] = false;
	faceArray[i].adherent_ID[edge] = -1;
	//faceArray[i].neighbor[edge] = NULL;
      } // for
      
      faceArray[i].drawn = false;
      faceArray[i].processed = false;
      faceArray[i].visited = false;
      faceArray[i].adherent = false;
    } // for
    
    cout << "trying offset: " << offset << endl;
    
    do {
      action = false;
      
      for(int index = offset; index < facecnt+offset; index++) {
	
	int i = index % facecnt;
	
	if(i < 0) i = facecnt + i;
	
	if(faceArray[i].drawn == false) {
	  
	  start_off++;
	  
	  calcExtend(facecnt, extend_x, extend_y);
	  
	  // cout << "EXTEND: " << extend_x << " " << extend_y << endl;
	  
	  if(extend_y < extend_x) {
	    last_offset = VektorR2(0.0, extend_y*1.05);
	  } // if
	  else{
	    last_offset = VektorR2(extend_x*1.05, 0.0);
	  } // else
	  
	  last_dir = VektorR2(1.0, 0.0);

	  // this is the priority queue for generateCutoutSheet

	  if(faceArray[i].generateCutoutSheet(DXFFile, faces,  last_offset, last_dir, 0, 
					      normal_vec, facecnt, extend3D, ALLOW_INTERSECTION, false) == true)
	    action = true;

	} // if
      } // for
    } while(action == true);
    
    cout << "start off " << start_off << endl;

    if(GENERATE_GLUE_TABS) {
      for(int i = 0; i < facecnt; i++) {
	for(int edge = 0; edge < faceArray[i].no_points; edge++)
	  faceArray[i].writeAdherent(DXFFile, extend3D, facecnt);
      } // for
    } // if
    
    if(GENERATE_GLUE_TABS)   
      cout << "suppressed glue tabs due to overlap: " << SUPPRESSED_ADHERENTS << endl;
    
    cout << "sum of enclosing angles (greater = better): " << paperFace::sumMinimalEnclosingAngles(facecnt) << endl;
  } // for

  closeDXFFile(DXFFile);

  delete[] gen_face_array;

  normalizeUvs(faces);

  root["faces"] = faces;
  std::cout << "JSON\n----\n" << root << std::endl;

}

static bool
convert_noext(dimeEntity::GeometryType type, dimeBSPTree &bsp,
	      dimeArray <int> &idx,
	      dimeArray <int> &lidx,
	      const dimeMatrix &matrix,
	      const int *indices,
	      const int numidx)
{
  int i;

  switch(type) {
  case dimeEntity::POLYGONS:
    for (i = 0; i < numidx; i++) {
      int index = indices[i];
      dimeVec3f tmp = vertexarray[index];
      matrix.multMatrixVec(tmp);
      idx.append(bsp.addPoint(dimeVec3f(tmp[0], tmp[1], tmp[2])));
    }
    idx.append(-1);
    break;
  case dimeEntity::LINES:
    for (i = 0; i < numidx; i++) {
      int index = indices[i];
      dimeVec3f tmp = vertexarray[index];
      matrix.multMatrixVec(tmp);
      lidx.append(bsp.addPoint(dimeVec3f(tmp[0], tmp[1], tmp[2])));
    }
    lidx.append(-1);
    break;
  default:
    break;
  }
  return true;
}

static bool
convert_ext(dimeEntity::GeometryType type, dimeBSPTree &bsp,
	    dimeArray <int> &idx,
	    const dimeVec3f &extrusion,
	    const dimeMatrix &matrix,
	    const int *indices,
	    const int numidx)
{
  int i,index;
  
  tmparray.setCount(0);

  switch(type) {
  case dimeEntity::POLYGONS:
    for (i = 0; i < numidx; i++) {
      dimeVec3f tmp = vertexarray[indices[i]];
      matrix.multMatrixVec(tmp);
      index = bsp.addPoint(dimeVec3f(tmp[0], tmp[1], tmp[2]));
      tmparray.append(index);
      idx.append(index);
    }
    idx.append(-1);
    for (i = 0; i < numidx; i++) {
      dimeVec3f tmp = vertexarray[indices[i]];
      tmp += extrusion;
      matrix.multMatrixVec(tmp);
      index = bsp.addPoint(dimeVec3f(tmp[0], tmp[1], tmp[2]));
      tmparray.append(index);
      idx.append(index);
    }
    idx.append(-1);
    
    for (i = 0; i < numidx-1; i++) {
      idx.append(tmparray[i]);
      idx.append(tmparray[i+1]);
      idx.append(tmparray[i+1+numidx]);
      idx.append(tmparray[i+numidx]);
      idx.append(-1);
    }
    idx.append(tmparray[i]);
    idx.append(tmparray[0]);
    idx.append(tmparray[numidx]);
    idx.append(tmparray[i+numidx]);
    idx.append(-1);
    break;
  case dimeEntity::LINES:
    for (i = 0; i < numidx; i++) {
      dimeVec3f tmp = vertexarray[indices[i]];
      matrix.multMatrixVec(tmp);
      index = bsp.addPoint(dimeVec3f(tmp[0], tmp[1], tmp[2]));
      tmparray.append(index);
    }
    for (i = 0; i < numidx; i++) {
      dimeVec3f tmp = vertexarray[indices[i]];
      tmp += extrusion;
      matrix.multMatrixVec(tmp);
      index = bsp.addPoint(dimeVec3f(tmp[0], tmp[1], tmp[2]));
      tmparray.append(index);
    }
    for (i = 0; i < numidx-1; i++) {
      idx.append(tmparray[i]);
      idx.append(tmparray[i+1]);
      idx.append(tmparray[i+numidx+1]);
      idx.append(tmparray[i+numidx]);
      idx.append(-1);
    }
    break;
  default:
    break;
  }
  return true;
}


static int
get_color_index(dimeEntity *entity)
{
  int colnum = entity->getColorNumber();
  if (colnum == 256) {
    const dimeLayer *layer = entity->getLayer();
    colnum = layer->getColorNumber();
  }
  return colnum;
}

static bool 
my_dime_callback(const dimeState * const state, 
		 dimeEntity *entity, void *)
{
  if (entity->isDeleted()) return true;

  static int blockcolidx = 7;
  
  if (entity->typeId() == dimeBase::dimeBlockType) {
    blockcolidx = entity->getColorNumber();
    if (blockcolidx == 256) {
      blockcolidx = entity->getLayer()->getColorNumber();
    }
  }

  dimeMatrix matrix = state->getMatrix();
  
  indexarray.setCount(0);
  vertexarray.setCount(0);

  float thickness = 0.0;
  dimeVec3f extrusionDir(0,0,1);

  dimeEntity::GeometryType type = 
    entity->extractGeometry(vertexarray, 
			    indexarray,
			    extrusionDir,
			    thickness);
  
  if (type == dimeEntity::NONE) return true;

  // normalize data before converting
  int n = vertexarray.count();
  int numidx = indexarray.count();
 
  if (numidx == 0) {
    if (n == 0) return true; // nothing to do
    for (int i = 0; i < n; i++) indexarray.append(i);
  }
  numidx = indexarray.count();
  if (numidx && indexarray[numidx-1] >= 0) 
    indexarray.append(-1);
  
  int id = entity->typeId();

  if (id == dimeBase::dimeTraceType ||
      id == dimeBase::dimeSolidType ||
      id == dimeBase::dimeCircleType) {
    if (extrusionDir != dimeVec3f(0,0,1)) {
      dimeMatrix m;
      dimeEntity::generateUCS(extrusionDir, m);
      matrix.multRight(m);
    }
  }
  
  else if ((type == dimeEntity::POLYGONS || 
	    id == dimeBase::dimeCircleType ||
	    id == dimeBase::dimeArcType ||
	    id == dimeBase::dimeEllipseType) && 
	   (thickness == 0.0 && extrusionDir != dimeVec3f(0,0,1))) {
    dimeMatrix m;
    dimeEntity::generateUCS(extrusionDir, m);
    matrix.multRight(m);
  }
  
  int colidx = get_color_index(entity);
  
  if (colidx == 0) { // BYBLOCK
    colidx = blockcolidx;
  }
  if (colidx < 0) colidx = -colidx;
  
  colidx--; // to use lookup table 
  
  if (colidx < 0 || colidx >= 255) { // safety check
    //    sim_trace("dime color changed from %d to 6.\n", colidx, 6);
    colidx = 6; // should never happen
  }
  
  dimeBSPTree &bsp = facesets[colidx].bspTree;
  dimeArray <int> &idx = facesets[colidx].indices;
  dimeArray <int> &lidx = facesets[colidx].lineIndices;
  //  dimeArray <dimeVec3f> &pts = facesets[colidx].points;

  // FIXME: set extrusionDir to 0,0,1 ?
  extrusionDir = dimeVec3f(0,0,1);
  extrusionDir *= thickness;
  int *indices = indexarray.arrayPointer();
  numidx = indexarray.count();
  int i = 0;
  while (i < numidx) {
    int start = i;
    while (indices[i] >= 0) i++; // search for next -1
    int cnt = i-start;
    if (cnt) {
      if (thickness == 0.0) {
	convert_noext(type, bsp, idx, lidx, matrix, 
		      &indices[start], cnt);
      }
      else {
	convert_ext(type, bsp, idx, extrusionDir, 
		    matrix, &indices[start], cnt);
      }
    }
    i++; // skip -1
  }
  return true;
}

int main(int argc, char **argv)
{
  int  c;
  int  digit_optind = 0;
  bool error = false;
  char* arg_ptr;


  STRATEGY = 5;
  MERGE_FACES = true;
  ALLOW_INTERSECTION = false;
  GENERATE_GLUE_TABS = true;
  CONNECT_FACES = true;
  NUMBER_FACES = false;
  FORCE_ADHERENTS = false;

  indexSplitFaces = 0; // index to faces that trigger 
                       // to generate the next continuous
                       // cutout sheet

  while (1)
    {
      int this_option_optind = optind ? optind : 1;
      int option_index = 0;

      static struct option long_options[] =
	{
	  {"nomerge", 1, 0, 'm'},
	  {"divide", 0, 0, 'd'},
	  {"overlap", 1, 0, 'o'},
	  {"hide", 0, 0, 'h'},
	  {"strategy", 1, 0, 's'},
	  {"number", 1, 0, 'n'},
	  {"force", 1, 0, 'f'},
	  {"split", 1, 0, 'p'},
	  {"help", 1, 0, '?'},
	  {0, 0, 0, 0}
	};

      c = getopt_long (argc, argv, "nfmdohs:p:?",
                       long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
        case 'm':
	  MERGE_FACES = false;
          break;

	case 's':
	  STRATEGY = atoi(optarg);
	  if((STRATEGY < 0) || (STRATEGY > 7))
	    error = true;
	  break;

        case 'd':
	  CONNECT_FACES = false;
          break;

        case 'f':
	  FORCE_ADHERENTS = true;
          break;

        case 'o':
	  ALLOW_INTERSECTION = true;
          break;

        case 'h':
	  GENERATE_GLUE_TABS = false;
          break;

        case 'n':
	  NUMBER_FACES = true;
          break;

        case 'p':
	  arg_ptr = optarg;

	  do {
	    arraySplitFaces[indexSplitFaces++] = atoi(arg_ptr);
	    if(indexSplitFaces > MAX_NO_SPLIT_FACES) {
	      cout << "ERROR: number of split faces exceeded. Please change MAX_NO_SPLIT_FACES and recompile" << endl;
	      exit(-1);
	    } // if
	    while((strlen(arg_ptr) > 0) && (*arg_ptr != ','))
	      arg_ptr++;
	    if(*arg_ptr == ',') arg_ptr++;
	  } while(strlen(arg_ptr) > 0);
          break;

        case '?':
	  error = true;
          break;

        default:
          printf ("?? getopt lieferte Zeichcode 0%o zurück ??\n", c);
	  error = true;
        }
    }
  /*
  if (optind < argc)
    {
      cout << "non option element of ARGV: ";
      while (optind < argc)
	printf ("%s ", argv[optind++]);
      printf ("\n");
    }
  */
  SHOW_DEBUG = false;

  debugFile = new ofstream("debug.dat");

  if ((argc < 3) || (error)) {
    cout << "Usage: " << argv[0] << " [options] infile3D.dxf outfile2D.dxf" << endl << endl;
    cout << "convert a polygonal 3D object into a 2D cut-out sheet for" << endl;
    cout << "producing a paper model of the object using glue and scissors" << endl;
    cout << "Copyright by Thomas Haenselmann <givenname@familyname.de>" << endl << endl;
    cout << "Options:" << endl;
    cout << endl;
    cout << " -m, --nomerge        no merging of faces into single polygon" << endl;
    cout << " -n, --number         print face numbers" << endl;
    cout << " -d, --divide         draw each face separate" << endl;
    cout << " -o, --overlap        allow overlapping faces in cut-out sheet" << endl;
    cout << " -h, --hide           hide glue tabs" << endl;
    cout << " -f, --force          force glue tabs, even if intersecting faces" << endl;
    cout << " -p, --split 8,17     face number 8 and 17 get disconnected from the rest" << endl;
    cout << "                      (use -n to see face numbers in 2D DXF file)" << endl;
    cout << " -s, --strategy 0..5  0: draw smallest polygon first / 1: draw largest first" << endl;
    cout << "                      2: as ordered in file / 3: keep adjacent faces continuous" << endl;
    cout << "                      4: stretch 2D layout wide / 5: keep layout dense" << endl;
    cout << " -?, --help           display this text" << endl;
  }
  else {
    init();
    if (!convertDxf2Wrl(argv[argc-2], argv[argc-1])) {
      fprintf(stderr,"Error while converting.\n");
    }
    end();
  }
}
