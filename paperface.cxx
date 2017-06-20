
#include <iomanip>
#include <stdlib.h>

#include "paperface.h"
#include "utilities.h"

#include "json/json.h"

extern double     epsilon;
extern ofstream*  DXFFile;
extern VektorR3*  pointArray;
extern paperFace* faceArray;
extern int        global_adherent_counter;
extern int        global_counter;
extern ofstream*  debugFile;
extern bool       SHOW_DEBUG;
extern int 	  SUPPRESSED_ADHERENTS;
extern int 	  STRATEGY;
extern bool       NUMBER_FACES;
extern bool       FORCE_ADHERENTS;
extern int        indexSplitFaces;
extern int        arraySplitFaces[];

bool paperFace::generateCutoutSheet(ofstream* dxfFile, Json::Value &faces_json, VektorR2 offset, VektorR2 dir, int start_with_index, 
				    VektorR2 normal, int no_faces, double extend, 
				    bool allow_intersection, bool allow_to_split = false)
{
  paperFace** neighborArray = new paperFace* [no_faces];
  int         neighbor_index = 0;

  double     FONT_SIZE = 50;
  paperFace* myFather, *optimal_face = this;
  bool       overlap;
  int        edge_index, neighbor_edge_index, new_neighbor;

  if(optimal_face->writeToDXF(dxfFile, faces_json, offset, dir, start_with_index, 
     normal, no_faces, allow_intersection, allow_to_split) == true) {
    
    if(NUMBER_FACES)
      writeNumber(dxfFile, optimal_face->getMidPoint()-VektorR2((log(global_counter+1)/log(10.0)+1)*extend/FONT_SIZE, 
      -0.5*extend/FONT_SIZE), VektorR2(1.0, 0.0), VektorR2(0.0, -1.0), extend/FONT_SIZE, global_counter++); 
  } // if

  // no face has been processed so far

  for(int i = 0; i < no_faces; i++)
    faceArray[i].visited = false;

  do {

    optimal_face->visited = true;    
    allow_to_split = true;

    // for all faces

    for(int j = 0; j < no_faces; j++) {

      // the considered face is neither me, nor has it been drawn so far
      
      if((&(faceArray[j]) != optimal_face) && (faceArray[j].drawn == false)) {

	// the considered face also has a common edge with me
	
	if(optimal_face->sharesSomeEdge(faceArray[j], edge_index, neighbor_edge_index) == true) {
	  
	  // chances are that this face is in the queue already. if so, do 
	  // not add it a second time but rather update its cost
	  
	  for(new_neighbor = 0; new_neighbor < neighbor_index; new_neighbor++)
	    if(neighborArray[new_neighbor] == &(faceArray[j])) break;
	  
	  if(new_neighbor == neighbor_index) {
    
	    // we did not find this face so far
	    // thus it will be stored
	    
	    neighborArray[neighbor_index] = &(faceArray[j]);
	    neighborArray[neighbor_index]->my_father = optimal_face;
	    neighbor_index++;
	  } // if

	  if(neighborArray[new_neighbor]->calcCostForCutoutSheet(STRATEGY, 
	     optimal_face->getMidPoint(), no_faces) == true) {
	    
	    // true means: the cost value improved. So store the new father-info
	    	    
	    neighborArray[new_neighbor]->my_father = optimal_face;
	    neighborArray[new_neighbor]->visited = false;
	          
	  } // if
	  
	} // if  
      } // if
    } // for
    
    double best_cost;
    int    best_index;

    // essentially, this loop will try to draw a
    // single neighbor from the queue before exiting
          
    do{

      best_index = -1;
      
      // determine the optimal cost neighbor
      
      for(int i = 0; i < neighbor_index; i++) {
	
	if((neighborArray[i]->drawn == false) && (neighborArray[i]->visited == false)) {
	  
	  if((best_index == -1) || (best_cost > neighborArray[i]->getCost())) {
	    best_cost = neighborArray[i]->getCost();
	    best_index = i;     
	  } // if      
	} // if
      } // for
      
      if(best_index == -1) {
	delete[] neighborArray;

	return true;
      } // if
      
      // PLEASE REMEMBER: Only the father's projection is
      //                  calculated so far. The optimal_face,
      //                  the son is to be rendered below.
      //                  Actually this os why the son has been
      //                  chosen above.

      myFather = neighborArray[best_index]->my_father;
      optimal_face = neighborArray[best_index];

      optimal_face->sharesSomeEdge(*myFather, edge_index, neighbor_edge_index);
      optimal_face->visited = true;
      offset = myFather->projection[neighbor_edge_index];      
      dir = myFather->projection[myFather->index_save(neighbor_edge_index+1)] - myFather->projection[neighbor_edge_index];

      start_with_index = edge_index;    
      normal = myFather->calcNormal(neighbor_edge_index);
           
      overlap = false;
      
      if(optimal_face->writeToDXF(dxfFile, faces_json, offset, dir, start_with_index, 
				  normal, no_faces, allow_intersection, allow_to_split) == true) {
	
	if(NUMBER_FACES)
	  writeNumber(dxfFile, optimal_face->getMidPoint()-VektorR2((log(global_counter+1)/log(10.0)+1)*
	     extend/FONT_SIZE, -0.5*extend/FONT_SIZE), VektorR2(1.0, 0.0), VektorR2(0.0, -1.0), 
	     extend/FONT_SIZE, global_counter++); 
      } // if
      else overlap = true;

    } while(overlap == true);
        
  } while(true);
  
} // paperFace::generateCutoutSheet

bool paperFace::calcCostForCutoutSheet(int strategy, VektorR2 my_midpoint, int no_faces)
{
  VektorR2   last_direction, temp_direction;
  double     temp_cost;

  // implement your strategy here
  
  switch(strategy) {
    
  case 0:
    // follow smallest polygon first
      temp_cost = estimateExtend();
    break;
    
  case 1:
    // follow largest polygon first
      temp_cost = estimateExtend();
    break;
    
  case 2:
    // random
    temp_cost = rand();
    break;
    
  case 3:
    // continue as straight as possible
    
    temp_direction = getMidPoint() - my_midpoint;
    temp_direction.Normalize();
    
    temp_cost = last_direction*temp_direction;
    break;
    
  case 4:
    // as extended as possible

      temp_cost = getCumulatedDistance(no_faces);
    break;
    
  case 5:
    // as close as possible

    temp_cost = getCumulatedDistance(no_faces);
    break;
    
  case 6:
    /*
    // long edges first
    
    length_shared_edge = (projection[edge_index[i]] -
    projection[index_save(edge_index[i]+1)]).Norm2();
    
    if(length_shared_edge < max_edge_length) {
    max_edge_length = length_shared_edge;
    best_neighbor = i;
    } // if
    */
    break;

  case 7:
    // face with largest minimal neighboring angle first

    temp_cost = 2.0*M_PI - getMinimalEnclosingAngle(no_faces);
    // temp_cost = getMinimalEnclosingAngle(no_faces);

    break;
    
  default:
    cout << "ERRORO: strategy not implemented" << endl;
    exit(-1);
    
  } // switch

  if(temp_cost < cost) {
    cost = temp_cost;
    return true;
  } // if
  return false;
} // paperFace::calcCostForCutoutSheet

double paperFace::getMinimalEnclosingAngle(int no_faces)
{
  double   min_angle = 100000;
  VektorR2 v1, v2, end_1, end_2;

  // for all faces which have been drawn so far

  for(int f = 0; f < no_faces; f++) {
    if((faceArray[f].drawn == true) && (&(faceArray[f]) != this)) {
      
      // for all edges of this face
      
      for(int i = 0; i < no_points; i++) {

	// for all edges of the other face

	for(int j = 0; j < faceArray[f].no_points; j++) {
	
	  // if two polygonal faces meet in one common vertex, then 
	  // there are four edges starting from that vertex, two from
	  // each face respectively. As a consequence, four pairs
	  // of edges can be build:

	  if((projection[i] - faceArray[f].projection[j]).Norm2() < epsilon) { // this is the same point

	    for(int c = 0; c < 4; c++) {
	      switch(c) {

	      case 0:
		end_1 = projection[index_save(i+1)];
		end_2 = faceArray[f].projection[faceArray[f].index_save(j+1)];
		break;

	      case 1:
		end_1 = projection[index_save(i+1)];
		end_2 = faceArray[f].projection[faceArray[f].index_save(j-1)];
		break;

	      case 2:
		end_1 = projection[index_save(i-1)];
		end_2 = faceArray[f].projection[faceArray[f].index_save(j+1)];
		break;

	      case 3:
		end_1 = projection[index_save(i-1)];
		end_2 = faceArray[f].projection[faceArray[f].index_save(j-1)];
		break;

	      } // switch

	      if((end_1 - end_2).Norm2() > epsilon) {

		v1 = end_1 - projection[index_save(i)];
		double length_1 = v1.Normalize();

		v2 = end_2 - projection[index_save(i)];
		double length_2 = v2.Normalize();
		
		double angle = fabs(acos(v1*v2));

		if(angle < min_angle) {
		  min_angle = angle;
		  //		  cout << "updated angle: " << angle*360.0/(2.0*M_PI) << endl;
		} // if
	      } // if
	    } // for
	  } // if

	} // for
      } // for
    } // if
  } // for

  //  cout << "returned angle: " << min_angle*360.0/(2.0*M_PI) << endl << endl;
  return min_angle;
} // paperFace::getMinimalEnclosingAngle

double paperFace::sumMinimalEnclosingAngles(int no_faces)
{
  double angle, enclosing_angle = 0.0;
  long   counter = 0;

  for(int f = 0; f < no_faces; f++) {
    angle = faceArray[f].getMinimalEnclosingAngle(no_faces);
    //    cout << "face #" << faceArray[f].ID << ": " << angle*360.0/(2.0*M_PI) << "°" << endl;
    if(angle <= 2*M_PI) {
      enclosing_angle += angle;
      counter++;
    } // if
  } // for

  return 360.0*(enclosing_angle/(double)counter)/(2*M_PI);
} // paperFace::sumMinimalEnclosingAngles

void paperFace::generateConnectedFaceGraph(int no_faces)
{
  paperFace* p_neighbor;
  int        neighbor_edge_index;

  // connect all faces to their neighbors by means of a depth-first traversal

  for(int i = 0; i < no_faces; i++) {

    if(&(faceArray[i]) != this) {

      for(int edge_index = 0; edge_index < no_points; edge_index++) {

	if(neighbor[edge_index] == NULL) { // we did not process this edge before

	  if(sharesThisEdge(faceArray[i], edge_index, neighbor_edge_index, true) == true) {
	    //cout << "graph: face.ID=" << this << " connected to face: " << &(faceArray[i]) << "In edge index " << edge_index << endl;
	    neighbor[edge_index] = &(faceArray[i]);
	    faceArray[i].generateConnectedFaceGraph(no_faces);
	  } // if       
	} // if
      } // for
    } // if
  } // for
} // paperFace::generateConnectedFaceGraph

bool paperFace::intersectsOtherFace(int no_faces)
{
  int sharesVertexCounter = 0;
    
  for(int i = 0; i < no_faces; i++) {
    
    if((&(faceArray[i]) != this) && (faceArray[i].drawn == true)) {
      
      for(int j = 0; j < faceArray[i].no_points; j++) {

	sharesVertexCounter = 0;
      
	for(int k = 0; k < no_points; k++) {
	  
	  if(IntersectLineLine(projection[k], projection[index_save(k+1)], faceArray[i].projection[j], 
			       faceArray[i].projection[faceArray[i].index_save(j+1)]) == true) {
	    
	    //	    cout << "   face is intersecting other faces" << endl;

	    return true;
	  } // if	  
	  
	  // no point on this face must be contained in another face
	  
	  VektorR2 inner_point, dir;

	  dir = VektorR2(0.30, 0.70);
	  dir.Norm2();

	  inner_point = projection[index_save(k)];
	  
	  if(faceArray[i].PointInside(inner_point, dir) == true) {
	    return true;
	  } // if
	  
	  /*
	  if(EdgeOnEdge(k, faceArray[i].projection[faceArray[i].index_save(j)],
			faceArray[i].projection[faceArray[i].index_save(j+1)]) == true) {
	    sharesVertexCounter++;
	  } // if
	  */
	} // for
	/*
	if(sharesVertexCounter > 2)
	  return true;
	*/
      } // for
      

      
    } // if
  } // for

  return false;
} // paperFace::intersectsOtherFace

bool paperFace::writeToDXF(ofstream* DXFFile, Json::Value & faces_json, VektorR2 offset, VektorR2 dir, int start_with_index, 
			   VektorR2 normal, int no_faces, bool allow_intersection, bool allow_to_split)
{
  VektorR3 v_old, v_new, temp;
  VektorR2 dir_rot, turtle;
  VektorR2 assumed_position;
  double   angle, length;
  double   dist_1, dist_2;

  Json::Value face_json;
  face_json["id"] = ID;
  Json::Value vtx(Json::arrayValue);
  get3DVertices(0, vtx);

  face_json["vtx"] = vtx;
  VektorR3 normalVec = calc3DNormal();
  normalVec.Normalize();

  Json::Value normalJson;
  normalJson["x"] = normalVec[0];
  normalJson["y"] = normalVec[1];
  normalJson["z"] = normalVec[2];

  face_json["normal"] = normalJson;

  processed = true; // this face has been considered

  dir.Normalize();
  turtle = offset;

  v_old = pointArray[point[index_save(start_with_index+1)]] - pointArray[point[index_save(start_with_index+0)]];
  length = v_old.Normalize();
  v_old.Normalize();

  for(int i = 1+start_with_index; i <= no_points+start_with_index; i++) {

    projection[index_save(i-1)] = turtle;

    turtle = turtle + length*dir;

    v_new = pointArray[point[index_save(i+1)]] - pointArray[point[index_save(i)]];
    length = v_new.Normalize();

    angle = acos(v_old*v_new);
    // cout << "winkel: " << angle*360.0/(2.0*M_PI) << endl;

    v_old = v_new;

    dir_rot[0] = cos(angle)*dir[0] - sin(angle)*dir[1];
    dir_rot[1] = sin(angle)*dir[0] + cos(angle)*dir[1];

    if(i == 1+start_with_index) {


      // we want to point into the opposite
      // direction of the previous face
      
      if((dir_rot * normal) > 0.0) {
	angle *= -1.0;
	
	dir_rot[0] = cos(angle)*dir[0] - sin(angle)*dir[1];
	dir_rot[1] = sin(angle)*dir[0] + cos(angle)*dir[1];
      } // if 
    } // if
    else {

      // in case of concave objects we may have rotated into the wrong direction
            
      assumed_position = turtle + length*dir_rot;
      
      dist_1 = (assumed_position - projection[index_save(i-2)]).Norm2();
      dist_2 = (pointArray[point[index_save(i+1)]] - pointArray[point[index_save(i-2)]]).Norm2();
      
      if(fabs(dist_1 - dist_2) > epsilon) {
	
	double save_dist = fabs(dist_1 - dist_2);
	/*
	  cout << ">>> distances do not match: i=" << i << endl;
	  cout << "   -> distance 2D: " << dist_1 << endl;
	  cout << "   -> distance 3D: " << dist_2 << endl;
	  
	  cout << "---- 2D coordinates" << endl;
	  printProjection(start_with_index);
	  cout << "assumed pos: " << assumed_position[0] << ", " << assumed_position[1] << endl;
	  cout << "turtle     : " << turtle[0] << ", " << turtle[1] << endl;
	  cout << "length     : " << length << endl;
	  cout << "dir        : " << dir_rot[0] << ", " << dir_rot[1] << endl;
	  cout << "---- 3D coordinates" << endl;
	  print3DVertices(start_with_index);
	*/
	angle *= -1.0;
	
	dir_rot[0] = cos(angle)*dir[0] - sin(angle)*dir[1];
	dir_rot[1] = sin(angle)*dir[0] + cos(angle)*dir[1];
	
	// test again!
	
	// cout << "    test again:" << endl;
	
	assumed_position = turtle + length*dir_rot;
	
	dist_1 = (assumed_position - projection[index_save(i-2)]).Norm2();
	dist_2 = (pointArray[point[index_save(i+1)]] - pointArray[point[index_save(i-2)]]).Norm2();
	
	if(fabs(dist_1 - dist_2) > epsilon) {
	  /*
	    cout << ">>> distances still do not match: i=" << i << endl;
	    cout << "   -> distance 2D: " << dist_1 << endl;
	    cout << "   -> distance 3D: " << dist_2 << endl;
	    cout << "assumed pos: " << assumed_position[0] << ", " << assumed_position[1] << endl;
	    cout << "turtle     : " << turtle[0] << ", " << turtle[1] << endl;
	    cout << "length     : " << length << endl;
	    cout << "dir        : " << dir_rot[0] << ", " << dir_rot[1] << endl;
	  */
	  if(save_dist < fabs(dist_1 - dist_2)) {
	    
	    // actually the first direction was still better
	    // so use old direction
	    
	    angle *= -1.0;
	    
	    dir_rot[0] = cos(angle)*dir[0] - sin(angle)*dir[1];
	    dir_rot[1] = sin(angle)*dir[0] + cos(angle)*dir[1];
	  } // if
	  
	} // if
	else {
	  // cout << "    everything ok!" << endl;
	} // else
      } // if
      
    } // else
    
    dir = dir_rot;
   
  } // for

  // check, whether the new polygon intersects a drawn one or whether
  // it has more than one common vertex with a previously drawn
  // face. This must not happen because then we can not add the
  // adherents anymore
  
  for(int i = 0; i < indexSplitFaces; i++) {

    if(arraySplitFaces[i] == ID) {

      if(allow_to_split == true) {
	clearConnectedFlags();
	drawn = false;
	return drawn;
      } // if
    } // if
  } // for
  

  if(((allow_intersection == false) && (intersectsOtherFace(no_faces) == true))) {
    clearConnectedFlags();
    drawn = false;
    return drawn;	    	    
  } // if

  (*DXFFile) << "  0\nPOLYLINE\n 70\n     1" << endl;

  Json::Value uvs(Json::arrayValue);

  for(int i = start_with_index; i < no_points+start_with_index; i++) {
    
    //    if(IsInnerEdge(index_save(i), index_save(i+1)) == true) {
    if(false) {

      // this is an inner edge

      writePoint(DXFFile, projection[index_save(i)]);         
      (*DXFFile) << "  0\nSEQEND" << endl;
      (*DXFFile) << "  0\nPOLYLINE\n 70\n     0" << endl;
    } // if
    else 
      writePoint(DXFFile, projection[index_save(i)]);   
      Json::Value uv;
      uv["u"] = Json::Value((int)projection[index_save(i)][0]);
      uv["v"] = Json::Value((int)projection[index_save(i)][1]);
      uvs.append(uv);
  } // for

  face_json["uvs"] = uvs;
    
  (*DXFFile) << "  0\nSEQEND" << endl;

  Json::Value neighbors_json(Json::arrayValue);

  // Neighbors
  for  ( int i = 0; i < MAX_POINTS_PER_FACE; i++ ) {
    paperFace* current_neighbor = neighbor[i];
    if (current_neighbor) {
        Json::Value neighbor_json;
        neighbor_json["id"] = current_neighbor->ID;
        Json::Value vertex_indexes_json(Json::arrayValue);
        vertex_indexes_json.append(i);
        vertex_indexes_json.append(i == 3 ? 0 : i + 1);
        neighbor_json["vertexIndex"] = vertex_indexes_json;
        neighbors_json.append(neighbor_json);
    }
  }

  face_json["neighbors"] = neighbors_json;



  faces_json.append(face_json);

  drawn = true;
  return drawn;
} // paperFace::writeToDXF

void paperFace::writeAdherent(ofstream* DXFFile, double extend, int& facecnt)
{
  VektorR2 turtle, dir, dir_rot, normal;
  VektorR2 writing_direction_x, write_offset;
  VektorR2 offset;
  double   length, angle, height;
  double   sign = -1.0;

  // return;

  for(int vertex = 0; vertex < no_points; vertex++) {
    //    if((adherent_connected[vertex] == false) && (IsInnerEdge(index_save(vertex), index_save(vertex+1)) == false)) {
    //    if(adherent_connected[vertex] == false) {
    if((adherent_connected[vertex] == false) && (sharesProjectedEdge(vertex, facecnt) == false) &&
       (IsInnerEdge(index_save(vertex), index_save(vertex+1)) == false)) {

      adherent_connected[vertex] = true;

      height = extend*ADHERENT_SIZE;
      int attempt = 0;
      
      do {
	faceArray[facecnt].initialize();
	
	normal = calcNormal(vertex);
	
	dir = projection[index_save(vertex+1)] - projection[index_save(vertex)];      
	writing_direction_x = dir;
	length = dir.Normalize();
	
	if(length < 2.0*sin(45.0/(360.0/(2.0*M_PI)))*height)
	  height = sin(45.0/(360.0/(2.0*M_PI)))*(length/1.0);
	
	turtle = projection[index_save(vertex+1)];     
	write_offset = turtle;
	
	// store adherent as face for later intersection test
	faceArray[facecnt].projection[0] = turtle;
	
	angle = sign*135.0/(360.0/(2.0*M_PI));
	
	dir_rot[0] = cos(angle)*dir[0] - sin(angle)*dir[1];
	dir_rot[1] = sin(angle)*dir[0] + cos(angle)*dir[1];
	
	if(dir_rot*normal > 0.0) {
	  sign *= -1;
	  angle *= -1;
	  
	  dir_rot[0] = cos(angle)*dir[0] - sin(angle)*dir[1];
	  dir_rot[1] = sin(angle)*dir[0] + cos(angle)*dir[1];
	  
	} // if
	
	dir = dir_rot;
	
	turtle = turtle + dir*height;      
	// store adherent as face for later intersection test
	faceArray[facecnt].projection[1] = turtle;
	
	angle = sign*45.0/(360.0/(2.0*M_PI));
	
	dir_rot[0] = cos(angle)*dir[0] - sin(angle)*dir[1];
	dir_rot[1] = sin(angle)*dir[0] + cos(angle)*dir[1];
	dir = dir_rot;
	
	turtle = turtle + dir*(length-2.0*(sin(45.0/(360.0/(2.0*M_PI)))*height));
	// store adherent as face for later intersection test
	faceArray[facecnt].projection[2] = turtle;
	
	angle = sign*45.0/(360.0/(2.0*M_PI));
	
	dir_rot[0] = cos(angle)*dir[0] - sin(angle)*dir[1];
	dir_rot[1] = sin(angle)*dir[0] + cos(angle)*dir[1];
	dir = dir_rot;
	
	turtle = turtle + dir*height;      
	// store adherent as face for later intersection test
	faceArray[facecnt].projection[3] = turtle;
	
	faceArray[facecnt].processed = true;
	faceArray[facecnt].drawn = false;
	faceArray[facecnt].adherent = true;
	faceArray[facecnt].no_points = 4;
	
	height /= 2.0;

      } while((FORCE_ADHERENTS == false) && (attempt++ < 3) && (faceArray[facecnt].intersectsOtherFace(facecnt) == true));

      if(attempt >= 3) {
	SUPPRESSED_ADHERENTS++;
	return;
      } // if

      /*
	
	for(int repeat = 0; repeat < 2; repeat++) {
	
	if(repeat == 0)
	angle = getEnclosedAngle(vertex, facecnt, offset, false);
	if(repeat == 1)
	angle = getEnclosedAngle(vertex+1, facecnt, offset, true);
	
	if(angle <= 360.0) {
	
	
	//  writeLine(DXFFile, projection[index_save(vertex+1,no_points)]+0.05*VektorR2(1, 1), this_vertex+0.05*VektorR2(1, 1));
	//  writeLine(DXFFile, neighbor_vertex+0.05*VektorR2(1, 1), this_vertex+0.05*VektorR2(1, 1));
	
	//	cout << "DEBUG: angle=" << dir_1*dir_2 << endl;
	
	double length_adherent = writing_direction_x.Normalize();
	
	if(normal*VektorR2(-1.0*writing_direction_x[1], writing_direction_x[0]) < 0.0) {
	normal = VektorR2(-1.0*writing_direction_x[1], writing_direction_x[0]);
	//	  offset = projection[vertex];
	writing_direction_x *= -1.0;
	} // if
	else {
	normal = -1.0*VektorR2(-1.0*writing_direction_x[1], writing_direction_x[0]);
	//	  offset = projection[vertex];
	} // else
	
	  //	  offset += 0.01*VektorR2(rand() % 100, rand() % 100);
	  
	  //	offset += (0.95-0.6)*height*normal + (0.4)*height*writing_direction_x;
	  
	  cout << "DEBUG: " << angle << endl;
	  
	  writeNumber(DXFFile, offset, writing_direction_x, normal, 0.2*height, angle);
	  return;
	  
	  if(fabs(angle) > 45) {
	  cout << "this adherent should not be drawn" << endl;
	  return;
	  } // if
	  
	  } // if
	  } // for
      */
      //      return;
      /*
	if(faceArray[facecnt].intersectsOtherFace(facecnt) == true) {
	cout << "this adherent should not be drawn" << endl;
	return;
	} // if
      */
      // draw the adherent
      
      (*DXFFile) << "  0\nPOLYLINE\n 70\n     1" << endl;
      
      for(int i = 0; i < faceArray[facecnt].no_points; i++) {
	writePoint(DXFFile, faceArray[facecnt].projection[i]);
      } // for
      
      (*DXFFile) << "  0\nSEQEND" << endl;
      
      // don't store adherent at the moment
      //      facecnt++;

      // allow user to identify adherent with its peer
      
      double length_adherent = writing_direction_x.Normalize();
      paperFace* p_connected_face = NULL;
      int        connected_edge;
      long       my_number = -1;
      
      my_number = getConnectedFace(vertex, facecnt);
      //      cout << "my_number = " << my_number << endl;

      if(my_number == -1) {
	my_number = adherent_ID[index_save(vertex)] = global_adherent_counter++ + 1; 
      } // if

      if(normal*VektorR2(-1.0*writing_direction_x[1], writing_direction_x[0]) < 0.0) {
	normal = VektorR2(-1.0*writing_direction_x[1], writing_direction_x[0]);
	offset = 0.5*(projection[index_save(vertex+1)]+projection[index_save(vertex)]);
	writing_direction_x *= -1.0;
      } // if
      else {
	normal = -1.0*VektorR2(-1.0*writing_direction_x[1], writing_direction_x[0]);
	offset = 0.5*(projection[index_save(vertex)]+projection[index_save(vertex+1)]);
      } // else

      offset += (0.95-0.6)*height*normal - (1.0)*height*writing_direction_x;
      
      while(length_adherent < height*0.2*0.4*1.5*4*7)
	height *= 0.95;
      	
      writeNumber(DXFFile, offset, writing_direction_x, normal, 0.4*height, my_number);
      
    } // if

  } // for

} // paperFace::writeAdherent

VektorR2 paperFace::getMidPoint()
{
  VektorR2 mid(0.0, 0.0);

  for(int i = 0; i < no_points; i++) {
    mid += projection[i];
  } // for

  if(no_points > 0)
    mid /= (double)no_points;

  return mid;
} // paperFace::getMidPoint

void paperFace::printProjection(int start_with_index)
{
  for(int i = start_with_index; i < no_points+start_with_index; i++) {

    double dist = (projection[index_save(i)] - projection[index_save(i-1)]).Norm2();

    cout << "(" << setprecision(5) << fixed << projection[index_save(i)][0] << "," << projection[index_save(i)][1] << ")" << " dist to prev.: " << dist << endl;
   
  } // for
  
} // paperFace::printProjection

void paperFace::plot3DVertices()
{
  for(int i = 0; i < no_points+1; i++) {

    // 3D version

    double dist = (pointArray[point[index_save(i)]] - pointArray[point[index_save(i-1)]]).Norm2();
    if(debugFile) (*debugFile) << pointArray[point[index_save(i)]][0] << " " << pointArray[point[index_save(i)]][1] << " " << pointArray[point[index_save(i)]][2] << endl;   
  } // for
  
} // paperFace::plot3DVertices

void paperFace::plotProjection()
{
  for(int i = 0; i < no_points+1; i++) {

    // 2D version

    double dist = (projection[index_save(i)] - projection[index_save(i-1)]).Norm2();
    if(debugFile) (*debugFile) << projection[index_save(i)][0] << " " << projection[index_save(i)][1] << endl;   
  } // for      
} // paperFace::plotProjection

void paperFace::get3DVertices(int start_with_index, Json::Value &vertices_json)
{
  for(int i = start_with_index; i < no_points+start_with_index; i++) {

    double dist = (pointArray[point[index_save(i)]] - pointArray[point[index_save(i-1)]]).Norm2();

    Json::Value vertex;
    vertex["x"] = Json::Value(pointArray[point[index_save(i)]][0]);
    vertex["y"] = Json::Value(pointArray[point[index_save(i)]][1]);
    vertex["z"] = Json::Value(pointArray[point[index_save(i)]][2]);
    vertices_json.append(vertex);
    
    //cout << "vertex #" << i << ": (" << setprecision(5) << fixed << pointArray[point[index_save(i)]][0] << "," << pointArray[point[index_save(i)]][1] << "," << pointArray[point[index_save(i)]][2] << ")" << " [" << point[index_save(i)] << "]" << " dist to prev.: " << dist << endl;
   

  } // for
  
} // paperFace::print3DVertices

bool paperFace::CheckColinear()
{
  VektorR3 dir, previous_dir;
  double    product = 1.0;

  if(no_points < 3) return true;

  for(int i = 0; i < no_points+1; i++) {

    dir = pointArray[point[index_save(i)]] - pointArray[point[index_save(i+1)]];
    dir.Normalize();

    if(i > 0) {

      product = fabs(dir*previous_dir);

      if(product < 0.99) return false;      
    } // if

    previous_dir = dir;
  } // for

  return true;
} // paperFace::CheckColinear


VektorR2 paperFace::calcNormal(int index)
{
  VektorR2 dir, normal;

  dir = projection[index_save(index+1)] - projection[index_save(index)];

  normal = VektorR2(-1.0*dir[1], dir[0]);
  normal.Normalize();

  // test whether normal points inward. Calc point
  // outside the polygon. Intersect with all vertices 
  // into one direction. Number of interections should
  // be uneven

  VektorR2 on_the_line = 0.356*(dir) + projection[index_save(index)];
  VektorR2 inner_point = on_the_line + normal*(dir.Norm2()/100.0);

  if(PointInside(inner_point, dir) == false) normal *= -1.0;


  return normal;
} // paperFace::calcNormal

bool paperFace::PointInside(VektorR2 inner_point, VektorR2 dir)
{
  int counter = 0;

  for(int i = 0; i < no_points; i++) {

    // if the point matches a vertex, it should be considered outside
    /*
    cout << "innter_point: "; inner_point.print(); cout << endl;
    cout << "projection[i]: "; projection[index_save(i)].print(); cout << endl;
    */
    if((inner_point-projection[index_save(i)]).Norm2() < epsilon) return false;

    // if the point is on an edge, it should also be considered outside
    /*
    if(PointBetweenLineEnds(inner_point, projection[index_save(i)], 
			    projection[index_save(i+1)]) == true) return false;
    */
    double p1, p2;

    if(IntersectLineLineParam(inner_point, inner_point+0.02*dir, projection[index_save(i)], 
			  projection[index_save(i+1)], p1, p2) == true) {


      if((p1 > 0.0) && (p2 > 0.0001) && (p2 < 0.9999)) {	
	counter++;
      } // if
    } // if
  } // for

  if((counter % 2) == 0) return false;
  else return true;
} // paperFace::PointInside

VektorR3 paperFace::calc3DNormal()
{
  VektorR3 p1, p2, p3;
  VektorR3 normal(0.0, 0.0, 0.0);
  double   dot_prod = 1000000.0;

  for(int i = 0; i < no_points; i++) {

    p1 = pointArray[point[index_save(i)]];
    p2 = pointArray[point[index_save(i+1)]];
    p3 = pointArray[point[index_save(i+2)]];
    
    VektorR3 d1, d2;

    d1 = p2-p1; d1.Normalize();
    d2 = p3-p1; d2.Normalize();

    double d_prod = fabs(d1*d2);

    if(d_prod < dot_prod) {
      dot_prod = d_prod;
      normal = VektorR3::CrossProduct(d1, d2);
    } // if
  } // for
  
  return normal;
} // calc3DNormal

// TODO: While merging the faces, keep track of the faces that were deleted. 
// While producing the projection, instead of adding the merged face, add the original faces with their respective UVs.
//

bool paperFace::mergeFace(paperFace& face)
{
  int      s, t;
  VektorR3 n1, n2;

  int tmp_points[MAX_POINTS_PER_FACE];

  if(sharesSomeEdge(face, s, t) == true) {

    if((point[s] != face.point[t]) || (point[index_save(s+1)] != face.point[face.index_save(t+1)])) {
      cout << "error: actually no common vertex" << endl;
      exit(-1);
    } // if

    // point s on this equals t on face, s+1 equals t+1

    // merging only takes place if both lie
    // within the same plane

    n1 = calc3DNormal();
    n1.Normalize();
    n2 = face.calc3DNormal();
    n2.Normalize();

    if(fabs(n1*n2) > 0.9999) {
      if((no_points + face.no_points) > MAX_POINTS_PER_FACE) {
	cerr << "error in paperFace::mergeFace(). Number of points exceed limit. Increase constant MAX_POINTS_PER_FACE" << endl;
	exit(-1);
      } // if

      int index = 0;

      // save all points from 0 to s-1 

      for(int i = s+1; i < s+1+no_points; i++) {
	tmp_points[index++] = point[index_save(i)];
      } // for

      // now copy all new points from face. The first and 
      // the last are duplicated and do not need to be copied

      for(int i = 1; i < face.no_points-1; i++) {
	tmp_points[index++] = face.point[face.index_save(t-i)];
      } // for

      // copy points from tmp buffer to this*

      for(int i = 0; i < index; i++) {
	point[i] = tmp_points[i];
      } // for

      no_points = index;

      return true;

    } // if
 } // if
  
  return false;
} // paperFace::mergeFace

void paperFace::copy(paperFace& copyFrom)
{
  for(int i = 0; i < copyFrom.no_points; i++) {
    point[i] = copyFrom.point[i];
    projection[i] = copyFrom.projection[i];
    adherent_connected[i] = copyFrom.adherent_connected[i];
    adherent_ID[i] = copyFrom.adherent_ID[i];
    neighbor[i] = copyFrom.neighbor[i];
  } // for

  no_points = copyFrom.no_points;
  processed = copyFrom.processed;
  drawn = copyFrom.drawn;
  adherent = copyFrom.adherent;
} // paperFace::copy

// sharesThisEdge test whether the edge s is shares

bool paperFace::sharesThisEdge(paperFace& face, int s, int& t, bool do_not_swap)
{
  for(t = 0; t < face.no_points; t++) {
    if((point[s] == face.point[t])) {
      
      if((point[index_save(s+1)] == face.point[face.index_save(t+1)])) {
	return true;
      } // if
      
      if((point[index_save(s+1)] == face.point[face.index_save(t-1)])) {
	
	if(do_not_swap == false) {

	  // turn around the order of the points

	  for(int w = 1; w <= face.no_points/2; w++) {
	    face.swapVertexData(t-w, t+w);
	    face.swapEdgeData(t+w-1, t-w);
	  } // for
	} // if

	return true;
      } // if     
    } // if
  } // for

  return false;
} // paperFace::sharesThisEdge

// sharesSomeEdge test whether some edge is shared. It is
// returned by the integer s, the other edge as integer t

bool paperFace::sharesSomeEdge(paperFace& face, int& s, int& t)
{
  for(s = 0; s < no_points; s++) {
    if(sharesThisEdge(face, s, t, false) == true) return true;
  } // for

  return false;
} // paperFace::sharesSomeEdge

void paperFace::swapVertexData(int pi, int pj)
{
  int i = index_save(pi);
  int j = index_save(pj);

  int        temp_point = point[i];
  VektorR2   temp_projection = projection[i];
	    
  point[i] = point[j];
  point[j] = temp_point;
  
  projection[i] = projection[j];
  projection[j] = temp_projection;
  
} // paperFace::swapVertexData

void paperFace::swapEdgeData(int pi, int pj)
{ 
  int i = index_save(pi);
  int j = index_save(pj);

  bool       temp_con = adherent_connected[i];
  int        temp_ID = adherent_ID[i];
  paperFace* temp_neighbor = neighbor[i];
	    
  adherent_connected[i] = adherent_connected[j];
  adherent_connected[j] = temp_con;
  
  adherent_ID[i] = adherent_ID[j];
  adherent_ID[j] = temp_ID;
  
  neighbor[i] = neighbor[j];
  neighbor[j] = temp_neighbor;
} // paperFace::swapEdgeData

int paperFace::getSharedEdgeID(paperFace& face, int edge_index)
{
  int s = edge_index;

  for(int t = 0; t < face.no_points; t++) {
    if((point[index_save(s)] == face.point[t])) {
      
      if((point[index_save(s+1)] == face.point[face.index_save(t+1)])) {
	
	return face.adherent_ID[t];
      } // if
      
      if((point[index_save(s+1)] == face.point[face.index_save(t-1)])) {

	return face.adherent_ID[face.index_save(t-1)];
      } // if
    } // if   
  } // for

  return -1;
} // paperFace::getSharedEdgeID

int paperFace::getConnectedFace(int edge_index, int facecnt)
{
  int id = -1;

  for(int i = 0; i < facecnt; i++) {

    if((&(faceArray[i]) != this) && (faceArray[i].adherent == false)) {
      id = getSharedEdgeID(faceArray[i], edge_index);
      
      if(id > -1) return id;    
    } // if   
  } // for

  //  cout << "WARNING: unconnected edge" << endl;

  return -1;
} // paperFace::getConnectedFace

void paperFace::removeColinearPoints()
{
  bool action;

  VektorR3 dir_1, dir_2;
  int      tmp_points[MAX_POINTS_PER_FACE];

  do {

    action = false;

    for(int i = 0; i < no_points; i++) {
      dir_1 = pointArray[point[index_save(i+1)]] - pointArray[point[index_save(i)]];
      dir_2 = pointArray[point[index_save(i+2)]] - pointArray[point[index_save(i)]];
      
      dir_1.Normalize();
      dir_2.Normalize();

      if(dir_1*dir_2 > 0.9999) {

	for(int j = i+2; j < i+2+no_points-1; j++) {
	  tmp_points[j-(i+2)] = point[index_save(j)];
	} // for

	no_points--;

	for(int j = 0; j < no_points; j++) {
	  point[index_save(j)] = tmp_points[j];
	} // for

	action = true;
      } // if
    } // for
  } while(action == true);
} // paperFace::removeColinearPoints

bool paperFace::SharesProjectedVertex(int index, VektorR2 p1, VektorR2 p2)
{
  if(((projection[index_save(index)]-p1).Norm2() < epsilon) || 
     ((projection[index_save(index)]-p2).Norm2() < epsilon)) {

    VektorR2 dir_1 = p2-p1; dir_1.Normalize();
    VektorR2 dir_2;

    dir_2 = projection[index_save(index)]-projection[index_save(index-1)];
    dir_2.Normalize();
    if(fabs(dir_1*dir_2) > 0.9999) return true;
    /*
    dir_2 = projection[index_save(index)]-projection[index_save(index+1)];
    dir_2.Normalize();
    if(fabs(dir_1*dir_2) > 0.9999) return true;
    */    
  } // if

  return false;  
} // paperFace::SharesProjectedVertex

bool paperFace::sharesProjectedEdge(int index, int no_faces)
{
  for(int i = 0; i < no_faces; i++) {

    if((&(faceArray[i]) != this) && (faceArray[i].drawn == true)) {

      for(int edge = 0; edge < faceArray[i].no_points; edge++) {

	if((projection[index] - faceArray[i].projection[edge]).Norm2() < epsilon) {
	   
	  if((projection[index_save(index+1)] - faceArray[i].projection[faceArray[i].index_save(edge+1)]).Norm2() < epsilon) return true;
	
	  if((projection[index_save(index+1)] - faceArray[i].projection[faceArray[i].index_save(edge-1)]).Norm2() < epsilon) return true;
	} // if
      } // for
    } // if
  } // for

    return false;
} // paperFace::sharesProjectedEdge

bool paperFace::IsInnerEdge(int p1, int p2)
{
  int counter = 0;

  for(int i = 0; i < no_points; i++) {
    if(SharesProjectedVertex(i, projection[index_save(p1)], projection[index_save(p2)]) == true) {
      counter++;
      if(counter == 2) return true;
    } // if
  } // for

  return false;
} // paperFace::IsInnerEdge

void paperFace::clearConnectedFlags()
{
  for(int i = 0; i < MAX_POINTS_PER_FACE; i++) {
    adherent_connected[i] = false;
  } // for
} // paperFace::clearConnectedFlags

bool paperFace::EdgeOnEdge(int index, VektorR2 x1, VektorR2 x2)
{
  // Lines do intersect if their end-points are equal

  VektorR2 y1, y2, dir;
  double    tmp, dist_1, dist_2;

  y1 = projection[index];
  y2 = projection[index_save(index+1)];
  /*
  cout << "    ********" << endl;
  cout << "    y1=(" << y1[0] << "," << y1[1] << ") y2=(" << y2[0] << "," << y2[1] << ")" << endl;
  cout << "    x1=(" << x1[0] << "," << x1[1] << ") x2=(" << x2[0] << "," << x2[1] << ")" << endl;
  cout << "    ********" << endl;
  */
  dist_1 = (y1 - x1).Norm2();
  tmp = (y1 - x2).Norm2();
  if(tmp < dist_1) dist_1 = tmp;

  dist_2 = (y2 - x1).Norm2();
  tmp = (y2 - x2).Norm2();
  if(tmp < dist_2) dist_2 = tmp;

  if((dist_1 + dist_2) < epsilon) {
    //    cout << "MATCH" << endl;
    return true;
  } // if

  // We concentrate on line (x1, x2). Either y1 or y2
  // might lie on the line
  
  dir = y2-y1;

  int pol_counter = 0;
  
  if(PointOnLine(y1, dir, x1, x2) == true) {
    //    cout << "MATCH" << endl;
    pol_counter++;
  } // if
  
  if(PointOnLine(y2, dir, x1, x2) == true) {
    //    cout << "MATCH" << endl;
    pol_counter++;
  } // if
  
  // We concentrate on line (y1, y2). Either x1 or x2
  // might lie on the line
  
  dir = x2-x1;
  
  if(PointOnLine(x1, dir, y1, y2) == true) {
    //    cout << "MATCH" << endl;
    pol_counter++;
  } // if

  if(PointOnLine(x2, dir, y1, y2) == true) {
    //    cout << "MATCH" << endl;
    pol_counter++;
  } // if

  if(pol_counter > 1) return true;
  
  return false;
} // paperFace::EdgeOnEdge

bool paperFace::EdgeOnEdge(int index, VektorR3 x1, VektorR3 x2)
{
  // Lines do intersect if their end-points are equal

  VektorR3 y1, y2, dir;
  double    tmp, dist_1, dist_2;

  y1 = pointArray[point[index]];
  y2 = pointArray[point[index_save(index+1)]];
  /*
  cout << "    ********" << endl;
  cout << "    y1=(" << y1[0] << "," << y1[1] << ") y2=(" << y2[0] << "," << y2[1] << ")" << endl;
  cout << "    x1=(" << x1[0] << "," << x1[1] << ") x2=(" << x2[0] << "," << x2[1] << ")" << endl;
  cout << "    ********" << endl;
  */
  dist_1 = (y1 - x1).Norm2();
  tmp = (y1 - x2).Norm2();
  if(tmp < dist_1) dist_1 = tmp;

  dist_2 = (y2 - x1).Norm2();
  tmp = (y2 - x2).Norm2();
  if(tmp < dist_2) dist_2 = tmp;

  if((dist_1 + dist_2) < epsilon) {
    //    cout << "MATCH" << endl;
    return true;
  } // if

  // We concentrate on line (x1, x2). Either y1 or y2
  // might lie on the line
  
  dir = y2-y1;

  int pol_counter = 0;
  
  if(PointOnLine(y1, dir, x1, x2) == true) {
    //    cout << "MATCH" << endl;
    pol_counter++;
  } // if
  
  if(PointOnLine(y2, dir, x1, x2) == true) {
    //    cout << "MATCH" << endl;
    pol_counter++;
  } // if
  
  // We concentrate on line (y1, y2). Either x1 or x2
  // might lie on the line
  
  dir = x2-x1;
  
  if(PointOnLine(x1, dir, y1, y2) == true) {
    //    cout << "MATCH" << endl;
    pol_counter++;
  } // if

  if(PointOnLine(x2, dir, y1, y2) == true) {
    //    cout << "MATCH" << endl;
    pol_counter++;
  } // if

  if(pol_counter > 1) return true;
  
  return false;
} // paperFace::EdgeOnEdge

void paperFace::printProjection()
{
  for(int i = 0; i < no_points; i++)
    cout << "(" << projection[i][0] << "," << projection[i][1] << ")" << endl;
} // paperFace::printProjection

void paperFace::initialize()
{
  for(int i = 0; i < MAX_POINTS_PER_FACE; i++) {
    point[i] = -1;
    adherent_connected[i] = false;
    adherent_ID[i] = -1;
    neighbor[i] = NULL;
  } // for

  no_points = 0;
  processed = false;
  drawn = false;
  adherent = false;
  visited = false;

  cost = 100000000;
} // paperFace::initialize

double paperFace::estimateExtend()
{
  double min_x, max_x, min_y, max_y;

  min_x = min_y = 1000000;
  max_x = max_y = -1000000;

  for(int i = 0; i < no_points; i++) {
    if(projection[i][0] < min_x) min_x = projection[i][0];
    if(projection[i][1] < min_y) min_y = projection[i][1];

    if(projection[i][0] > max_x) max_x = projection[i][0];
    if(projection[i][1] > max_y) max_y = projection[i][1];
  } // for
  return (max_x - min_x)*(max_y - min_y);
} // paperFace::estimateExtend

double paperFace::getCumulatedDistance(int no_faces)
{
  //  int counter = 0;
  double dist;

  for(int i = 0; i < no_faces; i++) {
    if((&(faceArray[i]) != this) && (faceArray[i].drawn == true)) {
      for(int vert = 0; vert < faceArray[i].no_points; vert++) {
	for(int comp = 0; comp < no_points; comp++) {
	  dist += (faceArray[i].projection[vert] - projection[comp]).Norm2();
	} // for
      } // for
    } // if
  } // for

  //  return counter;
  return dist;
} // paperFace::getCumulatedDistance

bool paperFace::getAdjacentEdge(int edge_index, int no_faces, VektorR2& mid_point, VektorR2& this_end, VektorR2& other_end, bool reverse)
{
  for(int i = 0; i < no_faces; i++) {
    if((&(faceArray[i]) != this) && (faceArray[i].drawn == true)) {
      for(int vertex = 0; vertex < faceArray[i].no_points; vertex++) {
	  
	if(reverse == false) {
	  if((projection[index_save(edge_index)]-faceArray[i].projection[vertex]).Norm2() < epsilon) {
	    
	    if(faceArray[i].sharesProjectedEdge(faceArray[i].index_save(vertex), no_faces) == false) {
	      mid_point = projection[index_save(edge_index)];
	      this_end = projection[index_save(edge_index+1)];
	      other_end = faceArray[i].projection[faceArray[i].index_save(vertex+1)];
	      return true;
	    } // if
	    
	    if(faceArray[i].sharesProjectedEdge(faceArray[i].index_save(vertex-1), no_faces) == false) {
	      mid_point = projection[index_save(edge_index)];
	      this_end = projection[index_save(edge_index+1)];
	      other_end = faceArray[i].projection[faceArray[i].index_save(vertex-1)];
	      return true;
	    } // if	  
	  } // if
	} // if
	else {
	  if((projection[index_save(edge_index)]-faceArray[i].projection[vertex]).Norm2() < epsilon) {
	    
	    if(faceArray[i].sharesProjectedEdge(faceArray[i].index_save(vertex), no_faces) == false) {
	      mid_point = projection[index_save(edge_index)];
	      this_end = projection[index_save(edge_index-1)];
	      other_end = faceArray[i].projection[faceArray[i].index_save(vertex+1)];
	      return true;
	    } // if
	    
	    if(faceArray[i].sharesProjectedEdge(faceArray[i].index_save(vertex-1), no_faces) == false) {
	      mid_point = projection[index_save(edge_index)];
	      this_end = projection[index_save(edge_index-1)];
	      other_end = faceArray[i].projection[faceArray[i].index_save(vertex-1)];
	      return true;
	    } // if	  
	  } // if
	} // else
	
      } // for
    } // if
  } // for

  return false;
} // paperFace::getAdjacentEdge

// return the angle of adjacent faces with commong point edge_index
// reverse means that we are at the end of the edge edge_index and
// go into the opposite direction

double paperFace::getEnclosedAngle(int edge_index, int no_faces, VektorR2& mid_point, bool reverse)
{
  VektorR2 this_vertex, neighbor_vertex;
  VektorR2 dir_1, dir_2;
  double   angle = 1000;

  if(getAdjacentEdge(edge_index, no_faces, mid_point, this_vertex, neighbor_vertex, reverse) == true) {
    
    dir_1 = this_vertex - mid_point;
    dir_1.Normalize();
    
    dir_2 = neighbor_vertex - mid_point;
    dir_2.Normalize();

    angle = acos(dir_1*dir_2)*360.0/(2.0*M_PI);
  } // if

  return angle;
} // paperFace::getEnclosedAngle
