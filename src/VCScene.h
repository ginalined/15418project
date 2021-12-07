#ifndef VCSCENE_H
#define VCSCENE_H

#include "image.h"
#include <cstring>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include "util.h"


struct SimplifyTri {
  int id;
  double p1[3], p2[3], p3[3];
};


class SimplfiedObj
{
  friend class VCScene;
  
private:
  
  int id;          //the id of the object.
  double trans[16]; //the current transformation matrix for the object.
  SimplifyTri *triList;   // list of triangles
  int triBufferSize;  // the size of the triangle buffer allocated
  int triSize;    // the real number of triangles

  void addTri(double v1[], double v2[], double v3[]);
};

class VCScene
{
private:
  Image* image;
  int      size;        //buffer size of the "vcobjects" array.
  SimplfiedObj **vc_objects;//array of pointers to VCObjects.

  int      next_id;     //next free id, since ids are generated by the program.
  int      current_id;  //the id of the object being worked on by "AddTri"

  void transTri(double *tri, double *trans);
 
public:
  VCScene(int size);
  ~VCScene();

  int NewObject(int *id); //create a new object in the database.

  int AddTri(double v1[], double v2[], double v3[]); //insert the geometry.

  int EndObject(void);    //tell VCollide that inserting the geometry is complete.
  
  int UpdateTrans(int id, double *trans); // after collision detection, insert the transformation matrix
 
  void allocateImage(int width, int height);
  
  Image *getImage();

  void clearImage();

  void render();

  void dumpTriangles(FILE *fp);
};
#endif