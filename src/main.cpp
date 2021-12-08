#include "VInternal.H"
#include <iostream>
#include <stdlib.h>
//#include "VCollide.H"
#include "VCScene.h"
using namespace std;

const int DATA_DUMP = 0;
const int RENDER_DUMP = 1;
void startRendererWithDisplay(VCScene* vs, int option, const std::string& frameFilename, int frame);

float toBW(int bytes, float sec) {
  return static_cast<float>(bytes) / (1024. * 1024. * 1024.) / sec;
}
int main(int argc, char *argv[]) {

  // if (argc != 1)
  //   {
  //     cerr<<argv[0]<<": USAGE: "<<argv[0]<<"\n";
  //     exit(1);
  //   }
  cout << "hello, mp" << endl;
  VCInternal vc(2, 100);
  VCScene vs(2);
  int id[2];

  int i;
  for (i = 0; i < 2; i++) // create both the objects.
  {
    vc.NewObject(&id[i]);
    vs.NewObject(&(id[i]));

    // the geometry is a unit cube with one vertex at the origin.
    double v1[3], v2[3], v3[3];

    v1[0] = 0;
    v1[1] = 0.5;
    v1[2] = 0.015;
    v2[0] = 0 +0.015;
    v2[1] = 0.5;
    v2[2] = 1.99019234502e-08;
    v3[0] = 0 +0.100000019902;
    v3[1] = 0.515;
    v3[2] = 1.99019234502e-08;
    vc.AddTri(v1, v2, v3);

    double p1[3], p2[3], p3[3];

        
        memcpy(&p1, &v1, sizeof(double) * 3);
        memcpy(&p2, &v2, sizeof(double) * 3);
        memcpy(&p3, &v3, sizeof(double) * 3);
        vs.AddTri(p1, p2, p3);

    v1[0] = 0;
    v1[1] = 0.5;
    v1[2] = 0.015;
    v2[0] = 0 +1.99019234501e-08;
    v2[1] = 0.515;
    v2[2] = 1.99019234502e-08;
    v3[0] = 0 -0.0149999999999;
    v3[1] = 0.5;
    v3[2] = 1.99019234502e-08;
    vc.AddTri(v1, v2, v3);

        
        memcpy(&p1, &v1, sizeof(double) * 3);
        memcpy(&p2, &v2, sizeof(double) * 3);
        memcpy(&p3, &v3, sizeof(double) * 3);
        vs.AddTri(p1, p2, p3);

    v1[0] = 0;
    v1[1] = 0.5;
    v1[2] = 0.015;
    v2[0] = 0 -0.0149999999999;
    v2[1] = 0.5;
    v2[2] = 1.99019234502e-08;
    v3[0] = 0 -5.9705770357e-08;
    v3[1] = 0.485;
    v3[2] = 1.99019234502e-08;
    vc.AddTri(v1, v2, v3);

        
        memcpy(&p1, &v1, sizeof(double) * 3);
        memcpy(&p2, &v2, sizeof(double) * 3);
        memcpy(&p3, &v3, sizeof(double) * 3);
        vs.AddTri(p1, p2, p3);

    v1[0] = 0;
    v1[1] = 0.5;
    v1[2] = 0.015;
    v2[0] = 0 -5.9705770357e-08;
    v2[1] = 0.485;
    v2[2] = 1.99019234502e-08;
    v3[0] = 0 +0.015;
    v3[1] = 0.5;
    v3[2] = 1.99019234502e-08;
    vc.AddTri(v1, v2, v3);

        
        memcpy(&p1, &v1, sizeof(double) * 3);
        memcpy(&p2, &v2, sizeof(double) * 3);
        memcpy(&p3, &v3, sizeof(double) * 3);
        vs.AddTri(p1, p2, p3);

    // lower

    v1[0] = 0;
    v1[1] = 0.5;
    v1[2] = -0.015;
    v2[0] = 0 +0.015;
    v2[1] = 0.5;
    v2[2] = 1.99019234502e-08;
    v3[0] = 0 +0.100000019902;
    v3[1] = 0.515;
    v3[2] = 1.99019234502e-08;
    vc.AddTri(v1, v2, v3);

        
        memcpy(&p1, &v1, sizeof(double) * 3);
        memcpy(&p2, &v2, sizeof(double) * 3);
        memcpy(&p3, &v3, sizeof(double) * 3);
        vs.AddTri(p1, p2, p3);

    v1[0] = 0;
    v1[1] = 0.5;
    v1[2] = -0.015;
    v2[0] = 0 +1.99019234501e-08;
    v2[1] = 0.515;
    v2[2] = 1.99019234502e-08;
    v3[0] = 0 -0.0149999999999;
    v3[1] = 0.5;
    v3[2] = 1.99019234502e-08;
    vc.AddTri(v1, v2, v3);

        
        memcpy(&p1, &v1, sizeof(double) * 3);
        memcpy(&p2, &v2, sizeof(double) * 3);
        memcpy(&p3, &v3, sizeof(double) * 3);
        vs.AddTri(p1, p2, p3);

    v1[0] = 0;
    v1[1] = 0.5;
    v1[2] = -0.015;
    v2[0] = 0 -0.0149999999999;
    v2[1] = 0.5;
    v2[2] = 1.99019234502e-08;
    v3[0] = 0 -5.9705770357e-08;
    v3[1] = 0.485;
    v3[2] = 1.99019234502e-08;
    vc.AddTri(v1, v2, v3);

        
        memcpy(&p1, &v1, sizeof(double) * 3);
        memcpy(&p2, &v2, sizeof(double) * 3);
        memcpy(&p3, &v3, sizeof(double) * 3);
        vs.AddTri(p1, p2, p3);

    v1[0] = 0;
    v1[1] = 0.5;
    v1[2] = -0.015;
    v2[0] = 0 -5.9705770357e-08;
    v2[1] = 0.485;
    v2[2] = 1.99019234502e-08;
    v3[0] = 0 +0.015;
    v3[1] = 0.5;
    v3[2] = 1.99019234502e-08;
    vc.AddTri(v1, v2, v3);

        
        memcpy(&p1, &v1, sizeof(double) * 3);
        memcpy(&p2, &v2, sizeof(double) * 3);
        memcpy(&p3, &v3, sizeof(double) * 3);
        vs.AddTri(p1, p2, p3);

    vc.EndObject();
    vs.EndObject();
  }

  double trans0[4][4], trans1[4][4]; // transformation matrices.
  double all_trans[4 * 4 * 2];

  // initialize the transformation matrices to identity.

  for (i = 0; i < 4; i++) {
    int j;
    for (j = 0; j < 4; j++)
      trans0[i][j] = trans1[i][j] = ((i == j) ? 1.0 : 0.0);
  }

  for (int simulation_step = -25; simulation_step <= 25; simulation_step++) // perform 51 frames of the simulation
  {
    cout << "Simulation step: " << simulation_step << "\n";

    // in successive frames of the simulation, the two objects
    // approach each other from far and finally collide and cross
    // each other.
    trans0[0][3] = 0.3 * simulation_step;  // we translate both the objects
    trans1[0][3] = -0.3 * simulation_step; // along the X-axis only.

    for (i = 0; i < 4; i++) {
      int j;
      for (j = 0; j < 4; j++) {
        all_trans[i * 4 + j] = trans0[i][j];
      }
    }
    for (i = 0; i < 4; i++) {
      int j;
      for (j = 0; j < 4; j++) {
        all_trans[16 + i * 4 + j] = trans1[i][j];
      }
    }

    vc.UpdateAllTrans(id, 2, all_trans);

    vc.Collide();

    startRendererWithDisplay(&vs, DATA_DUMP, "./output/nbody", i);
  }


  cout << " Finish Detected collision between objects\n";
  return 0;
}
