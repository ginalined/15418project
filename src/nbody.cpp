#include "CycleTimer.h"
#include "VInternal.H"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
//#include "VCollide.H"
#include "VCScene.h"
using namespace std;

const int NO_OF_OBJECTS = 32;    // number of instances
const int SIMULATION_STEPS = 10; // number of steps in the simulation.
const int SCREEN_SIZE = 100;

const int DATA_DUMP = 0;
const int RENDER_DUMP = 1;
void startRendererWithDisplay(VCScene *vs, int option,
                              const std::string &frameFilename, int frame, int num_tri);

int main(int argc, char *argv[]) {

  double startTime = CycleTimer::currentSeconds();
  double dumpTime = 0;

  if (argc != 3) {
    cerr << argv[0] << ": USAGE: " << argv[0]
         << " <input-file> <transformation-file>\n";
    exit(1);
  }

  int num_tri;
  VCInternal vc(NO_OF_OBJECTS, SCREEN_SIZE);
  VCScene vs(NO_OF_OBJECTS);
  int id[NO_OF_OBJECTS];

  int i;

  FILE *fp = fopen(argv[1], "r");
  fscanf(fp, "%d", &num_tri);
  cout << "num_tri = " << num_tri << endl;
  for (i = 0; i < NO_OF_OBJECTS; i++) // add the objects to the library.
  {
    // cout<<"Reading object "<<i<<"\n";

    vc.NewObject(&(id[i]));
    vs.NewObject(&(id[i]));
    // cout<<"Adding triangles\n";

    for (int j = 1; j <= num_tri; j++) {
      double v1[3], v2[3], v3[3];
      fscanf(fp, "%lf %lf %lf", &(v1[0]), &(v1[1]), &(v1[2]));
      fscanf(fp, "%lf %lf %lf", &(v2[0]), &(v2[1]), &(v2[2]));
      fscanf(fp, "%lf %lf %lf", &(v3[0]), &(v3[1]), &(v3[2]));
      // cout << v3[0] << " " << v3[1] << " "<< v3[2] << endl;

      vc.AddTri(v1, v2, v3);

      double p1[3], p2[3], p3[3];
      memcpy(&p1, &v1, sizeof(double) * 3);
      memcpy(&p2, &v2, sizeof(double) * 3);
      memcpy(&p3, &v3, sizeof(double) * 3);
      vs.AddTri(p1, p2, p3);
    }
    // std::cout<<"closing files\n";

    // cout<<"Calling finish_object\n";
    vc.EndObject();
    // vs.EndObject();

    cout << "Inserted object " << i << "\n";
  }

  fclose(fp);

  // vc.EndAllObjects();
  // FILE *fp = fopen(argv[2], "r");

  for (i = 1; i <= SIMULATION_STEPS; i++) // perform the simulation.
  {
    cout << "Simulation step : " << i << "\n";
    int j;
    double *all_trans = new double[NO_OF_OBJECTS * 16];

    //     for (j=0; j<NO_OF_OBJECTS; j++)
    // {
    //   for (int j1=0; j1<16; j1++){
    //     fscanf(fp, "%lf", &(all_trans[j*16+j1]));
    //   }
    // }

    for (j = 0; j < NO_OF_OBJECTS; j++) {
      all_trans[j * 16] = 1;
      all_trans[j * 16 + 5] = 1;
      all_trans[j * 16 + 10] = 1;
      all_trans[j * 16 + 15] = 1;

      all_trans[j * 16 + 1] = 0;
      all_trans[j * 16 + 2] = 0;
      all_trans[j * 16 + 3] = 0;
      all_trans[j * 16 + 4] = 0;
      all_trans[j * 16 + 6] = 0;
      all_trans[j * 16 + 7] = 0;
      all_trans[j * 16 + 8] = 0;
      all_trans[j * 16 + 9] = 0;
      all_trans[j * 16 + 11] = 0;
      all_trans[j * 16 + 12] = 0;
      all_trans[j * 16 + 13] = 0;
      all_trans[j * 16 + 14] = 0;
    }

    // (along y-axis) half move down, half move up
    for (j = 0; j < NO_OF_OBJECTS / 2; j++) {
      all_trans[j * 16 + 7] = -0.5 * (i - 3);
    }
    for (; j < NO_OF_OBJECTS; j++) {
      all_trans[j * 16 + 7] = 0.5 * (i - 3);
    }

    vc.UpdateAllTrans(id, NO_OF_OBJECTS, all_trans);

    vc.Collide();
    
    double dumpStartTime = CycleTimer::currentSeconds();
    // output trans matrix
    for (j = 0; j < NO_OF_OBJECTS; j++) {
      double *per_trans = new double[16];
      for (int jtrans = j * 16; jtrans < (j + 1) * 16; jtrans++) {
        per_trans[jtrans - (j * 16)] = all_trans[jtrans];
      }
      vs.UpdateTrans(id[j], per_trans);
    }
    startRendererWithDisplay(&vs, DATA_DUMP, "./output/nbody", i, num_tri);
    dumpTime += CycleTimer::currentSeconds() - dumpStartTime;
  }

  double endTime = CycleTimer::currentSeconds();
  std::cout << endTime - startTime - dumpTime << std::endl;
  // double seconds = difftime(endtime, now);
  // printf ("%.f running time\n", seconds);
  // cout<<" Finish Detected collision between objects\n";

  return 0;
}
