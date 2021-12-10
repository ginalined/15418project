#include "CycleTimer.h"
#include "VInternal.H"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
//#include "VCollide.H"

using namespace std;

const int NO_OF_OBJECTS = 32;    // number of instances
const int SIMULATION_STEPS = 10; // number of steps in the simulation.
const int SCREEN_SIZE = 100;

const int DATA_DUMP = 0;
const int RENDER_DUMP = 1;


int main(int argc, char *argv[]) {

 

  if (argc != 3) {
    cerr << argv[0] << ": USAGE: " << argv[0]
         << " <input-file> <transformation-file>\n";
    exit(1);
  }

  int num_tri = 12;
  VCInternal vc(NO_OF_OBJECTS, SCREEN_SIZE);

  int id[NO_OF_OBJECTS];

  int i;


  for (i = 0; i < NO_OF_OBJECTS; i++) // add the objects to the library.
  {
    // cout<<"Reading object "<<i<<"\n";

    vc.NewObject(&(id[i]));
  
    // cout<<"Adding triangles\n";

    // for (int j = 1; j <= num_tri; j++) {
    //   double v1[3], v2[3], v3[3];
    //   fscanf(fp, "%lf %lf %lf", &(v1[0]), &(v1[1]), &(v1[2]));
    //   fscanf(fp, "%lf %lf %lf", &(v2[0]), &(v2[1]), &(v2[2]));
    //   fscanf(fp, "%lf %lf %lf", &(v3[0]), &(v3[1]), &(v3[2]));
    //   // cout << v3[0] << " " << v3[1] << " "<< v3[2] << endl;

    //   vc.AddTri(v1, v2, v3);

    //   double p1[3], p2[3], p3[3];
    //   memcpy(&p1, &v1, sizeof(double) * 3);
    //   memcpy(&p2, &v2, sizeof(double) * 3);
    //   memcpy(&p3, &v3, sizeof(double) * 3);
     
    // }
    double v1[3], v2[3], v3[3];
      
      v1[0] = 0.0 + i*10; v1[1] = 0.0; v1[2] = 0.0;
      v2[0] = 1.0 + i*10; v2[1] = 0.0; v2[2] = 0.0;
      v3[0] =  1.0 + i*10; v3[1] = 0.0; v3[2] = 1.0;
      vc.AddTri(v1, v2, v3);
      
      v1[0] = 0.0 + i*10; v1[1] = 0.0; v1[2] = 0.0;
      v2[0] =  i*10; v2[1] = 0.0; v2[2] = 1.0;
      v3[0] =  1.0 + i*10; v3[1] = 0.0; v3[2] = 1.0;
      vc.AddTri(v1, v2, v3);
      
      v1[0] = 0.0 + i*10; v1[1] = 1.0; v1[2] = 0.0;
      v2[0] = 1.0 + i*10; v2[1] = 1.0; v2[2] = 0.0;
      v3[0] =  1.0 + i*10; v3[1] = 1.0; v3[2] = 1.0;
      vc.AddTri(v1, v2, v3);
      
      v1[0] = 0.0 + i*10; v1[1] = 1.0; v1[2] = 0.0;
      v2[0] =  i*10; v2[1] = 1.0; v2[2] = 1.0;
      v3[0] =  1.0 + i*10; v3[1] = 1.0; v3[2] = 1.0;
      vc.AddTri(v1, v2, v3);
      
      v1[0] =  1.0 + i*10; v1[1] = 0.0; v1[2] = 0.0;
      v2[0] = 1.0 + i*10; v2[1] = 1.0; v2[2] = 0.0;
      v3[0] =  1.0 + i*10; v3[1] = 1.0; v3[2] = 1.0;
      vc.AddTri(v1, v2, v3);

      v1[0] =  1.0 + i*10; v1[1] = 0.0; v1[2] = 0.0;
      v2[0] = 1.0 + i*10; v2[1] = 0.0; v2[2] = 1.0;
      v3[0] =  1.0 + i*10; v3[1] = 1.0; v3[2] = 1.0;
      vc.AddTri(v1, v2, v3);
      
      v1[0] = 0.0 + i*10; v1[1] = 0.0; v1[2] = 0.0;
      v2[0] =  i*10; v2[1] = 1.0; v2[2] = 0.0;
      v3[0] =  i*10; v3[1] = 1.0; v3[2] = 1.0;
      vc.AddTri(v1, v2, v3);
      
      v1[0] = 0.0 + i*10; v1[1] = 0.0; v1[2] = 0.0;
      v2[0] =  i*10; v2[1] = 0.0; v2[2] = 1.0;
      v3[0] =  i*10; v3[1] = 1.0; v3[2] = 1.0;
      vc.AddTri(v1, v2, v3);

      v1[0] =  1.0 + i*10; v1[1] = 0.0; v1[2] = 1.0;
      v2[0] = 1.0 + i*10; v2[1] = 1.0; v2[2] = 1.0;
      v3[0] =  i*10; v3[1] = 1.0; v3[2] = 1.0;
      vc.AddTri(v1, v2, v3);

      v1[0] =  1.0 + i*10; v1[1] = 0.0; v1[2] = 1.0;
      v2[0] =  i*10; v2[1] = 0.0; v2[2] = 1.0;
      v3[0] =  i*10; v3[1] = 1.0; v3[2] = 1.0;
      vc.AddTri(v1, v2, v3);

      v1[0] =  1.0 + i*10; v1[1] = 0.0; v1[2] = 0.0;
      v2[0] = 1.0 + i*10; v2[1] = 1.0; v2[2] = 0.0;
      v3[0] =  i*10; v3[1] = 1.0; v3[2] = 0.0;
      vc.AddTri(v1, v2, v3);

      v1[0] =  1.0 + i*10; v1[1] = 0.0; v1[2] = 0.0;
      v2[0] =  i*10; v2[1] = 0.0; v2[2] = 0.0;
      v3[0] =  i*10; v3[1] = 1.0; v3[2] = 0.0;
      vc.AddTri(v1, v2, v3);
    // std::cout<<"closing files\n";

    // cout<<"Calling finish_object\n";
    vc.EndObject();
    // vs.EndObject();

    //cout << "Inserted object " << i << "\n";
  }



// FILE *fp = fopen(argv[1], "r");
//   fscanf(fp, "%d", &num_tri);
//   for (i = 0; i < NO_OF_OBJECTS; i++) // add the objects to the library.
//   {
//     // cout<<"Reading object "<<i<<"\n";

//     vc.NewObject(&(id[i]));

//     // cout<<"Adding triangles\n";

//     for (int j = 1; j <= num_tri; j++) {
//       double v1[3], v2[3], v3[3];
//       fscanf(fp, "%lf %lf %lf", &(v1[0]), &(v1[1]), &(v1[2]));
//       fscanf(fp, "%lf %lf %lf", &(v2[0]), &(v2[1]), &(v2[2]));
//       fscanf(fp, "%lf %lf %lf", &(v3[0]), &(v3[1]), &(v3[2]));

//       vc.AddTri(v1, v2, v3);


//     }
//     // std::cout<<"closing files\n";

//     vc.EndObject();

//     // cout << "Inserted object " << i << "\n";
//   }

//   fclose(fp);

//   fp = fopen(argv[2], "r");

  double *collision_pos = new double[NO_OF_OBJECTS * 16];





double dumpTime = 0;

  // vc.EndAllObjects();
  // FILE *fp = fopen(argv[2], "r");
  double *all_trans = new double[NO_OF_OBJECTS * 16];
  int j;
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
    for (j = 0; j < NO_OF_OBJECTS / 2; j++) {
      all_trans[j * 16 + 3] = 5;
    }
    for (; j < NO_OF_OBJECTS; j++) {
      all_trans[j * 16 + 3] = -5;
    }
  for (i = 1; i <= SIMULATION_STEPS; i++) // perform the simulation.
  {
        cout << "Simulation step : " << i << "\n";
    int j;
    // double *all_trans = new double[NO_OF_OBJECTS * 16];

    // for (j = 0; j < NO_OF_OBJECTS; j++) {
    //   for (int j1 = 0; j1 < 16; j1++) {
    //     fscanf(fp, "%lf\n", &(all_trans[j * 16 + j1]));
    //   }
    // }

    
double startTime = CycleTimer::currentSeconds();
    //     for (j=0; j<NO_OF_OBJECTS; j++)
    // {
    //   for (int j1=0; j1<16; j1++){
    //     fscanf(fp, "%lf", &(all_trans[j*16+j1]));
    //   }
    // }

    // (along y-axis) half move down, half move up


    vc.UpdateAllTrans(id, NO_OF_OBJECTS, all_trans);
printf ("%.5f running time for updating \n",  CycleTimer::currentSeconds() - startTime);
double startTime2 = CycleTimer::currentSeconds();
    std::vector<int>result = vc.all_Collide(all_trans);

    dumpTime += CycleTimer::currentSeconds() - startTime;
    printf ("%.5f running time for updating \n",  CycleTimer::currentSeconds() - startTime2);
    for (int k = 0;k<result.size();k++){
      all_trans[result[k]*16+3] = -all_trans[result[k]*16+3];
    }
    
    
    // output trans matrix
   // printf("the time spend this iteration is %)
    
  }

printf ("%.5f running time for all iterations \n",  dumpTime );
  //double endTime = CycleTimer::currentSeconds();
  //std::cout << endTime - startTime - dumpTime << std::endl;
  // double seconds = difftime(endtime, now);
  // printf ("%.f running time\n", seconds);
  // cout<<" Finish Detected collision between objects\n";

  return 0;
}