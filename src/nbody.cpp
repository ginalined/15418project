#include <iostream>
#include <stdlib.h>
#include "VInternal.H"
#include <stdio.h>
#include <time.h>  
#include "CycleTimer.h"
//#include "VCollide.H"
using namespace std;

const int NO_OF_OBJECTS=32;      //number of instances
const int SIMULATION_STEPS=1;  //number of steps in the simulation.
const int SCREEN_SIZE = 100;
int main(int argc, char *argv[])
{

    double startTime = CycleTimer::currentSeconds();

    
  if (argc != 3)
    {
      cerr<<argv[0]<<": USAGE: "<<argv[0]<<" <input-file> <transformation-file>\n";
      exit(1);
    }
  
  int num_tri;
  VCInternal vc(NO_OF_OBJECTS, SCREEN_SIZE);
  int id[NO_OF_OBJECTS];
  
  int i;

  
  for (i=0; i<NO_OF_OBJECTS; i++)  //add the objects to the library.
  {
      //cout<<"Reading object "<<i<<"\n";

      vc.NewObject(&(id[i]));
      //cout<<"Adding triangles\n";
      FILE *fp = fopen(argv[1], "r");
      fscanf(fp, "%d", &num_tri);
      
      for (int j=1; j<=num_tri; j++)
      {
        double v1[3], v2[3], v3[3];
        fscanf(fp, "%lf %lf %lf", &(v1[0]), &(v1[1]), &(v1[2]));
        fscanf(fp, "%lf %lf %lf", &(v2[0]), &(v2[1]), &(v2[2]));
        fscanf(fp, "%lf %lf %lf", &(v3[0]), &(v3[1]), &(v3[2]));
        
        vc.AddTri(v1, v2, v3);
      }
  //std::cout<<"closing files\n";
      
      fclose(fp);
      
      //cout<<"Calling finish_object\n";
      vc.EndObject();
      
      
      cout<<"Inserted object "<<i<<"\n";
    }

  

  //vc.EndAllObjects();
  FILE *fp = fopen(argv[2], "r");

  for (i=1; i<=SIMULATION_STEPS; i++)  //perform the simulation.
  {
      cout<<"Simulation step : "<<i<<"\n";
      int j;
      double* all_trans = new double[NO_OF_OBJECTS*16];

      for (j=0; j<NO_OF_OBJECTS; j++)
	{
    for (int j1=0; j1<16; j1++){
      fscanf(fp, "%lf", &(all_trans[j*16+j1]));
    }
	}

  vc.UpdateAllTrans(id, NO_OF_OBJECTS, all_trans);
  vc.Collide();  
      
    
    }
 
  double endTime = CycleTimer::currentSeconds();
  std::cout << endTime << ' '<< startTime << std::endl;
  // double seconds = difftime(endtime, now);
  // printf ("%.f running time\n", seconds);
  // cout<<" Finish Detected collision between objects\n";

   
    return 0;
}
