/************************************************************************\

  Copyright 1997 The University of North Carolina at Chapel Hill.
  All Rights Reserved.

  Permission to use, copy, modify and distribute this software
  and its documentation for educational, research and non-profit
  purposes, without fee, and without a written agreement is
  hereby granted, provided that the above copyright notice and
  the following three paragraphs appear in all copies.

  IN NO EVENT SHALL THE UNIVERSITY OF NORTH CAROLINA AT CHAPEL
  HILL BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL,
  INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
  ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION,
  EVEN IF THE UNIVERSITY OF NORTH CAROLINA HAVE BEEN ADVISED OF
  THE POSSIBILITY OF SUCH DAMAGES.


  Permission to use, copy, modify and distribute this software
  and its documentation for educational, research and non-profit
  purposes, without fee, and without a written agreement is
  hereby granted, provided that the above copyright notice and
  the following three paragraphs appear in all copies.

  THE UNIVERSITY OF NORTH CAROLINA SPECIFICALLY DISCLAIM ANY
  WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
  PURPOSE.  THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS"
  BASIS, AND THE UNIVERSITY OF NORTH CAROLINA HAS NO OBLIGATION
  TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
  MODIFICATIONS.


   --------------------------------- 
  |Please send all BUG REPORTS to:  |
  |                                 |
  |   geom@cs.unc.edu               |
  |                                 |
   ---------------------------------
  
     
  The authors may be contacted via:

  US Mail:  A. Pattekar/J. Cohen/T. Hudson/S. Gottschalk/M. Lin/D. Manocha
            Department of Computer Science
            Sitterson Hall, CB #3175
            University of N. Carolina
            Chapel Hill, NC 27599-3175
	    
  Phone:    (919)962-1749
	    
  EMail:    geom@cs.unc.edu

\************************************************************************/
#include <iostream>
#include <stdlib.h>
#include "VInternal.H"
#include <stdio.h>

//#include "VCollide.H"
using namespace std;

const int NO_OF_OBJECTS=20;      //number of instances
const int SIMULATION_STEPS=99;  //number of steps in the simulation.

int main(int argc, char *argv[])
{

  if (argc != 3)
    {
      cerr<<argv[0]<<": USAGE: "<<argv[0]<<" <input-file> <transformation-file>\n";
      exit(1);
    }
  
  int num_tri;
  VCInternal vc(NO_OF_OBJECTS + 20);
  int id[NO_OF_OBJECTS+ 20];
  
  int i;
  for (i=0; i<NO_OF_OBJECTS; i++)  //add the objects to the library.
    {
      cout<<"Reading object "<<i<<"\n";

      vc.NewObject(&(id[i]));
      cout<<"Adding triangles\n";
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
  std::cout<<"closing files\n";
      
      fclose(fp);
      
      cout<<"Calling finish_object\n";
      vc.EndObject();
      
      
      cout<<"Inserted object "<<i<<"\n";
    }
  
  
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
    // for (j=0; j<NO_OF_OBJECTS; j++){
    //   vc.UpdateTrans(id[j], all_trans[j]);
    // }

  //vc.UpdateTrans(id, NO_OF_OBJECTS, all_trans);
      
      vc.Collide();  //perform collision test.
      
      //report the results.
      int VCReportSize=10;
      VCReportType *vcrep = new VCReportType[VCReportSize];
      
      int no_of_colliding_pairs = vc.Report(VCReportSize, vcrep);
      
      if (no_of_colliding_pairs > VCReportSize)
	{
	  VCReportSize=no_of_colliding_pairs;
	  delete vcrep;
	  vcrep = new VCReportType[VCReportSize];
	  no_of_colliding_pairs = vc.Report(VCReportSize, vcrep);
	}
    
      for (j=0; j<no_of_colliding_pairs; j++)
	      cout<<"Detected collision between objects "<<vcrep[j].id1<<" and "<<vcrep[j].id2<<"\n";

    }
    cout<<" Finish Detected collision between objects\n";

   
    return 0;
}
