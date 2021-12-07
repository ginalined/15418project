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

#include "VInternal.H"
#include <iostream>
#include <stdlib.h>
//#include "VCollide.H"
using namespace std;

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
  int id[2];

  int i;
  for (i = 0; i < 2; i++) // create both the objects.
  {
    vc.NewObject(&id[i]);

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

    vc.EndObject();
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
  }


  cout << " Finish Detected collision between objects\n";
  return 0;
}
