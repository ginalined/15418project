
   
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


/************************************************************************\
Filename: NBody.H
--
Description: This file declares the class NBody.
\************************************************************************/

#ifndef NBODY_H
#define NBODY_H

#include "objects.h"
/************************************************************************\
Class: NBody
--
Description: The class NBody forms the core of VCollide. It
             maintains a set of three linked lists and keeps
             updating them using the sweep and prune
             algorithm. At any point in execution, the class
             NBody maintains a set of pairs of objects whose
             AABBs overlap. VCollide uses this set while making
             collision calls to RAPID. In general, this set would
             be very small compared to the n(n-1)/2 possible
             pairs.
\************************************************************************/
struct AABB;  //defined in file AABB.H
struct EndPoint;
struct cudaPoint;
const int MIN=1;
const int MAX=2;

struct Elem
{

  int   id;
  Elem  *next;
};

struct PairData
{

  int    size;   //the size of array pointed to by "arr"
  Elem   **arr;  //arr points to a dynamic array. Each element of
};


struct AABB
{

  int          id;             
  double       center[3];      // center of the AABB
  double       radius;         // radius of the AABB
  EndPoint     *lo;            // the (min_x, min_y, min_z) and  
  EndPoint     *hi;            // (max_x, max_y, max_z) corners of
  cudaPoint    *cuda_lo;
  cudaPoint    *cuda_hi;

   
};   

struct EndPoint
{

  char        minmax;     //whether the instance represents a
                          //"min" or a "max" end-point.
  AABB        *aabb;      //back pointer to the parent AABB.
  
  double      val[3];     //the coordinates of the EndPoint.
  EndPoint    *prev[3];   //for maintaining the three linked
  EndPoint    *next[3];   //lists.


};
struct cudaPoint
{
  char        minmax;     
  AABB        *aabb;      
  double      val[3];     
};

struct NBody
{
  EndPoint *head[3];  //pointers to the three linked lists, one for
                       //each axis.
  int      size;       //the size of dynamic array - AABB_arr
  AABB     **AABB_arr; //AABB_arr points to a dynamic array of pointers.
                       //Element at index "i" in this array points to an
                       //AABB corresponding to the object with id = i. If
                       //this object doesn't exist, then the pointer is NULL.
  PairData overlapping_pairs; //at any instance during the execution of the
                              //program, overlapping_pairs contains the set
                              //of pairs of ids whose AABBs overlap.
EndPoint *cudaHead[3];

};





//void updateTrans(int id, double trans[][4], NBody * obj);

#endif  /* NBODY_H */
