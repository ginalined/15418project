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
Filename: NBody.C
--
Description: This file implements the member functions of the class NBody.

\************************************************************************/

#include <iostream>
#include <math.h>
#include "objects.h"
#include "NBody.H"
#include <math.h>  
const int REVERSE  = 1;//whether the direction of movement of the interval
const int FORWARD  = 2;//along the list is forward (FORWARD), reverse (REVERSE)
const int NOCHANGE = 3;//or there is no movement at all (NOCHANGE).

const int DEFAULT_SIZE=10; //default number of AABBs. Arbitrarily set here to
                           //10.

inline double GT(double a, double b)
{
  return (( (a) > (b) ) ? (a) : (b));
}

NBody::NBody()  //constructor.
{
  elist[0] = NULL;
  elist[1] = NULL;
  elist[2] = NULL;
  size = DEFAULT_SIZE;
  AABB_arr = new AABB*[size];  //allocate the dynamic array and initialize
  int i;
  for (i=0; i<size; i++)      //all its elements to NULL.
    AABB_arr[i] = NULL;
  
}


NBody::~NBody()  //destructor.
{
  int i;
  for (i=0; i<size; i++)
    {
      if (AABB_arr[i])
	{
	  delete AABB_arr[i]->hi;
	  delete AABB_arr[i]->lo;
	  delete AABB_arr[i];
	  AABB_arr[i] = NULL;
	}
    }
}





// void NBody::AddObject(int id, Object *b) //add a new object
// {
//   AABB *curr = new AABB;
  
//   curr->id = id; //set the id to the given value.
  
//   //The centroid of the object is computed and this is taken to be the
//   //center of the AABB. 找到AABB的中心 
//   curr->center[0] = curr->center[1] = curr->center[2] = 0.0;
  
//   int i;
//   for (i=0; i<(b->num_tris); i++)
//     {
//       curr->center[0] += b->tris[i].p1[0] + b->tris[i].p2[0] + b->tris[i].p3[0];
//       curr->center[1] += b->tris[i].p1[1] + b->tris[i].p2[1] + b->tris[i].p3[1];
//       curr->center[2] += b->tris[i].p1[2] + b->tris[i].p2[2] + b->tris[i].p3[2];
//     }
  
//   curr->center[0] /= (3*b->num_tris); 
//   curr->center[1] /= (3*b->num_tris);
//   curr->center[2] /= (3*b->num_tris);
//   //------------------
  

//   //The "radius" of the AABB is computed as the maximum distance of the AABB
//   //center from any of the vertices of the object.
//   curr->radius = 0.0;

//   for (i=0; i<(b->num_tris); i++)
//     {
//       double cur_rad1_sq = 0;
//       double cur_rad2_sq = 0;
//       double cur_rad3_sq = 0;
//       for (int w=0; w<3; w++)
//       {
//         double my_num1 = curr->center[w] - b->tris[i].p1[w];
//         double my_num2 = curr->center[w] - b->tris[i].p2[w];
//         double my_num3 = curr->center[w] - b->tris[i].p3[w];
//         cur_rad1_sq += pow(my_num1, 2);
//         cur_rad2_sq += pow(my_num2, 2);
//         cur_rad3_sq += pow(my_num3, 2);
//       }
               
      
//       double max_rad_sq = GT(cur_rad1_sq, GT(cur_rad2_sq,cur_rad3_sq));
      
//       curr->radius = GT(max_rad_sq, curr->radius);
      
//     }

//   curr->radius = sqrt(curr->radius);
//   curr->radius *= 1.0001;  //add a 0.01% buffer.
//   curr->lo = new EndPoint;
//   curr->hi = new EndPoint;
//   curr->lo->minmax = MIN;
//   curr->hi->minmax = MAX;
//   curr->lo->aabb = curr;
//   curr->hi->aabb = curr;
//   double min[3], max[3];
  
//   for (int w=0; w<3; w++){
//   min[w] = curr->center[w] - curr->radius; 
//   max[w] = curr->center[w] + curr->radius;
//   curr->lo->val[w] = min[w];
//   curr->hi->val[w] = max[w];
//   }


//   for (i=0; i<size; i++)      //Now, check the overlap of this AABB with 
//     {                         //with all other AABBs and add the pair to
//       if (AABB_arr[i])        //the set of overlapping pairs if reqd.
// 	if (overlaps(curr, AABB_arr[i]))
// 	  add_overlap_pair(curr->id, i);
//     }

//   if (id >= size)    //increase the size of the dynamic array if necessary.
//       {
// 	int newsize = (id >= 2*size) ? (id+1) : 2*size;

// 	AABB **temp = new AABB*[newsize];
// 	int i;
// 	for (i=0; i<size; i++)
// 	  temp[i] = AABB_arr[i];
// 	for (i=size; i<newsize; i++)
// 	  temp[i] = NULL;
// 	delete [] AABB_arr;
// 	AABB_arr = temp;
// 	size = newsize;
//       }
  
//   AABB_arr[id] = curr;  //finally, insert the AABB in AABB_arr.
  

//   //Now, for each of the three co-ordinates, insert the interval
//   //in the correspoding list. 
//   int coord;
//   for (coord=0; coord <3; coord++)
//     {
//       EndPoint *current = elist[coord];
      
//       //first insert the "hi" endpoint.
//       if (current == NULL)    //if the list is empty, insert in front.
// 	{
// 	  elist[coord] = curr->hi;
// 	  curr->hi->prev[coord] = curr->hi->next[coord] = NULL;
// 	}
//       else  //otherwise, find the correct location in the list and
// 	{   //insert there. Note: the list is sorted.
// 	  while ( (current->next[coord] != NULL) && (current->val[coord] < curr->hi->val[coord]) )
// 	    current = current->next[coord];
	  
	  
// 	  if (current->val[coord] >= curr->hi->val[coord])
// 	    {
// 	      curr->hi->prev[coord] = current->prev[coord];
// 	      curr->hi->next[coord] = current;
// 	      if (current->prev[coord] == NULL)
// 		elist[coord] = curr->hi;
// 	      else
// 		current->prev[coord]->next[coord] = curr->hi;
	      
// 	      current->prev[coord] = curr->hi;
// 	    }
// 	  else
// 	    {
// 	      curr->hi->prev[coord] = current;
// 	      curr->hi->next[coord] = NULL;
// 	      current->next[coord] = curr->hi;
// 	    }
// 	}
      
//       //now, insert the "lo" endpoint.
//       current = elist[coord];
      
//       //at this point, the list cannot be empty since we have already 
//       //inserted the "hi" endpoint. So, we straightaway look for the 
//       //correct location in the non-empty list and insert at that location.
//       while ( (current->next[coord] != NULL) && (current->val[coord] < curr->lo->val[coord]) )
// 	current = current->next[coord];
      
//       if (current->val[coord] >= curr->lo->val[coord])
// 	{
// 	  curr->lo->prev[coord] = current->prev[coord];
// 	  curr->lo->next[coord] = current;
// 	  if (current->prev[coord] == NULL)
// 	    elist[coord] = curr->lo;
// 	  else
// 	    current->prev[coord]->next[coord] = curr->lo;
	  
// 	  current->prev[coord] = curr->lo;
// 	}
//       else
// 	{
// 	  curr->lo->prev[coord] = current;
// 	  curr->lo->next[coord] = NULL;
// 	  current->next[coord] = curr->lo;
// 	}
      
//     }
  
// }



// void NBody::DeleteObject(int id) //deleting an AABB with given id.
// {
//   if (id >= size)
//     {
//       //cerr<<"Should not get here since VCollide should send only valid ids\n";
//       return;
//     }
  
//   if (AABB_arr[id] == NULL)
//     {
//       //cerr<<"Should not get here since VCollide should send only valid ids\n";
//       return;
//     }
  
//   AABB *curr = AABB_arr[id];  //this is the AABB to be deleted.
//   AABB_arr[id] = NULL;        //remove it from the AABB array.
  
//   //first, we delete all the three intervals from the corresponding lists.
//   int coord;
//   for (coord=0; coord<3; coord++)
//     {
//       //first delete the "lo" endpoint of the interval.
//       if (curr->lo->prev[coord] == NULL)
// 	elist[coord] = curr->lo->next[coord];
//       else
// 	curr->lo->prev[coord]->next[coord] = curr->lo->next[coord];
      
//       curr->lo->next[coord]->prev[coord] = curr->lo->prev[coord];
      
//       //then, delete the "hi" endpoint.
//       if (curr->hi->prev[coord] == NULL)
// 	elist[coord] = curr->hi->next[coord];
//       else
// 	curr->hi->prev[coord]->next[coord] = curr->hi->next[coord];
      
//       if (curr->hi->next[coord] != NULL)
// 	curr->hi->next[coord]->prev[coord] = curr->hi->prev[coord];
      
//     }
  
//   //delete all entries involving this id from the set of 
//   //overlapping pairs.
//   overlapping_pairs.DelPairsInvolvingId(id);
  
//   //de-allocate the memory
//   delete curr->lo;
//   delete curr->hi;
//   delete curr;
// }
