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
Filename: VInternal.C
--
Description: This file implements the member functions of the class vinternal.c

\************************************************************************/



#include <iostream>
#include <string.h>     //for memset and memcpy.
#include "VInternal.H"
#include <cuda.h>
#include <cuda_runtime.h>
#include <driver_functions.h>
#include "CycleTimer.h"
#include "objects.h"
#include "NBody.H"
#include <math.h>  
const int DEFAULT_SIZE=10; //some arbitrary default size for "vc_objects" array.

  void add_overlap_pair(int id1, int id2, NBody * obj) //add a pair to the set of
    {                                     //overlapping pairs.
      if (id1 != id2)
	obj->overlapping_pairs.AddPair(id1, id2);
    }
  
  void del_overlap_pair(int id1, int id2, NBody * obj) //delete a pair from the set.
    {
      if (id1 != id2)
	obj->overlapping_pairs.DelPair(id1, id2);
    }

int overlaps(AABB *obj1, AABB *obj2) //to check it the two AABBs overlap.
{
  int coord;
  for (coord=0; coord<3; coord++)
    {
      if (obj1->lo->val[coord] < obj2->lo->val[coord])
	{
	  if (obj2->lo->val[coord] > obj1->hi->val[coord])
	    return 0;
	}
      else
	{
	  if (obj1->lo->val[coord] > obj2->hi->val[coord])
	    return 0;
	}
    }
  
  return 1;
}
void add_node(EndPoint* node, int dim, EndPoint* prevNode ) {
        EndPoint* temp = prevNode->next[dim];
        node->next[dim] = temp;
        node->prev[dim] = prevNode;
        prevNode->next[dim] = node;
        if (temp != NULL)
          temp->prev[dim] = node;
    }
    
void delete_node(EndPoint* delnode, int dim) {
        EndPoint* delprev = delnode->prev[dim];
        EndPoint* delnext = delnode->next[dim];
        if (delprev !=NULL)
          delprev->next[dim] = delnext;
        if (delnext != NULL)
          delnext->prev[dim] = delprev;
    }
void updatetempTrans(int id, double trans[][4], NBody * obj){
  
  AABB *current = obj->AABB_arr[id]; //the given object exists !
  
  //compute the new position of the AABB center.
  double new_center[3], min[3], max[3];
  AABB dummy;       //we need these so that we can use the same function
  EndPoint lo, hi; 
  dummy.lo = &lo;
  dummy.hi = &hi;
  lo.minmax = MIN;
  hi.minmax = MAX;
  lo.aabb = &dummy;
  hi.aabb = &dummy;

  for (int dim = 0; dim < 3; dim ++){
    new_center[dim] = current->center[0] * trans[dim][0] + current->center[1] * trans[dim][1] + current->center[2] * trans[dim][2] + trans[dim][3];
    min[dim] = lo.val[dim] = new_center[dim] - current->radius;
    max[dim] =  hi.val[dim] = new_center[dim] + current->radius;
    
  }
  
  //update all the three lists by moving the endpoint to correct position.
  int coord;
  for (coord=0; coord<3; coord++)
    {
      int direction;
      EndPoint *temp;
      
      //set the direction of motion of the endpoint along the list.
      if (current->lo->val[coord] > min[coord])
	        direction = REVERSE;
      else if (current->lo->val[coord] <min[coord])
	        direction = FORWARD;
      else
	        direction = NOCHANGE;
      
      if (direction == REVERSE) //backward motion....
	{
	  //first update the "lo" endpoint of the interval
	  if (current->lo->prev[coord] != NULL)
	    {
	      temp = current->lo;
	      while ((temp != NULL) && (temp->val[coord] > min[coord]))
		{
		  if (temp->minmax == MAX){
		    if (overlaps(temp->aabb, &dummy))
		      add_overlap_pair(temp->aabb->id, current->id, obj);
      } 
		  temp = temp->prev[coord];
		}

    delete_node(current->lo,coord); 

  // if it is the smallest
	  if (temp == NULL)
		{
		  current->lo->prev[coord] = NULL;
		  current->lo->next[coord] = obj->elist[coord];
		  obj->elist[coord]->prev[coord] = current->lo;
		  obj->elist[coord] = current->lo;
		}
	  else
		{
      add_node(current->lo, coord, temp);
		} 
	    }
	  
	  current->lo->val[coord] = min[coord];
	  
	  //then update the "hi" endpoint of the interval.
	  if (current->hi->val[coord] != max[coord])
	    {
	      temp = current->hi;
	      
	  while (temp->val[coord] > max[coord])
		{
		if ( (temp->minmax == MIN) && (overlaps(temp->aabb, current)) )
		    del_overlap_pair(temp->aabb->id, current->id, obj);
		  temp = temp->prev[coord];
		}
	  
    delete_node(current->hi, coord);
    add_node(current->hi, coord, temp);

	      current->hi->val[coord] = max[coord];
	    }
	}
      else if (direction == FORWARD) //forward motion....
	{
	  //here, we first update the "hi" endpoint.
	  if (current->hi->next[coord] != NULL)
	    {
	      temp = current->hi;
	      while ( (temp->next[coord] != NULL) && (temp->val[coord] < max[coord]) )
		{
		  if (temp->minmax == MIN)
		    if (overlaps(temp->aabb, &dummy))
		      add_overlap_pair(temp->aabb->id, current->id, obj);
		  
		  temp = temp->next[coord];
		}
	  
	  if (temp->val[coord] < max[coord])
		{
      delete_node(current->hi, coord);
		  current->hi->prev[coord] = temp;
		  current->hi->next[coord] = NULL;
		  temp->next[coord] = current->hi;
		}
	  else if (current->hi->val[coord] != max[coord])
		{
      delete_node(current->hi, coord);
      add_node(current->hi, coord, temp->prev[coord]);
		}
	    }
	  current->hi->val[coord] = max[coord];
	  
	  //then, update the "lo" endpoint of the interval.
	  temp = current->lo;
	  
	  while (temp->val[coord] < min[coord])
	    {
	      if ( (temp->minmax == MAX) && (overlaps(temp->aabb, current)) )
		del_overlap_pair(temp->aabb->id, current->id, obj);
	      
	      temp = temp->next[coord];
	    }
	  
	  if (current->lo->prev[coord] != NULL)
	    current->lo->prev[coord]->next[coord] = current->lo->next[coord];
	  else
	    obj->elist[coord] = current->lo->next[coord];
	  current->lo->next[coord]->prev[coord] = current->lo->prev[coord];
	  current->lo->prev[coord] = temp->prev[coord];
	  current->lo->next[coord] = temp;
	  if (temp->prev[coord] != NULL)
	    temp->prev[coord]->next[coord] = current->lo;
	  else
	    obj->elist[coord] = current->lo;
	  temp->prev[coord] = current->lo;
	  current->lo->val[coord] = min[coord];
	}   
    }
}

inline double GT(double a, double b)
{
  return (( (a) > (b) ) ? (a) : (b));
}

void AddObject(int id, Object *b, NBody * obj) //add a new object
{
  AABB *curr = new AABB;
  
  curr->id = id; //set the id to the given value.
  
  //The centroid of the object is computed and this is taken to be the
  //center of the AABB. 找到AABB的中心 
  curr->center[0] = curr->center[1] = curr->center[2] = 0.0;
  
  int i;
  for (i=0; i<(b->num_tris); i++)
    {
      curr->center[0] += b->tris[i].p1[0] + b->tris[i].p2[0] + b->tris[i].p3[0];
      curr->center[1] += b->tris[i].p1[1] + b->tris[i].p2[1] + b->tris[i].p3[1];
      curr->center[2] += b->tris[i].p1[2] + b->tris[i].p2[2] + b->tris[i].p3[2];
    }
  
  curr->center[0] /= (3*b->num_tris); 
  curr->center[1] /= (3*b->num_tris);
  curr->center[2] /= (3*b->num_tris);
  //------------------
  

  //The "radius" of the AABB is computed as the maximum distance of the AABB
  //center from any of the vertices of the object.
  curr->radius = 0.0;

  for (i=0; i<(b->num_tris); i++)
    {
      double cur_rad1_sq = 0;
      double cur_rad2_sq = 0;
      double cur_rad3_sq = 0;
      for (int w=0; w<3; w++)
      {
        double my_num1 = curr->center[w] - b->tris[i].p1[w];
        double my_num2 = curr->center[w] - b->tris[i].p2[w];
        double my_num3 = curr->center[w] - b->tris[i].p3[w];
        cur_rad1_sq += pow(my_num1, 2);
        cur_rad2_sq += pow(my_num2, 2);
        cur_rad3_sq += pow(my_num3, 2);
      }
               
      
      double max_rad_sq = GT(cur_rad1_sq, GT(cur_rad2_sq,cur_rad3_sq));
      
      curr->radius = GT(max_rad_sq, curr->radius);
      
    }

  curr->radius = sqrt(curr->radius);
  curr->radius *= 1.0001;  //add a 0.01% buffer.
  curr->lo = new EndPoint;
  curr->hi = new EndPoint;
  curr->lo->minmax = MIN;
  curr->hi->minmax = MAX;
  curr->lo->aabb = curr;
  curr->hi->aabb = curr;
  double min[3], max[3];
  
  for (int w=0; w<3; w++){
  min[w] = curr->center[w] - curr->radius; 
  max[w] = curr->center[w] + curr->radius;
  curr->lo->val[w] = min[w];
  curr->hi->val[w] = max[w];
  }


  for (i=0; i<obj->size; i++)      //Now, check the overlap of this AABB with 
    {                         //with all other AABBs and add the pair to
      if (obj->AABB_arr[i])        //the set of overlapping pairs if reqd.
	if (overlaps(curr, obj->AABB_arr[i]))
	  add_overlap_pair(curr->id, i, obj);
    }

  if (id >= obj->size)    //increase the size of the dynamic array if necessary.
      {
	int newsize = (id >= 2*obj->size) ? (id+1) : 2*obj->size;

	AABB **temp = new AABB*[newsize];
	int i;
	for (i=0; i<obj->size; i++)
	  temp[i] = obj->AABB_arr[i];
	for (i=obj->size; i<newsize; i++)
	  temp[i] = NULL;
	delete [] obj->AABB_arr;
	obj->AABB_arr = temp;
	obj->size = newsize;
      }
  
  obj->AABB_arr[id] = curr;  //finally, insert the AABB in AABB_arr.
  

  //Now, for each of the three co-ordinates, insert the interval
  //in the correspoding list. 
  int coord;
  for (coord=0; coord <3; coord++)
    {
      EndPoint *current = obj->elist[coord];
      
      //first insert the "hi" endpoint.
      if (current == NULL)    //if the list is empty, insert in front.
	{
	 obj->elist[coord] = curr->hi;
	  curr->hi->prev[coord] = curr->hi->next[coord] = NULL;
	}
      else  //otherwise, find the correct location in the list and
	{   //insert there. Note: the list is sorted.
	  while ( (current->next[coord] != NULL) && (current->val[coord] < curr->hi->val[coord]) )
	    current = current->next[coord];
	  
	  
	  if (current->val[coord] >= curr->hi->val[coord])
	    {
	      curr->hi->prev[coord] = current->prev[coord];
	      curr->hi->next[coord] = current;
	      if (current->prev[coord] == NULL)
		obj->elist[coord] = curr->hi;
	      else
		current->prev[coord]->next[coord] = curr->hi;
	      
	      current->prev[coord] = curr->hi;
	    }
	  else
	    {
	      curr->hi->prev[coord] = current;
	      curr->hi->next[coord] = NULL;
	      current->next[coord] = curr->hi;
	    }
	}
      
      //now, insert the "lo" endpoint.
      current = obj->elist[coord];
      
      //at this point, the list cannot be empty since we have already 
      //inserted the "hi" endpoint. So, we straightaway look for the 
      //correct location in the non-empty list and insert at that location.
      while ( (current->next[coord] != NULL) && (current->val[coord] < curr->lo->val[coord]) )
	current = current->next[coord];
      
      if (current->val[coord] >= curr->lo->val[coord])
	{
	  curr->lo->prev[coord] = current->prev[coord];
	  curr->lo->next[coord] = current;
	  if (current->prev[coord] == NULL)
	    obj->elist[coord] = curr->lo;
	  else
	    current->prev[coord]->next[coord] = curr->lo;
	  
	  current->prev[coord] = curr->lo;
	}
      else
	{
	  curr->lo->prev[coord] = current;
	  curr->lo->next[coord] = NULL;
	  current->next[coord] = curr->lo;
	}
      
    }
  
}



void deleteObjects(int id, NBody * obj) //deleting an AABB with given id.
{
  if (id >= obj->size)
    {
      //cerr<<"Should not get here since VCollide should send only valid ids\n";
      return;
    }
  
  if (obj->AABB_arr[id] == NULL)
    {
      //cerr<<"Should not get here since VCollide should send only valid ids\n";
      return;
    }
  
  AABB *curr = obj->AABB_arr[id];  //this is the AABB to be deleted.
  obj->AABB_arr[id] = NULL;        //remove it from the AABB array.
  
  //first, we delete all the three intervals from the corresponding lists.
  int coord;
  for (coord=0; coord<3; coord++)
    {
      //first delete the "lo" endpoint of the interval.
      if (curr->lo->prev[coord] == NULL)
	obj->elist[coord] = curr->lo->next[coord];
      else
	curr->lo->prev[coord]->next[coord] = curr->lo->next[coord];
      
      curr->lo->next[coord]->prev[coord] = curr->lo->prev[coord];
      
      //then, delete the "hi" endpoint.
      if (curr->hi->prev[coord] == NULL)
	obj->elist[coord] = curr->hi->next[coord];
      else
	curr->hi->prev[coord]->next[coord] = curr->hi->next[coord];
      
      if (curr->hi->next[coord] != NULL)
	curr->hi->next[coord]->prev[coord] = curr->hi->prev[coord];
      
    }
  
  //delete all entries involving this id from the set of 
  //overlapping pairs.
  obj->overlapping_pairs.DelPairsInvolvingId(id);
  
  //de-allocate the memory
  delete curr->lo;
  delete curr->hi;
  delete curr;
}




VCInternal::VCInternal()
{
  state = VCstate_default;
  next_id = 0;
  
  size = DEFAULT_SIZE;                   //set the size of the array.
  vc_objects = new VCObject*[size]; //allocate the array.
  int i;
  for (i=0; i<size; i++)
    vc_objects[i] = NULL;
  
  disabled.Clear();  //to begin with, no pairs are disabled.
}


VCInternal::~VCInternal()
{

  //deallocate the memory.
  int i;
  for (i=0; i<size; i++)
    {
      if (vc_objects[i])
	{
	  delete vc_objects[i]->b;
	  delete vc_objects[i];
	}
    }
  delete [] vc_objects;
}


//1. check if the size fit in
//2 assign the object an id and activate the object
int VCInternal::NewObject(int *id) //create a new object in the database.
{

  //increase the size of the "vc_objects" array if required.
  if (next_id >= size) 
    {
      int newsize = (next_id >= 2*size) ? (next_id+1) : 2*size;
      VCObject **temp = new VCObject*[newsize];
      int i;
      for (i=0; i<size; i++)
	temp[i] = vc_objects[i];
      for (i=size; i<newsize; i++)
	temp[i] = NULL;
      delete [] vc_objects;
      vc_objects = temp;
      size = newsize;
      
    }
  
  //allocate a new object.
  vc_objects[next_id] = new VCObject;
  
  *id = next_id;  //for returning the id generated by VCollide.
  current_id = next_id;
  vc_objects[next_id]->id = next_id;
  vc_objects[next_id]->b = new Object;
  vc_objects[next_id]->b->BeginModel();
  //_state = 1;//default the object is activate
  next_id++; 
  
  return 0;
}

int VCInternal::AddTri(double v1[], double v2[], double v3[]) 
{                     

  vc_objects[current_id]->b->AddTri(v1, v2, v3, 0);  //add triangle.
  return 0;
}

// 1. add current object to n body
// 2. have RAPID build the OBB tree.
// 3. initialize trans
int VCInternal::EndObject(void)
{  

  AddObject(current_id, vc_objects[current_id]->b, &nbody);
  
  vc_objects[current_id]->b->EndModel();

  memset( ( (void *)vc_objects[current_id]->trans), 0, 16*sizeof(double) );
  vc_objects[current_id]->trans[0][0] = 1.0;
  vc_objects[current_id]->trans[1][1] = 1.0;
  vc_objects[current_id]->trans[2][2] = 1.0;
  vc_objects[current_id]->trans[3][3] = 1.0;
  
  return 0;
  
}


int VCInternal::UpdateTrans(int id, double t[][4])
{           

  VCObject *current = vc_objects[id];
  
  //update the private copy of the transformation matrix.
  memcpy((void *)current->trans, (void *)t, 16*sizeof(double));
  
  //have the nbody database update itself appropriately.
  //updateteTrans(current->id, t, &nbody);
  updatetempTrans(current->id, t, &nbody);
  

  return 0;
  
}



int VCInternal::ActivatePair(int id1, int id2)
{
      disabled.DelPair(id1, id2);
      return 0;
}

int VCInternal::DeactivatePair(int id1, int id2)
{

      if (id1!=id2)
	disabled.AddPair(id1, id2);
      
    return 0;
}


int VCInternal::DeleteObject(int id) //delete an object from the database.
{


      delete vc_objects[id]->b; //delete the RAPID box.
      delete vc_objects[id];    //delete the object.
      vc_objects[id] = NULL; 
      
      disabled.DelPairsInvolvingId(id);

      deleteObjects(id, &nbody); //delete the object from the nbody database.
      return 0;

  
}



int VCInternal::Collide(void)  //perform collision detection.
{

  
  //Clear the results from earlier collision tests.
  report_data.Clear();
  
  //Simultaneously traverse the "overlapping_pairs" database and the 
  //"disabled_pairs" database and make calls to the RAPID collision
  //detection routine where required.
  int i;
   //std::cout<< nbody.overlapping_pairs.size<< std::endl;
  for (i=0; i<nbody.overlapping_pairs.size; i++)
    {
      Elem *curr_ovrlp = nbody.overlapping_pairs.arr[i];
      
      Elem *curr_disabled = i<disabled.size ? disabled.arr[i]: NULL;
      
      
	  while (curr_ovrlp != NULL)
		{

		  while (curr_disabled && curr_disabled->id <= curr_ovrlp->id)
			      curr_disabled = curr_disabled->next;

      double R1[3][3], T1[3], R2[3][3], T2[3];
      for (int index = 0; index < 9; index++){
        int x = index/3;
        int y = index%3;
        R1[x][y] = vc_objects[i]->trans[x][y];
        R2[x][y] = vc_objects[curr_ovrlp->id]->trans[x][y];
      }

      for (int index = 0; index < 3; index++){
        T1[index] = vc_objects[i]->trans[index][3];
        T2[index] = vc_objects[curr_ovrlp->id]->trans[index][3];
      }

			      //call the RAPID collision detection routine.
			::Collide(R1, T1, vc_objects[i]->b, R2, T2, vc_objects[curr_ovrlp->id]->b, FIRST_CONTACT);
			      
			      //if there is a collision, then add the pair to the
			      //collision report database.
			if (Object_num_contacts != 0)
			      report_data.AddPair(i, vc_objects[curr_ovrlp->id]->id);
			      
		  curr_ovrlp = curr_ovrlp->next;
		}
    }
  return 0;
}


//report the results of collision detection.
//sz is the size of the array pointed to by vcrep. If sz is less than
//the number of collision pairs, then fill the array with first sz number
//of collision pairs.
//Returns the total number of collision pairs.
int VCInternal::Report(int sz, VCReportType *vcrep)
{
  int no_of_colliding_pairs=0;
  int vc_rep_count = 0;
  
  int i;
  for (i=0; i<report_data.size; i++)
    {
      Elem *current;
      for (current=report_data.arr[i]; current != NULL; current=current->next)
	{
	  no_of_colliding_pairs++;
	  if (vc_rep_count <sz) //if the array is not full yet, then 
	    {                   //fill the data in it.
	      vcrep[vc_rep_count].id1 = i;
	      vcrep[vc_rep_count].id2 = current->id;
	      vc_rep_count++;
	    }
	}
      
    }
  return no_of_colliding_pairs;
}


