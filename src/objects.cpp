

// to silence the compiler's complaints about unreferenced identifiers.
//static void r1(char *f){  r1(f);  r1(Objecttag_data);  r1(Object_version);}


#include <stdio.h>
#include <stdlib.h>

#include "objects.h"
#include "moments.H"

extern int Object_initialized;
void Object_initialize();

static moment *Object_moment = 0;
static tri *Object_tri = 0;
static box *Object_boxes = 0;
static int Object_boxes_inited = 0;

double Object_mR[3][3];
double Object_mT[3];
double Object_ms;

int Object_first_contact;
int Object_num_box_tests;
int Object_num_tri_tests;
int Object_num_contacts;

int Object_num_cols_alloced = 0;
collision_pair *Object_contact = 0;

int add_collision(int id1, int id2);


Object::Object() 
{
  if (!Object_initialized) Object_initialize();

  b = 0;
  num_boxes_alloced = 0;

  tris = 0;
  num_tris = 0;
  num_tris_alloced = 0;
}

Object::~Object()
{
  if (!Object_initialized) Object_initialize();

  // the boxes pointed to should be deleted.
  delete [] b;

  // the triangles pointed to should be deleted.
  delete [] tris;
}

int Object_initialized = 0;

void
Object_initialize()
{
  Object_num_box_tests = 0;
  Object_num_contacts = 0;
  Object_contact = 0;

  Object_initialized = 1;
}             



int
Object::BeginModel()
{
  int bs = build_state;
  
  if (!Object_initialized) Object_initialize();
  
  // free whatever storage we had.  Remember, it's okay to delete null
  // pointers in C++, so we don't have to check them first.
  delete [] b;   b = 0;
  num_boxes_alloced = 0;
  delete [] tris;   tris = 0;
  num_tris = 0;
  num_tris_alloced = 0;

  return 0;
}

int
Object::EndModel()
{
  int myrc = build_hierarchy();
  return 0;
}

int
Object::AddTri(const double *p1, const double *p2, const double *p3, int id)
{

  // first make sure that we haven't filled up our allocation.
  // if we have, allocate a new array of twice the size, and copy
  // the old data to it.

  if (num_tris == num_tris_alloced)
    {
      // decide on new size -- accounting for first time, where none are 
      // allocated
      int n = num_tris_alloced*2;
      if (n == 0) n = 1;

      // make new array, and copy the old one to it
      tri *t = new tri[n];

      int i;
      for(i=0; i<num_tris; i++) {
        t[i] = tris[i]; 
      }

      // free the old array and reassign.  
      delete [] tris;
      tris = t;
      
      // update the allocation counter.
      num_tris_alloced = n;
    }

  // now copy the new tri into the array
  tris[num_tris].p1[0] = p1[0];
  tris[num_tris].p1[1] = p1[1];
  tris[num_tris].p1[2] = p1[2];
  tris[num_tris].p2[0] = p2[0];
  tris[num_tris].p2[1] = p2[1];
  tris[num_tris].p2[2] = p2[2];
  tris[num_tris].p3[0] = p3[0];
  tris[num_tris].p3[1] = p3[1];
  tris[num_tris].p3[2] = p3[2];
  tris[num_tris].id = id;

  // update the counter
  num_tris++;

  return 0;
}



int
Object::build_hierarchy()
{
  // allocate the boxes and set the box list globals

  num_boxes_alloced = num_tris * 2;
  b = new box[num_boxes_alloced];
  Object_boxes = b;
  Object_boxes_inited = 1;   // we are in process of initializing b[0].
  
  // Determine initial orientation, mean point, and splitting axis.

  int i; 
  accum M;
  
  //  double F1[3];
  //  double S1[6];
  double C[3][3];
  
  Object_moment = new moment[num_tris]; 

  compute_moments(Object_moment, tris, num_tris);

  clear_accum(M);  
  for(i=0; i<num_tris; i++)
    accum_moment(M, Object_moment[i]);
  
  mean_from_accum(b[0].pT, M);
  covariance_from_accum(C, M);
  
  eigen_and_sort1(b[0].pR, C);

  // create the index list
  int *t = new int[num_tris];

  for(i=0; i<num_tris; i++) t[i] = i;

  // set the tri pointer
  Object_tri = tris;
  
  // do the build
  int rc = b[0].split_recurse(t, num_tris);

  
  // free the moment list
  delete [] Object_moment;  Object_moment = 0;

  // null the tri pointer
  Object_tri = 0;
  
  // free the index list
  delete [] t;

  return 0;
}


inline
void
reaccum_moments(accum &A, int *t, int n)
{
  clear_accum(A);
  for(int i=0; i<n; i++)
      accum_moment(A, Object_moment[t[i]]);
}


int
box::split_recurse(int *t, int n)
{
  // The orientation for the parent box is already assigned to this->pR.
  // The axis along which to split will be column 0 of this->pR.
  // The mean point is passed in on this->pT.

  // When this routine completes, the position and orientation in model
  // space will be established, as well as its dimensions.  Child boxes
  // will be constructed and placed in the parent's CS.

  if (n == 1)
    {
      return split_recurse(t);
    }
  
  // walk along the tris for the box, and do the following:
  //   1. collect the max and min of the vertices along the axes of <or>.
  //   2. decide which group the triangle goes in, performing appropriate swap.
  //   3. accumulate the mean point and covariance data for that triangle.

  accum M1, M2;
  double C[3][3];
  double c[3];
  double minval[3], maxval[3];

  int rc;   // for return code on procedure calls.
  int in;
  tri *ptr;
  int i;
  double axdmp;
  int n1 = 0;  // The number of tris in group 1.  
  // Group 2 will have n - n1 tris.

  // project approximate mean point onto splitting axis, and get coord.
  axdmp = (pR[0][0] * pT[0] + pR[1][0] * pT[1] + pR[2][0] * pT[2]);

  clear_accum(M1);
  clear_accum(M2);

  MTxV(c, pR, Object_tri[t[0]].p1);
  minval[0] = maxval[0] = c[0];
  minval[1] = maxval[1] = c[1];
  minval[2] = maxval[2] = c[2];
  for(i=0; i<n; i++)
    {
      in = t[i];
      ptr = Object_tri + in;
      
      MTxV(c, pR, ptr->p1);
      minmax(minval[0], maxval[0], c[0]);
      minmax(minval[1], maxval[1], c[1]);
      minmax(minval[2], maxval[2], c[2]);

      MTxV(c, pR, ptr->p2);
      minmax(minval[0], maxval[0], c[0]);
      minmax(minval[1], maxval[1], c[1]);
      minmax(minval[2], maxval[2], c[2]);

      MTxV(c, pR, ptr->p3);
      minmax(minval[0], maxval[0], c[0]);
      minmax(minval[1], maxval[1], c[1]);
      minmax(minval[2], maxval[2], c[2]);

      // grab the mean point of the in'th triangle, project
      // it onto the splitting axis (1st column of pR) and
      // see where it lies with respect to axdmp.
      mean_from_moment(c, Object_moment[in]);
      
      if (((pR[0][0]*c[0] + pR[1][0]*c[1] + pR[2][0]*c[2]) < axdmp)
	  && ((n!=2)) || ((n==2) && (i==0)))    
	{
	  // accumulate first and second order moments for group 1
	  accum_moment(M1, Object_moment[in]);

	  // put it in group 1 by swapping t[i] with t[n1]
	  int temp = t[i];
	  t[i] = t[n1];
	  t[n1] = temp;
	  n1++;
	}
      else
	{
	  // accumulate first and second order moments for group 2
	  accum_moment(M2, Object_moment[in]);

	  // leave it in group 2
	  // do nothing...it happens by default
	}
    }

  // done using this->pT as a mean point.


  // error check!
  if ((n1 == 0) || (n1 == n))
    {
      // our partitioning has failed: all the triangles fell into just
      // one of the groups.  So, we arbitrarily partition them into
      // equal parts, and proceed.

      n1 = n/2;
      
      // now recompute accumulated stuff
      reaccum_moments(M1, t, n1);
      reaccum_moments(M2, t + n1, n - n1);
    }

  // With the max and min data, determine the center point and dimensions
  // of the parent box.

  c[0] = (minval[0] + maxval[0])*0.5;
  c[1] = (minval[1] + maxval[1])*0.5;
  c[2] = (minval[2] + maxval[2])*0.5;

  pT[0] = c[0] * pR[0][0] + c[1] * pR[0][1] + c[2] * pR[0][2];
  pT[1] = c[0] * pR[1][0] + c[1] * pR[1][1] + c[2] * pR[1][2];
  pT[2] = c[0] * pR[2][0] + c[1] * pR[2][1] + c[2] * pR[2][2];
  d[0] = (maxval[0] - minval[0])*0.5;
  d[1] = (maxval[1] - minval[1])*0.5;
  d[2] = (maxval[2] - minval[2])*0.5;

  // allocate new boxes
  P = Object_boxes + Object_boxes_inited++;
  N = Object_boxes + Object_boxes_inited++;

  // Compute the orienations for the child boxes (eigenvectors of
  // covariance matrix).  Select the direction of maximum spread to be
  // the split axis for each child.
  
  double tR[3][3];
  
  if (n1 > 1)
    {
      mean_from_accum(P->pT, M1);
      covariance_from_accum(C, M1);

      if (eigen_and_sort1(tR, C) > 30)
	{
	  // unable to find an orientation.  We'll just pick identity.
	  Midentity(tR);
	}

      McM(P->pR, tR);
      if ((rc = P->split_recurse(t, n1)) != 0) return rc;
    }
  else
    {
      if ((rc = P->split_recurse(t)) != 0) return rc;
    }
  McM(C, P->pR);  MTxM(P->pR, pR, C);   // and F1
  VmV(c, P->pT, pT);  MTxV(P->pT, pR, c);
  
  if ((n-n1) > 1)
    {      
      mean_from_accum(N->pT, M2);
      covariance_from_accum (C, M2);

      if (eigen_and_sort1(tR, C) > 30)
	{
	  // unable to find an orientation.  We'll just pick identity.
	  Midentity(tR);
	}
      
      McM(N->pR, tR);
      if ((rc = N->split_recurse(t + n1, n - n1)) != 0) return rc;
    }
  else
    {
      if ((rc = N->split_recurse(t+n1)) != 0) return rc;
    }
  McM(C, N->pR); MTxM(N->pR, pR, C);
  VmV(c, N->pT, pT);  MTxV(N->pT, pR, c);  

  return 0;
}

int
box::split_recurse(int *t)
{
  // For a single triangle, orientation is easily determined.
  // The major axis is parallel to the longest edge.
  // The minor axis is normal to the triangle.
  // The in-between axis is determine by these two.

  // this->pR, this->d, and this->pT are set herein.

  P = N = 0;
  tri *ptr = Object_tri + t[0];

  // Find the major axis: parallel to the longest edge.
  double u12[3], u23[3], u31[3];

  // First compute the squared-lengths of each edge
  VmV(u12, ptr->p1, ptr->p2);  
  double d12 = VdotV(u12,u12);
  VmV(u23, ptr->p2, ptr->p3);  
  double d23 = VdotV(u23,u23);
  VmV(u31, ptr->p3, ptr->p1);  
  double d31 = VdotV(u31,u31);

  // Find the edge of longest squared-length, normalize it to
  // unit length, and put result into a0.
  double a0[3];
  double l;  
  if (d12 > d23)
    {
      if (d12 > d31)
	{
	  l = 1.0 / sqrt(d12); 
	  a0[0] = u12[0] * l; 
	  a0[1] = u12[1] * l;
	  a0[2] = u12[2] * l;
	}
      else 
	{
	  l = 1.0 / sqrt(d31);
	  a0[0] = u31[0] * l;
	  a0[1] = u31[1] * l;
	  a0[2] = u31[2] * l;
	}
    }
  else 
    {
      if (d23 > d31)
	{
	  l = 1.0 / sqrt(d23);
	  a0[0] = u23[0] * l;
	  a0[1] = u23[1] * l;
	  a0[2] = u23[2] * l;
	}
      else
	{
	  l = 1.0 / sqrt(d31);
	  a0[0] = u31[0] * l;
	  a0[1] = u31[1] * l;
	  a0[2] = u31[2] * l;
	}
    }

  // Now compute unit normal to triangle, and put into a2.
  double a2[3];
  VcrossV(a2, u12, u23);
  l = 1.0 / Vlength(a2);  a2[0] *= l;  a2[1] *= l;  a2[2] *= l;

  // a1 is a2 cross a0.
  double a1[3];
  VcrossV(a1, a2, a0);

  // Now make the columns of this->pR the vectors a0, a1, and a2.
  pR[0][0] = a0[0];  pR[0][1] = a1[0];  pR[0][2] = a2[0];
  pR[1][0] = a0[1];  pR[1][1] = a1[1];  pR[1][2] = a2[1];
  pR[2][0] = a0[2];  pR[2][1] = a1[2];  pR[2][2] = a2[2];
  
  // Now compute the maximum and minimum extents of each vertex 
  // along each of the box axes.  From this we will compute the 
  // box center and box dimensions.
  double minval[3], maxval[3];
  double c[3];
  
  MTxV(c, pR, ptr->p1);
  minval[0] = maxval[0] = c[0];
  minval[1] = maxval[1] = c[1];
  minval[2] = maxval[2] = c[2];

  MTxV(c, pR, ptr->p2);
  minmax(minval[0], maxval[0], c[0]);
  minmax(minval[1], maxval[1], c[1]);
  minmax(minval[2], maxval[2], c[2]);
  
  MTxV(c, pR, ptr->p3);
  minmax(minval[0], maxval[0], c[0]);
  minmax(minval[1], maxval[1], c[1]);
  minmax(minval[2], maxval[2], c[2]);
  
  // With the max and min data, determine the center point and dimensions
  // of the box
  c[0] = (minval[0] + maxval[0])*0.5;
  c[1] = (minval[1] + maxval[1])*0.5;
  c[2] = (minval[2] + maxval[2])*0.5;

  pT[0] = c[0] * pR[0][0] + c[1] * pR[0][1] + c[2] * pR[0][2];
  pT[1] = c[0] * pR[1][0] + c[1] * pR[1][1] + c[2] * pR[1][2];
  pT[2] = c[0] * pR[2][0] + c[1] * pR[2][1] + c[2] * pR[2][2];

  d[0] = (maxval[0] - minval[0])*0.5;
  d[1] = (maxval[1] - minval[1])*0.5;
  d[2] = (maxval[2] - minval[2])*0.5;
  
  // Assign the one triangle to this box
  trp = ptr;

  return 0;
}



int
tri_contact(box *b1, box *b2)
{
  // assume just one triangle in each box.

  // the vertices of the tri in b2 is in model1 C.S.  The vertices of
  // the other triangle is in model2 CS.  Use Object_mR, Object_mT, and
  // Object_ms to transform into model2 CS.

  double i1[3];
  double i2[3];
  double i3[3];
  int rc;  // return code
  
  sMxVpV(i1, Object_ms, Object_mR, b1->trp->p1, Object_mT);
  sMxVpV(i2, Object_ms, Object_mR, b1->trp->p2, Object_mT);
  sMxVpV(i3, Object_ms, Object_mR, b1->trp->p3, Object_mT);

  Object_num_tri_tests++;

  int f = tri_contact(i1, i2, i3, b2->trp->p1,b2->trp->p2, b2->trp->p3);

  if (f) 
    {
      // add_collision may be unable to allocate enough memory,
      // so be prepared to pass along an OUT_OF_MEMORY return code.
      if ((rc = add_collision(b1->trp->id, b2->trp->id)) != 0)
	return rc;
    }
  
  return 0;
}



int 
collide_recursive(box *b1, box *b2, double R[3][3], double T[3], double s)
{
  double d[3]; // temp storage for scaled dimensions of box b2.
  int rc;      // return codes
  
  if (1)
    {

      
      if (Object_first_contact && (Object_num_contacts > 0)) return 0;

      // test top level

      Object_num_box_tests++;
  
      int f1;
  
      d[0] = s * b2->d[0];
      d[1] = s * b2->d[1];
      d[2] = s * b2->d[2];
      f1 = obb_disjoint(R, T, b1->d, d);



      if (f1 != 0) 
	{
	  return 0;  // stop processing this test, go to top of loop
	}

      // contact between boxes
      if (b1->leaf() && b2->leaf()) 
	{
	  // it is a leaf pair - compare the polygons therein
          // tri_contact uses the model-to-model transforms stored in
	  // Object_mR, Object_mT, Object_ms.

	  // this will pass along any OUT_OF_MEMORY return codes which
	  // may be generated.
	  return tri_contact(b1, b2);
	}

      double U[3];

      double cR[3][3], cT[3], cs;
      
      // Currently, the transform from model 2 to model 1 space is
      // given by [B T s], where y = [B T s].x = s.B.x + T.

      if (b2->leaf() || (!b1->leaf() && (b1->size() > b2->size())))
	{
	  // here we descend to children of b1.  The transform from
	  // a child of b1 to b1 is stored in [b1->N->pR,b1->N->pT],
	  // but we will denote it [B1 T1 1] for short.  Notice that
	  // children boxes always have same scaling as parent, so s=1
	  // for such nested transforms.

	  // Here, we compute [B1 T1 1]'[B T s] = [B1'B B1'(T-T1) s]
	  // for each child, and store the transform into the collision
	  // test queue.

	  MTxM(cR, b1->N->pR, R); 
	  VmV(U, T, b1->N->pT); MTxV(cT, b1->N->pR, U);
	  cs = s;

	  if ((rc = collide_recursive(b1->N, b2, cR, cT, cs)) != 0)
	    return rc;
	  
	  MTxM(cR, b1->P->pR, R); 
	  VmV(U, T, b1->P->pT); MTxV(cT, b1->P->pR, U);
	  cs = s;

	  if ((rc = collide_recursive(b1->P, b2, cR, cT, cs)) != 0)
	    return rc;
	  
	  return 0;
	}
      else 
	{
	  // here we descend to the children of b2.  See comments for
	  // other 'if' clause for explanation.

	  MxM(cR, R, b2->N->pR); sMxVpV(cT, s, R, b2->N->pT, T);
	  cs = s;
	  
	  if ((rc = collide_recursive(b1, b2->N, cR, cT, cs)) != 0)
	    return rc;
	  
	  MxM(cR, R, b2->P->pR); sMxVpV(cT, s, R, b2->P->pT, T);
	  cs = s;

	  if ((rc = collide_recursive(b1, b2->P, cR, cT, cs)) != 0)
	    return rc;
	  
	  return 0; 
	}
  
    }
  
  return 0;
} 
  
int 
Collide(double R1[3][3], double T1[3], Object *Object_model1,
	double R2[3][3], double T2[3], Object *Object_model2,
	int flag)
{
  return Collide(R1, T1, 1.0, Object_model1, R2, T2, 1.0, Object_model2, flag);
}


int 
Collide(double R1[3][3], double T1[3], double s1, Object *Object_model1,
	double R2[3][3], double T2[3], double s2, Object *Object_model2,
	int flag)
{

  box *b1 = Object_model1->b;
  box *b2 = Object_model2->b;
  
  Object_first_contact = 0; 
  if (flag == FIRST_CONTACT) Object_first_contact = 1;
  
  double tR1[3][3], tR2[3][3], R[3][3];
  double tT1[3], tT2[3], T[3], U[3];
  double s;
  
  // [R1,T1,s1] and [R2,T2,s2] are how the two triangle sets
  // (i.e. models) are positioned in world space.  [tR1,tT1,s1] and
  // [tR2,tT2,s2] are how the top level boxes are positioned in world
  // space
  
  MxM(tR1, R1, b1->pR);                  // tR1 = R1 * b1->pR;
  sMxVpV(tT1, s1, R1, b1->pT, T1);       // tT1 = s1 * R1 * b1->pT + T1;
  MxM(tR2, R2, b2->pR);                  // tR2 = R2 * b2->pR;
  sMxVpV(tT2, s2, R2, b2->pT, T2);       // tT2 = s2 * R2 * b2->pT + T2;
  
  // (R,T,s) is the placement of b2's top level box within
  // the coordinate system of b1's top level box.

  MTxM(R, tR1, tR2);                            // R = tR1.T()*tR2;
  VmV(U, tT2, tT1);  sMTxV(T, 1.0/s1, tR1, U);  // T = tR1.T()*(tT2-tT1)/s1;
  
  s = s2/s1;

  // To transform tri's from model1's CS to model2's CS use this:
  //    x2 = ms . mR . x1 + mT

  {
    MTxM(Object_mR, R2, R1);
    VmV(U, T1, T2);  sMTxV(Object_mT, 1.0/s2, R2, U);
    Object_ms = s1/s2;
  }
  

  // reset the report fields
  Object_num_box_tests = 0;
  Object_num_tri_tests = 0;
  Object_num_contacts = 0;

  // make the call
  return collide_recursive(b1, b2, R, T, s);
}

int
add_collision(int id1, int id2)
{
  if (!Object_contact)
    {
      Object_contact = new collision_pair[10];

      Object_num_cols_alloced = 10;
      Object_num_contacts = 0;
    }
  
  if (Object_num_contacts == Object_num_cols_alloced)
    {
      collision_pair *t = new collision_pair[Object_num_cols_alloced*2];

      Object_num_cols_alloced *= 2;
      
      for(int i=0; i<Object_num_contacts; i++) t[i] = Object_contact[i];
      delete [] Object_contact;
      Object_contact = t;
    }
  
  Object_contact[Object_num_contacts].id1 = id1;
  Object_contact[Object_num_contacts].id2 = id2;
  Object_num_contacts++;

  return 0;
}




inline
double
max(double a, double b, double c)
{
  double t = a;
  if (b > t) t = b;
  if (c > t) t = c;
  return t;
}

inline
double
min(double a, double b, double c)
{
  double t = a;
  if (b < t) t = b;
  if (c < t) t = c;
  return t;
}


int
project6(double *ax, 
	 double *p1, double *p2, double *p3, 
	 double *q1, double *q2, double *q3)
{
  double P1 = VdotV(ax, p1);
  double P2 = VdotV(ax, p2);
  double P3 = VdotV(ax, p3);
  double Q1 = VdotV(ax, q1);
  double Q2 = VdotV(ax, q2);
  double Q3 = VdotV(ax, q3);
  
  double mx1 = max(P1, P2, P3);
  double mn1 = min(P1, P2, P3);
  double mx2 = max(Q1, Q2, Q3);
  double mn2 = min(Q1, Q2, Q3);

  if (mn1 > mx2) return 0;
  if (mn2 > mx1) return 0;
  return 1;
}


// very robust triangle intersection test
// uses no divisions
// works on coplanar triangles

int 
tri_contact (double *P1, double *P2, double *P3,
		    double *Q1, double *Q2, double *Q3) 
{

  /*
     One triangle is (p1,p2,p3).  Other is (q1,q2,q3).
     Edges are (e1,e2,e3) and (f1,f2,f3).
     Normals are n1 and m1
     Outwards are (g1,g2,g3) and (h1,h2,h3).

     We assume that the triangle vertices are in the same coordinate system.

     First thing we do is establish a new c.s. so that p1 is at (0,0,0).

     */

  double p1[3], p2[3], p3[3];
  double q1[3], q2[3], q3[3];
  double e1[3], e2[3], e3[3];
  double f1[3], f2[3], f3[3];
  double g1[3], g2[3], g3[3];
  double h1[3], h2[3], h3[3];
  double n1[3], m1[3];
  double z[3];

  double ef11[3], ef12[3], ef13[3];
  double ef21[3], ef22[3], ef23[3];
  double ef31[3], ef32[3], ef33[3];
  
  z[0] = 0.0;  z[1] = 0.0;  z[2] = 0.0;
  
  p1[0] = P1[0] - P1[0];  p1[1] = P1[1] - P1[1];  p1[2] = P1[2] - P1[2];
  p2[0] = P2[0] - P1[0];  p2[1] = P2[1] - P1[1];  p2[2] = P2[2] - P1[2];
  p3[0] = P3[0] - P1[0];  p3[1] = P3[1] - P1[1];  p3[2] = P3[2] - P1[2];
  
  q1[0] = Q1[0] - P1[0];  q1[1] = Q1[1] - P1[1];  q1[2] = Q1[2] - P1[2];
  q2[0] = Q2[0] - P1[0];  q2[1] = Q2[1] - P1[1];  q2[2] = Q2[2] - P1[2];
  q3[0] = Q3[0] - P1[0];  q3[1] = Q3[1] - P1[1];  q3[2] = Q3[2] - P1[2];
  
  e1[0] = p2[0] - p1[0];  e1[1] = p2[1] - p1[1];  e1[2] = p2[2] - p1[2];
  e2[0] = p3[0] - p2[0];  e2[1] = p3[1] - p2[1];  e2[2] = p3[2] - p2[2];
  e3[0] = p1[0] - p3[0];  e3[1] = p1[1] - p3[1];  e3[2] = p1[2] - p3[2];

  f1[0] = q2[0] - q1[0];  f1[1] = q2[1] - q1[1];  f1[2] = q2[2] - q1[2];
  f2[0] = q3[0] - q2[0];  f2[1] = q3[1] - q2[1];  f2[2] = q3[2] - q2[2];
  f3[0] = q1[0] - q3[0];  f3[1] = q1[1] - q3[1];  f3[2] = q1[2] - q3[2];
  
  VcrossV(n1, e1, e2);
  VcrossV(m1, f1, f2);

  VcrossV(g1, e1, n1);
  VcrossV(g2, e2, n1);
  VcrossV(g3, e3, n1);
  VcrossV(h1, f1, m1);
  VcrossV(h2, f2, m1);
  VcrossV(h3, f3, m1);

  VcrossV(ef11, e1, f1);
  VcrossV(ef12, e1, f2);
  VcrossV(ef13, e1, f3);
  VcrossV(ef21, e2, f1);
  VcrossV(ef22, e2, f2);
  VcrossV(ef23, e2, f3);
  VcrossV(ef31, e3, f1);
  VcrossV(ef32, e3, f2);
  VcrossV(ef33, e3, f3);
  
  // now begin the series of tests

  if (!project6(n1, p1, p2, p3, q1, q2, q3)) return 0;
  if (!project6(m1, p1, p2, p3, q1, q2, q3)) return 0;
  
  if (!project6(ef11, p1, p2, p3, q1, q2, q3)) return 0;
  if (!project6(ef12, p1, p2, p3, q1, q2, q3)) return 0;
  if (!project6(ef13, p1, p2, p3, q1, q2, q3)) return 0;
  if (!project6(ef21, p1, p2, p3, q1, q2, q3)) return 0;
  if (!project6(ef22, p1, p2, p3, q1, q2, q3)) return 0;
  if (!project6(ef23, p1, p2, p3, q1, q2, q3)) return 0;
  if (!project6(ef31, p1, p2, p3, q1, q2, q3)) return 0;
  if (!project6(ef32, p1, p2, p3, q1, q2, q3)) return 0;
  if (!project6(ef33, p1, p2, p3, q1, q2, q3)) return 0;

  if (!project6(g1, p1, p2, p3, q1, q2, q3)) return 0;
  if (!project6(g2, p1, p2, p3, q1, q2, q3)) return 0;
  if (!project6(g3, p1, p2, p3, q1, q2, q3)) return 0;
  if (!project6(h1, p1, p2, p3, q1, q2, q3)) return 0;
  if (!project6(h2, p1, p2, p3, q1, q2, q3)) return 0;
  if (!project6(h3, p1, p2, p3, q1, q2, q3)) return 0;

  return 1;
}



/*

int
obb_disjoint(double B[3][3], double T[3], double a[3], double b[3]);

This is a test between two boxes, box A and box B.  It is assumed that
the coordinate system is aligned and centered on box A.  The 3x3
matrix B specifies box B's orientation with respect to box A.
Specifically, the columns of B are the basis vectors (axis vectors) of
box B.  The center of box B is located at the vector T.  The
dimensions of box B are given in the array b.  The orientation and
placement of box A, in this coordinate system, are the identity matrix
and zero vector, respectively, so they need not be specified.  The
dimensions of box A are given in array a.

This test operates in two modes, depending on how the library is
compiled.  It indicates whether the two boxes are overlapping, by
returning a boolean.  

The second version of the routine will return a conservative bounds on
the distance between the polygon sets which the boxes enclose.  It is
used when Object is being used to estimate the distance between two
models.

*/


int
obb_disjoint(double B[3][3], double T[3], double a[3], double b[3])
{
  register double t, s;
  register int r;
  double Bf[3][3];
  const double reps = 1e-6;
  
  // Bf = fabs(B)
  Bf[0][0] = myfabs(B[0][0]);  Bf[0][0] += reps;
  Bf[0][1] = myfabs(B[0][1]);  Bf[0][1] += reps;
  Bf[0][2] = myfabs(B[0][2]);  Bf[0][2] += reps;
  Bf[1][0] = myfabs(B[1][0]);  Bf[1][0] += reps;
  Bf[1][1] = myfabs(B[1][1]);  Bf[1][1] += reps;
  Bf[1][2] = myfabs(B[1][2]);  Bf[1][2] += reps;
  Bf[2][0] = myfabs(B[2][0]);  Bf[2][0] += reps;
  Bf[2][1] = myfabs(B[2][1]);  Bf[2][1] += reps;
  Bf[2][2] = myfabs(B[2][2]);  Bf[2][2] += reps;

  
#if TRACE1
  printf("Box test: Bf[3][3], B[3][3], T[3], a[3], b[3]\n");
  Mprintg(Bf);
  Mprintg(B);
  Vprintg(T);
  Vprintg(a);
  Vprintg(b);
#endif
  
  // if any of these tests are one-sided, then the polyhedra are disjoint
  r = 1;

  // A1 x A2 = A0
  t = myfabs(T[0]);
  
  r &= (t <= 
	  (a[0] + b[0] * Bf[0][0] + b[1] * Bf[0][1] + b[2] * Bf[0][2]));
  if (!r) return 1;
  
  // B1 x B2 = B0
  s = T[0]*B[0][0] + T[1]*B[1][0] + T[2]*B[2][0];
  t = myfabs(s);

  r &= ( t <=
	  (b[0] + a[0] * Bf[0][0] + a[1] * Bf[1][0] + a[2] * Bf[2][0]));
  if (!r) return 2;
    
  // A2 x A0 = A1
  t = myfabs(T[1]);
  
  r &= ( t <= 
	  (a[1] + b[0] * Bf[1][0] + b[1] * Bf[1][1] + b[2] * Bf[1][2]));
  if (!r) return 3;

  // A0 x A1 = A2
  t = myfabs(T[2]);

  r &= ( t <= 
	  (a[2] + b[0] * Bf[2][0] + b[1] * Bf[2][1] + b[2] * Bf[2][2]));
  if (!r) return 4;

  // B2 x B0 = B1
  s = T[0]*B[0][1] + T[1]*B[1][1] + T[2]*B[2][1];
  t = myfabs(s);

  r &= ( t <=
	  (b[1] + a[0] * Bf[0][1] + a[1] * Bf[1][1] + a[2] * Bf[2][1]));
  if (!r) return 5;

  // B0 x B1 = B2
  s = T[0]*B[0][2] + T[1]*B[1][2] + T[2]*B[2][2];
  t = myfabs(s);

  r &= ( t <=
	  (b[2] + a[0] * Bf[0][2] + a[1] * Bf[1][2] + a[2] * Bf[2][2]));
  if (!r) return 6;

  // A0 x B0
  s = T[2] * B[1][0] - T[1] * B[2][0];
  t = myfabs(s);
  
  r &= ( t <= 
	(a[1] * Bf[2][0] + a[2] * Bf[1][0] +
	 b[1] * Bf[0][2] + b[2] * Bf[0][1]));
  if (!r) return 7;
  
  // A0 x B1
  s = T[2] * B[1][1] - T[1] * B[2][1];
  t = myfabs(s);

  r &= ( t <=
	(a[1] * Bf[2][1] + a[2] * Bf[1][1] +
	 b[0] * Bf[0][2] + b[2] * Bf[0][0]));
  if (!r) return 8;

  // A0 x B2
  s = T[2] * B[1][2] - T[1] * B[2][2];
  t = myfabs(s);

  r &= ( t <=
	  (a[1] * Bf[2][2] + a[2] * Bf[1][2] +
	   b[0] * Bf[0][1] + b[1] * Bf[0][0]));
  if (!r) return 9;

  // A1 x B0
  s = T[0] * B[2][0] - T[2] * B[0][0];
  t = myfabs(s);

  r &= ( t <=
	  (a[0] * Bf[2][0] + a[2] * Bf[0][0] +
	   b[1] * Bf[1][2] + b[2] * Bf[1][1]));
  if (!r) return 10;

  // A1 x B1
  s = T[0] * B[2][1] - T[2] * B[0][1];
  t = myfabs(s);

  r &= ( t <=
	  (a[0] * Bf[2][1] + a[2] * Bf[0][1] +
	   b[0] * Bf[1][2] + b[2] * Bf[1][0]));
  if (!r) return 11;

  // A1 x B2
  s = T[0] * B[2][2] - T[2] * B[0][2];
  t = myfabs(s);

  r &= (t <=
	  (a[0] * Bf[2][2] + a[2] * Bf[0][2] +
	   b[0] * Bf[1][1] + b[1] * Bf[1][0]));
  if (!r) return 12;

  // A2 x B0
  s = T[1] * B[0][0] - T[0] * B[1][0];
  t = myfabs(s);

  r &= (t <=
	  (a[0] * Bf[1][0] + a[1] * Bf[0][0] +
	   b[1] * Bf[2][2] + b[2] * Bf[2][1]));
  if (!r) return 13;

  // A2 x B1
  s = T[1] * B[0][1] - T[0] * B[1][1];
  t = myfabs(s);

  r &= ( t <=
	  (a[0] * Bf[1][1] + a[1] * Bf[0][1] +
	   b[0] * Bf[2][2] + b[2] * Bf[2][0]));
  if (!r) return 14;

  // A2 x B2
  s = T[1] * B[0][2] - T[0] * B[1][2];
  t = myfabs(s);

  r &= ( t <=
	  (a[0] * Bf[1][2] + a[1] * Bf[0][2] +
	   b[0] * Bf[2][1] + b[1] * Bf[2][0]));
  if (!r) return 15;

  return 0;  // should equal 0
}
