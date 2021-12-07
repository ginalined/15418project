#include "CycleTimer.h"
#include "VInternal.H"
#include <cuda.h>
#include <cuda_runtime.h>
#include <driver_functions.h>
#include <iostream>
#include <string.h> //for memset and memcpy.

#include <math.h>

#include "matvecc.cu_inl"

static moment *Object_moment = 0;
static tri *Object_tri = 0;
static box *Object_boxes = 0;
static int Object_boxes_inited = 0;
static tri *cuda_tris = 0;
double Object_mR[3][3];
double Object_mT[3];
double Object_ms;

int Object_first_contact;
int Object_num_box_tests = 0;
int Object_num_tri_tests;
int Object_num_contacts = 0;
int Object_num_cols_alloced = 0;
collision_pair *Object_contact = 0;

__device__ double cuda_mR[3][3];
__device__  double cuda_mT[3];
__device__  double cuda_ms;
__device__ int cuda_first_contact;
__device__ int cuda_num_box_tests = 0;
__device__ int cuda_num_tri_tests;
__device__ int cuda_num_contacts = 0;
__device__ int cuda_num_cols_alloced = 0;
// collision_pair *Object_contact = 0;

int BLOCK_SIZE = 128;
int add_collision(int id1, int id2);

Object::Object() {
  b = 0;
  num_boxes_alloced = 0;
  tris = 0;
  num_tris = 0;
  num_tris_alloced = 0;
}

Object::~Object() {}

int Object::BeginModel() {
  delete[] b;
  b = 0;
  num_boxes_alloced = 0;
  delete[] tris;
  tris = 0;
  num_tris = 0;
  num_tris_alloced = 0;
  return 0;
}

__global__ void test_tris(tri *t) { printf("the p %f", t[0].p1); }

int Object::EndModel() {
  cudaMalloc(&cuda_tris, sizeof(tri) * num_tris);
  cudaMemcpy(cuda_tris, tris, sizeof(tri) * num_tris, cudaMemcpyHostToDevice);
  // test_tris<<<1,1>>>(cuda_tris);

  int myrc = build_hierarchy();

  return 0;
}

int Object::AddTri(const double *p1, const double *p2, const double *p3,
                   int id) {

  // first make sure that we haven't filled up our allocation.
  // if we have, allocate a new array of twice the size, and copy
  // the old data to it.

  if (num_tris == num_tris_alloced) {
    // decide on new size -- accounting for first time, where none are
    // allocated
    int n = num_tris_alloced * 2;
    if (n == 0)
      n = 1;

    // make new array, and copy the old one to it
    tri *t = new tri[n];

    int i;
    for (i = 0; i < num_tris; i++) {
      t[i] = tris[i];
    }

    // free the old array and reassign.
    delete[] tris;
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

int Object::build_hierarchy() {
  // allocate the boxes and set the box list globals

  num_boxes_alloced = num_tris * 2;
  b = new box[num_boxes_alloced];
  Object_boxes = b;
  Object_boxes_inited = 1; // we are in process of initializing b[0].

  int i;
  accum M;

  double C[3][3];
  // num_tris = 16;
  Object_moment = new moment[num_tris];

  // every tris has its moment
  compute_moments(Object_moment, tris, num_tris);

  clear_accum(M);
  // added all moments up
  for (i = 0; i < num_tris; i++)
    accum_moment(M, Object_moment[i]);

  // calculate mean and covariance
  mean_from_accum(b[0].pT, M);
  covariance_from_accum(C, M);
  eigen_and_sort1(b[0].pR, C);

  // create the index list
  int *t = new int[num_tris];

  for (i = 0; i < num_tris; i++)
    t[i] = i;

  // set the tri pointer
  Object_tri = tris;

  // do the build
  int rc = b[0].split_recurse(t, num_tris);

  // free the moment list
  delete[] Object_moment;
  Object_moment = 0;

  // null the tri pointer
  Object_tri = 0;

  // free the index list
  delete[] t;

  return 0;
}



inline void reaccum_moments(accum &A, int *t, int n) {
  clear_accum(A);
  for (int i = 0; i < n; i++)
    accum_moment(A, Object_moment[t[i]]);
}

__global__ void split_cuda(double *all_box, int *t, tri *cuda_tris, int N,
                           double *pR, double *minval) {

  int j = blockIdx.x * blockDim.x + threadIdx.x;
  if (j >= N)
    return;
  int in = t[j];
  int maxThread =
      N < (blockIdx.x * (blockDim.x + 1)) ? N : (blockIdx.x * (blockDim.x + 1));
  __shared__ double temp[8 * 3];
  // printf("ptr %f", ptr->p1[0]);
  cuda_MTxV(&all_box[j * 9], pR, cuda_tris[t[j]].p1);
  cuda_MTxV(&all_box[j * 9 + 3], pR, cuda_tris[t[j]].p2);
  cuda_MTxV(&all_box[j * 9 + 6], pR, cuda_tris[t[j]].p3);
  // temp[threadIdx.x*3]=temp[threadIdx.x*3+1]=temp[threadIdx.x*3+2]=1<<29;

  for (int k = 3; k < 9; k += 1) {
    if (all_box[j * 9 + k % 3] > all_box[j * 9 + k]) {
      // temp[threadId.x*3 + k%3] = all_box[j*9+k];
      double temp = all_box[j * 9 + k % 3];
      all_box[j * 9 + k % 3] = all_box[j * 9 + k];
      all_box[j * 9 + k] = temp;
    }
  }
  for (int k = 3; k < 6; k += 1) {
    if (all_box[j * 9 + 3 + k] < all_box[j * 9 + k]) {
      all_box[j * 9 + 3 + k] = all_box[j * 9 + k];
    }
  }

  // printf("j is %d, %d, %d\n", j, blockIdx.x , blockDim.x);
  // __syncthreads();
  // if (threadIdx.x == 0){
  //   all_box[j*9+6]
  // for (int k = 0; k < 24; k++){
  //   if (temp[k%3] > temp[k]){
  //     temp[k%3] = temp[k];
  //   }
  // }
}

int box::split_recurse(int *t, int n) {
  // The orientation for the parent box is already assigned to this->pR.
  // The axis along which to split will be column 0 of this->pR.
  // The mean point is passed in on this->pT.

  // When this routine completes, the position and orientation in model
  // space will be established, as well as its dimensions.  Child boxes
  // will be constructed and placed in the parent's CS.

  if (n == 1) {
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

  int rc; // for return code on procedure calls.
  int in;
  tri *ptr;
  int i;
  double axdmp;
  int n1 = 0; // The number of tris in group 1.
  // Group 2 will have n - n1 tris.

  // project approximate mean point onto splitting axis, and get coord.
  axdmp = (pR[0][0] * pT[0] + pR[1][0] * pT[1] + pR[2][0] * pT[2]);

  clear_accum(M1);
  clear_accum(M2);

  MTxV(c, pR, Object_tri[t[0]].p1);
  minval[0] = maxval[0] = c[0];
  minval[1] = maxval[1] = c[1];
  minval[2] = maxval[2] = c[2];
  double *all_box = new double[9 * n];
  double *cuda_box;
  // tri * cuda_tris;
  int *cuda_key;
  // double * cuda_moment;
  double *cuda_pr;
  double *cuda_min;
  cudaMalloc(&cuda_box, sizeof(double) * 9 * n);
  // cudaMalloc(&cuda_moment, sizeof(double) * 13 * n);
  cudaMalloc(&cuda_key, sizeof(int) * n);
  cudaMalloc(&cuda_min, sizeof(double) * 3);
  cudaMalloc(&cuda_pr, sizeof(double) * 9);
  // cudaMalloc(&cuda_tris, sizeof(tri) * n);
  //  cudaMemcpy(cuda_tris, Object_tri, n * sizeof(tri),
  //  cudaMemcpyHostToDevice);
  cudaMemcpy(cuda_key, t, sizeof(int) * n, cudaMemcpyHostToDevice);
  double *check_pr = new double[9];

  for (int i = 0; i < 9; i++) {
    check_pr[i] = pR[i / 3][i % 3];
  }
  for (int i = 0; i < 3; i++) {
    minval[i] = minval[i];
  }

  cudaMemcpy(cuda_pr, check_pr, sizeof(double) * 9, cudaMemcpyHostToDevice);
  cudaMemcpy(cuda_min, minval, sizeof(double) * 3, cudaMemcpyHostToDevice);

  // cudaMalloc(&cuda_moment, sizeof(moment) * n);
  // cudaMalloc(&output, sizeof(int) * n);

  double checkMin[3];
  split_cuda<<<8, 8>>>(cuda_box, cuda_key, cuda_tris, n, cuda_pr, cuda_min);
  cudaMemcpy(all_box, cuda_box, sizeof(double) * n * 9, cudaMemcpyDeviceToHost);
  cudaMemcpy(checkMin, cuda_min, sizeof(double) * 3, cudaMemcpyDeviceToHost);

  //  for (int k = 0;k < 3;k++){
  //   printf("check difference %f, and %f \n", minval[0], checkMin[0]);
  // }

  for (int j = 0; j < 9 * n; j += 9) {
    for (int k = 0; k < 3; k++) {
      minval[k] = std::min(minval[k], all_box[j + k]);
      maxval[k] = std::max(maxval[k], all_box[j + k + 6]);
    }
  }

  for (i = 0; i < n; i++) {
    in = t[i];
    ptr = Object_tri + in;
    mean_from_moment(c, Object_moment[in]);

    if (((pR[0][0] * c[0] + pR[1][0] * c[1] + pR[2][0] * c[2]) < axdmp) &&
            ((n != 2)) ||
        ((n == 2) && (i == 0))) {
      // accumulate first and second order moments for group 1
      accum_moment(M1, Object_moment[in]);

      // put it in group 1 by swapping t[i] with t[n1]
      int temp = t[i];
      t[i] = t[n1];
      t[n1] = temp;
      n1++;
    } else {
      // accumulate first and second order moments for group 2
      accum_moment(M2, Object_moment[in]);

      // leave it in group 2
      // do nothing...it happens by default
    }
  }

  // error check!
  if ((n1 == 0) || (n1 == n)) {
    // our partitioning has failed: all the triangles fell into just
    // one of the groups.  So, we arbitrarily partition them into
    // equal parts, and proceed.

    n1 = n / 2;

    // now recompute accumulated stuff
    reaccum_moments(M1, t, n1);
    reaccum_moments(M2, t + n1, n - n1);
  }

  // With the max and min data, determine the center point and dimensions
  // of the parent box.

  c[0] = (minval[0] + maxval[0]) * 0.5;
  c[1] = (minval[1] + maxval[1]) * 0.5;
  c[2] = (minval[2] + maxval[2]) * 0.5;

  pT[0] = c[0] * pR[0][0] + c[1] * pR[0][1] + c[2] * pR[0][2];
  pT[1] = c[0] * pR[1][0] + c[1] * pR[1][1] + c[2] * pR[1][2];
  pT[2] = c[0] * pR[2][0] + c[1] * pR[2][1] + c[2] * pR[2][2];
  d[0] = (maxval[0] - minval[0]) * 0.5;
  d[1] = (maxval[1] - minval[1]) * 0.5;
  d[2] = (maxval[2] - minval[2]) * 0.5;

  // allocate new boxes
  prev_index = Object_boxes_inited++;
  next_index = Object_boxes_inited++;
  P = &Object_boxes[prev_index];
  N = &Object_boxes[next_index];

  // Compute the orienations for the child boxes (eigenvectors of
  // covariance matrix).  Select the direction of maximum spread to be
  // the split axis for each child.

  double tR[3][3];

  if (n1 > 1) {
    mean_from_accum(P->pT, M1);
    covariance_from_accum(C, M1);

    if (eigen_and_sort1(tR, C) > 30) {
      // unable to find an orientation.  We'll just pick identity.
      Midentity(tR);
    }

    McM(P->pR, tR);
    if ((rc = P->split_recurse(t, n1)) != 0)
      return rc;
  } else {
    if ((rc = P->split_recurse(t)) != 0)
      return rc;
  }
  McM(C, P->pR);
  MTxM(P->pR, pR, C); // and F1
  VmV(c, P->pT, pT);
  MTxV(P->pT, pR, c);

  if ((n - n1) > 1) {
    mean_from_accum(N->pT, M2);
    covariance_from_accum(C, M2);

    if (eigen_and_sort1(tR, C) > 30) {
      // unable to find an orientation.  We'll just pick identity.
      Midentity(tR);
    }

    McM(N->pR, tR);
    if ((rc = N->split_recurse(t + n1, n - n1)) != 0)
      return rc;
  } else {
    if ((rc = N->split_recurse(t + n1)) != 0)
      return rc;
  }
  McM(C, N->pR);
  MTxM(N->pR, pR, C);
  VmV(c, N->pT, pT);
  MTxV(N->pT, pR, c);

  return 0;
}

int box::split_recurse(int *t) {

  P = N = 0;
  tri *ptr = Object_tri + t[0];

  // Find the major axis: parallel to the longest edge.
  double u12[3], u23[3], u31[3];

  // First compute the squared-lengths of each edge
  VmV(u12, ptr->p1, ptr->p2);
  double d12 = VdotV(u12, u12);
  VmV(u23, ptr->p2, ptr->p3);
  double d23 = VdotV(u23, u23);
  VmV(u31, ptr->p3, ptr->p1);
  double d31 = VdotV(u31, u31);

  // Find the edge of longest squared-length, normalize it to
  // unit length, and put result into a0.
  double a0[3];
  double l;
  if (d12 > d23) {
    if (d12 > d31) {
      l = 1.0 / sqrt(d12);
      a0[0] = u12[0] * l;
      a0[1] = u12[1] * l;
      a0[2] = u12[2] * l;
    } else {
      l = 1.0 / sqrt(d31);
      a0[0] = u31[0] * l;
      a0[1] = u31[1] * l;
      a0[2] = u31[2] * l;
    }
  } else {
    if (d23 > d31) {
      l = 1.0 / sqrt(d23);
      a0[0] = u23[0] * l;
      a0[1] = u23[1] * l;
      a0[2] = u23[2] * l;
    } else {
      l = 1.0 / sqrt(d31);
      a0[0] = u31[0] * l;
      a0[1] = u31[1] * l;
      a0[2] = u31[2] * l;
    }
  }

  // Now compute unit normal to triangle, and put into a2.
  double a2[3];
  VcrossV(a2, u12, u23);
  l = 1.0 / Vlength(a2);
  a2[0] *= l;
  a2[1] *= l;
  a2[2] *= l;

  // a1 is a2 cross a0.
  double a1[3];
  VcrossV(a1, a2, a0);

  // Now make the columns of this->pR the vectors a0, a1, and a2.
  pR[0][0] = a0[0];
  pR[0][1] = a1[0];
  pR[0][2] = a2[0];
  pR[1][0] = a0[1];
  pR[1][1] = a1[1];
  pR[1][2] = a2[1];
  pR[2][0] = a0[2];
  pR[2][1] = a1[2];
  pR[2][2] = a2[2];

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
  c[0] = (minval[0] + maxval[0]) * 0.5;
  c[1] = (minval[1] + maxval[1]) * 0.5;
  c[2] = (minval[2] + maxval[2]) * 0.5;

  pT[0] = c[0] * pR[0][0] + c[1] * pR[0][1] + c[2] * pR[0][2];
  pT[1] = c[0] * pR[1][0] + c[1] * pR[1][1] + c[2] * pR[1][2];
  pT[2] = c[0] * pR[2][0] + c[1] * pR[2][1] + c[2] * pR[2][2];

  d[0] = (maxval[0] - minval[0]) * 0.5;
  d[1] = (maxval[1] - minval[1]) * 0.5;
  d[2] = (maxval[2] - minval[2]) * 0.5;

  // Assign the one triangle to this box
  trp = ptr;

  return 0;
}

int tri_contact(box *b1, box *b2) {
  // assume just one triangle in each box.

  // the vertices of the tri in b2 is in model1 C.S.  The vertices of
  // the other triangle is in model2 CS.  Use Object_mR, Object_mT, and
  // Object_ms to transform into model2 CS.

  double i1[3];
  double i2[3];
  double i3[3];
  int rc; // return code

  sMxVpV(i1, Object_ms, Object_mR, b1->trp->p1, Object_mT);
  sMxVpV(i2, Object_ms, Object_mR, b1->trp->p2, Object_mT);
  sMxVpV(i3, Object_ms, Object_mR, b1->trp->p3, Object_mT);

  Object_num_tri_tests++;

  int f = tri_contact(i1, i2, i3, b2->trp->p1, b2->trp->p2, b2->trp->p3);

  if (f) {
    // add_collision may be unable to allocate enough memory,
    // so be prepared to pass along an OUT_OF_MEMORY return code.
    if ((rc = add_collision(b1->trp->id, b2->trp->id)) != 0)
      return rc;
  }

  return 0;
}

__device__ int cuda_tri_contact(box *b1, box *b2, int *collisions ) {
  // assume just one triangle in each box.

  // the vertices of the tri in b2 is in model1 C.S.  The vertices of
  // the other triangle is in model2 CS.  Use Object_mR, Object_mT, and
  // Object_ms to transform into model2 CS.

  double i1[3];
  double i2[3];
  double i3[3];
  int rc; // return code

  sMxVpV(i1, cuda_ms, cuda_mR, b1->trp->p1, cuda_mT);
  sMxVpV(i2, cuda_ms, cuda_mR, b1->trp->p2, cuda_mT);
  sMxVpV(i3, cuda_ms, cuda_mR, b1->trp->p3, cuda_mT);

  cuda_num_tri_tests++;

  int f = tri_contact(i1, i2, i3, b2->trp->p1, b2->trp->p2, b2->trp->p3);


  if (f) {
    printf("hello, I see a collision! \n");
  }

  return 0;
}




int collide_recursive(box *b1, box *b2, double R[3][3], double T[3], double s) {
  double d[3]; // temp storage for scaled dimensions of box b2.
  int rc;      // return codes

  if (Object_first_contact && (Object_num_contacts > 0))
    return 0;

  // test top level

  Object_num_box_tests++;

  int f1;

  d[0] = s * b2->d[0];
  d[1] = s * b2->d[1];
  d[2] = s * b2->d[2];
  f1 = obb_disjoint(R, T, b1->d, d);

  if (f1 != 0) {
    return 0; // stop processing this test, go to top of loop
  }

  // contact between boxes
  if (b1->leaf() && b2->leaf()) {
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

  if (b2->leaf() || (!b1->leaf() && (b1->size() > b2->size()))) {
    // here we descend to children of b1.  The transform from
    // a child of b1 to b1 is stored in [b1->N->pR,b1->N->pT],
    // but we will denote it [B1 T1 1] for short.  Notice that
    // children boxes always have same scaling as parent, so s=1
    // for such nested transforms.

    // Here, we compute [B1 T1 1]'[B T s] = [B1'B B1'(T-T1) s]
    // for each child, and store the transform into the collision
    // test queue.

    MTxM(cR, b1->N->pR, R);
    VmV(U, T, b1->N->pT);
    MTxV(cT, b1->N->pR, U);
    cs = s;

    if ((rc = collide_recursive(b1->N, b2, cR, cT, cs)) != 0)
      return rc;

    MTxM(cR, b1->P->pR, R);
    VmV(U, T, b1->P->pT);
    MTxV(cT, b1->P->pR, U);
    cs = s;

    if ((rc = collide_recursive(b1->P, b2, cR, cT, cs)) != 0)
      return rc;

    return 0;
  } else {
    // here we descend to the children of b2.  See comments for
    // other 'if' clause for explanation.

    MxM(cR, R, b2->N->pR);
    sMxVpV(cT, s, R, b2->N->pT, T);
    cs = s;

    if ((rc = collide_recursive(b1, b2->N, cR, cT, cs)) != 0)
      return rc;

    MxM(cR, R, b2->P->pR);
    sMxVpV(cT, s, R, b2->P->pT, T);
    cs = s;

    if ((rc = collide_recursive(b1, b2->P, cR, cT, cs)) != 0)
      return rc;

    return 0;
  }

  return 0;
}

int Collide(double R1[3][3], double T1[3], Object *Object_model1,
            double R2[3][3], double T2[3], Object *Object_model2, int flag) {
  return Collide(R1, T1, 1.0, Object_model1, R2, T2, 1.0, Object_model2, flag);
}

__device__ int cuda_collide_recursive(box *b1, box *b2, double R[3][3],
                                      double T[3], double s) {
  double d[3]; // temp storage for scaled dimensions of box b2.
  int rc;

  if (cuda_first_contact && (cuda_num_contacts > 0))
    return 0;

  cuda_num_box_tests++;

  int f1;

  d[0] = s * b2->d[0];
  d[1] = s * b2->d[1];
  d[2] = s * b2->d[2];
  f1 = obb_disjoint(R, T, b1->d, d);

  if (f1 != 0) {
    return 0; // stop processing this test, go to top of loop
  }

  // contact between boxes
  if (b1->leaf() && b2->leaf()) {
    return 0;
    //return cuda_tri_contact(b1, b2);
  }

  double U[3];

  double cR[3][3], cT[3], cs;

  // Currently, the transform from model 2 to model 1 space is
  // given by [B T s], where y = [B T s].x = s.B.x + T.

  if (b2->leaf() || (!b1->leaf() && (b1->size() > b2->size()))) {

    MTxM(cR, b1->N->pR, R);
    VmV(U, T, b1->N->pT);
    MTxV(cT, b1->N->pR, U);
    cs = s;

    if ((rc = cuda_collide_recursive(b1->N, b2, cR, cT, cs)) != 0)
      return rc;

    MTxM(cR, b1->P->pR, R);
    VmV(U, T, b1->P->pT);
    MTxV(cT, b1->P->pR, U);
    cs = s;

    if ((rc = cuda_collide_recursive(b1->P, b2, cR, cT, cs)) != 0)
      return rc;

    return 0;
  } else {

    MxM(cR, R, b2->N->pR);
    sMxVpV(cT, s, R, b2->N->pT, T);
    cs = s;

    if ((rc = cuda_collide_recursive(b1, b2->N, cR, cT, cs)) != 0)
      return rc;

    MxM(cR, R, b2->P->pR);
    sMxVpV(cT, s, R, b2->P->pT, T);
    cs = s;

    if ((rc = cuda_collide_recursive(b1, b2->P, cR, cT, cs)) != 0)
      return rc;
    return 0;
  }

  return 0;
}

__device__ int test_Collide(double R1[3][3], double T1[3], double s1,
                            box * b1, double R2[3][3],
                            double T2[3], double s2, box *b2) {


  cuda_first_contact = 1;

  double tR1[3][3], tR2[3][3], R[3][3];
  double tT1[3], tT2[3], T[3], U[3];
  double s;

  MxM(tR1, R1, b1->pR);            // tR1 = R1 * b1->pR;
  sMxVpV(tT1, s1, R1, b1->pT, T1); // tT1 = s1 * R1 * b1->pT + T1;
  MxM(tR2, R2, b2->pR);            // tR2 = R2 * b2->pR;
  sMxVpV(tT2, s2, R2, b2->pT, T2); // tT2 = s2 * R2 * b2->pT + T2;
  MTxM(R, tR1, tR2);               // R = tR1.T()*tR2;
  VmV(U, tT2, tT1);
  sMTxV(T, 1.0 / s1, tR1, U); // T = tR1.T()*(tT2-tT1)/s1;

  s = s2 / s1;

  MTxM(cuda_mR, R2, R1);
  VmV(U, T1, T2);
  sMTxV(cuda_mT, 1.0 / s2, R2, U);
  cuda_ms = s1 / s2;

  // reset the report fields
  cuda_num_box_tests = 0;
  cuda_num_tri_tests = 0;
  cuda_num_contacts = 0;

  // make the call
  return cuda_collide_recursive(b1, b2, R, T, s);
}

int Collide(double R1[3][3], double T1[3], double s1, Object *Object_model1,
            double R2[3][3], double T2[3], double s2, Object *Object_model2,
            int flag) {

  box *b1 = Object_model1->b;
  box *b2 = Object_model2->b;

  Object_first_contact = 1;

  double tR1[3][3], tR2[3][3], R[3][3];
  double tT1[3], tT2[3], T[3], U[3];
  double s;

  // [R1,T1,s1] and [R2,T2,s2] are how the two triangle sets
  // (i.e. models) are positioned in world space.  [tR1,tT1,s1] and
  // [tR2,tT2,s2] are how the top level boxes are positioned in world
  // space

  MxM(tR1, R1, b1->pR);            // tR1 = R1 * b1->pR;
  sMxVpV(tT1, s1, R1, b1->pT, T1); // tT1 = s1 * R1 * b1->pT + T1;
  MxM(tR2, R2, b2->pR);            // tR2 = R2 * b2->pR;
  sMxVpV(tT2, s2, R2, b2->pT, T2); // tT2 = s2 * R2 * b2->pT + T2;

  // (R,T,s) is the placement of b2's top level box within
  // the coordinate system of b1's top level box.

  MTxM(R, tR1, tR2); // R = tR1.T()*tR2;
  VmV(U, tT2, tT1);
  sMTxV(T, 1.0 / s1, tR1, U); // T = tR1.T()*(tT2-tT1)/s1;

  s = s2 / s1;

  // To transform tri's from model1's CS to model2's CS use this:
  //    x2 = ms . mR . x1 + mT

  {
    MTxM(Object_mR, R2, R1);
    VmV(U, T1, T2);
    sMTxV(Object_mT, 1.0 / s2, R2, U);
    Object_ms = s1 / s2;
  }

  // reset the report fields
  Object_num_box_tests = 0;
  Object_num_tri_tests = 0;
  Object_num_contacts = 0;

  // make the call
  return collide_recursive(b1, b2, R, T, s);
}

int add_collision(int id1, int id2) {
  if (!Object_contact) {
    Object_contact = new collision_pair[10];

    Object_num_cols_alloced = 10;
    Object_num_contacts = 0;
  }

  if (Object_num_contacts == Object_num_cols_alloced) {
    collision_pair *t = new collision_pair[Object_num_cols_alloced * 2];

    Object_num_cols_alloced *= 2;

    for (int i = 0; i < Object_num_contacts; i++)
      t[i] = Object_contact[i];
    delete[] Object_contact;
    Object_contact = t;
  }

  Object_contact[Object_num_contacts].id1 = id1;
  Object_contact[Object_num_contacts].id2 = id2;
  Object_num_contacts++;
  ;

  return 0;
}

int get_info(AABB *input) {

  AABB *temp = new AABB[32];
  cudaMemcpy(temp, input, 32 * sizeof(AABB), cudaMemcpyDeviceToHost);
  printf("no way!\n");
  for (int i = 0; i < 32; i++)
    printf("the idea is that %f, %d\n", temp[i].center[0], temp[i].id);
  return 0;
}
int get_info1(int *input) {

  int *temp = new int[32];
  cudaMemcpy(temp, input, 32 * sizeof(int), cudaMemcpyDeviceToHost);
  printf("no way!\n");
  for (int i = 0; i < 32; i++)
    printf("the idea is that %d\n", temp[i]);
  return 0;
}

__global__ void MergeSort(AABB *input, int N, AABB *output, int total,
                          int dim) {

  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int start = index * N;
  int end = (index + 1) * N;
  if (start >= total)
    return;

  int j = start;
  int k = start + (end - start) / 2;
  for (int i = start; i < end; i++) {
    if (j >= start + (end - start) / 2) {
      output[i] = input[k];
      k++;
    } else if (k >= end) {
      output[i] = input[j];
      j++;
    } else if (input[j].lo.val[dim] <= input[k].lo.val[dim]) {
      output[i] = input[j];
      j++;
    } else {
      output[i] = input[k];
      k++;
    }
  }

  memcpy(&(input[start]), &(output[start]), sizeof(AABB) * N);
}

__global__ void findOverlap(AABB *input, int batchSize, int *overlap, int total,
                            int dim) {

  int index = blockIdx.x * blockDim.x + threadIdx.x;
  if (index > total) {
    return;
  }
  // printf("check %f \n", input[index].lo.val[2]);
  int startIndex = batchSize * index;
  int endIndex = batchSize * (index + 1);
  for (int i = startIndex; i < endIndex; i++) {
    for (int j = i + 1; j < total; j++) {

      if (input[i].hi.val[dim] < input[j].lo.val[dim])
        break;
      overlap[input[i].id * total + input[j].id] = 1;
      // printf("index is %d, %d, %d\n", input[i].id ,input[j].id, input[i].id
      // *total + input[j].id);
    }

    for (int j = i - 1; j >= 0; j--) {

      if (input[i].lo.val[dim] > input[j].hi.val[dim])
        break;
      overlap[input[j].id * total + input[i].id] = 1;
    }
  }
}

int sort_AABB(AABB *res, int N, int *overlap) {
  int original_size = N;
  int sort_block = 2;

  // N = nextPow2(N);

  AABB *output;
  cudaMalloc(&output, sizeof(AABB) * N);
  // get_info(res);

  for (int dim = 0; dim < 3; dim++) {
    while (sort_block <= N) {
      MergeSort<<<BLOCK_SIZE, (N / BLOCK_SIZE) + 1>>>(res, sort_block, output,
                                                      N, dim);
      cudaDeviceSynchronize();
      sort_block *= 2;
    }
    // get_info(output);
    // int* dev = new int[1024];

    // cudaMemcpy(dev, overlap, sizeof(int)*1024, cudaMemcpyDeviceToHost);
    // for (int i = 0; i < 1024;i++){
    //   printf("%d ", dev[i]);
    // }
    AABB *dev = new AABB[32];

    cudaMemcpy(dev, res, sizeof(AABB) * 32, cudaMemcpyDeviceToHost);
    // for (int i = 0; i < 32;i++){
    //   printf("%f %d, \n", dev[i].lo.val[dim], dev[i].id);
    // }
    // printf("\n\n");

    findOverlap<<<BLOCK_SIZE, (N / BLOCK_SIZE) + 1>>>(res, 1, overlap, N, dim);
    // get_info(res);
    cudaDeviceSynchronize();

    // printf("\n");

    // printf("shall print something %d\n", value);
  }
  int value = find_peaks(32 * 32, overlap);

  // get_info1(overlap);
  return value;
}

inline double GT(double a, double b) { return (((a) > (b)) ? (a) : (b)); }

double findRadius(AABB *curr, Object *b) {
  double val = 0.0;

  for (int i = 0; i < (b->num_tris); i++) {
    double cur_rad1_sq = 0;
    double cur_rad2_sq = 0;
    double cur_rad3_sq = 0;
    for (int w = 0; w < 3; w++) {
      double my_num1 = curr->center[w] - b->tris[i].p1[w];
      double my_num2 = curr->center[w] - b->tris[i].p2[w];
      double my_num3 = curr->center[w] - b->tris[i].p3[w];
      cur_rad1_sq += pow(my_num1, 2);
      cur_rad2_sq += pow(my_num2, 2);
      cur_rad3_sq += pow(my_num3, 2);
    }

    double max_rad_sq = GT(cur_rad1_sq, GT(cur_rad2_sq, cur_rad3_sq));

    val = GT(max_rad_sq, val);
  }
  return sqrt(val) * 1.0001;
}

void findCenter(AABB *curr, Object *b) {
  curr->center[0] = curr->center[1] = curr->center[2] = 0.0;

  for (int dim = 0; dim < 3; dim++) {
    for (int i = 0; i < (b->num_tris); i++) {

      curr->center[dim] +=
          b->tris[i].p1[dim] + b->tris[i].p2[dim] + b->tris[i].p3[dim];
    }

    curr->center[dim] /= (3 * b->num_tris);
  }
}

VCInternal::VCInternal(int mySize, int ss) {
  state = VCstate_default;
  next_id = 0;
  screen_size = ss;
  vc_objects = new VCObject *[mySize]; // allocate the array.
  size = mySize;

  cudaMalloc(&overlaps, sizeof(int) * mySize * mySize + 2);

  int i;
  for (i = 0; i < mySize; i++)
    vc_objects[i] = NULL;
  boxes = new AABB[mySize];

  cudaMalloc(&cuda_boxes, mySize * sizeof(AABB));
  cudaMemcpy(cuda_boxes, boxes, mySize * sizeof(AABB), cudaMemcpyHostToDevice);
}

VCInternal::~VCInternal() {}

void VCInternal::AddObject(int id, Object *b) // add a new object
{
  AABB *curr = &boxes[id];
  curr->id = id; // set the id to the given value.

  findCenter(curr, b);
  curr->radius = findRadius(curr, b);
  EndPoint lo = (EndPoint){.minmax = MIN};
  EndPoint hi = (EndPoint){.minmax = MAX};
  curr->lo = lo;
  curr->hi = hi;
  for (int w = 0; w < 3; w++) {
    curr->lo.val[w] = curr->center[w] - curr->radius;
    curr->hi.val[w] = curr->center[w] + curr->radius;
  }

  // cudaMalloc(&cudacurr, sizeof(AABB));
  cudaMemcpy(&cuda_boxes[id], curr, sizeof(AABB), cudaMemcpyHostToDevice);

  // print_kernel<<<1, 1>>>(cudacurr);
  // cudaDeviceSynchronize();
}
// 1. check if the size fit in
// 2 assign the object an id and activate the object
int VCInternal::NewObject(int *id) // create a new object in the database.
{

  // allocate a new object.
  vc_objects[next_id] = new VCObject;

  *id = next_id; // for returning the id generated by VCollide.
  current_id = next_id;
  vc_objects[next_id]->id = next_id;
  vc_objects[next_id]->b = new Object;
  vc_objects[next_id]->b->BeginModel();
  //_state = 1;//default the object is activate
  next_id++;

  return 0;
}

int VCInternal::AddTri(double v1[], double v2[], double v3[]) {

  vc_objects[current_id]->b->AddTri(v1, v2, v3, 0); // add triangle.
  return 0;
}

// 1. add current object to n body
// 2. have RAPID build the OBB tree.
// 3. initialize trans
int VCInternal::EndObject(void) {

  vc_objects[current_id]->b->EndModel();
  AddObject(current_id, vc_objects[current_id]->b);

  memset(((void *)vc_objects[current_id]->trans), 0, 16 * sizeof(double));
  vc_objects[current_id]->trans[0] = 1.0;
  vc_objects[current_id]->trans[4 + 1] = 1.0;
  vc_objects[current_id]->trans[8 + 2] = 1.0;
  vc_objects[current_id]->trans[12 + 3] = 1.0;
  return 0;
}

// __global__ void checking(){
//   printf("hello from cuda!\n");
// }
int VCInternal::EndAllObjects(void) {

  AABB *temp = new AABB[size];
  cudaMemcpy(temp, cuda_boxes, size * sizeof(AABB), cudaMemcpyDeviceToHost);
  printf("no way!\n");
  for (int i = 0; i < 20; i++)
    // printf("the idea is that %d\n",temp[i].id );
    return 0;
}

__global__ void cuda_update_trans(int id_max, double *trans, AABB *cuda_boxes) {
  int id = blockIdx.x * blockDim.x + threadIdx.x;
  if (id >= id_max)
    return;
  // printf("%d, %d\n", id, id_max);
  //   for (int i = id*16; i < (id+1)*16;i++){
  // printf("%f\n", trans[i]);
  //   }
  AABB *current = &cuda_boxes[id];

  // for (int dim = 0; dim < 3; dim ++){
  //   printf("%d, %f\n", id*16+dim*4, current->center[0] * trans[id*16+dim*4] +
  //   current->center[1] * trans[id*16+dim*4+1] + current->center[2] *
  //   trans[id*16+dim*4+2] + trans[id*16+dim*4+3]);
  // }

  double new_center[3];
  for (int dim = 0; dim < 3; dim++) {
    new_center[dim] = current->center[0] * trans[id * 16 + dim * 4] +
                      current->center[1] * trans[id * 16 + dim * 4 + 1] +
                      current->center[2] * trans[id * 16 + dim * 4 + 2] +
                      trans[id * 16 + dim * 4 + 3];
    current->lo.val[dim] = new_center[dim] - current->radius;
    current->hi.val[dim] = new_center[dim] + current->radius;
  }
  for (int dim = 0; dim < 3; dim++) {
    current->center[dim] = new_center[dim];
  }
  current->sorting_center = new_center[0] * id + new_center[1] + new_center[2];
}

int VCInternal::UpdateAllTrans(int id[], int total, double *trans) {

  // EndAllObjects();
  // for (int i = 0; i < 16*total;i++){
  //   printf("%f ", trans[i]);
  // }

  for (int j = 0; j < total; j++) {
    VCObject *current = vc_objects[id[j]];
    memcpy((void *)current->trans, (void *)(&trans[16 * id[j]]),
           16 * sizeof(double));
  }
  double *temp;
  cudaMalloc(&temp, sizeof(double) * total * 17);

  cudaMemcpy(temp, trans, sizeof(double) * total * 17, cudaMemcpyHostToDevice);

  cuda_update_trans<<<BLOCK_SIZE, (total / BLOCK_SIZE) + 1>>>(total, temp,
                                                              cuda_boxes);
  cudaDeviceSynchronize();

  // EndAllObjects();

  cudaMemset(overlaps, 0, (size * size) * sizeof(int));
  // printf("size is %d", size);

  overlap_count = sort_AABB(cuda_boxes, size, overlaps);

  // cuda_get_collision<<<32, 32>>>(&total, cuda_nbody);
  // cudaDeviceSynchronize();
  return 0;
}

__device__ void cuda_test(double R1[3][3], double T1[3], double s1, Object *o1,
                          double R2[3][3], double T2[3], double s2,
                          Object *o2) {}
__global__ void cuda_collide(int N, int *overlaps, double *trans) {
  int index = blockIdx.x + blockDim.x + threadIdx.x;
  if (index > N)
    return;
  int val = overlaps[index];
  int i = val / N;
  int j = val % N;

  double R1[3][3], T1[3], R2[3][3], T2[3];
  for (int index = 0; index < 9; index++) {
    int x = index / 3;
    int y = index % 3;
    R1[x][y] = trans[index * 16 + x * 4 + y];
    R2[x][y] = trans[index * 16 + x * 4 + y];
  }

  for (int x = 0; x < 3; x++) {
    T1[x] = trans[index * 16 + x * 4 + 3];
    T2[x] = trans[index * 16 + x * 4 + 3];
  }
}
int VCInternal::all_Collide(void) { // perform collision detection.

  return 1;
}

int VCInternal::Collide(void) // perform collision detection.
{

  // std::cout<< nbody.overlapping_pairs.size<< std::endl;
  int *dev = new int[overlap_count];
  cudaMemcpy(dev, overlaps, sizeof(int) * overlap_count,
             cudaMemcpyDeviceToHost);
  // printf("overlapCount %d \n", overlap_count);
  for (int k = 0; k < overlap_count; k++) {

    // for(int i = 0; i< 32;i++){
    //   for(int j = i+1; j< 32;j++){
    int val = dev[k];
    int i = val / size;
    int j = val % size;
    if (i == j)
      continue;
    double R1[3][3], T1[3], R2[3][3], T2[3];
    for (int index = 0; index < 9; index++) {
      int x = index / 3;
      int y = index % 3;
      R1[x][y] = vc_objects[i]->trans[x * 4 + y];
      R2[x][y] = vc_objects[j]->trans[x * 4 + y];
    }

    for (int index = 0; index < 3; index++) {
      T1[index] = vc_objects[i]->trans[index * 4 + 3];
      T2[index] = vc_objects[j]->trans[index * 4 + 3];
    }

    // call the RAPID collision detection routine.
    ::Collide(R1, T1, vc_objects[i]->b, R2, T2, vc_objects[j]->b,
              FIRST_CONTACT);

    if (Object_num_contacts != 0) {
      printf("collision between object %d, and object %d!\n", i, j);
    }
  }

  return 0;
}