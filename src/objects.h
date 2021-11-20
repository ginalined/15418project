
#ifndef OBJECTS_H_
#define OBJECTS_H_


//definition of tri 
struct tri
{
  int id;
  double p1[3], p2[3], p3[3];
};


class box
{
public:

  // placement in parent's space
  // box to parent space: x_m = pR*x_b + pT
  // parent to box space: x_b = pR.T()*(x_m - pT)
  double pR[3][3];
  double pT[3];
  
  // dimensions
  double d[3];        // this is "radius", that is, 
                      // half the measure of a side length

  box *P;  // points to but does not "own".  
  box *N;

  tri *trp;

  int leaf() { return (!P && !N); } 
  double size() { return d[0]; } 

  int split_recurse(int *t, int n);
  int split_recurse(int *t);               // specialized for leaf nodes
};


int 
tri_contact (double *P1, double *P2, double *P3,
	     double *Q1, double *Q2, double *Q3);

int
obb_disjoint(double B[3][3], double T[3], double a[3], double b[3]);


class Object
{
public:
  // these are only for internal use

  box *b;
  int num_boxes_alloced;

  tri *tris;
  int num_tris;
  int num_tris_alloced;

  int build_state;
  
  int build_hierarchy();
  
  int friend Collide(double R1[3][3], double T1[3], 
		 double s1, Object *RAPID_model1,
		 double R2[3][3], double T2[3], 
		 double s2, Object *RAPID_model2,
		 int flag);
public:

  // these are for the client

  Object();
  ~Object();
  
  int BeginModel();
  int AddTri(const double *p1, const double *p2, const double *p3, int id);
  int EndModel();
  
};

/****************************************************************************/

// these are for the client

const int ALL_CONTACTS = 1;    // Find all pairwise intersecting triangles

const int FIRST_CONTACT = 2;   // Just report one intersecting triangle pair
                               //   if there are any.

// this is the collision query invocation.  It assumes that the 
// models are not being scaled up or down, but have their native
// dimensions.
int 
Collide(double R1[3][3], double T1[3], Object *o1,
	double R2[3][3], double T2[3], Object *o2,
	int flag = ALL_CONTACTS);

// this collision query permits the models to each be scaled by
// some nonnegative factor.
int 
Collide(double R1[3][3], double T1[3], double s1, Object *o1,
	double R2[3][3], double T2[3], double s2, Object *o2,
	int flag = ALL_CONTACTS);

// this is for the client
struct collision_pair
{
  int id1;
  int id2;
};

/****************************************************************************/

extern int Object_first_contact;
extern  int Object_num_box_tests;
extern  int Object_num_tri_tests;
extern  int Object_num_contacts;
extern  struct collision_pair *Object_contact;

#endif





