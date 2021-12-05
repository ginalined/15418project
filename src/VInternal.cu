#include <iostream>
#include <string.h>     //for memset and memcpy.
#include "VInternal.H"
#include <cuda.h>
#include <cuda_runtime.h>
#include <driver_functions.h>
#include "CycleTimer.h"
#include "objects.h"

#include <math.h>  
#include <thrust/sort.h>
#include <thrust/execution_policy.h>


static inline int nextPow2(int n)
{
    n--;
    n |= n >> 1;
    n |= n >> 2;
    n |= n >> 4;
    n |= n >> 8;
    n |= n >> 16;
    n++;
    return n;
}

int get_info(AABB * input){
  
  AABB *temp = new AABB[16];
  cudaMemcpy(temp, input, 16 * sizeof(AABB), cudaMemcpyDeviceToHost);
  printf("no way!\n");
  for (int i = 0; i < 16;i++)
    printf("the idea is that %f, %d\n",temp[i].center[0] , temp[i].id);
  return 0;
}

__global__ 
void MergeSort(AABB * input, int N, AABB * output, int total, int dim)
{
    
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int start = index * N;
    int end = (index+1) * N ;
    if (start >= total)
        return;
    
    int j = start;
    int k = start + (end - start)/2;
    for (int i = start; i < end; i++){
        if (j >= start + (end - start)/2){
            output[i] = input[k];
            k++;
        }
        else if (k >= end){
            output[i] = input[j];
            j++;
        }
        else if (input[j].center[dim] <= input[k].center[dim]){
            output[i] = input[j];
            j++;
        }else{
            output[i] = input[k];
            k++;
        }
    }
    
    memcpy(&(input[start]), &(output[start]),sizeof(AABB) * N);
}

__global__ void findOverlap(AABB * input, int batchSize, int * overlap, int total, int dim)
{
  //printf("this and next is %f, %f\n",input[i].hi.val[dim], input[j].lo.val[dim] );
  printf("shall print something \n");
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int startIndex = batchSize * index;
  int endIndex = batchSize * (index+1);
  for(int i = startIndex; i < endIndex; i++){
    for (int j = i+1; j < total; i++){
      
      if (input[i].hi.val[dim] < input[j].lo.val[dim])
        break;
      overlap[input[i].id * total + input[j].id] = 1;
    }
  }      
}

void sort_AABB(AABB * res, int N, int * overlap){
  int original_size = N;
  int sort_block = 2;
  
  N = nextPow2(N);

  AABB * output;
  cudaMalloc(&output, sizeof(AABB) * N);
  //get_info(res);  
  for (int dim = 0; dim < 1; dim++){
  while (sort_block <= N)
	{
		MergeSort<<<4,4>>>(res, sort_block, output, N,  dim);
		cudaDeviceSynchronize();        
		sort_block *= 2; 
	}

  findOverlap<<<4,4>>>(res, 1, overlap, N, dim);
  printf("shall print something \n");

}
}





inline double GT(double a, double b)
{
  return (( (a) > (b) ) ? (a) : (b));
}


double  findRadius(AABB *curr, Object *b){
  double val = 0.0;

  for (int i=0; i<(b->num_tris); i++)
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
      
      val = GT(max_rad_sq, val);
    }
    return sqrt(val) * 1.0001;
}

void findCenter(AABB *curr, Object *b){
  curr->center[0] = curr->center[1] = curr->center[2] = 0.0;
  
  for (int dim = 0; dim < 3; dim++){
  for (int i=0; i<(b->num_tris); i++){
    
    curr->center[dim] += b->tris[i].p1[dim] + b->tris[i].p2[dim] + b->tris[i].p3[dim];
  }
      
  curr->center[dim] /= (3*b->num_tris); 
  
  }
}







// __global__ void initAABB(int mySize, AABB * box_pointer){
//   box_pointer = new AABB[mySize];
//   box_pointer[0].id = 0;
//   printf("hello everyone\n");
//   printf("hello  %d, \n", box_pointer[0].id);
// }


VCInternal::VCInternal(int mySize, int ss)
{
  state = VCstate_default;
  next_id = 0;
  screen_size = ss;
  vc_objects = new VCObject*[mySize+20]; //allocate the array.
  size = mySize;
  cudaMalloc(&overlaps, sizeof(int) * mySize*mySize+2);
  
  int i;
  for (i=0; i<mySize; i++)
    vc_objects[i] = NULL;
  boxes = new AABB[mySize];
  
  cudaMalloc(&cuda_boxes, mySize * sizeof(AABB));
  cudaMemcpy(cuda_boxes, boxes, mySize * sizeof(AABB), cudaMemcpyHostToDevice);
}


VCInternal::~VCInternal()
{}

void VCInternal::AddObject(int id, Object *b) //add a new object
{
  AABB *curr = &boxes[id];
  curr->id = id; //set the id to the given value.
  
  findCenter(curr, b);
  curr->radius = findRadius(curr, b);
  EndPoint lo = (EndPoint){.minmax = MIN};
  EndPoint hi = (EndPoint){.minmax = MAX};
  curr->lo = lo;
  curr->hi = hi;
  for (int w=0; w<3; w++){
  curr->lo.val[w] = curr->center[w] - curr->radius; 
  curr->hi.val[w] = curr->center[w] + curr->radius;
  }

  // cudaMalloc(&cudacurr, sizeof(AABB));
  cudaMemcpy(&cuda_boxes[id], curr, sizeof(AABB), cudaMemcpyHostToDevice);

  // print_kernel<<<1, 1>>>(cudacurr);
  // cudaDeviceSynchronize();
  

}
//1. check if the size fit in
//2 assign the object an id and activate the object
int VCInternal::NewObject(int *id) //create a new object in the database.
{
 

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

 
  vc_objects[current_id]->b->EndModel();
  AddObject(current_id, vc_objects[current_id]->b);

  memset( ( (void *)vc_objects[current_id]->trans), 0, 16*sizeof(double) );
  vc_objects[current_id]->trans[0] = 1.0;
  vc_objects[current_id]->trans[4+1] = 1.0;
  vc_objects[current_id]->trans[8+2] = 1.0;
  vc_objects[current_id]->trans[12+3] = 1.0;
  return 0;
  
}

// __global__ void checking(){
//   printf("hello from cuda!\n");
// }
int VCInternal::EndAllObjects(void){
  
  AABB *temp = new AABB[size];
  cudaMemcpy(temp, cuda_boxes, size * sizeof(AABB), cudaMemcpyDeviceToHost);
  printf("no way!\n");
  for (int i = 0; i < 16;i++)
    printf("the idea is that %f\n",temp[i].center[0] );
  return 0;
}

__global__ void cuda_update_trans(int id_max, double * trans,  AABB * cuda_boxes){
  int id = blockIdx.x * blockDim.x + threadIdx.x;
  if (id >= id_max)
    return;
  //printf("%d, %d\n", id, id_max);
//   for (int i = id*16; i < (id+1)*16;i++){
// printf("%f\n", trans[i]);
//   }
 AABB *current = &cuda_boxes[id];
 
// for (int dim = 0; dim < 3; dim ++){
//   printf("%d, %f\n", id*16+dim*4, current->center[0] * trans[id*16+dim*4] + current->center[1] * trans[id*16+dim*4+1] + current->center[2] * trans[id*16+dim*4+2] + trans[id*16+dim*4+3]);
// }
 

  double new_center[3];
  for (int dim = 0; dim < 3; dim ++){
    new_center[dim] = current->center[0] * trans[id*16+dim*4] + current->center[1] * trans[id*16+dim*4+1] + current->center[2] * trans[id*16+dim*4+2] + trans[id*16+dim*4+3];
    current->lo.val[dim] = new_center[dim] - current->radius;
    current->hi.val[dim] = new_center[dim] + current->radius;
  }
  for (int dim = 0; dim < 3; dim ++){
    current->center[dim] = new_center[dim];
  }
  current->sorting_center = new_center[0] * id + new_center[1]+ new_center[2];
  

}




int VCInternal::UpdateAllTrans(int id[], int total, double * trans)
{     
  
  //EndAllObjects();
  // for (int i = 0; i < 16*total;i++){
  //   printf("%f ", trans[i]);
  // }
  //printf("\n");
  for (int j=0; j<total; j++){
    
   VCObject *current = vc_objects[id[j]];
  memcpy((void *)current->trans, (void *)(&trans[16*j]), 16*sizeof(double));
       
    }
    double *temp;
    cudaMalloc(&temp, sizeof(double) * total*17);

    cudaMemcpy(temp, trans, sizeof(double) * total*17, cudaMemcpyHostToDevice);


    cuda_update_trans<<<4, 4>>>(total, temp, cuda_boxes);
    
    cudaDeviceSynchronize();
    //EndAllObjects();
    cudaMemset(overlaps, 0, size * size);
    sort_AABB(cuda_boxes, size, overlaps);

    //cuda_get_collision<<<32, 32>>>(&total, cuda_nbody);
    // cudaDeviceSynchronize();
    return 0;
}


int VCInternal::Collide(void)  //perform collision detection.
{

  

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
        R1[x][y] = vc_objects[i]->trans[x*4+y];
        R2[x][y] = vc_objects[curr_ovrlp->id]->trans[x*4+y];
      }

      for (int index = 0; index < 3; index++){
        T1[index] = vc_objects[i]->trans[index*4+3];
        T2[index] = vc_objects[curr_ovrlp->id]->trans[index*4+3];
      }

			      //call the RAPID collision detection routine.
			::Collide(R1, T1, vc_objects[i]->b, R2, T2, vc_objects[curr_ovrlp->id]->b, FIRST_CONTACT);
			      

			if (Object_num_contacts != 0)
          {
            printf("hello!\n");
          }
			      
		  curr_ovrlp = curr_ovrlp->next;
		}
    }
  return 0;
}




