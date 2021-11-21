
#include <stdio.h>
#include <iostream>
#include <cuda.h>
#include <cuda_runtime.h>
#include <driver_functions.h>
#include "CycleTimer.h"

extern float toBW(int bytes, float sec);
__global__ void
saxpy_kernel(int N, float alpha, float* x, float* y, float* result) {

    // compute overall index from position of thread in current block,
    // and given the block we are in
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    if (index < N)
       result[index] = alpha * x[index] + y[index];
}

void
saxpyCuda(int N, float alpha, float* xarray, float* yarray, float* resultarray) {

    int totalBytes = sizeof(float) * 3 * N;

    // compute number of blocks and threads per block
    const int threadsPerBlock = 512;
    const int blocks = (N + threadsPerBlock - 1) / threadsPerBlock;

    float* device_x;
    float* device_y;
    float* device_result;

    //
    // allocate device memory buffers on the GPU using cudaMalloc
    //
    cudaMalloc(&device_x, N * sizeof(float));
    cudaMalloc(&device_y, N * sizeof(float));
    cudaMalloc(&device_result, N * sizeof(float));


    // start timing after allocation of device memory
    double startTime = CycleTimer::currentSeconds();

    //
    // copy input arrays to the GPU using cudaMemcpy
    //
    double kernelCopyStartTime = CycleTimer::currentSeconds();
    cudaMemcpy(device_x, xarray, N, cudaMemcpyHostToDevice);
    cudaMemcpy(device_y, yarray, N, cudaMemcpyHostToDevice);
    double kernelCopyEndTime = CycleTimer::currentSeconds();


    // run kernel
    double kernelStartTime = CycleTimer::currentSeconds();
    saxpy_kernel<<<blocks, threadsPerBlock>>>(N, alpha, device_x, device_y, device_result);
    cudaThreadSynchronize();
    double kernelEndTime = CycleTimer::currentSeconds();

    //
    // copy result from GPU using cudaMemcpy
    //
    cudaMemcpy(resultarray, device_result, N, cudaMemcpyDeviceToHost);

    // end timing after result has been copied back into host memory
    double endTime = CycleTimer::currentSeconds();

    cudaError_t errCode = cudaPeekAtLastError();
    if (errCode != cudaSuccess) {
        fprintf(stderr, "WARNING: A CUDA error occured: code=%d, %s\n", errCode, cudaGetErrorString(errCode));
    }

    double overallDuration = endTime - startTime;
    double kernelRuntime = kernelEndTime - kernelStartTime;
    printf("Overall: %.3f ms\t\t[%.3f GB/s]\nKernel:%.3f ms\t\t[%.3f GB/s]\n\n", 
            1000.f * overallDuration, toBW(totalBytes, overallDuration), 
            1000.f * kernelRuntime, toBW(totalBytes, kernelRuntime));

    // free memory buffers on the GPU
    cudaFree(device_x);
    cudaFree(device_y);
    cudaFree(device_result);

}