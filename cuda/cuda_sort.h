#ifndef CUDA_SORT_H
#define CUDA_SORT_H
#include <chrono>
#include <stack>

void thrustSort(float *V, int *K, int N);
extern std::stack<clock_t> gtictoc_stack;
void gtic();
void gtoc();


#endif