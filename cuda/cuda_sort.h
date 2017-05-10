#ifndef CUDA_SORT_H
#define CUDA_SORT_H
#include <chrono>
#include <stack>
#include "LCMS.h"

void thrustSort(float *V, int *K, int N);
extern std::stack<clock_t> gtictoc_stack;
void gtic();
void gtoc();

void processLCMS(const LCMS & lcms);


#endif