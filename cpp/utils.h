#ifndef _UTILS_H_
#define _UTILS_H_

#include <Eigen/Dense>
#include <iostream>    
#include <algorithm>    
#include <vector> 
#include <stack>
#include <ctime>
#include "pymass_export.h"

template <typename T>
T slice(const T& full, const Eigen::VectorXi & ind)
{
    int num_indices = static_cast<int>(ind.innerSize());
    T target(num_indices);
    for (int i = 0; i < num_indices; i++)
    {
        target[i] = full[ind[i]];
    }
    return target;
}

void PYMASS_EXPORT clip(Eigen::VectorXi & idx, int s, int e);
Eigen::VectorXi PYMASS_EXPORT searchsorted(const Eigen::VectorXf& v, const Eigen::VectorXf& t);
Eigen::VectorXi PYMASS_EXPORT findclosest(const Eigen::VectorXf& v, const Eigen::VectorXf& t);
float PYMASS_EXPORT ReverseFloat(const float inFloat);
void PYMASS_EXPORT finishProcess();

extern std::stack<clock_t> tictoc_stack;
void tic();
void toc();

#endif