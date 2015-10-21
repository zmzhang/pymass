#ifndef _UTILS_H_
#define _UTILS_H_

#include <Eigen/Dense>
#include <iostream>    
#include <algorithm>    
#include <vector> 

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

void clip(Eigen::VectorXi & idx, int s, int e);
Eigen::VectorXi searchsorted(const Eigen::VectorXf& v, const Eigen::VectorXf& t);
Eigen::VectorXi findclosest(const Eigen::VectorXf& v, const Eigen::VectorXf& t);

#endif