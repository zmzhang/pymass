#ifndef CUDA_SORT_H
#define CUDA_SORT_H
#include <chrono>
#include <stack>
#include <set>
#include "LCMS.h"

extern std::stack<clock_t> gtictoc_stack;
void gtic();
void gtoc();
void sort_by_col(Eigen::MatrixXf & m, int col);
std::set<Eigen::VectorXf> pic_seed(const Eigen::MatrixXf & m,float mz_tol);
void processLCMS(LCMS & lcms);


#endif