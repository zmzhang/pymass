#ifndef CUDA_SORT_H
#define CUDA_SORT_H
#include <chrono>
#include <stack>
#include <set>
#include <functional>
#include "LCMS.h"

extern std::stack<clock_t> gtictoc_stack;
void gtic();
void gtoc();
void printVV(const std::vector<Eigen::Vector3f> & vec, int n = -1);
void sort_by_col(Eigen::MatrixXf & m, int col);

typedef std::function<bool(const Eigen::VectorXf &lhs, const Eigen::VectorXf &rhs)> mz_comp;
std::vector<Eigen::Vector3f> pic_seeds(const Eigen::MatrixXf & m,float mz_tol, int num_seed);
void processLCMS(LCMS & lcms);


#endif