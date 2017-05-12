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
void sort_by_col(Eigen::MatrixXf & m, int col);

typedef std::function<bool(const Eigen::VectorXf &lhs, const Eigen::VectorXf &rhs)> mz_comp;
std::set<Eigen::VectorXf, mz_comp> pic_seed(const Eigen::MatrixXf & m,float mz_tol, int num_seed);
void processLCMS(LCMS & lcms);


#endif