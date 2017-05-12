#include <iostream>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/sort.h>
#include <thrust/generate.h>
#include <thrust/copy.h>
#include <algorithm>
#include <cstdlib>
#include "cuda_sort.h"


using namespace std;



std::stack<clock_t> gtictoc_stack;
void gtic() {
	gtictoc_stack.push(clock());
}

void gtoc() {
	std::cout << "Time elapsed: "
		<< ((double)(clock() - gtictoc_stack.top())) / CLOCKS_PER_SEC
		<< std::endl;
	gtictoc_stack.pop();
}

void sort_by_col(Eigen::MatrixXf & m, int col)
{
	for (int i = 0; i < m.cols(); i++)
	{
		thrust::device_vector<float> d_K(m.col(col).data(), m.col(col).data() + m.rows());
		thrust::device_vector<float> d_V(m.col(i).data(), m.col(i).data() + m.rows());
		thrust::sort_by_key(d_K.begin(), d_K.end(), d_V.begin(), thrust::greater<float>());
		if (i!=col)
		{
			thrust::copy(d_V.begin(), d_V.end(), m.col(i).data());
		}
		if (i == m.cols()-1)
		{
			thrust::copy(d_K.begin(), d_K.end(), m.col(col).data());
		}
	}
}

std::set<Eigen::VectorXf> pic_seed(const Eigen::MatrixXf & m, float mz_tol)
{
	std::set<Eigen::VectorXf> ret;
	return ret;
}

void processLCMS(LCMS & lcms)
{
	cudaFree(0);
	cout << "using lcms object in CUDA, and its scan size is: " << lcms.m_massScans.size() << endl;

	Eigen::MatrixXf rmv = lcms.getAll();

	gtic();
	sort_by_col(rmv, 2);
	gtoc();

	gtic();
	pic_seed(rmv, 0.5);
	gtoc();

	cout << rmv.topRows(10) << endl;
	cout << "##########################"<<endl;
}

