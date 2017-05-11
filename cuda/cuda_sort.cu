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


void thrustSort(float *V, int *K, int N)
{
	thrust::device_vector<float> d_V(V, V+N);
	thrust::device_vector<int> d_K(N);
	thrust::sequence(d_K.begin(), d_K.end(), 0, 1);
	thrust::sort_by_key(d_V.begin(), d_V.end(), d_K.begin(), thrust::greater<float>());
	thrust::copy(d_K.begin(), d_K.end(), K);
	thrust::copy(d_V.begin(), d_V.end(), V);
}

void processLCMS(LCMS & lcms)
{
	cudaFree(0);
	cout << "using lcms object in CUDA, and its scan size is: " << lcms.m_massScans.size() << endl;

	Eigen::MatrixXf rmv = lcms.getAll();
	Eigen::VectorXi ids(rmv.rows());
	gtic();
	thrustSort(rmv.col(2).data(), ids.data(), rmv.rows());
	gtoc();

	//cout << lcms.m_massScans[0].mz << endl;
	cout << rmv.col(2).head(20) << endl;
	cout << "##########################"<<endl;
	cout << ids.head(20) << endl;
}

