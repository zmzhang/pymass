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

int main() {
	cudaFree(0);
	int size = 100*1024*1024;
	thrust::host_vector<float> h_values(size);
	thrust::host_vector<int> h_keys(size);
	
	std::generate(h_values.begin(), h_values.end(), rand);

	return 0;
}