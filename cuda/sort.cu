#include <iostream>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/sort.h>
#include <thrust/generate.h>
#include <thrust/copy.h>
#include <algorithm>
#include <cstdlib>


using namespace std;


#include <chrono>
#include <stack>
std::stack<clock_t> tictoc_stack;
void tic() {
	tictoc_stack.push(clock());
}

void toc() {
	std::cout << "Time elapsed: "
		<< ((double)(clock() - tictoc_stack.top())) / CLOCKS_PER_SEC
		<< std::endl;
	tictoc_stack.pop();
}


void thrustSort(float *V, int *K, int N)
{
	thrust::device_vector<float> d_V(V, V+N);
	thrust::device_vector<int> d_K(N);
	thrust::sequence(d_K.begin(), d_K.end(), 0, 1);
	thrust::sort_by_key(d_V.begin(), d_V.end(), d_K.begin());
	thrust::copy(d_K.begin(), d_K.end(), K);
	thrust::copy(d_V.begin(), d_V.end(), V);
}

int main() {
	cudaFree(0);
	int size = 100*1024*1024;
	thrust::host_vector<float> h_values(size);
	thrust::host_vector<int> h_keys(size);
	
	std::generate(h_values.begin(), h_values.end(), rand);

	thrustSort(h_values.data(), h_keys.data(), h_values.size());

	if (size < 100)
	{
		for (int i = 0; i < size; i++)
		{
			cout << "(" << h_values[i] << ", " << h_keys[i] << ")";
		}
	}

	return 0;
}