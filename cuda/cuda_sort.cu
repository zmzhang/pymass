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



namespace Kernel
{

	static void HandleError(cudaError_t err, const char *file, int line)
	{
		if (err != cudaSuccess)
		{
			printf("%s in %s at line %d\n", cudaGetErrorString(err), file, line);
			exit(EXIT_FAILURE);
		}
	}

    #define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))

	__global__ void cu_dot(Eigen::Vector3d *v1, Eigen::Vector3d *v2, double *out, size_t N)
	{
		int idx = blockIdx.x * blockDim.x + threadIdx.x;
		if (idx < N)
		{
			out[idx] = v1[idx].dot(v2[idx]);
		}
		return;
	}
	
	double dot(const std::vector<Eigen::Vector3d> & v1, const std::vector<Eigen::Vector3d> & v2)
	{
		int n = v1.size();
		double *ret = new double[n];

		Eigen::Vector3d *dev_v1, *dev_v2;
		HANDLE_ERROR(cudaMalloc((void **)&dev_v1, sizeof(Eigen::Vector3d)*n));
		HANDLE_ERROR(cudaMalloc((void **)&dev_v2, sizeof(Eigen::Vector3d)*n));
		double* dev_ret;
		HANDLE_ERROR(cudaMalloc((void **)&dev_ret, sizeof(double)*n));

		HANDLE_ERROR(cudaMemcpy(dev_v1, v1.data(), sizeof(Eigen::Vector3d)*n, cudaMemcpyHostToDevice));
		HANDLE_ERROR(cudaMemcpy(dev_v2, v2.data(), sizeof(Eigen::Vector3d)*n, cudaMemcpyHostToDevice));

		cu_dot << <(n + 1023) / 1024, 1024 >> > (dev_v1, dev_v2, dev_ret, n);

		HANDLE_ERROR(cudaMemcpy(ret, dev_ret, sizeof(double)*n, cudaMemcpyDeviceToHost));

		for (int i = 1; i < n; ++i)
		{
			ret[0] += ret[i];
		}

		return ret[0];
	}
}



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


std::set<Eigen::VectorXf, mz_comp> pic_seed(const Eigen::MatrixXf & m, float mz_tol, int num_seed)
{
	auto comp = [](const Eigen::VectorXf& lhs, const Eigen::VectorXf& rhs) -> bool {
		return lhs[1] < rhs[1]; 
	};
	std::set<Eigen::VectorXf, mz_comp> ret(comp);

	for (int i =0; i< m.rows(); i++)
	{
		auto it = ret.lower_bound(m.row(i));
		if (ret.size()==0)
		{
			ret.insert(m.row(i));
		}
		else
		{
			if (it == ret.end())
			{
				if (m.row(i)[1] - (*std::prev(ret.end()))[1] > mz_tol)
				{
					ret.insert(m.row(i));
				}
			}
			else if (it == ret.begin())
			{
				if ((*ret.begin())[1] - m.row(i)[1] > mz_tol)
				{
					ret.insert(m.row(i));
				}
			}
			else
			{
				if (  ((*it)[1] - m.row(i)[1] > mz_tol ) && 
					  ( m.row(i)[1] - (*std::prev(it))[1]> mz_tol))
				{
					ret.insert(m.row(i));
				}
			}
		}

		if (ret.size()==num_seed)
		{
			break;
		}
	}


	Eigen::MatrixXf seed(ret.size(), 3);
	int i = 0;
	std::for_each(ret.begin(), ret.end(), [&seed, &i](const Eigen::VectorXf & v) {
		seed.row(i) = v;
		i++; });

	cout << seed << endl;

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
	pic_seed(rmv, 0.05, 50);
	gtoc();

	//cout << rmv.topRows(10) << endl;


	std::vector<Eigen::Vector3d> v1(4000, Eigen::Vector3d{ 1.0, 1.0, 1.0 });
	std::vector<Eigen::Vector3d> v2(4000, Eigen::Vector3d{ -1.0, 1.0, 1.0 });
	gtic();
	double x = Kernel::dot(v1, v2);
	gtoc();
	cout << "Dot calculated by CUDA kernel: " << x << endl;

	cout << "##########################"<<endl;
}