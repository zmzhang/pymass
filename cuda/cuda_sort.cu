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

	__global__ void find_pics_k(Eigen::Vector3f * seeds_dev, Eigen::Vector3f * regions_dev, int *ids_dev, int n)
	{
		int idx = blockIdx.x * blockDim.x + threadIdx.x;
		if (idx < n)
		{
			Eigen::Vector3f seed = seeds_dev[idx];
			Eigen::Vector3f region = regions_dev[ids_dev[idx]];

		}
		return;
	}
	
	double find_pics(const std::vector<Eigen::Vector3f> & seeds, const std::vector<Eigen::Vector3f> & regions, const std::vector<int> & ids)
	{
		int sz_seed = seeds.size();
		int sz_rg   = regions.size();
		Eigen::Vector3f *seeds_dev;
		HANDLE_ERROR(cudaMalloc((void **)&seeds_dev, sizeof(Eigen::Vector3f)*sz_seed));
		HANDLE_ERROR(cudaMemcpy(seeds_dev, seeds.data(), sizeof(Eigen::Vector3f)*sz_seed, cudaMemcpyHostToDevice));


		Eigen::Vector3f * regions_dev;
		HANDLE_ERROR(cudaMalloc((void **)&regions_dev, sizeof(Eigen::Vector3f)*sz_rg));
		HANDLE_ERROR(cudaMemcpy(regions_dev, regions.data(), sizeof(Eigen::Vector3f)*sz_rg, cudaMemcpyHostToDevice));

		int * ids_dev;
		HANDLE_ERROR(cudaMalloc((void **)&ids_dev, sizeof(int)*(sz_seed+1)));
		HANDLE_ERROR(cudaMemcpy(ids_dev, ids.data(), sizeof(int)*(sz_seed + 1), cudaMemcpyHostToDevice));


		find_pics_k << <(sz_seed + 1023) / 1024, 1024 >> >(seeds_dev, regions_dev, ids_dev, sz_seed);

		cudaFree(ids_dev);
		cudaFree(regions_dev);
		cudaFree(seeds_dev);

		return 0.0;
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

void printVV(const std::vector<Eigen::Vector3f> & vec, int n)
{
	Eigen::MatrixXf m(vec.size(), 3);
	int i = 0;
	std::for_each(vec.begin(), vec.end(), [&m, &i](const Eigen::Vector3f & v) {
		m.row(i) = v;
		i++; });
	if (n<= vec.size() && n>0)
	{
		cout << m.topRows(n) << endl;
	}
	else
	{
		cout << m << endl;
	}
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

#include <Eigen/Core>
#include <algorithm>
#include <vector>

template <typename Scalar, int rows, int cols, int options, int maxRows, int maxCols>
Eigen::Matrix<Scalar, rows, cols, options, maxRows, maxCols> sortMatrix(
	Eigen::Matrix<Scalar, rows, cols, options, maxRows, maxCols> target, int col
){
	std::vector<Eigen::Matrix<Scalar, 1, cols>> matrixRows;
	for (unsigned int i = 0; i < target.rows(); i++)
		matrixRows.push_back(target.row(i));
	std::sort(
		matrixRows.begin(),
		matrixRows.end(),
		[&col](const Eigen::Matrix<Scalar, 1, cols> & a,const Eigen::Matrix<Scalar, 1, cols> & b)->bool
		{
			return a(0, col) < b(0, col);
		}
	);

	Eigen::Matrix<Scalar, rows, cols, options, maxRows, maxCols> sorted;
	for (unsigned int i = 0; i < matrixRows.size(); i++)
		sorted.row(i) = matrixRows[i];
	return sorted;
}



std::vector<Eigen::Vector3f> pic_seeds(const Eigen::MatrixXf & m, float mz_tol, int num_seed)
{
	auto comp = [](const Eigen::VectorXf& lhs, const Eigen::VectorXf& rhs) -> bool {
		return lhs[1] < rhs[1]; 
	};
	std::set<Eigen::VectorXf, mz_comp> seed_set(comp);

	for (int i =0; i< m.rows(); i++)
	{
		auto it = seed_set.lower_bound(m.row(i));
		if (seed_set.size()==0)
		{
			seed_set.insert(m.row(i));
		}
		else
		{
			if (it == seed_set.end())
			{
				if (m.row(i)[1] - (*std::prev(seed_set.end()))[1] > mz_tol)
				{
					seed_set.insert(m.row(i));
				}
			}
			else if (it == seed_set.begin())
			{
				if ((*seed_set.begin())[1] - m.row(i)[1] > mz_tol)
				{
					seed_set.insert(m.row(i));
				}
			}
			else
			{
				if (  ((*it)[1] - m.row(i)[1] > mz_tol ) && 
					  ( m.row(i)[1] - (*std::prev(it))[1]> mz_tol))
				{
					seed_set.insert(m.row(i));
				}
			}
		}

		if (seed_set.size()==num_seed)
		{
			break;
		}
	}

	std::vector<Eigen::Vector3f> ret(seed_set.size());

	int i = 0;
	std::for_each(seed_set.begin(), seed_set.end(), [&ret, &i](const Eigen::VectorXf & v) {
		ret[i] = v;
		i++; });

	return ret;
}

std::tuple<std::vector<Eigen::Vector3f>, std::vector<int> > regions_of_seeds(LCMS & lcms, const std::vector<Eigen::Vector3f> & seeds, float peak_width, float mz_tol)
{
	std::vector<std::vector<Eigen::Vector3f>> region_vec;
	std::vector<int>                          ids(seeds.size()+1);

	int rows = 0;
	ids[0] = 0;
	for (int i=0; i<seeds.size(); i++)
	{
		Eigen::Vector3f seed = seeds[i];
		std::vector<Eigen::Vector3f> region = lcms.getRegion(seed[0] - peak_width, seed[0] + peak_width, seed[1] - mz_tol, seed[1] + mz_tol);
		region_vec.push_back(region);
		rows += region.size();
		ids[i + 1] = rows;
	}
	std::vector<Eigen::Vector3f>              regions(rows);

	for (int i = 0; i < seeds.size(); i++)
	{
		int sz = ids[i + 1] - ids[i];
		std::copy_n(region_vec[i].begin(), sz, &regions[ids[i]]);
	}

	return std::make_tuple(regions, ids);
}

void processLCMS(LCMS & lcms)
{
	cudaFree(0);
	cout << "using lcms object in NVCC, and its scan size is: " << lcms.m_massScans.size() << endl;

	Eigen::MatrixXf rmv = lcms.getAll();

	gtic();
	sort_by_col(rmv, 2);
	gtoc();

	gtic();
	std::vector<Eigen::Vector3f> seeds = pic_seeds(rmv, 0.05f, 8000);
	gtoc();

	gtic();
	std::vector<Eigen::Vector3f> regions;
	std::vector<int>             ids;
	std::tie(regions, ids) = regions_of_seeds(lcms, seeds, 50.0f, 0.05f);
	gtoc();

	gtic();
	double x = Kernel::find_pics(seeds, regions, ids);
	gtoc();
	cout << "Calculated by CUDA kernel: " << x << endl;

	cout << "##########################"<<endl;
}