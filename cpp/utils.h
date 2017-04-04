#ifndef _UTILS_H_
#define _UTILS_H_

#include <Eigen/Dense>
#include <iostream>    
#include <algorithm>    
#include <vector> 
#include <stack>
#include <ctime>
#include "pymass_export.h"

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

void PYMASS_EXPORT clip(Eigen::VectorXi & idx, int s, int e);
Eigen::VectorXi PYMASS_EXPORT searchsorted(const Eigen::VectorXf& v, const Eigen::VectorXf& t);
Eigen::VectorXi PYMASS_EXPORT findclosest(const Eigen::VectorXf& v, const Eigen::VectorXf& t);
float PYMASS_EXPORT ReverseFloat(const float inFloat);
void PYMASS_EXPORT finishProcess();

PYMASS_EXPORT extern std::stack<clock_t> tictoc_stack;
void PYMASS_EXPORT tic();
void PYMASS_EXPORT toc();

#include <immintrin.h>
static inline void byteswap_avx2(uint32_t* intsToSwap, uint64_t num_points) {

	uint64_t number;

	const unsigned int nPerSet = 8;
	const uint64_t     nSets = num_points / nPerSet;

	uint32_t* inputPtr = intsToSwap;

	const uint8_t shuffleVector[32] = { 3, 2, 1, 0, 7, 6, 5, 4, 11, 10, 9, 8, 15, 14, 13, 12, 19, 18, 17, 16, 23, 22, 21, 20, 27, 26, 25, 24, 31, 30, 29, 28 };

	const __m256i myShuffle = _mm256_loadu_si256((__m256i*) &shuffleVector);

	for (number = 0; number < nSets; number++) {

		// Load the 32t values, increment inputPtr later since we're doing it in-place.
		const __m256i input = _mm256_loadu_si256((__m256i*)inputPtr);
		const __m256i output = _mm256_shuffle_epi8(input, myShuffle);

		// Store the results
		_mm256_storeu_si256((__m256i*)inputPtr, output);
		inputPtr += nPerSet;
	}
	_mm256_zeroupper();

	// Byteswap any remaining points:
	for (number = nSets * nPerSet; number < num_points; number++) {
		uint32_t outputVal = *inputPtr;
		outputVal = (((outputVal >> 24) & 0xff) | ((outputVal >> 8) & 0x0000ff00) | ((outputVal << 8) & 0x00ff0000) | ((outputVal << 24) & 0xff000000));
		*inputPtr = outputVal;
		inputPtr++;
	}
}

#endif