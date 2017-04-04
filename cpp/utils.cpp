#include "utils.h"
#include <chrono>
#include <thread>
#include <iostream>

using namespace std;
using namespace Eigen; 

void clip(Eigen::VectorXi & idx, int s, int e)
{
	for_each(idx.data(), idx.data()+idx.size(), [&]( int & n){ 
		if(n<s) n=s;
		if(n>e) n=e;
	});
}

Eigen::VectorXi searchsorted(const Eigen::VectorXf& v, const Eigen::VectorXf& t)
{
	VectorXi out(t.size());
	for_each(t.data(), t.data()+t.size(), [&](const float & n){ 
	  auto lower=lower_bound (v.data(), v.data()+v.size(),n);          
	  out[&n-t.data()]= static_cast<int>(lower- v.data());
	});
	return out;
} 

Eigen::VectorXi findclosest(const Eigen::VectorXf& v, const Eigen::VectorXf& t)
{
	if (v.size()==1 || v.size() == 0) return VectorXi::Zero(t.size());
	VectorXi idx = searchsorted(v,t);
	clip(idx,1, static_cast<int>(v.size()-1));
	Eigen::VectorXf right=slice(v,idx);
	Eigen::VectorXf left=slice(v,(idx.array()-1));
	idx.array()-=(t.array() - left.array() < right.array() - t.array()).select(1,VectorXi::Zero(t.size())).array();
	return idx;
}

void finishProcess()
{
	cout << "finish Process"<<endl;
	std::this_thread::sleep_for(std::chrono::milliseconds(10));
}


float ReverseFloat(const float inFloat)
{
	float retVal;
	char *floatToConvert = (char*)& inFloat;
	char *returnFloat = (char*)& retVal;

	// swap the bytes into a temporary buffer
	returnFloat[0] = floatToConvert[3];
	returnFloat[1] = floatToConvert[2];
	returnFloat[2] = floatToConvert[1];
	returnFloat[3] = floatToConvert[0];

	return retVal;
}


PYMASS_EXPORT std::stack<clock_t> tictoc_stack;
void tic() {
	tictoc_stack.push(clock());
}

void toc() {
	std::cout << "Time elapsed: "
		<< ((double)(clock() - tictoc_stack.top())) / CLOCKS_PER_SEC
		<< std::endl;
	tictoc_stack.pop();
}
