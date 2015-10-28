#include "utils.h"

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