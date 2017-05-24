#pragma once

#include <list>
#include <set>
#include "utils.h"
#include "LCMS.h"

#include <numeric>
using namespace std;


void LCMS::update()
{
	if (m_vecBIC.size() != m_massScans.size() || m_vecRT.size() != m_massScans.size() || m_vecTIC.size() != m_massScans.size())
	{
		m_vecBIC.resize(m_massScans.size());
		m_vecTIC.resize(m_massScans.size());
		m_vecRT.resize(m_massScans.size());
		for (int i = 0; i < m_massScans.size(); i++) {
			m_vecBIC[i]= m_massScans[i].BIC;
			m_vecTIC[i] = m_massScans[i].TIC;
			m_vecRT[i] = m_massScans[i].RT;
		}
	}
}

Eigen::MatrixXf LCMS::getAll()
{
	update();
	int rows = std::accumulate(m_massScans.begin(), m_massScans.end(), 0, [](int s, const MassScan & mss) {
		return s + mss.mz.size();
	});
	Eigen::MatrixXf ret;
	ret.resize(rows, 3);
	int s = 0;
	for (int i = 0; i < m_massScans.size(); i++)
	{
		int step = m_massScans[i].mz.size();
		ret.col(0).segment(s, step) = Eigen::VectorXf::Constant(step, m_massScans[i].RT);
		ret.col(1).segment(s, step) = m_massScans[i].mz;
		ret.col(2).segment(s, step) = m_massScans[i].val;
		s += step;
	}
	return ret;
}

Eigen::VectorXf LCMS::getBIC()
{
	update();
	return m_vecBIC;
}

Eigen::VectorXf LCMS::getRT()
{
	update();
	return m_vecRT;
}

Eigen::VectorXf LCMS::getTIC()
{
	update();
	return m_vecTIC;
}


Eigen::VectorXf LCMS::getMS(int i, int level)
{
	update();
	if (level == 1)
	{
		return m_massScans[i].mz;
	}
	else if (level == 2)
	{
		return m_massScans[i].childs[0]->mz;
	}
	return Eigen::VectorXf(0, 0);
}



Eigen::VectorXf LCMS::getVal(int i, int level)
{
	update();
	if (level == 1)
	{
		return m_massScans[i].val;
	}
	else if (level == 2)
	{
		return m_massScans[i].childs[0]->val;
	}
	return Eigen::VectorXf(0, 0);
}

std::vector<Eigen::Vector4f> LCMS::getRegion(float rt_begin, float rt_end, float mz_begin, float mz_end)
{
	update();
	std::vector<Eigen::VectorXi> mzi_vec;
	Eigen::VectorXf rts(2); rts << rt_begin, rt_end;
	Eigen::VectorXf mzs(2); mzs << mz_begin, mz_end;
	Eigen::VectorXi rti = findclosest(m_vecRT, rts);
	for (int i = rti[0]; i<=rti[1]; i++)
	{
		Eigen::VectorXi mzi = findclosest(m_massScans[i].mz, mzs);
		mzi_vec.push_back(mzi);
	}

	int rows = std::accumulate(mzi_vec.begin(), mzi_vec.end(), 0, [](int s, const Eigen::VectorXi & mzi) {
		return s + (mzi[1] - mzi[0] + 1);
	}); 

	std::vector<Eigen::Vector4f> ret;
	ret.resize(rows);
	int s = 0;
	for (int i = 0; i<mzi_vec.size(); i++)
	{
		for (int j = mzi_vec[i][0]; j<=mzi_vec[i][1]; j++)
		{
			ret[s][0] = m_vecRT[i + rti[0]];
			ret[s][1] = m_massScans[i + rti[0]].mz[j];
			ret[s][2] = m_massScans[i + rti[0]].val[j];
			ret[s][3] = m_massScans[i + rti[0]].id[j];
			s = s + 1;
		}
	}

	return ret;
}

std::vector<Eigen::Vector3f> pic_seeds(const Eigen::MatrixXf & m, const int & idx , const Eigen::VectorXi  & b_inc, float mz_tol)
{
	auto comp = [](const Eigen::VectorXf& lhs, const Eigen::VectorXf& rhs) -> bool {
		return lhs[1] < rhs[1];
	};
	std::set<Eigen::VectorXf, mz_comp> seed_set(comp);

	for (int i = 0; i < idx; i++)
	{
		if (b_inc[i] == 1)
		{
			continue;
		}
		auto it = seed_set.lower_bound(m.row(i));
		if (seed_set.size() == 0)
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
				if (((*it)[1] - m.row(i)[1] >= mz_tol) &&
					(m.row(i)[1] - (*std::prev(it))[1] >= mz_tol))
				{
					seed_set.insert(m.row(i));
				}
			}
		}
	}

	std::vector<Eigen::Vector3f> ret(seed_set.size());

	int i = 0;
	std::for_each(seed_set.begin(), seed_set.end(), [&ret, &i](const Eigen::VectorXf & v) {
		ret[i][0] = v[0];
		ret[i][1] = v[1];
		ret[i][2] = v[2];
		i++; });

	//std::stable_sort(ret.begin(), ret.end(), [](const Eigen::Vector3f & v1, const Eigen::Vector3f & v2)
	//{
	//	return v1[2] > v2[2];
	//});

	return ret;
}

inline int find_closest(const std::vector<Eigen::Vector4f>& rg, const int & l, const int & u, const float & v, const int & col)
{
	auto lower = std::lower_bound(rg.begin() + l, rg.begin() + u, v, [&col](const Eigen::Vector4f & v1, const float& v2) -> bool {
		return v1[col] < v2; });

	if (lower == rg.begin() + l)
	{
		return l;
	}
	if (lower == rg.begin() + u)
	{
		return u - 1;
	}

	if (abs((*lower)[col] - v) < abs((*std::prev(lower))[col] - v))
	{
		return lower - rg.begin();
	}
	else
	{
		return lower - rg.begin() - 1;
	}
}

inline int find_idx(const std::vector<Eigen::Vector4f>& rg, const float & rt, const float & mz, float threshold)
{
	float _rt = rg[find_closest(rg, 0, (int)rg.size(), rt, 0)][0];

	int lower = std::lower_bound(rg.begin() , rg.end(), _rt, [](const Eigen::Vector4f & v1, const float& v2) -> bool {
		return v1[0] < v2; }) -rg.begin();

	int upper = std::upper_bound(rg.begin(), rg.end(), _rt, [](const float& v2, const Eigen::Vector4f & v1) -> bool {
		return v1[0] > v2; }) - rg.begin();
	
	int idx = find_closest(rg, lower, upper, mz, 1);

	if (rg[idx][1] > mz - threshold && rg[idx][1] < mz + threshold)
	{
		return idx;
	}
	else
	{
		return -1;
	}
}

Eigen::MatrixXf FPIC(LCMS & lcms, const Eigen::Vector3f & seed, float rt_width, float mz_width)
{
	std::vector<Eigen::Vector4f> rg = lcms.getRegion(seed[0] - rt_width, seed[0] + rt_width, seed[1] - mz_width, seed[1] + mz_width);
	Eigen::VectorXf rts = lcms.getRT();
	float rt_gap = (lcms.getRT().segment(1, rts.size() - 1) - lcms.getRT().segment(0, rts.size() - 1)).mean();


	std::list<int>  pic_ids;
	bool b_left = true, b_right = true;
	pic_ids.push_back(find_idx(rg, seed[0], seed[1], std::numeric_limits<float>::max()));

	float threshold = std::numeric_limits<float>::max();
	for (int i = 0; i< int(rt_width); i++)
	{
		if (pic_ids.size()==5)
		{
			Eigen::VectorXf mzs(pic_ids.size());
			std::transform(pic_ids.begin(), pic_ids.end(), mzs.data(), [&rg](const int & i) {
				return rg[i][1];
			});
			threshold =10 * sqrt((mzs.array() - mzs.mean()).pow(2).sum()/mzs.size());
		}

		if (b_left)
		{
			float rt_left = rg[pic_ids.front()][0] - rt_gap;
			int idx_left = find_idx(rg, rt_left, seed[1], threshold);
			if (idx_left != -1)
			{
				pic_ids.push_front(idx_left);
			}
			else
			{
				b_left = false;
			}
		}

		if (b_right)
		{
			float rt_right = rg[pic_ids.back()][0] + rt_gap;
			int idx_right = find_idx(rg, rt_right, seed[1], threshold);
			if (idx_right != -1)
			{
				pic_ids.push_back(idx_right);
			}
			else
			{
				b_right = false;
			}
		}
		if (!b_left && !b_right)
		{
			break;
		}
	}

	Eigen::MatrixXf ret(pic_ids.size(), 4);


	std::accumulate(pic_ids.begin(), pic_ids.end(), 0, [&ret, &rg](int idx, const int & i) -> int {
		ret.row(idx) = rg[i];
		return idx + 1;
	});
	return ret;
}

#include "parallel_stable_sort.h"
Eigen::MatrixXf sort_by_col(const Eigen::MatrixXf & target, int col)
{
	std::vector<int> ids(target.rows());
	std::iota(ids.begin(), ids.end(), 0);

	pss::parallel_stable_sort(
		ids.begin(),
		ids.end(),
		[&col, &target](const int & i1, const int & i2)->bool
	{
		return target(i1,col) > target(i2,col);
	}
	);

	std::vector<int> rids(target.rows());
	std::iota(rids.begin(), rids.end(), 0);

	pss::parallel_stable_sort(rids.begin(), rids.end(), [&ids](const int & i1, const int & i2) {
		return ids[i1] < ids[i2];
	});

	Eigen::MatrixXf sorted;
	sorted.resize(target.rows(), target.cols() + 2);
	for (unsigned int i = 0; i < target.rows(); i++)
	{
		sorted.row(i).segment(0, target.cols()) = target.row(ids[i]);
		sorted.row(i)[target.cols()] = ids[i];
		sorted.row(i)[target.cols() + 1] = rids[i];
	}

	return sorted;
}

#include <omp.h>
std::vector<Eigen::MatrixXf> FPICs(LCMS & lcms, float min_peak, float rt_width, float mz_width)
{

	Eigen::MatrixXf rmv = lcms.getAll();
	Eigen::MatrixXf rmv_sort = sort_by_col(rmv, 2);
	int idx = std::upper_bound(rmv_sort.col(2).data(), rmv_sort.col(2).data() + rmv_sort.rows(), min_peak,
		[](const float & a, const float & b)->bool
		{
			return a > b;
		}) - rmv_sort.col(2).data();
	Eigen::VectorXi b_inc(rmv_sort.rows()); b_inc.fill(0);

	std::vector<Eigen::MatrixXf> pics;


	int total = 0;
	while (std::any_of(b_inc.data(), b_inc.data() + idx, [&b_inc](const int & i) {
		return i == 0;
	}))
	{
		std::vector<Eigen::Vector3f> seeds = pic_seeds(rmv_sort, idx, b_inc, mz_width);

		
		#pragma omp parallel for 
		for (int i=0; i< seeds.size(); i++)
		{
			Eigen::MatrixXf pic = FPIC(lcms, seeds[i], rt_width, mz_width);
			#pragma omp critical (result)
			{
				pics.push_back(pic);
				std::for_each(pic.col(3).data(), pic.col(3).data() + pic.rows(), [&b_inc, &rmv_sort](const float & v) {
					b_inc[(int)rmv_sort(int(v), 4)] = 1;
				});
			}
			total = total + 1;
		}
	}
	return pics;
}