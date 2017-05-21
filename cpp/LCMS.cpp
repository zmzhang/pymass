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
	tic();
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
		ret[i] = v;
		i++; });
	toc();
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

	if (abs((*lower)[col] - v) < abs((*(--lower))[col] - v))
	{
		return lower - rg.begin() + 1;
	}
	else
	{
		return lower - rg.begin();
	}
}

inline int find_idx(const std::vector<Eigen::Vector4f>& rg, const float & rt, const float & mz, float threshold)
{
	return -1;
}

Eigen::MatrixXf FPIC(LCMS & lcms, const Eigen::Vector3f & seed, float rt_width, float mz_width)
{
	Eigen::MatrixXf ret;
	std::vector<Eigen::Vector4f> rg = lcms.getRegion(seed[0] - rt_width, seed[0] + rt_width, seed[1] - mz_width, seed[1] + mz_width);
	Eigen::VectorXf rts = lcms.getRT();
	float rt_gap = (lcms.getRT().segment(1, rts.size() - 1) - lcms.getRT().segment(0, rts.size() - 1)).mean();
	

	return ret;
}
