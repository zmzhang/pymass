#pragma once

#include <list>
#include "base64.h"
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

std::vector<Eigen::Vector3f> LCMS::getRegion(float rt_begin, float rt_end, float mz_begin, float mz_end)
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

	std::vector<Eigen::Vector3f> ret;
	ret.resize(rows);
	int s = 0;
	for (int i = 0; i<mzi_vec.size(); i++)
	{
		for (int j = mzi_vec[i][0]; j<=mzi_vec[i][1]; j++)
		{
			ret[s][0] = m_vecRT[i + rti[0]];
			ret[s][1] = m_massScans[i + rti[0]].mz[j];
			ret[s][2] = m_massScans[i + rti[0]].val[j];
			s = s + 1;
		}
	}

	return ret;
}