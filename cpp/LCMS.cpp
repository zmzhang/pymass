#pragma once

#include <list>
#include "base64.h"
#include "utils.h"
#include "LCMS.h"

using namespace std;


void LCMS::update()
{
	if(m_vecBIC.size()!=m_massScans.size() || m_vecRT.size()!=m_massScans.size() || m_vecTIC.size()!= m_massScans.size())
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

Eigen::MatrixXf LCMS::getRegion(float rt_begin, float rt_end, float mz_begin, float mz_end)
{
	update();
	Eigen::MatrixXf ret;
	std::list<Eigen::VectorXi> mzi_vec;
	Eigen::VectorXf rts(2); rts << rt_begin, rt_end;
	Eigen::VectorXf mzs(2); mzs << mz_begin, mz_end;
	Eigen::VectorXi rti = findclosest(m_vecRT, rts);
	for (int i = rti[0]; i<rti[1]; i++)
	{
		Eigen::VectorXi mzi = findclosest(m_massScans[i].mz, mzs);
		mzi_vec.push_back(mzi);
	}

	//todo:
	//	std::accumulate(mzi_vec)
	//	resize ret matrix
	//	slice(m_massScans[i].mz, mzi);
	//	slice(m_massScans[i].val, mzi);

	return ret;
}