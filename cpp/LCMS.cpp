#pragma once

#include "base64.h"
#include "utils.h"
#include "LCMS.h"

using namespace std;


Eigen::VectorXf LCMS::getBIC()
{
	return Eigen::VectorXf::Map(m_vecBIC.data(), m_vecBIC.size());
}

Eigen::VectorXf LCMS::getRT()
{
	return Eigen::VectorXf::Map(m_vecRT.data(), m_vecRT.size());
}

Eigen::VectorXf LCMS::getTIC()
{
	return Eigen::VectorXf::Map(m_vecTIC.data(), m_vecTIC.size());
}


Eigen::VectorXf LCMS::getMS(int i, int level)
{
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
	Eigen::MatrixXf ret;
	return ret;
}