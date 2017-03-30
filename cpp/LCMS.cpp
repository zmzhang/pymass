#pragma once

#include "base64.h"
#include "utils.h"
#include "LCMS.h"

using namespace std;


Eigen::VectorXd LCMS::getBIC()
{
	return Eigen::VectorXd::Map(m_vecBIC.data(), m_vecBIC.size());
}

Eigen::VectorXd LCMS::getRT()
{
	return Eigen::VectorXd::Map(m_vecRT.data(), m_vecRT.size());
}

Eigen::VectorXd LCMS::getTIC()
{
	return Eigen::VectorXd::Map(m_vecTIC.data(), m_vecTIC.size());
}


Eigen::VectorXd LCMS::getMS(int i, int level)
{
	if (level == 1)
	{
		return m_massScans[i].mz;
	}
	else if (level == 2)
	{
		return m_massScans[i].childs[0]->mz;
	}
	return Eigen::VectorXd(0, 0);
}



Eigen::VectorXd LCMS::getVal(int i, int level)
{
	if (level == 1)
	{
		return m_massScans[i].val;
	}
	else if (level == 2)
	{
		return m_massScans[i].childs[0]->val;
	}
	return Eigen::VectorXd(0, 0);
}

Eigen::MatrixXd LCMS::getRegion(double rt_begin, double rt_end, double mz_begin, double mz_end)
{
	Eigen::MatrixXd ret;
	return ret;
}