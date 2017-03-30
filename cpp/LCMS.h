#pragma once

#include <iostream>
#include <map>
#include <functional>
#include <vector>
#include <memory>
#include <string>
#include <expat.h>
#include <Eigen/Core>
#include "pymass_export.h"


struct PYMASS_EXPORT MassScan
{
public:
	Eigen::VectorXd mz;
	Eigen::VectorXd val;
	double precursor_mz;
	double RT;
	double BIC;
	double TIC;
	std::vector<std::shared_ptr<MassScan> > childs;
};


class PYMASS_EXPORT LCMS {

public:
	LCMS(){}
	~LCMS(){}
	void push_back(const MassScan& scan) { m_massScans.push_back(scan); }
	std::vector<MassScan> m_massScans;

	std::vector<double> m_vecBIC;
	std::vector<double> m_vecRT;
	std::vector<double> m_vecTIC;

	Eigen::VectorXd getBIC();
	Eigen::VectorXd getRT();
	Eigen::VectorXd getTIC();
	Eigen::VectorXd getMS(int i, int level = 1);
	Eigen::VectorXd getVal(int i, int level = 1);
	Eigen::MatrixXd getRegion(double rt_begin, double rt_end, double mz_begin, double mz_end);

};

