#pragma once

#include <iostream>
#include <map>
#include <functional>
#include <vector>
#include <memory>
#include <string>
#include <numeric>
#include <expat.h>
#include <Eigen/Core>
#include "pymass_export.h"

struct PYMASS_EXPORT MassScan
{
public:
	Eigen::VectorXf mz;
	Eigen::VectorXf val;
	Eigen::VectorXf id;
	float precursor_mz;
	float RT;
	float BIC;
	float TIC;
	std::vector<std::shared_ptr<MassScan> > childs;
};


class PYMASS_EXPORT LCMS {

public:
	LCMS(){}
	~LCMS(){}
	void push_back(const MassScan& scan) { m_massScans.push_back(scan); }
	std::vector<MassScan> m_massScans;

	Eigen::VectorXf m_vecBIC;
	Eigen::VectorXf m_vecRT;
	Eigen::VectorXf m_vecTIC;

	Eigen::MatrixXf getAll();
	Eigen::VectorXf getBIC();
	Eigen::VectorXf getRT();
	Eigen::VectorXf getTIC();
	Eigen::VectorXf getMS(int i, int level = 1);
	Eigen::VectorXf getVal(int i, int level = 1);
	std::vector<Eigen::Vector4f> getRegion(float rt_begin, float rt_end, float mz_begin, float mz_end);


private:
	void update();

};


typedef std::function<bool(const Eigen::VectorXf &lhs, const Eigen::VectorXf &rhs)> mz_comp;
std::vector<Eigen::Vector3f> PYMASS_EXPORT pic_seeds(const Eigen::MatrixXf & m, const int & idx, const Eigen::VectorXi & b_inc, float mz_tol);
Eigen::MatrixXf PYMASS_EXPORT FPIC(LCMS & lcms, const Eigen::Vector3f & seed, float rt_width, float mz_width);
Eigen::MatrixXf PYMASS_EXPORT sort_by_col(const Eigen::MatrixXf & target, int col);

