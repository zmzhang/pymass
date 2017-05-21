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


#include <Eigen/Core>
#include <algorithm>
#include <vector>

template <typename Scalar, int rows, int cols, int options, int maxRows, int maxCols>
Eigen::Matrix<Scalar, rows, cols, options, maxRows, maxCols> sort_by_col(
	Eigen::Matrix<Scalar, rows, cols, options, maxRows, maxCols> target, int col
) {
	std::vector<Eigen::Matrix<Scalar, 1, cols>> matrixRows;
	for (unsigned int i = 0; i < target.rows(); i++)
		matrixRows.push_back(target.row(i));
	std::sort(
		matrixRows.begin(),
		matrixRows.end(),
		[&col](const Eigen::Matrix<Scalar, 1, cols> & a, const Eigen::Matrix<Scalar, 1, cols> & b)->bool
	{
		return a(0, col) > b(0, col);
	}
	);

	Eigen::Matrix<Scalar, rows, cols, options, maxRows, maxCols> sorted;
	sorted.resize(target.rows(), target.cols());
	for (unsigned int i = 0; i < matrixRows.size(); i++)
		sorted.row(i) = matrixRows[i];
	return sorted;
}

typedef std::function<bool(const Eigen::VectorXf &lhs, const Eigen::VectorXf &rhs)> mz_comp;
std::vector<Eigen::Vector3f> PYMASS_EXPORT pic_seeds(const Eigen::MatrixXf & m, const int & idx, const Eigen::VectorXi & b_inc, float mz_tol);
Eigen::MatrixXf PYMASS_EXPORT FPIC(LCMS & lcms, const Eigen::Vector3f & seed, float rt_width, float mz_width);

