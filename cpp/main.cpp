#include <iostream>
#include <fstream>
#include "mzXMLParser.h"
#include "utils.h"

using namespace std;
using namespace Eigen;

void testMZXML() {

	for (int i = 0; i < 1; i++)
	{
		mzXMLParser e;
		tic();
		//LCMS lcms = e.parseFile("D:/workspace/pymass/python/标2-方法5-正负离子_Seg1Ev1.mzXML");
		//LCMS lcms = e.parseFile("../../python/mixture_bsa300fmol_n3.mzXML");
		//LCMS lcms = e.parseFile("../../python/detnoise_sigma3_PICKED.mzXML");
		LCMS lcms = e.parseFile("../../python/MM14_20um.mzxml");

		toc();
		//cout << lcms.getTIC().transpose() << endl;
		//cout << lcms.getRT().transpose() << endl;
		//cout << lcms.getRegion(740,760,301.1,301.15) << endl;
		Eigen::MatrixXf rmv = lcms.getAll();

		tic();
		Eigen::MatrixXf rmv_sort = sort_by_col(rmv, 2);
		toc();

		FPIC(lcms, rmv_sort.row(0), 100, 0.5);

		cout << rmv_sort.topRows(10) << endl;
		cout << rmv_sort.bottomRows(10) << endl;


	}
}

void testFindClosest()
{
	VectorXf v(10);
	v << 1, 2, 3, 4, 5, 5, 6, 7, 8, 9;
	VectorXf t(3);
	t << 1.1f, 6.5f, 9.5f;


	cout << searchsorted(v, t) << endl;
	cout << "##########" << endl;
	cout << findclosest(v, t) << endl;
}


int main(int argc, const char * argv[]) {
	testMZXML();  
	//testFindClosest();

    return 0;
}
