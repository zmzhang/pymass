#include <iostream>
#include <fstream>
#include "mzXMLParser.h"
#include "utils.h"
#include "cuda_sort.h"

using namespace std;
using namespace Eigen;

void testMZXML() {

	for (int i = 0; i < 1; i++)
	{
		mzXMLParser e;
		tic();
		//LCMS lcms = e.parseFile("D:/workspace/pymass/python/标2-方法5-正负离子_Seg1Ev1.mzXML");
		LCMS lcms = e.parseFile("../../python/mixture_bsa300fmol_n3.mzXML");
		//LCMS lcms = e.parseFile("../../python/detnoise_sigma3_PICKED.mzXML");
		toc();
		//cout << lcms.getTIC().transpose() << endl;
		//cout << lcms.getRT().transpose() << endl;
		cout << lcms.getRegion(1000,1100,500,501) << endl;
		Eigen::MatrixXf rmv = lcms.getAll();

		processLCMS(lcms);

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

	VectorXf t(4);
	t << 1.1f, 6.5f, 9.5f, 3.0f;
	VectorXi ids(4);

	thrustSort(t.data(),ids.data(), 4);

	cout << t <<endl;
	cout << ids;


    return 0;
}
