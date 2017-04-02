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
		LCMS lcms = e.parseFile("C:/workspace/pymass/python/mixture_bsa300fmol_n3.mzXML");
		//cout << lcms.getMS(3, 2).transpose() << endl;
		//cout << lcms.getVal(3, 2).transpose() << endl;
		toc();
	}
}

void testFindClosest()
{
	VectorXf v(10);
	v << 1, 2, 3, 4, 5, 5, 6, 7, 8, 9;
	VectorXf t(3);
	t << 1.1, 6.5, 9.5;


	cout << searchsorted(v, t) << endl;
	cout << "##########" << endl;
	cout << findclosest(v, t) << endl;
}


int main(int argc, const char * argv[]) {
    
	testMZXML();  
	testFindClosest();
    return 0;
}
