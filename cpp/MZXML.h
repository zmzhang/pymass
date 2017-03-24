#pragma once

#include <iostream>
#include <map>
#include <functional>
#include <vector>
#include <memory>
#include <expat.h>
#include <Eigen/Core>
#include "pymass_export.h"

using namespace std;



struct PYMASS_EXPORT MassScan
{
public:
	Eigen::VectorXd mz;
	Eigen::VectorXd val;
	double precursor_mz;
	double RT;
	double BIC;
	double TIC;
	vector<shared_ptr<MassScan> > childs;
};


class PYMASS_EXPORT LCMS {

public:
	LCMS(){}
	~LCMS(){}
	void push_back(const MassScan& scan) { m_massScans.push_back(scan); }
	vector<MassScan> m_massScans;

};

class PYMASS_EXPORT MZXML {

public:
	typedef function< void(MZXML& a) > Handler;
	void InitHandlers();
	void initParser();
	MZXML();
	~MZXML();
	long line();
	long column();
	void parseString(string s);
	void parseFile(const std::string& filename);
	void AddStartElementHandler(string elementName, Handler call);
	void AddEndElementHandler(string elementName, Handler call);
	void AddValueHandler(string elementName, Handler call);
	string getCurrentElement();
	void setCurrentElement(string str);
	void setCurrentText(string str);
	map<string, string> attributes();
	const string& getCurrentText();


	LCMS m_LCMS;
	vector<map<string, string> > m_scanAttributes;

	
	vector<double> m_vecBIC;
	vector<double> m_vecRT;
	vector<double> m_vecTIC;

	Eigen::VectorXd getBIC();
	Eigen::VectorXd getRT();
	Eigen::VectorXd getTIC();
	Eigen::VectorXd getMS(int i, int level=1);
	Eigen::VectorXd getVal(int i, int level=1);

private:
	XML_Parser parser;
	static void startElement(void *data, const char *name, const char **atts);
	static void endElement(void *data, const char *name);
	static void characterDataHandler(void *data, char const *d, int len);
	map<string, Handler> startElementHandlers;
	map<string, Handler> endElementHandlers;
	map<string, Handler> valueHandlers;
	map<string, string> currentAttributes_;
	string currentElement_;
	string currentText_;
};



