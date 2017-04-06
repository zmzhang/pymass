#pragma once

#include <iostream>
#include <map>
#include <functional>
#include <vector>
#include <memory>
#include <expat.h>
#include <Eigen/Core>
#include "pymass_export.h"
#include "LCMS.h"

class PYMASS_EXPORT mzXMLParser {

public:
	typedef std::function< void(mzXMLParser& a) > Handler;
	void InitHandlers();
	void initParser();
	mzXMLParser();
	~mzXMLParser();
	long line();
	long column();
	void parseString(std::string s);
	LCMS parseFile(const std::string& filename);
	void AddStartElementHandler(std::string elementName, Handler call);
	void AddEndElementHandler(std::string elementName, Handler call);
	void AddValueHandler(std::string elementName, Handler call);
	std::string getCurrentElement();
	void setCurrentElement(std::string str);
	void setCurrentText(std::string str);
	std::map<std::string, std::string> attributes();
	const std::string& getCurrentText();

	LCMS m_LCMS;
	std::vector<std::map<std::string, std::string> > m_scanAttributes;


private:
	XML_Parser parser;
	static void startElement(void *data, const char *name, const char **atts);
	static void endElement(void *data, const char *name);
	static void characterDataHandler(void *data, char const *d, int len);
	std::map<std::string, Handler> startElementHandlers;
	std::map<std::string, Handler> endElementHandlers;
	std::map<std::string, Handler> valueHandlers;
	std::map<std::string, std::string> currentAttributes_;
	std::string currentElement_;
	std::string currentText_;
};



