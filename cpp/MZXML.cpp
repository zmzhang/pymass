#pragma once

#include <iostream>
#include <map>
#include <functional>
#include <string>
#include "MZXML.h"
#include "base64.h"

using namespace std;

float ReverseFloat(const float inFloat)
{
	float retVal;
	char *floatToConvert = (char*)& inFloat;
	char *returnFloat = (char*)& retVal;

	// swap the bytes into a temporary buffer
	returnFloat[0] = floatToConvert[3];
	returnFloat[1] = floatToConvert[2];
	returnFloat[2] = floatToConvert[1];
	returnFloat[3] = floatToConvert[0];

	return retVal;
}


void MZXML::InitHandlers() {
	AddStartElementHandler("scan", [](MZXML& a) -> void {
		map<string, string> atts = a.attributes();
		a.m_scanAttributes.push_back(atts);

	});

	AddEndElementHandler("scan", [](MZXML& a) -> void {
		a.m_scanAttributes.pop_back();
	});

	AddStartElementHandler("peaks", [](MZXML& a) -> void {
		a.setCurrentText("");
		if (a.m_scanAttributes.size() == 1)
		{
			MassScan scan;
			a.m_LCMS.push_back(scan);
		}
		if (a.m_scanAttributes.size() == 2)
		{
			a.m_LCMS.m_massScans.back().childs.push_back(make_shared<MassScan>());
		}
	});

	AddEndElementHandler("peaks", [](MZXML& a) -> void {

		if (a.m_scanAttributes.size() == 1)
		{
			
			string str = a.getCurrentText();
			vector<BYTE> raw = base64_decode(a.getCurrentText());
			int nNum = raw.size() / (2 * sizeof(float));
			float* floatArray = reinterpret_cast<float*>(raw.data());
			MassScan& scan = a.m_LCMS.m_massScans.back();
			scan.mz.resize(nNum);
			scan.val.resize(nNum);
			scan.precursor_mz = -1;

			map<string, string> atts = a.m_scanAttributes.back();
			size_t nLen = atts["retentionTime"].length();
			scan.RT = stod(atts["retentionTime"].substr(2, nLen - 3));
			scan.BIC = stod(atts["basePeakIntensity"]);
			for (int i = 0; i < nNum; i++)
			{
				scan.mz[i] = ReverseFloat(floatArray[2 * i]);
				scan.val[i] = ReverseFloat(floatArray[2 * i + 1]);
			}
			scan.TIC = scan.val.sum();
		}

		if (a.m_scanAttributes.size() == 2)
		{

			string str = a.getCurrentText();
			vector<BYTE> raw = base64_decode(a.getCurrentText());
			int nNum = raw.size() / (2 * sizeof(float));
			float* floatArray = reinterpret_cast<float*>(raw.data());
			shared_ptr<MassScan> pscan = a.m_LCMS.m_massScans.back().childs.back();
			pscan->mz.resize(nNum);
			pscan->val.resize(nNum);
			

			map<string, string> atts = a.m_scanAttributes.back();
			size_t nLen = atts["retentionTime"].length();
			pscan->RT = stod(atts["retentionTime"].substr(2, nLen - 3));
			pscan->BIC = stod(atts["basePeakIntensity"]);
			pscan->precursor_mz = stod(atts["basePeakMz"]);;
			for (int i = 0; i < nNum; i++)
			{
				pscan->mz[i] = ReverseFloat(floatArray[2 * i]);
				pscan->val[i] = ReverseFloat(floatArray[2 * i + 1]);
			}
			pscan->TIC = pscan->val.sum();
		}

		a.setCurrentText("");
	});



}


void MZXML::initParser()
{
	if (parser)
	{
		XML_ParserFree(parser);
	}
	parser = XML_ParserCreate(NULL);
	if (!parser) {
		throw std::bad_alloc();
	}

	// pass ref from Expat(class) to the expat parser
	XML_SetUserData(parser, this);

	// callbacks for the elements handlers
	XML_SetElementHandler(parser, &MZXML::startElement, &MZXML::endElement);

	XML_SetCharacterDataHandler(parser, &MZXML::characterDataHandler);

	InitHandlers();
}


MZXML::MZXML() {
	parser = NULL;
}

MZXML::~MZXML() {
	if (parser)
	{
		XML_ParserFree(parser);
	}
}

long MZXML::line() {
	return XML_GetCurrentLineNumber(parser);
}

long MZXML::column() {
	return XML_GetCurrentColumnNumber(parser);
}

void MZXML::parseString(string s) {
	int len = static_cast<int>(strlen(s.c_str()));
	XML_Parse(parser, s.c_str(), len, 1);
}

#include <boost/progress.hpp>

void MZXML::parseFile(const std::string& filename) {
	
	initParser();
	try
	{
		FILE* fd = fopen(filename.c_str(), "r");

		fseek(fd, 0L, SEEK_END);
		long sz = ftell(fd); rewind(fd);

		//cout << "filename: " << filename.c_str() << endl;
		//cout << "fd: " << fd << endl;
		if (fd == NULL) {
			throw std::runtime_error("File does not exist");
		}
		const size_t BUFFER_SIZE = 1*1024*1024;
		int i = 0;
		boost::progress_display show_progress(sz);
		for (;;) {
			show_progress+= BUFFER_SIZE;
			void *buffer = XML_GetBuffer(parser, BUFFER_SIZE);
			if (buffer == NULL) {
				throw std::runtime_error("out of memory");
			}

			int bytes_read = static_cast<int>(fread(buffer, 1, BUFFER_SIZE, fd));
			if (bytes_read < 0) {
				throw std::runtime_error("error reading file");
			}

			if (!XML_ParseBuffer(parser, bytes_read, bytes_read == 0)) {
				throw std::runtime_error("could not parse buffer");
			}

			if (bytes_read == 0) {
				break;
			}
			i++;
		}
		fclose(fd);

		for (auto &i : m_LCMS.m_massScans) {
			m_vecBIC.push_back(i.BIC);
			m_vecTIC.push_back(i.TIC);
			m_vecRT.push_back(i.RT);
		}
		cout << endl;

	}
	catch (std::exception const& e)
	{
		std::cout << "Exception: " << e.what() << "\n";
	}
}



void MZXML::AddStartElementHandler(string elementName, Handler call) {
	startElementHandlers[elementName] = call;
}

void MZXML::AddEndElementHandler(string elementName, Handler call) {
	endElementHandlers[elementName] = call;
}

void MZXML::AddValueHandler(string elementName, Handler call) {
	valueHandlers[elementName] = call;
}


string  MZXML::getCurrentElement() {
	return currentElement_;
}

void MZXML::setCurrentElement(string str)
{
	currentElement_ = str;
}

void MZXML::setCurrentText(string str)
{
	currentText_ = str;
}

map<string, string> MZXML::attributes() {
	return currentAttributes_;
}

const string& MZXML::getCurrentText() {
	return currentText_;
}


Eigen::VectorXd MZXML::getBIC()
{
	return Eigen::VectorXd::Map(m_vecBIC.data(),m_vecBIC.size());
}

Eigen::VectorXd MZXML::getRT()
{
	return Eigen::VectorXd::Map(m_vecRT.data(), m_vecRT.size());
}

Eigen::VectorXd MZXML::getTIC()
{
	return Eigen::VectorXd::Map(m_vecTIC.data(), m_vecTIC.size());
}


Eigen::VectorXd MZXML::getMS(int i, int level)
{
	if (level == 1)
	{
		return m_LCMS.m_massScans[i].mz;
	}
	else if (level==2)
	{
		return m_LCMS.m_massScans[i].childs[0]->mz;
	}
	return Eigen::VectorXd(0, 0);
}



Eigen::VectorXd MZXML::getVal(int i, int level)
{
	if (level==1)
	{
		return m_LCMS.m_massScans[i].val;
	}
	else if (level == 2)
	{
		return m_LCMS.m_massScans[i].childs[0]->val;
	}
	return Eigen::VectorXd(0,0);
}

void MZXML::startElement(void *data, const char *name, const char **atts) {

	MZXML* e = static_cast<MZXML *>(data);

	e->currentElement_ = name;

	if (e->startElementHandlers.find(name) != e->startElementHandlers.end()) {

		for (int i = 0; atts[i]; i += 2) {
			e->currentAttributes_[atts[i]] = atts[i + 1];
		}
		e->startElementHandlers.find(name)->second(*e);
	}
}

void MZXML::endElement(void *data, const char *name) {

	MZXML* e = static_cast<MZXML *>(data);
	if (e->endElementHandlers.find(name) != e->endElementHandlers.end()) {
		e->endElementHandlers.find(name)->second(*e);
	}
}

void MZXML::characterDataHandler(void *data, char const *d, int len) {
	

	MZXML* e = static_cast<MZXML *>(data);
	e->currentText_ += std::string(d, len);
	if (e->valueHandlers.find(e->currentElement_) != e->valueHandlers.end()) {
		e->valueHandlers.find(e->currentElement_)->second(*e);
	}
}
      