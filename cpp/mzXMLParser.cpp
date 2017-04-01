#pragma once

//#include "base64.h"
#include "utils.h"
#include "mzXMLParser.h"

#include <libbase64.h>

using namespace std;

void mzXMLParser::InitHandlers() {
	AddStartElementHandler("scan", [](mzXMLParser& a) -> void {
		map<string, string> atts = a.attributes();
		a.m_scanAttributes.push_back(atts);

	});

	AddEndElementHandler("scan", [](mzXMLParser& a) -> void {
		a.m_scanAttributes.pop_back();
	});

	AddStartElementHandler("peaks", [](mzXMLParser& a) -> void {
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
	
	AddEndElementHandler("peaks", [](mzXMLParser& a) -> void {

		if (a.m_scanAttributes.size() == 1)
		{
			string str = a.getCurrentText();
			char * raw = new char[str.size()];
			size_t outlen;
			base64_decode(str.data(), str.size(), raw, &outlen, BASE64_FORCE_AVX2);
			size_t nNum = outlen / (2 * sizeof(float));
			float* floatArray = reinterpret_cast<float*>(raw);
			MassScan& scan = a.m_LCMS.m_massScans.back();
			scan.mz.resize(nNum);
			scan.val.resize(nNum);
			scan.precursor_mz = -1;

			map<string, string> atts = a.m_scanAttributes.back();
			size_t nLen = atts["retentionTime"].length();
			scan.RT = stod(atts["retentionTime"].substr(2, nLen - 3));
			if (atts.find("basePeakIntensity")!=atts.end())
			{
				scan.BIC = stod(atts["basePeakIntensity"]);
			}			
			for (int i = 0; i < nNum; i++)
			{
				scan.mz[i] = ReverseFloat(floatArray[2 * i]);
				scan.val[i] = ReverseFloat(floatArray[2 * i + 1]);
			}
			scan.TIC = scan.val.sum();
			delete[] raw;
		}


		if (a.m_scanAttributes.size() == 2)
		{

			string str = a.getCurrentText();
			char * raw = new char[str.size()];
			size_t outlen;
			base64_decode(str.data(), str.size(), raw, &outlen, BASE64_FORCE_AVX2);
			size_t nNum = outlen / (2 * sizeof(float));
			float* floatArray = reinterpret_cast<float*>(raw);

			shared_ptr<MassScan> pscan = a.m_LCMS.m_massScans.back().childs.back();
			pscan->mz.resize(nNum);
			pscan->val.resize(nNum);


			map<string, string> atts = a.m_scanAttributes.back();
			size_t nLen = atts["retentionTime"].length();
			pscan->RT = stod(atts["retentionTime"].substr(2, nLen - 3));
			if (atts.find("basePeakIntensity") != atts.end())
			{
				pscan->BIC = stod(atts["basePeakIntensity"]);
			}
			pscan->precursor_mz = stod(atts["basePeakMz"]);;
			for (int i = 0; i < nNum; i++)
			{
				pscan->mz[i] = ReverseFloat(floatArray[2 * i]);
				pscan->val[i] = ReverseFloat(floatArray[2 * i + 1]);
			}
			pscan->TIC = pscan->val.sum();
			delete[] raw;
		}

		a.setCurrentText("");
	});

	

}


void mzXMLParser::initParser()
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
	XML_SetElementHandler(parser, &mzXMLParser::startElement, &mzXMLParser::endElement);

	XML_SetCharacterDataHandler(parser, &mzXMLParser::characterDataHandler);

	InitHandlers();
}


mzXMLParser::mzXMLParser() {
	parser = NULL;
}

mzXMLParser::~mzXMLParser() {
	if (parser)
	{
		XML_ParserFree(parser);
	}
}

long mzXMLParser::line() {
	return XML_GetCurrentLineNumber(parser);
}

long mzXMLParser::column() {
	return XML_GetCurrentColumnNumber(parser);
}

void mzXMLParser::parseString(string s) {
	int len = static_cast<int>(strlen(s.c_str()));
	XML_Parse(parser, s.c_str(), len, 1);
}

#include <boost/progress.hpp>

LCMS mzXMLParser::parseFile(const std::string& filename) {
	cout << "Parsing " + filename << endl;
	initParser();
	m_LCMS = LCMS();
	try
	{
		FILE* fd = fopen(filename.c_str(), "r");

		if (fd == NULL) {
			throw std::runtime_error("file doesn't exist");
		}

		fseek(fd, 0L, SEEK_END);
		long sz = ftell(fd); rewind(fd);

		const size_t BUFFER_SIZE = std::min(long(50*1024*1024), sz/50);
		boost::progress_display show_progress(sz);
		for (;;) {
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
			else
			{
				show_progress += bytes_read;
			}
		}
		fclose(fd);

		for (auto &i : m_LCMS.m_massScans) {
			m_LCMS.m_vecBIC.push_back(i.BIC);
			m_LCMS.m_vecTIC.push_back(i.TIC);
			m_LCMS.m_vecRT.push_back(i.RT);
		}
		cout << endl;

	}
	catch (std::exception const& e)
	{
		std::cout << "Exception: " << e.what() << "\n";
	}
	return m_LCMS;
}



void mzXMLParser::AddStartElementHandler(string elementName, Handler call) {
	startElementHandlers[elementName] = call;
}

void mzXMLParser::AddEndElementHandler(string elementName, Handler call) {
	endElementHandlers[elementName] = call;
}

void mzXMLParser::AddValueHandler(string elementName, Handler call) {
	valueHandlers[elementName] = call;
}


string  mzXMLParser::getCurrentElement() {
	return currentElement_;
}

void mzXMLParser::setCurrentElement(string str)
{
	currentElement_ = str;
}

void mzXMLParser::setCurrentText(string str)
{
	currentText_ = str;
}

map<string, string> mzXMLParser::attributes() {
	return currentAttributes_;
}

const string& mzXMLParser::getCurrentText() {
	return currentText_;
}

void mzXMLParser::startElement(void *data, const char *name, const char **atts) {

	mzXMLParser* e = static_cast<mzXMLParser *>(data);

	e->currentElement_ = name;

	if (e->startElementHandlers.find(name) != e->startElementHandlers.end()) {

		for (int i = 0; atts[i]; i += 2) {
			e->currentAttributes_[atts[i]] = atts[i + 1];
		}
		e->startElementHandlers.find(name)->second(*e);
	}
}

void mzXMLParser::endElement(void *data, const char *name) {

	mzXMLParser* e = static_cast<mzXMLParser *>(data);
	if (e->endElementHandlers.find(name) != e->endElementHandlers.end()) {
		e->endElementHandlers.find(name)->second(*e);
	}
}

void mzXMLParser::characterDataHandler(void *data, char const *d, int len) {
	

	mzXMLParser* e = static_cast<mzXMLParser *>(data);
	e->currentText_ += std::string(d, len);
	if (e->valueHandlers.find(e->currentElement_) != e->valueHandlers.end()) {
		e->valueHandlers.find(e->currentElement_)->second(*e);
	}
}
      