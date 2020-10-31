#ifndef Parser_H
#define Parser_H

#include <iostream>
#include <fstream>
#include <cstring>
#include <string> 
#include <vector>
#include <map>
#include <cctype>

#include "PGSolver.hpp"

using namespace PowerGridSolver;

class Parser : public PGSolver
{
	
public :

	Parser();
	
	Parser(const std::string& file_path);
	
	void ParseInput(const std::string &file_path, const std::string &cktname);
	
protected : 	
	
	void ParseLine(const std::string& line, const std::string& seperators, std::vector<std::string>& line_element);
	
	PowerGridModel::Resistor *res;
	
	PowerGridModel::VoltageSource *volt;
	
	PowerGridModel::CurrentSource *curr;
	
	PowerGridModel::Node *node;
	
};

#endif 
