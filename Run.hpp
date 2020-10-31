#ifndef Run_H
#define Run_H

#include <iostream>
#include <string>
#include <ctime>
#include <memory>

#include "Parser.hpp"

using namespace PowerGridSolver;

class Run : public Parser
{

public : 

	Run();

	~Run();
	
	void Run_DC_Simulation();

	void Run_DC_Verification(int mode);

	void ClearVerificationDatas();
	
	void DestructPowerGridDatas();

};


#endif
