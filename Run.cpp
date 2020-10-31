#include "Run.hpp"

Run::Run() 
{}

Run::~Run()
{
	DestructPowerGridDatas();
}

void Run::Run_DC_Simulation()
{
	clock_t start, end;
	start = clock();
	SolvePowerGridPlane(mTotalNodes);
	SolvePowerGridVdrop(mTotalNodes);
	end = clock();
	std::cout << "Solving elapsed time = " << (double)(end-start)/CLOCKS_PER_SEC << std::endl;
}


void Run::DestructPowerGridDatas()
{
	for(auto& elements : mTotalNodes)
	{
		delete elements.second;
	}
	mTotalNodes.clear();
	
	for(auto& elements : mtotal_supply_Vsource)
	{
		delete elements.first;
	}
	mtotal_supply_Vsource.clear();
	
	for(auto& res : mtotal_unshorted_Res)
	{
		delete res;
	}
	mtotal_unshorted_Res.clear();
	
	for(auto& curr : mtotal_Isource)
	{
		delete curr;
	}
	mtotal_Isource.clear();
}


void Run::ClearVerificationDatas()
{
	mPG_OriginConductanceMatrix.clear();
	mStampNodeList.clear();
	mPG_ThermalMatrix.clear();	
	mEachRowPGInverseMatrix.clear();	
	mPG_EigenConductanceSubMatrix.clear();
	mLocalCurrentConstraint.clear();
	mGlobalCurrentConstraint.clear();
	mTotalThermalGrids.clear();
	mCurrentSourceIdx.clear();
	mCurrentSourceIdxMapping.clear();	
}


void Run::Run_DC_Verification(int mode)
{
	clock_t start;
	switch(mode)
	{
		case 1 : 	// Direct QP  (Don't run this, use Efficient QP !) 
			start = clock();
			ConstructThermalAwareVerificationData_QP(0);
			mConstructDataTime =  (double)(clock()-start) / (CLOCKS_PER_SEC) ;
			PGVerificationQP();
			break;
			
		case 2 : 	// Efficient QP without IR drop recalculation (in thesis)
			start = clock();
			ConstructThermalAwareVerificationData_QP(1);
			mConstructDataTime =  (double)(clock()-start) / (CLOCKS_PER_SEC) ;
			PGVerificationQP_Efficient();
			break;

		case 3 : 	// Direct Expansion Point Method without IR drop recalculation (not in thesis) 
			ConstructThermalAwareVerificationData_ExpandPoint(1);
			PGVerificationExpandPoint_EfficientExpand();
			break;

		case 4 : 	// Expansion Point Method without IR drop recalculation
			start = clock();
			ConstructThermalAwareVerificationData_ExpandPoint(1);
			mConstructDataTime =  (double)(clock()-start) / (CLOCKS_PER_SEC) ;
			PGVerificationExpandPoint_Efficient();
			break;

		case 5 : 	// Thermal-less LP
			start = clock();
			ConstructThermalLessVerificationData(0);
			mConstructDataTime =  (double)(clock()-start) / (CLOCKS_PER_SEC) ;
			PGVerificationLP();
			break;

		case 6 : 	// Efficient Thermal-less LP 
			start = clock();
			ConstructThermalLessVerificationData(1);
			mConstructDataTime =  (double)(clock()-start) / (CLOCKS_PER_SEC) ;
			PGVerificationLP_Efficient();
			break;

		case 7 : 	// Efficient Expand Point (Method2) Direct (not in thesis)
			start = clock();
			ConstructThermalAwareVerificationData_ExpandPoint_Method2(1);
			mConstructDataTime =  (double)(clock()-start) / (CLOCKS_PER_SEC) ;
			PGVerificationExpandPoint_Efficient_Method2(1);
			break;	

		case 8 : 	// Efficient Expand Point (Method2) Iterative (not in thesis)
			start = clock();
			ConstructThermalAwareVerificationData_ExpandPoint_Method2(1);
			mConstructDataTime =  (double)(clock()-start) / (CLOCKS_PER_SEC) ;
			PGVerificationExpandPoint_Efficient_Method2(2);
			break;	

		case 9 :	// Efficient QP with IR drop recalculation (Combined with Efficient QP without IR drop recalculation) (in thesis)
			start = clock();
			ConstructThermalAwareVerificationData_QP(1);
			mConstructDataTime =  (double)(clock()-start) / (CLOCKS_PER_SEC) ;
			PGVerificationQP_Efficient_Simulation();
			break;

		case 10 :	// Direct Expansion Point Method with IR drop recalculation (in thesis)
			start = clock();
			ConstructThermalAwareVerificationData_ExpandPoint_Method2(1);
			mConstructDataTime =  (double)(clock()-start) / (CLOCKS_PER_SEC) ;
			PGVerificationExpandPoint_Efficient_Method2(3);
			break;	

		case 11 : 	// Expansion Point Method with IR drop recalculation (in thesis)
			start = clock();
			ConstructThermalAwareVerificationData_ExpandPoint(1);
			mConstructDataTime =  (double)(clock()-start) / (CLOCKS_PER_SEC) ;
			PGVerificationExpandPoint_Efficient_Simulation();
			break;

		case 12 :	// 做兩次QP,第二次QP會在第一次QP找到的solution current pattern所造成的溫度之下展開 (not in thesis)
			start = clock();
			ConstructThermalAwareVerificationData_QP(1);
			mConstructDataTime =  (double)(clock()-start) / (CLOCKS_PER_SEC) ;
			PGVerificationQP_Efficient_Simulation_Verify();
			break;
		case 13 :  // Same approach as case 12 but using iterative method (not in thesis)
			start = clock();
			ConstructThermalAwareVerificationData_QP(1);
			mConstructDataTime =  (double)(clock()-start) / (CLOCKS_PER_SEC) ;
			PGVerificationQP_Efficient_Simulation_Verify_Iterative();
			break;	


		case 14 : 	// Expansion Point Method with IR drop recalculation using Matlab
			start = clock();
			ConstructThermalAwareVerificationData_ExpandPoint_Matlab(1);
			mConstructDataTime =  (double)(clock()-start) / (CLOCKS_PER_SEC) ;
			PGVerificationExpandPoint_Efficient_Simulation_Matlab();
			break;	
		case 15 : 	// Expansion Point Method without IR drop recalculation using Matlab
			start = clock();
			ConstructThermalAwareVerificationData_ExpandPoint_Matlab(1);
			mConstructDataTime =  (double)(clock()-start) / (CLOCKS_PER_SEC) ;
			PGVerificationExpandPoint_Efficient_Matlab();
			break;		

		/*
		case 14 :   // MonteCarlo (uncomplete)
			start = clock();
			ConstructMonteCarloData();
			mConstructDataTime =  (double)(clock()-start) / (CLOCKS_PER_SEC) ;
			MonteCarlo();
			break;
		*/	
		default : 
			break;
	}
	ClearVerificationDatas();
}

