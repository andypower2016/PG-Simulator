#ifndef PG_SOLVER_H
#define PG_SOLVER_H

#include "PowerGridModel.hpp"
#include "PowerGrid.hpp"
#include "SparseMatrix.hpp"
#include "mosek.h"		
#include "engine.h"	

void MSKAPI printstr(void *handle, const char str[]);

namespace PowerGridSolver {
	
class PGSolver : public PowerGridPlane::PowerGrid
{

public : 

	PGSolver();

	void SolvePowerGridPlane(PG_PLANE & Plane);

	void SolvePowerGridVdrop(PG_PLANE &Plane);

	void SolveThermalPlane(THERMAL_PLANE & thermal_plane);
	
	void Cal_Total_Current_Power(); 
	
	void OnePointThermalAnalysis
	(unsigned int NodeIndex, double *CurrentSolution, 
		EigenDenseMatrix & ThermalSystemMatrix, int type);

	void PGVerificationQP();

	void PGVerificationQP_Efficient();  

	void PGVerificationQP_Efficient_Simulation();

	void PGVerificationQP_Efficient_Simulation_Verify();

	void PGVerificationQP_Efficient_Simulation_Verify_Iterative();

	void PGVerificationExpandPoint_Efficient();	  

	void PGVerificationExpandPoint_Efficient_Simulation();

	template<typename RetType>
	Eigen::SparseMatrix<RetType> 
	SimulationForConductanceMatrix
	(const Eigen::MatrixXd & AkMatrix, 
 	 const double *xx,	
 	 unsigned int thermalgridsize,
 	 const STAMPLIST & stampList,
 	 EigenDenseVector & CurrentSourceVector);		// returns the stamped conductance matrix at Texp induced by current vector

	
	void PGVerificationExpandPoint_Efficient_Method2(int Method);	

	void PGVerificationExpandPoint_EfficientExpand();	
	
	void PGVerificationLP();

	void PGVerificationLP_Efficient();

	void MonteCarlo();

	// random double function M~N
	auto randMToN(double M, double N) 
		-> decltype( M + (rand() / ( RAND_MAX / (N-M) ) ) )
	{
    	return M + (rand() / ( RAND_MAX / (N-M) ) );  
	}

	double mConstructDataTime;


	// Added
	// added 2018.8.18 ~
	void PGVerificationExpandPoint_Efficient_Matlab();
	void PGVerificationExpandPoint_Efficient_Simulation_Matlab();
	
	
};

}
#endif
