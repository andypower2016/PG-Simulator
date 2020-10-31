#ifndef PowerGrid_H
#define PowerGrid_H

#include <fstream>
#include "PowerGridModel.hpp"
#include "ThermalModel.hpp"

using namespace PowerGridModel;
using namespace ThermalGridModel;

namespace PowerGridPlane {


class PowerGrid :
		public PowerGridModel::PowerGridData,
		public PowerGridModel::PowerGridSpecs,
		public PowerGridModel::PowerGridGeometry,
		public PowerGridModel::PowerGridConstraints,
		public ThermalGridModel::ThermalPlane
{

public:

		PowerGrid();

		void ConstructPowerGridPlane();
		
		void ConstructThermalPlane();

		void Output_Element_Solution();

		void OutputVdrop();
	
		void Output_DC_Solution();
	
		void OutputThermalSolution();

		void OutuptDistributionMaps();

		void OutputDesignInfo();

		void OutputThermalGrids();

		void OutputPointReport();

		void ConstructCurrentConstraint(int VERIFICATION_MODE, int METHOD_MODE);

		void ConstructThermalGrids();

		void StampPGConductanceMatrix();

		void StampThermalAwareMatrix();

		void StampThermalAwareMatrix_Method2();

		void ComputeInversePGConductanceMatrix();

		void ComputeInverseThermalConductanceMatrix();

		void ConstructThermalLessVerificationData(int METHOD);

		void ConstructThermalAwareVerificationData_QP(int METHOD);

		void ConstructThermalAwareVerificationData_ExpandPoint(int METHOD);

		void ConstructThermalAwareVerificationData_ExpandPoint_Method2(int METHOD);

		void ConstructMonteCarloData();

		// added 2018.8.18 ~
		// Expand Point using matlab
		void ConstructThermalAwareVerificationData_ExpandPoint_Matlab(int METHOD);

};

}
#endif
