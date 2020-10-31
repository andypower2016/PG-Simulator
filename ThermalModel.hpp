#ifndef ThermalModel_H
#define ThermalModel_H

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <cmath>
#include <utility> 
#include <vector>
#include <map>

namespace ThermalGridModel {


class ThermalGrid;

using THERMAL_PLANE = std::vector<std::vector<std::vector<ThermalGridModel::ThermalGrid>>> ;

class ThermalGrid {

public : 
	
	ThermalGrid() {}
	
	ThermalGrid
	(double dx, double dy, double dz, double ux, double uy, double uz, 
		double k, double heat_trans_coeff);

	double getXmidpoint() { return fabs( dx_loc_idx + ux_loc_idx ) / 2 ; }

	double getYmidpoint() { return fabs( dy_loc_idx + uy_loc_idx ) / 2 ; }

	double getZmidpoint() { return fabs( dz_loc_idx + uz_loc_idx ) / 2 ; }


	double dx_loc_idx, dy_loc_idx, dz_loc_idx,
		   ux_loc_idx, uy_loc_idx, uz_loc_idx,
		   x_res, y_res, z_res,
		   z_bound_res;
	
	double mTemperature;
	
	double mPower;
	
	unsigned int mGridIndex;	

};


class PowerMap
{

public : 
	
	PowerMap() : mTotalPower(0) {}
	
	double **m2DPowerMap;

	double mTotalPower;
};


class ThermalPlane : public PowerMap
{
	
public : 
	
	ThermalPlane();

	THERMAL_PLANE & getTotalThermalGrids() { return mTotalThermalGrids; }

	THERMAL_PLANE mTotalThermalGrids;
  	
	std::vector<double> mPowerVector;

	std::vector<double> mThermalSolution;
 
    double mXCut, mYCut, mZlayer;

    double mZ_metalCut, mZ_subCut;

};	
	
	
}

#endif