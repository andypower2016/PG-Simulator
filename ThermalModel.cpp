#include "ThermalModel.hpp"

namespace ThermalGridModel {

	
ThermalGrid::ThermalGrid
(double dx, double dy, double dz, double ux, double uy, double uz, double k, 
	double heat_trans_coeff)
 : dx_loc_idx(dx), dy_loc_idx(dy), dz_loc_idx(dz),
   ux_loc_idx(ux), uy_loc_idx(uy), uz_loc_idx(uz)
 {
    mPower = 0;
    mTemperature = 300;	// 300K = 27C
    
	double x_length = fabs(dx_loc_idx - ux_loc_idx);
	double y_length = fabs(dy_loc_idx - uy_loc_idx);
	double z_length = fabs(dz_loc_idx - uz_loc_idx);
	
	x_res = (0.5 * x_length) / (k * y_length * z_length);
	y_res = (0.5 * y_length) / (k * x_length * z_length);
	z_res = (0.5 * z_length) / (k * x_length * y_length);
	z_bound_res = 1 / (heat_trans_coeff * x_length * y_length);
 }

 
ThermalPlane::ThermalPlane() 
 {} 
 

}
