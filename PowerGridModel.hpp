#ifndef PowerGridModel_H
#define PowerGridModel_H

#include <iostream>
#include <cstdio>
#include <algorithm>
#include <fstream>
#include <string>
#include <utility> 
#include <vector>
#include <unordered_map>
#include <map>
#include <queue>
#include <cmath>
#include <cassert>
#include <iomanip>
#include <cstdlib>
#include <thread>
#include <ctime>

#include "DebugMacro.hpp"
#include "Eigen"
#include "SparseCholesky"
#include "SparseLU"
#include "IterativeLinearSolvers"
#include "Dense"
#include "Eigenvalues"	

namespace PowerGridModel {

using triplet 			= Eigen::Triplet<double>;
using EigenSparseMatrix = Eigen::SparseMatrix<double>;
using EigenDenseMatrix  = Eigen::MatrixXd;
using EigenDenseVector  = Eigen::VectorXd;

struct CurrentSource;
struct VoltageSource;
struct Resistor;
struct Node;

using MATRIX    = std::vector< std::unordered_map<unsigned int, double> >;
using PG_PLANE  = std::unordered_map<std::string, Node*>;
using STAMPLIST = std::unordered_map<Node*, unsigned int>;

struct Resistor {
	
	Resistor();
	
	Resistor(std::string n, std::string pos, std::string neg, double value);
	
	void Add_Node(Node* n1, Node* n2);
	
	std::string res_name;
	
	std::string pos_node, neg_node; 
	
	double res_value;
	
	double curr_value;
	
	double power_value;
	
	double Vdrop;
	
	std::pair<Node*, Node*> shorted_node;	

	bool is_in_submatrix;

	unsigned int thermal_grid_index;

};


struct VoltageSource {
	
	VoltageSource();
	
	VoltageSource(std::string n, std::string pos, std::string neg, double v);

	std::string volt_name;
	
	std::string pos_node, neg_node; 
	
	double volt_value;
	
	double curr_value;
	
	double power_value;

};


struct CurrentSource {
	
	CurrentSource();
	
	CurrentSource(std::string n, std::string pos, std::string neg, double c);
	
	void ParseNodeInfo(const std::string &line, const std::string &seperators); 
	
	std::string curr_name;
	
	std::string pos_node, neg_node; 
	
	double curr_value;
	
	double power;
	
	unsigned int x_loc_idx, y_loc_idx;
	
};


struct Node {

	Node();
	
	Node(std::string n);

	~Node()
	{
		neighbor_NodeRes.clear();
		connected_Vsource.clear();
		connected_Isource.clear();
	}
	
	void Add_Neighbor_Node(Node* nei, Resistor* res);		
	
	void Add_Vsource(std::string n, VoltageSource* v);
	
	void Add_Isource(std::string n, CurrentSource* c);
	
	std::string node_name;
	
	double node_voltage;
	
	bool is_connect_V;	// if node connects to Vdd
	
	bool is_connect_I;	// if node connects to current source		
	
	bool is_grounded;	// if node connects to gnd	
	
	bool is_stamp;		
	
	unsigned int x_loc_idx, y_loc_idx;

	unsigned int thermal_grid_index;	// node at which thermal grid index
	
	void ParseNodeInfo(const std::string &line, const std::string &seperators);
	
	std::unordered_map<Node*, std::vector<Resistor*> > neighbor_NodeRes;	
	
	std::unordered_map<std::string, std::vector<VoltageSource*> > connected_Vsource;		
	
	std::unordered_map<std::string, std::vector<CurrentSource*> > connected_Isource;  
 	
};


class PowerGridData {
		
public : 

	PowerGridData();

	// Power Grid elements
	
	std::unordered_map<std::string, Node*> mTotalNodes;

	std::unordered_map<VoltageSource*, Node*> mtotal_supply_Vsource;	

	std::vector<Resistor*> mtotal_unshorted_Res;
	
	std::vector<CurrentSource*> mtotal_Isource;

	std::unordered_map<std::string, double> mVdropList;

	// Verification datas

	STAMPLIST mStampNodeList;	// is used to stamp PG conductance and subconductance matrix

	unsigned int mPGEquivalentNodeSize;

	std::vector<Node*> mPGEquivalentNodes;
	
	MATRIX mPG_OriginConductanceMatrix;
	
	MATRIX mPG_ThermalMatrix;

 	std::vector<EigenDenseVector> mEachRowPGInverseMatrix;	// stores each row

	EigenDenseMatrix mDenseInverseThermalMatrix;

	std::unordered_map<unsigned int, MATRIX> mPG_ConductanceSubMatrix;	  // thermal grid index to each sub matrix due to the resistors in each thermal grid

	std::unordered_map<unsigned int, EigenSparseMatrix> mPG_EigenConductanceSubMatrix; 
    
	EigenSparseMatrix mElectricalGridtoThermalGridIncidenceMatrix;  // E to T incidence matrix

	unsigned int mCurrentSourceSize;

	std::vector<unsigned int> mCurrentSourceIdx; 	// efficient index -> origin index

	std::unordered_map<unsigned int, unsigned int> mCurrentSourceIdxMapping;     // origin index -> efficient index

	std::vector<unsigned int> mThermalGridIndexConstainer;

};


class PowerGridSpecs
{

public : 

	PowerGridSpecs()
	{
		mSupplyVoltage  = 0; 
		mtotal_current  = 0; 
		mR_temp_coeff   = 4e-3;
		mViolations     = 0;
	}

	std::string mCktName;

	double mSupplyVoltage;		// Vdd voltage 

	unsigned int mViolations;	// voltage drop violations 

	double mtotal_current;	    // PG total peak current 

	// coefficients 

	double mR_temp_coeff;		// resistor temperature coefficient 

	double mwire_seg_Xlen;		// wire segment length 

	double mwire_seg_Ylen;

};


class PowerGridGeometry
{

public : 

	PowerGridGeometry() 
	{
		mX_loc_min = 1;
		mY_loc_min = 1;
	} 

	double mPowerPlaneXLength, mPowerPlaneYLength;
	
	double mX_loc_min, mX_loc_max, mY_loc_min, mY_loc_max;
	
	double mZ_metal_len, mZ_powerWidth, mZ_silicon_len;   // mZ_powerWidth corresponds to the active layer height which generates power
	
};


class PowerGridConstraints
{

public : 

	PowerGridConstraints() 
	{
		mXBlockDim  = 0;
		mYBlockDim  = 0;
		mGlobalCCSF = 0;
	}

	// Current constraint partition datas 

	double mXBlockDim, mYBlockDim;

	double mGlobalCCSF;

	unsigned int mGlobalConstraintSize;

	// Current constraint datas

	std::vector<double> mLocalCurrentConstraint;  

	std::vector<double> mGlobalCurrentConstraint;   

	EigenSparseMatrix mGlobalCurrentConstraintIncidenceMatrix; 

	EigenSparseMatrix mMatlabAmatrix;

};



}
#endif

