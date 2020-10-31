// PGSolver.cpp
// 1. Power grid DC simulation solver
// 2. Power grid verification solvers

#include "PGSolver.hpp"

void MSKAPI printstr(void *handle, const char str[])
{
	printf("%s",str);
} 

namespace PowerGridSolver {

PGSolver::PGSolver() 
{}



void PGSolver::SolvePowerGridPlane(PG_PLANE &Plane)
{
	
	std::cout << "Stamping Plane ..." << std::endl;
	
	STAMPLIST stampNodeList;

	unsigned int mapping_idx = 0;
	for(auto& node : Plane)
	{
		if(node.second->is_stamp)
		{
			stampNodeList.insert(std::pair<Node*, unsigned int>(node.second, mapping_idx));
			mapping_idx++;
		}
	}
	
	MATRIX MNA;
	
	double *RHS = NULL;
	
	unsigned int MatrixDim(0);
	
	MatrixDim = stampNodeList.size();
	
	MNA.resize(MatrixDim);
	
	RHS = new double[MatrixDim]();

	// Stamp Conductance Matrix
	for (auto& stamp_node : stampNodeList)
	{
		for (auto& nei_node : stamp_node.first->neighbor_NodeRes)    
		{
			if (nei_node.first->is_grounded == true)
			{
				for (auto& each_res : nei_node.second)
				{
					if(each_res->res_value !=0)
						MNA[stamp_node.second][stamp_node.second] += 1 / each_res->res_value;
				}
			}
			else if (nei_node.first->is_connect_V == true && nei_node.first->is_grounded == false)
			{
				for (auto& each_res : nei_node.second)
				{				
					if (each_res->res_value != 0) {
						MNA[stamp_node.second][stamp_node.second] += 1 / each_res->res_value;
						RHS[stamp_node.second]
							+= nei_node.first->node_voltage / each_res->res_value;
					}
				}
			}
			else if (nei_node.first->is_connect_V == false && nei_node.first->is_grounded == false && nei_node.first->is_stamp == true)
			{
				for (auto& each_res : nei_node.second)
				{
					if (each_res->res_value != 0) {
						MNA[stamp_node.second][stamp_node.second] += 1 / each_res->res_value;
						MNA[stamp_node.second][stampNodeList[nei_node.first]]
							+= -1 / each_res->res_value;
					}
				}
			}
		
		}
		for (auto& curr : stamp_node.first->connected_Isource)	
		{
			if (curr.first == "pos")
			{
				for (auto& each_curr : curr.second)
					RHS[stamp_node.second] += -each_curr->curr_value;
			}
			else if (curr.first == "neg")
			{
				for (auto& each_curr : curr.second)
					RHS[stamp_node.second] += +each_curr->curr_value;
			}
		}
	}
	std::cout << "Stammping Plane Complete !" << std::endl;
	
	
	std::cout << "Solving Plane ..." << std::endl;
	
	// Eigen matrix and vector 
	std::vector<triplet> tripletList;
	
	Eigen::SparseMatrix<double> Mat(MatrixDim, MatrixDim);
	
	Eigen::VectorXd b(MatrixDim), x;
	
	// Load matrix and vector elements to Eigen
	for(unsigned int row_idx = 0; row_idx < MNA.size(); ++row_idx)
    {
        for(auto& row_data : MNA[row_idx])
        {           
            tripletList.push_back(triplet(row_idx, row_data.first, row_data.second));   
        }
    }
	Mat.setFromTriplets(tripletList.begin(), tripletList.end());
	
	for(unsigned int i=0 ; i < MatrixDim ; ++i)
	{
		b[i] = RHS[i];
	}


	// Eigen Sparse Matrix Solver
	Eigen::SimplicialLLT <Eigen::SparseMatrix<double> > solver;
	
	solver.compute(Mat);
	if(solver.info()!= Eigen::Success) {
		std::cerr << "Matrix decomposition failed !" << std::endl;
		return ;
	}
	x = solver.solve(b);
	if(solver.info()!= Eigen::Success) {
		std::cerr << "Solving failed !" << std::endl;
		return ;
	}
	
	// assign voltage solution
	for (auto& node : stampNodeList)
	{
		node.first->node_voltage = x[node.second];
	}

	delete [] RHS;
}



// Solve PG voltage drop
void PGSolver::SolvePowerGridVdrop(PG_PLANE &Plane)
{
	
	std::cout << "Stamping Plane ..." << std::endl;
	
	STAMPLIST stampNodeList;

	unsigned int mapping_idx = 0;
	for (auto& node : Plane)
	{
		if(node.second->is_stamp)
		{
			stampNodeList.insert(std::pair<Node*, unsigned int>(node.second, mapping_idx));
			mapping_idx++;
		}
	}
	
	MATRIX MNA;
	
	double *RHS = NULL;
	
	unsigned int MatrixDim(0);
	
	MatrixDim = stampNodeList.size();
	
	MNA.resize(MatrixDim);
	
	RHS = new double[MatrixDim]();

	// Stamp Conductance Matrix
	for (auto& stamp_node : stampNodeList)
	{
		for (auto& nei_node : stamp_node.first->neighbor_NodeRes)    
		{
			if (nei_node.first->is_grounded == true)
			{
				for (auto& each_res : nei_node.second)
				{
					if(each_res->res_value !=0)
						MNA[stamp_node.second][stamp_node.second] += 1 / each_res->res_value;
				}
			}
			else if (nei_node.first->is_connect_V == true && nei_node.first->is_grounded == false)
			{
				for (auto& each_res : nei_node.second)
				{				
					if (each_res->res_value != 0) {
						
						MNA[stamp_node.second][stamp_node.second] += 1 / each_res->res_value;
						
						/*	Short voltage sources */ 
						// RHS[stamp_node.second] += 0 / each_res->res_value; 
					}
				}
			}
			else if (nei_node.first->is_connect_V == false && nei_node.first->is_grounded == false && nei_node.first->is_stamp == true)
			{
				for (auto& each_res : nei_node.second)
				{
					if (each_res->res_value != 0) {
						MNA[stamp_node.second][stamp_node.second] += 1 / each_res->res_value;
						MNA[stamp_node.second][stampNodeList[nei_node.first]] += -1 / each_res->res_value;
					}
				}
			}
		
		}
		for (auto& curr : stamp_node.first->connected_Isource)	
		{
			if (curr.first == "pos")	// Reverse each current source
			{
				for (auto& each_curr : curr.second)
					RHS[stamp_node.second] += each_curr->curr_value;
			}
			else if (curr.first == "neg")
			{
				for (auto& each_curr : curr.second)
					RHS[stamp_node.second] += -each_curr->curr_value;
			}
		}
	}
	std::cout << "Stammping Plane Complete !" << std::endl;
	
	
	std::cout << "Solving Plane ..." << std::endl;
	
	// Eigen matrix and vector 
	std::vector<triplet> tripletList;
	
	Eigen::SparseMatrix<double> Mat(MatrixDim, MatrixDim);
	
	Eigen::VectorXd b(MatrixDim), x;
	
	// Load matrix and vector elements to Eigen
	for(unsigned int row_idx = 0; row_idx < MNA.size(); ++row_idx)
    {
        for(auto& row_data : MNA[row_idx])
        {           
            tripletList.push_back(triplet(row_idx, row_data.first, row_data.second));   
        }
    }
	Mat.setFromTriplets(tripletList.begin(), tripletList.end());
	
	for(unsigned int i=0 ; i < MatrixDim ; ++i)
	{
		b[i] = RHS[i];
	}


	// Eigen Sparse Matrix Solver
	Eigen::SimplicialLLT <Eigen::SparseMatrix<double> > solver;
	
	solver.compute(Mat);
	if(solver.info()!= Eigen::Success) {
		std::cerr << "Matrix decomposition failed !" << std::endl;
		return ;
	}
	x = solver.solve(b);
	if(solver.info()!= Eigen::Success) {
		std::cerr << "Solving failed !" << std::endl;
		return ;
	}
	
	// assign voltage drop solution
	for (auto& node : stampNodeList)
	{
		mVdropList[node.first->node_name] = x[node.second];
	}

	delete [] RHS;
}


void PGSolver::Cal_Total_Current_Power()
{
		
		std::cout << "Calculating current and power dissipation ..." << std::endl;
		
		for (auto& supply_node : mtotal_supply_Vsource)		
		{
			double total_current = 0;
			for (auto& neighbor : supply_node.second->neighbor_NodeRes)
			{
				if(neighbor.first->is_stamp)
					for (auto& nei_res : neighbor.second)
					{
						total_current +=
							(supply_node.first->volt_value - neighbor.first->node_voltage) / nei_res->res_value;
					}
			}
			supply_node.first->curr_value = -1 * total_current;
			supply_node.first->power_value = supply_node.first->volt_value * fabs(supply_node.first->curr_value);
		}
		
		
		for (auto& each_res : mtotal_unshorted_Res)	
		{	
			each_res->Vdrop = each_res->shorted_node.first->node_voltage - each_res->shorted_node.second->node_voltage;
			each_res->curr_value = each_res->Vdrop / each_res->res_value;
			each_res->power_value = pow(each_res->curr_value, 2) * each_res->res_value;	
		}
	
		std::cout << "Current and power dissipation Calculation complete !" << std::endl;
}


void PGSolver::SolveThermalPlane(THERMAL_PLANE & thermal_plane)
{
	unsigned int MatrixDim = (unsigned int)(mZlayer*mYCut*mXCut);
	
	MATRIX ThermalMat;
	
	mPowerVector.resize(MatrixDim, 0); 
	
	ThermalMat.resize(MatrixDim);

	std::cout << "Stamping Thermal Plane ..." << std::endl;

	for(int z = 0 ; z < mZlayer ; ++z)
		for(int y = 0; y < mYCut ; ++y)
			for(int x = 0; x < mXCut ; ++x)
	{
				
				if(x-1 < 0)   // at x plane min bound
				{
					ThermalMat[thermal_plane[z][y][x].mGridIndex][thermal_plane[z][y][x].mGridIndex]
					+= 1/(thermal_plane[z][y][x].x_res + thermal_plane[z][y][x+1].x_res);
					
					ThermalMat[thermal_plane[z][y][x].mGridIndex][thermal_plane[z][y][x+1].mGridIndex]
					+= -1/(thermal_plane[z][y][x].x_res + thermal_plane[z][y][x+1].x_res);
				}
				else if(x-1 >= 0 && x+1 < mXCut) // at x plane between
				{
					ThermalMat[thermal_plane[z][y][x].mGridIndex][thermal_plane[z][y][x].mGridIndex]
					+= 1/(thermal_plane[z][y][x].x_res + thermal_plane[z][y][x-1].x_res)+
					   1/(thermal_plane[z][y][x].x_res + thermal_plane[z][y][x+1].x_res);
					
					ThermalMat[thermal_plane[z][y][x].mGridIndex][thermal_plane[z][y][x+1].mGridIndex]
					+= -1/(thermal_plane[z][y][x].x_res + thermal_plane[z][y][x+1].x_res);
					
					ThermalMat[thermal_plane[z][y][x].mGridIndex][thermal_plane[z][y][x-1].mGridIndex]
					+= -1/(thermal_plane[z][y][x].x_res + thermal_plane[z][y][x-1].x_res);
				}
				else if(x+1 >= mXCut) // at x plane max bound
				{
					ThermalMat[thermal_plane[z][y][x].mGridIndex][thermal_plane[z][y][x].mGridIndex]
					+= 1/(thermal_plane[z][y][x].x_res + thermal_plane[z][y][x-1].x_res);
					
					ThermalMat[thermal_plane[z][y][x].mGridIndex][thermal_plane[z][y][x-1].mGridIndex]
					+= -1/(thermal_plane[z][y][x].x_res + thermal_plane[z][y][x-1].x_res);	
				}


				
				if(y-1 < 0)   
				{
					ThermalMat[thermal_plane[z][y][x].mGridIndex][thermal_plane[z][y][x].mGridIndex]
					+= 1/(thermal_plane[z][y][x].y_res + thermal_plane[z][y+1][x].y_res);
					
					ThermalMat[thermal_plane[z][y][x].mGridIndex][thermal_plane[z][y+1][x].mGridIndex]
					+= -1/(thermal_plane[z][y][x].y_res + thermal_plane[z][y+1][x].y_res);
				}
				else if(y-1 >= 0 && y+1 < mYCut) 
				{
					ThermalMat[thermal_plane[z][y][x].mGridIndex][thermal_plane[z][y][x].mGridIndex]
					+= 1/(thermal_plane[z][y][x].y_res + thermal_plane[z][y-1][x].y_res)+
					   1/(thermal_plane[z][y][x].y_res + thermal_plane[z][y+1][x].y_res);
					
					ThermalMat[thermal_plane[z][y][x].mGridIndex][thermal_plane[z][y+1][x].mGridIndex]
					+= -1/(thermal_plane[z][y][x].y_res + thermal_plane[z][y+1][x].y_res);
					
					ThermalMat[thermal_plane[z][y][x].mGridIndex][thermal_plane[z][y-1][x].mGridIndex]
					+= -1/(thermal_plane[z][y][x].y_res + thermal_plane[z][y-1][x].y_res);	
				}
				else if(y+1 >= mYCut) 
				{
					ThermalMat[thermal_plane[z][y][x].mGridIndex][thermal_plane[z][y][x].mGridIndex]
					+= 1/(thermal_plane[z][y][x].y_res + thermal_plane[z][y-1][x].y_res);
					
					ThermalMat[thermal_plane[z][y][x].mGridIndex][thermal_plane[z][y-1][x].mGridIndex]
					+= -1/(thermal_plane[z][y][x].y_res + thermal_plane[z][y-1][x].y_res);		
				}
					
				if(z-1 < 0)   
				{
					ThermalMat[thermal_plane[z][y][x].mGridIndex][thermal_plane[z][y][x].mGridIndex]
					+= 1/(thermal_plane[z][y][x].z_res + thermal_plane[z][y][x].z_bound_res)+
					   1/(thermal_plane[z][y][x].z_res + thermal_plane[z+1][y][x].z_res);
					
					ThermalMat[thermal_plane[z][y][x].mGridIndex][thermal_plane[z+1][y][x].mGridIndex]
					+= -1/(thermal_plane[z][y][x].z_res + thermal_plane[z+1][y][x].z_res);
				}
				else if(z-1 >= 0 && z+1 < mZlayer) 
				{
					ThermalMat[thermal_plane[z][y][x].mGridIndex][thermal_plane[z][y][x].mGridIndex]
					+= 1/(thermal_plane[z][y][x].z_res + thermal_plane[z-1][y][x].z_res)+
					   1/(thermal_plane[z][y][x].z_res + thermal_plane[z+1][y][x].z_res);
					
					ThermalMat[thermal_plane[z][y][x].mGridIndex][thermal_plane[z+1][y][x].mGridIndex]
					+= -1/(thermal_plane[z][y][x].z_res + thermal_plane[z+1][y][x].z_res);
					
					ThermalMat[thermal_plane[z][y][x].mGridIndex][thermal_plane[z-1][y][x].mGridIndex]
					+= -1/(thermal_plane[z][y][x].z_res + thermal_plane[z-1][y][x].z_res);	
				}
				else if(z+1 >= mZlayer) 
				{
					ThermalMat[thermal_plane[z][y][x].mGridIndex][thermal_plane[z][y][x].mGridIndex]
					+= 1/(thermal_plane[z][y][x].z_res + thermal_plane[z][y][x].z_bound_res)+
					   1/(thermal_plane[z][y][x].z_res + thermal_plane[z-1][y][x].z_res);
					
					ThermalMat[thermal_plane[z][y][x].mGridIndex][thermal_plane[z-1][y][x].mGridIndex]
					+= -1/(thermal_plane[z][y][x].z_res + thermal_plane[z-1][y][x].z_res);	
				}
				
				// assign power value at active layer
				if(z == mZ_subCut)
					mPowerVector[thermal_plane[z][y][x].mGridIndex] += thermal_plane[z][y][x].mPower;	
	
	}


	std::cout << "Solving Plane ..." << std::endl;

	std::vector<triplet> tripletList;
	
	Eigen::SparseMatrix<double> Mat(MatrixDim, MatrixDim);
	
	Eigen::VectorXd b(MatrixDim), x;

	for(unsigned int row_idx = 0; row_idx < ThermalMat.size(); ++row_idx)
    {
        for(auto& row_data : ThermalMat[row_idx])
        {           
            tripletList.push_back(triplet(row_idx, row_data.first, row_data.second));   
        }
    }
	Mat.setFromTriplets(tripletList.begin(), tripletList.end());
	
	for(unsigned int i=0 ; i < MatrixDim ; ++i)
	{
		b[i] = mPowerVector[i];
	}

	Eigen::SimplicialLLT < Eigen::SparseMatrix<double> > solver;
	
	solver.compute(Mat);
	if(solver.info()!= Eigen::Success) {
		std::cerr << "Matrix decomposition failed !" << std::endl;
		return ;
	}
	x = solver.solve(b);
	if(solver.info()!= Eigen::Success) {
		std::cerr << "Solving failed !" << std::endl;
		return ;
	}
	
	mThermalSolution.resize(MatrixDim, 0); 
	
	for (unsigned int i=0 ; i < MatrixDim ; ++i)
	{
		mThermalSolution[i] = x[i];
	}

	mPowerVector.clear();
	ThermalMat.clear();
	
}


// --------------------------------------------------- For verification ---------------------------------------------------

//2017.7.31 ~ 
// Efficient Method Thermal Analysis
// types
// 0  QP solution
// 1  Expand point 
// 2  Expand point solution
void PGSolver::OnePointThermalAnalysis
(unsigned int NodeIndex, double *CurrentSolution, EigenDenseMatrix & ThermalSystemMatrix, int type)
{
	
	std::cout << "Node x location = " << mPGEquivalentNodes[NodeIndex]->x_loc_idx << std::endl;
	std::cout << "Node y location = " << mPGEquivalentNodes[NodeIndex]->y_loc_idx << std::endl;
	double *xx = CurrentSolution;
	EigenDenseVector CurrentPattern(mPGEquivalentNodeSize), TemperatureSolution;
	for(unsigned int i = 0 ; i < mPGEquivalentNodeSize ; ++i) 
	{ 
		CurrentPattern[i] = 0; 
	}
	for(unsigned int i = 0 ; i < mCurrentSourceSize ; ++i) 
	{ 
		CurrentPattern[mCurrentSourceIdx[i]] = xx[i]; 
	}
	TemperatureSolution = mSupplyVoltage * (ThermalSystemMatrix * CurrentPattern);
	
	std::string MTMap;
	std::string CurrentMap;

	if( type == 0 )
	{
		MTMap = "_EfficientQPMTMap_";
		CurrentMap = "_EfficientQPCurrentMap_";
	}
	else if( type == 1 )     //  Efficient Expand Point
	{
		MTMap = "_EfficientExpandPointMTMap_";
		CurrentMap = "_EfficientExpandPointCurrentMap_";
	}
	else if( type == 2 )     //  Efficient Expand Point
	{
		MTMap = "_EfficientExpandPointSolutionMTMap_";
		CurrentMap = "_EfficientExpandPointSolutionCurrentMap_";
	}

	// Output Temperature Map
	std::fstream fout("../PointMap/" + mCktName + MTMap + std::to_string(NodeIndex), std::ios::out);
	for(unsigned int y = 0; y < mYCut ; ++y)
	{
		for(unsigned int x = 0; x < mXCut ; ++x)
		{
			fout << x
		         << " " 
			     << y
			     << " " 
			     << 23 + TemperatureSolution[mTotalThermalGrids[mZlayer-1][y][x].mGridIndex] << std::endl;	
		}
	}
	fout.close();
	// Output Current Map
	fout.open("../PointMap/" + mCktName + CurrentMap + std::to_string(NodeIndex), std::ios::out);
	for(auto& PNode : mPGEquivalentNodes)
	{
		unsigned int idx((PNode->x_loc_idx - 1) * mY_loc_max + (PNode->y_loc_idx - 1));
		
		if(PNode->is_connect_I)
			fout << PNode->x_loc_idx
		  		 << " " 
			 	 << PNode->y_loc_idx
			     << " " 
			     << xx[mCurrentSourceIdxMapping[idx]] << std::endl;	
		else {
			fout << PNode->x_loc_idx
		         << " " 
			 	 << PNode->y_loc_idx
			 	 << " " 
			 	 << 0 << std::endl;	
		}
	}	
	fout.close();


}


// Direct QP 
void PGSolver::PGVerificationQP()
{

	double Dmatrixtime(0);

	clock_t start;
	start = clock();

	EigenDenseMatrix InversePGConductanceMatrix(mPGEquivalentNodeSize, mPGEquivalentNodeSize);
	for(unsigned int i = 0 ; i < mPGEquivalentNodeSize ; ++i)
	{
		InversePGConductanceMatrix.row(i) = mEachRowPGInverseMatrix[i];
	}
	EigenSparseMatrix sInversePGConductanceMatrix(mPGEquivalentNodeSize, mPGEquivalentNodeSize);
	sInversePGConductanceMatrix = InversePGConductanceMatrix.sparseView();

	std::cout << "Construct pg inverse matrix complete !" << std::endl;
	
	// Store D matrix
	// k is the thermal grid index
	std::map<unsigned int, EigenSparseMatrix> mDMatrix, mPreDMatrix; // thermal grid index to pg subconductance matrix
	
	for(auto& each_conductance_submatrix : mPG_EigenConductanceSubMatrix)
	{
		mPreDMatrix[each_conductance_submatrix.first] = sInversePGConductanceMatrix * (each_conductance_submatrix.second);
	}
	// std::cout << "Pre mat complete" << std::endl;
	for(auto& each_conductance_submatrix : mPreDMatrix)
	{
		mDMatrix[each_conductance_submatrix.first] = each_conductance_submatrix.second * sInversePGConductanceMatrix;
	}
	mPreDMatrix.clear();

	// Construct INVGT * M_EtoT, find the grid index to calculate the temperature, ak
	EigenDenseMatrix AkMatrix = mDenseInverseThermalMatrix * Eigen::MatrixXd(mElectricalGridtoThermalGridIncidenceMatrix);
	
	// load AkMatrix to sparse vector data structure
	unsigned int thermalgridsize = mXCut*mYCut*mZlayer;
	unsigned int ArraySize = thermalgridsize*mPGEquivalentNodeSize;
	double *cAkMatrix = new double[ArraySize]();
	Eigen::Map<Eigen::MatrixXd>( cAkMatrix, AkMatrix.rows(), AkMatrix.cols() ) = AkMatrix;
	MATRIX SparseAkMatrix;
	SparseAkMatrix.resize(thermalgridsize);
	for(unsigned int row = 0 ; row < thermalgridsize ; ++row)
	{
		for(unsigned int col = 0 ; col < mPGEquivalentNodeSize ; ++col)
		{
			if(cAkMatrix[row+col*thermalgridsize] != 0)
			{
				SparseAkMatrix[row][col] = cAkMatrix[row+col*thermalgridsize];
			}
		}
	}

	Dmatrixtime = (double)(clock()-start) / (CLOCKS_PER_SEC);
	std::cout << "Construct D matrix time = " << Dmatrixtime << std::endl;


	// start matlab engine

	Engine *ep;

	if (!(ep = engOpen(""))) {
		fprintf(stderr, "\nCan't start MATLAB engine\n");
		return ;
	}
	else {
		std::cout << "Started matlab engine !" << std::endl;
	}

	// construct matlab QP A matrix
	//-----------------------------------------------------------------------------------------------------------
	// load mMatlabAmatrix & lower bound constraint to Amatrix  => 0 < Ui < Ig
	unsigned int Amatrixsize = 2 * mGlobalConstraintSize * mPGEquivalentNodeSize;
	double *Amatrix = new double[Amatrixsize];
	Eigen::Map<Eigen::MatrixXd>( Amatrix, mMatlabAmatrix.rows(), mMatlabAmatrix.cols() ) = Eigen::MatrixXd(mMatlabAmatrix);
    
  	// load Amatrix to matAmatrix
	mxArray *matAmatrix = NULL;
	matAmatrix = mxCreateDoubleMatrix(2 * mGlobalConstraintSize, mPGEquivalentNodeSize, mxREAL);
	memcpy((void *)mxGetPr(matAmatrix), (void *)Amatrix, sizeof(double) * Amatrixsize);
	engPutVariable(ep, "matAmatrix", matAmatrix);	// put variable to matlab workspace 
	//-----------------------------------------------------------------------------------------------------------
	
	// construct matlab QP b vector
	double *b = new double[2*mGlobalConstraintSize]();
	for(unsigned int i=0 ; i < mGlobalConstraintSize ; ++i)
	{
		b[i] = mGlobalCurrentConstraint[i];
	}
	mxArray *matb = NULL;
	matb = mxCreateDoubleMatrix(1, 2 * mGlobalConstraintSize, mxREAL);
	memcpy((void *)mxGetPr(matb), (void *)b, sizeof(double) * 2 * mGlobalConstraintSize);
	engPutVariable(ep, "matb", matb);

	// construct matlab QP local upper and lower variable constraints => 0 < i < Il
	//-----------------------------------------------------------------------------------------------------------
	double *clower = new double[mPGEquivalentNodeSize]();
	double *cupper = new double[mPGEquivalentNodeSize];
	for(unsigned int i=0 ; i < mPGEquivalentNodeSize ; ++i)
	{
		cupper[i] = mLocalCurrentConstraint[i];
	}
	mxArray *matlower = NULL;
	mxArray *matupper = NULL;
	matlower = mxCreateDoubleMatrix(1, mPGEquivalentNodeSize, mxREAL);
	matupper = mxCreateDoubleMatrix(1, mPGEquivalentNodeSize, mxREAL);
	memcpy((void *)mxGetPr(matlower), (void *)clower, sizeof(double) * mPGEquivalentNodeSize);
	memcpy((void *)mxGetPr(matupper), (void *)cupper, sizeof(double) * mPGEquivalentNodeSize);
	engPutVariable(ep, "matlower", matlower);
	engPutVariable(ep, "matupper", matupper);
	//-----------------------------------------------------------------------------------------------------------
	
	//-----------------------------------------------------------------------------------------------------------
	// construct matlab linear terms 
	double *cf = new double[mPGEquivalentNodeSize];
	mxArray *matf = NULL;
	matf = mxCreateDoubleMatrix(1, mPGEquivalentNodeSize, mxREAL);

	// H matrix
	mxArray *matHmatrix = NULL;
	matHmatrix = mxCreateDoubleMatrix(mPGEquivalentNodeSize, mPGEquivalentNodeSize, mxREAL);
	double *cHbarmatrix = new double[mPGEquivalentNodeSize*mPGEquivalentNodeSize]();
	double *cHmatrix = new double[mPGEquivalentNodeSize*mPGEquivalentNodeSize]();
	

	// Factor
	double sFactor = -1 * mR_temp_coeff * mSupplyVoltage;	

	// solution
	std::vector<double> VdropSolution;

	// cpu time data
	engEvalString(ep, "qptime = 0;");	// for calculating matlab quadprog time 
	clock_t matrixStart; 
	double QPtime(0), MatrixTime(0);
	


for(unsigned int each_node_idx = 0 ; each_node_idx < mPGEquivalentNodeSize ; ++each_node_idx)
{
	
	// construct matlab f vector
	matrixStart = clock();

	for(unsigned int i=0 ; i<mPGEquivalentNodeSize ; ++i)
	{
		cf[i] = -1 * mEachRowPGInverseMatrix[each_node_idx][i];
	}
	memcpy((void *)mxGetPr(matf), (void *)cf, sizeof(double) * mPGEquivalentNodeSize);
	engPutVariable(ep, "matf", matf);	
	// std::cout << "row time = " << std::setprecision(15) << (double)(clock()-start) / (CLOCKS_PER_SEC) << std::endl;
		
	// Construct cHbarmatrix
	for(unsigned int i=0 ; i < mPGEquivalentNodeSize*mPGEquivalentNodeSize ; ++i)
	{
		cHbarmatrix[i] = 0;
	}

	for(auto& each_grid : mDMatrix)		
	{
		EigenDenseVector DkRow = each_grid.second.row(each_node_idx);

		for(auto each_row : SparseAkMatrix[each_grid.first])
		{
			for(unsigned int col = 0 ; col < mPGEquivalentNodeSize ; ++col)
			{
				if(DkRow[col] != 0)
				{
					unsigned int idx = each_row.first + col * mPGEquivalentNodeSize;
					cHbarmatrix[idx] += each_row.second * DkRow[col];	
				}
			}
		}
	}
	
	// load sFactor * (cHbarmatrix + cHbarmatrix transpose) to cHmatrix
	for(unsigned int row = 0 ; row < mPGEquivalentNodeSize ; ++row)
	{
		for(unsigned int col = 0 ; col < mPGEquivalentNodeSize ; ++col)
		{
			unsigned int idx     =  row+col*mPGEquivalentNodeSize;
			unsigned int tranidx =  col+row*mPGEquivalentNodeSize;
			cHmatrix[idx] = sFactor * (cHbarmatrix[idx] + cHbarmatrix[tranidx]);
		}
	}
	memcpy((void *)mxGetPr(matHmatrix), (void *)cHmatrix, sizeof(double) * mPGEquivalentNodeSize * mPGEquivalentNodeSize);
	engPutVariable(ep, "matHmatrix", matHmatrix);	// put variable to matlab workspace 

	MatrixTime += (double)(clock() - matrixStart) / (CLOCKS_PER_SEC);


	// obtain QP solution & CPU time
	engEvalString(ep, "tStart = tic;");
	engEvalString(ep, "[x,fval,exitflag] = quadprog(matHmatrix,matf,matAmatrix,matb,[],[],matlower,matupper);");
	engEvalString(ep, "qptime = toc(tStart);");	
	QPtime += *(double *)mxGetPr(engGetVariable(ep, "qptime"));
	

	// store voltage drop solution
	VdropSolution.emplace_back( -1 * *(double *)mxGetPr(engGetVariable(ep, "fval")) );


	// thermal analysis
	//---------------------------------------------------------------------------------------------------------
	/*
	std::cout << "Node x location = " << mPGEquivalentNodes[each_node_idx]->x_loc_idx << std::endl;
	std::cout << "Node y location = " << mPGEquivalentNodes[each_node_idx]->y_loc_idx << std::endl;

	double *xx = (double *)mxGetPr(engGetVariable(ep, "x"));
	EigenDenseVector CurrentPattern(mPGEquivalentNodeSize), TemperatureSolution;
	for(unsigned int i = 0 ; i < mPGEquivalentNodeSize ; ++i) { CurrentPattern[i] = xx[i];}
	TemperatureSolution = mSupplyVoltage * (AkMatrix * CurrentPattern);
	// Output Temperature Map
	std::fstream fout("../result/" + mCktName + "_QPMTMap_" + std::to_string(each_node_idx), std::ios::out);
	for(unsigned int y = 0; y < mYCut ; ++y)\
	{
		for(unsigned int x = 0; x < mXCut ; ++x)
		{
			fout << x
		         << " " 
			     << y
			     << " " 
			     << 23 + TemperatureSolution[mTotalThermalGrids[mZlayer-1][y][x].mGridIndex] << std::endl;	
		}
	}
	fout.close();
	// Output Current Map
	fout.open("../result/" + mCktName + "_QPCurrentMap_" + std::to_string(each_node_idx), std::ios::out);
	for(auto& PNode : mPGEquivalentNodes)
	{
		unsigned int idx((PNode->x_loc_idx - 1) * mY_loc_max + (PNode->y_loc_idx - 1));
		fout << PNode->x_loc_idx
		     << " " 
			 << PNode->y_loc_idx
			 << " " 
			 << xx[idx] << std::endl;	
	}	
	fout.close();
	*/
	//---------------------------------------------------------------------------------------------------------

}
	
	std::cout << "Matrix Time = " << std::setprecision(15) << MatrixTime << std::endl;

	std::cout << "Matlab QP Time = " << std::setprecision(15) << QPtime << std::endl;

	std::cout << "Total time (Matrix + Matlab) = " << std::setprecision(15) << Dmatrixtime + MatrixTime + QPtime << std::endl;

	// Output solution 
	
	double sum(0);
	for(auto& vdrop : VdropSolution) { sum +=vdrop; }

	std::cout << "Average voltage drop = " << std::setprecision(15) << sum / VdropSolution.size() << std::endl;
	std::cout << "Maximum voltage drop = " << std::setprecision(15) << *std::max_element(VdropSolution.begin(), VdropSolution.end()) << std::endl;
	
	std::ofstream fout("../LPresult/" + mCktName + "_QPVdropMap");
	for(unsigned int i=0;i<mPGEquivalentNodeSize;++i)
	{
		fout << mPGEquivalentNodes[i]->x_loc_idx << " " << mPGEquivalentNodes[i]->y_loc_idx << " " << std::setprecision(15) << VdropSolution[i] << std::endl;
	}
	fout.close();
	
	// Delete memory

	mxDestroyArray(matHmatrix);
	mxDestroyArray(matf);
	mxDestroyArray(matAmatrix);
	mxDestroyArray(matb);
	mxDestroyArray(matlower);
	mxDestroyArray(matupper);
	engClose(ep);

	delete [] Amatrix;
	delete [] b;
	delete [] clower;
	delete [] cupper;
	delete [] cf;
	delete [] cHbarmatrix;
	delete [] cHmatrix;
	
}



// Efficient QP
// using Matlab for QP
void PGSolver::PGVerificationQP_Efficient()
{

    double Dmatrixtime(0);

	clock_t start;
	
	start = clock();

	std::map<unsigned int, MySparseMatrixClass::SparseMatrix> MyPGConductanceSubMatrix;	 

	for(auto& each_conductance_submatrix : mPG_EigenConductanceSubMatrix)
	{
		(each_conductance_submatrix.second).makeCompressed();	

		int    *aptrb    = (each_conductance_submatrix.second).outerIndexPtr();	  
		int    *asub     = (each_conductance_submatrix.second).innerIndexPtr();	 
		double *aval     = (each_conductance_submatrix.second).valuePtr();		
		unsigned int NonZerosIndexAdder(0);

		MySparseMatrixClass::SparseMatrix tmpSparseMatrix;
		for(unsigned int Col = 0 ; Col < mPGEquivalentNodeSize ; ++Col)
		{
			MySparseMatrixClass::SparseVector tmpSparseVector((aptrb[Col+1] - aptrb[Col]), asub, aval, &NonZerosIndexAdder);
			tmpSparseMatrix.PushBack(tmpSparseVector);
		}
		MyPGConductanceSubMatrix[each_conductance_submatrix.first] = tmpSparseMatrix;
	}
	mPG_EigenConductanceSubMatrix.clear();


	double *cG0INV_Q = new double[mPGEquivalentNodeSize * mCurrentSourceSize]();

	for(unsigned int each_row = 0 ; each_row < mPGEquivalentNodeSize ; ++each_row)
	{
		for(unsigned int each_col = 0 ; each_col < mCurrentSourceSize ; ++each_col)
		{
			double value = mEachRowPGInverseMatrix[each_row][mCurrentSourceIdx[each_col]];
			if(value != 0)
			{
				cG0INV_Q[each_row + each_col*mPGEquivalentNodeSize] = value;
			}
		}
	}

	// Construct INVGT * M_EtoT, find the grid index to calculate the temperature, ak
	EigenDenseMatrix AkMatrix = mDenseInverseThermalMatrix * Eigen::MatrixXd(mElectricalGridtoThermalGridIncidenceMatrix);
	
	// load AkMatrix to sparse vector data structure
	unsigned int thermalgridsize = mXCut * mYCut * mZlayer;
	
	unsigned int ArraySize = thermalgridsize * mPGEquivalentNodeSize;
	
	double *cAkMatrix = new double[ArraySize]();
	
	Eigen::Map<Eigen::MatrixXd>( cAkMatrix, AkMatrix.rows(), AkMatrix.cols() ) = AkMatrix;
	

	std::vector< std::vector<double> > DenseAkMatrix(thermalgridsize);	// Q * Ak
	
	for(unsigned int i=0 ; i<thermalgridsize ; ++i )
		DenseAkMatrix[i].resize(mCurrentSourceSize);

	// Factor
	double sFactor = -1 * mR_temp_coeff * mSupplyVoltage;	// -1 * alpha * Vdd 	
	
	for(unsigned int row = 0 ; row < thermalgridsize ; ++row)
	{
		for(unsigned int col = 0 ; col < mPGEquivalentNodeSize ; ++col)
		{
			unsigned int index = row + col * thermalgridsize;
			if(cAkMatrix[index] != 0)
			{
				DenseAkMatrix[row][mCurrentSourceIdxMapping[col]] = sFactor * cAkMatrix[index];
			}
		}
	}
	
	Dmatrixtime = (double)(clock()-start) / (CLOCKS_PER_SEC);
	
	std::cout << "Construct D matrix time = " << Dmatrixtime << std::endl;

	// start matlab engine

	Engine *ep;

	if (!(ep = engOpen(""))) {
		fprintf(stderr, "\nCan't start MATLAB engine\n");
		return ;
	}
	else {
		std::cout << "Started matlab engine !" << std::endl;
	}

	// engEvalString(ep, "LASTN = maxNumCompThreads(1);");

	// construct matlab QP A matrix
	//-----------------------------------------------------------------------------------------------------------
	// load mMatlabAmatrix & lower bound constraint to Amatrix  => 0 < Ui < Ig
	unsigned int Amatrixsize = 2 * mGlobalConstraintSize * mCurrentSourceSize;
	double *Amatrix = new double[Amatrixsize];
	Eigen::Map<Eigen::MatrixXd>( Amatrix, mMatlabAmatrix.rows(), mMatlabAmatrix.cols() ) = Eigen::MatrixXd(mMatlabAmatrix);
    
  	// load Amatrix to matAmatrix
	mxArray *matAmatrix = NULL;
	matAmatrix = mxCreateDoubleMatrix(2 * mGlobalConstraintSize, mCurrentSourceSize, mxREAL);
	memcpy((void *)mxGetPr(matAmatrix), (void *)Amatrix, sizeof(double) * Amatrixsize);
	engPutVariable(ep, "matAmatrix", matAmatrix);	// put variable to matlab workspace 
	//-----------------------------------------------------------------------------------------------------------
	
	// construct matlab QP b vector
	double *b = new double[2*mGlobalConstraintSize]();
	for(unsigned int i=0 ; i < mGlobalConstraintSize ; ++i)
	{
		b[i] = mGlobalCurrentConstraint[i];
	}
	mxArray *matb = NULL;
	matb = mxCreateDoubleMatrix(1, 2 * mGlobalConstraintSize, mxREAL);
	memcpy((void *)mxGetPr(matb), (void *)b, sizeof(double) * 2 * mGlobalConstraintSize);
	engPutVariable(ep, "matb", matb);

	// construct matlab QP local upper and lower variable constraints => 0 < i < Il
	//-----------------------------------------------------------------------------------------------------------
	double *clower = new double[mCurrentSourceSize]();
	double *cupper = new double[mCurrentSourceSize];
	for(unsigned int i=0 ; i < mCurrentSourceSize ; ++i)
	{
		cupper[i] = mLocalCurrentConstraint[i];
	}
	mxArray *matlower = NULL;
	mxArray *matupper = NULL;
	matlower = mxCreateDoubleMatrix(1, mCurrentSourceSize, mxREAL);
	matupper = mxCreateDoubleMatrix(1, mCurrentSourceSize, mxREAL);
	memcpy((void *)mxGetPr(matlower), (void *)clower, sizeof(double) * mCurrentSourceSize);
	memcpy((void *)mxGetPr(matupper), (void *)cupper, sizeof(double) * mCurrentSourceSize);
	engPutVariable(ep, "matlower", matlower);
	engPutVariable(ep, "matupper", matupper);
	//-----------------------------------------------------------------------------------------------------------
	
	//-----------------------------------------------------------------------------------------------------------
	// construct matlab linear terms 
	double *cf = new double[mCurrentSourceSize];
	mxArray *matf = NULL;
	matf = mxCreateDoubleMatrix(1, mCurrentSourceSize, mxREAL);

	// H matrix
	mxArray *matHmatrix = NULL;
	matHmatrix = mxCreateDoubleMatrix(mCurrentSourceSize, mCurrentSourceSize, mxREAL);
	double *cHbarmatrix = new double[mCurrentSourceSize*mCurrentSourceSize]();
	double *cHmatrix = new double[mCurrentSourceSize*mCurrentSourceSize]();
	
	// solution
	std::vector<double> VdropSolution;

	
	engEvalString(ep, "qptime = 0;");	// for calculating matlab quadprog time 
	clock_t matrixStart; 
	double QPtime(0), MatrixTime(0);

	EigenDenseVector DkRow(mCurrentSourceSize);

    //int debug_idx(0);
for(unsigned int each_node_idx = 0 ; each_node_idx < mPGEquivalentNodeSize ; ++each_node_idx)
//for(unsigned int each_node_idx = 5050 ; each_node_idx < 5050+1 ; ++each_node_idx)
{
	
if(mPGEquivalentNodes[each_node_idx]->is_connect_I)
{

	// construct matlab f vector
	matrixStart = clock();

	for(unsigned int i=0 ; i<mCurrentSourceSize ; ++i)
	{
		cf[i] = -1 * mEachRowPGInverseMatrix[each_node_idx][mCurrentSourceIdx[i]];
	}
	memcpy((void *)mxGetPr(matf), (void *)cf, sizeof(double) * mCurrentSourceSize);
	engPutVariable(ep, "matf", matf);	
	//std::cout << "row time = " << std::setprecision(15) << (double)(clock()-start) / (CLOCKS_PER_SEC) << std::endl;
	
	// Construct cHbarmatrix
	for(unsigned int i = 0 ; i < mCurrentSourceSize * mCurrentSourceSize ; ++i)
	{
		cHbarmatrix[i] = 0;
	}

	for(auto& each_grid : MyPGConductanceSubMatrix)		
	{
		MySparseMatrixClass::SparseVector MyPreRow; 	// G0INV * G0,k

		for(unsigned int each_col = 0 ; each_col < mPGEquivalentNodeSize ; ++each_col)
		{
			double adjointValue(0);
			for(unsigned int i = 0 ; i < each_grid.second.getSparseVector(each_col).getVector().size() ; ++i)
			{
				adjointValue += each_grid.second.getSparseVector(each_col).getVector()[i] * mEachRowPGInverseMatrix[each_node_idx][each_grid.second.getSparseVector(each_col).getLocation()[i]];
			}
			if(adjointValue != 0)
			{
				MyPreRow.PushBackElement(each_col, adjointValue);
			}
		}

		for(unsigned int each_col = 0 ; each_col < mCurrentSourceSize ; ++each_col)	  //   Dk = (G0INV * G0,k) * (G0INV * Q)
		{
			double adjointValue(0);
			for(unsigned int i = 0 ; i < MyPreRow.getVector().size() ; ++i)
			{
				unsigned int index = MyPreRow.getLocation()[i] + each_col * mPGEquivalentNodeSize;
				if( cG0INV_Q[index] != 0)
				{
					adjointValue += MyPreRow.getVector()[i] * cG0INV_Q[index];
				}
			}
			DkRow[each_col] = adjointValue;
		}
		
		for(unsigned int row = 0 ; row < mCurrentSourceSize ; ++row)  // ak * Dk
		{
			for(unsigned int col = 0 ; col < mCurrentSourceSize ; ++col)
			{
				cHbarmatrix[row + col * mCurrentSourceSize] += DenseAkMatrix[each_grid.first][row] * DkRow[col];
			}
		}

	}
	
	// load (cHbarmatrix + cHbarmatrix transpose) to cHmatrix
	for(unsigned int row = 0 ; row < mCurrentSourceSize ; ++row)
	{
		for(unsigned int col = 0 ; col < mCurrentSourceSize ; ++col)
		{
			unsigned int idx     =  row + col * mCurrentSourceSize;
			unsigned int tranidx =  col + row * mCurrentSourceSize;
			cHmatrix[idx] = (cHbarmatrix[idx] + cHbarmatrix[tranidx]);
		}
	}
	memcpy((void *)mxGetPr(matHmatrix), (void *)cHmatrix, sizeof(double) * mCurrentSourceSize * mCurrentSourceSize);
	engPutVariable(ep, "matHmatrix", matHmatrix);	
	MatrixTime += (double)(clock() - matrixStart) / (CLOCKS_PER_SEC);

	// Run Matlab Quadratic Program
	// engEvalString(ep, "tStart = tic;");
	engEvalString(ep, "tStart = cputime;");
	engEvalString(ep, "[x,fval,exitflag] = quadprog(matHmatrix,matf,matAmatrix,matb,[],[],matlower,matupper);");
	engEvalString(ep, "qptime = cputime - tStart;");	
	// engEvalString(ep, "qptime = toc(tStart);");	
	
	QPtime += *(double *)mxGetPr(engGetVariable(ep, "qptime"));

	//std::cout << -1 * *(double *)mxGetPr(engGetVariable(ep, "fval")) << std::endl;
	//std::cout <<  *(double *)mxGetPr(engGetVariable(ep, "exitflag")) << std::endl;

	// Store Solution
	VdropSolution.emplace_back( -1 * *(double *)mxGetPr(engGetVariable(ep, "fval")) );

	//OnePointThermalAnalysis(each_node_idx, (double *)mxGetPr(engGetVariable(ep, "x")), AkMatrix, 0);

}

}

	std::cout << "Complete Efficient QP Verification !" << std::endl;

	// Output solution 
	unsigned int violations(0);
	
	double sum(0);
	for(auto& vdrop : VdropSolution) { sum +=vdrop; }

	std::ofstream fout("../LPresult/" + mCktName + "_EfficientQPVdropMap");
	
	int sol_idx(0);
	
	for(unsigned int i = 0 ; i < mPGEquivalentNodeSize ; ++i)
	{
		if(mPGEquivalentNodes[i]->is_connect_I)
		{
			fout << mPGEquivalentNodes[i]->x_loc_idx << " " << mPGEquivalentNodes[i]->y_loc_idx 
				<< " " << std::setprecision(15) << VdropSolution[sol_idx] << std::endl;
	
			if(VdropSolution[sol_idx] > mSupplyVoltage * 0.1)
			{
				++violations;
			}
			++sol_idx;
		}
	}
	fout.close();

	fout.open("../LPresult/" + mCktName + "_EfficientQP_TimeInfo");
	fout << "Matrix Time = " << std::setprecision(15) << MatrixTime << " (S) " << std::endl;
	fout << "Matlab QP Time = " << std::setprecision(15) << QPtime << " (S) " << std::endl;
	fout << "Total time (Matrix + Matlab) = " << std::setprecision(15) << MatrixTime + QPtime << " (S) " << std::endl;
	fout << "Total verification time = " << mConstructDataTime + Dmatrixtime + MatrixTime + QPtime << " (S) " << std::endl;
	fout << "Average voltage drop = " << std::setprecision(15) << sum / VdropSolution.size() << " (V) " << std::endl;
	fout << "Maximum voltage drop = " << std::setprecision(15) << *std::max_element(VdropSolution.begin(), VdropSolution.end()) << " (V) " << std::endl;
	fout << "Total violations = " << violations << std::endl;
	fout.close();

	// Delete memory
	mxDestroyArray(matHmatrix);
	mxDestroyArray(matf);
	mxDestroyArray(matAmatrix);
	mxDestroyArray(matb);
	mxDestroyArray(matlower);
	mxDestroyArray(matupper);
	engClose(ep);

	delete [] cG0INV_Q;
	delete [] Amatrix;
	delete [] b;
	delete [] clower;
	delete [] cupper;
	delete [] cf;
	delete [] cHbarmatrix;
	delete [] cHmatrix;

}



// Efficient QP with IR drop recalculation (Combined with Efficient QP without IR drop recalculation)
void PGSolver::PGVerificationQP_Efficient_Simulation()
{

    double Dmatrixtime(0);

	clock_t start;
	
	start = clock();

	// in order to perform sparse vector dot product
	std::map<unsigned int, MySparseMatrixClass::SparseMatrix> MyPGConductanceSubMatrix;	 
	for(auto& each_conductance_submatrix : mPG_EigenConductanceSubMatrix)
	{
		(each_conductance_submatrix.second).makeCompressed();	

		int    *aptrb    = (each_conductance_submatrix.second).outerIndexPtr();	  
		int    *asub     = (each_conductance_submatrix.second).innerIndexPtr();	 
		double *aval     = (each_conductance_submatrix.second).valuePtr();		
		unsigned int NonZerosIndexAdder(0);

		MySparseMatrixClass::SparseMatrix tmpSparseMatrix;
		for(unsigned int Col = 0 ; Col < mPGEquivalentNodeSize ; ++Col)
		{
			MySparseMatrixClass::SparseVector tmpSparseVector((aptrb[Col+1] - aptrb[Col]), asub, aval, &NonZerosIndexAdder);
			tmpSparseMatrix.PushBack(tmpSparseVector);
		}
		MyPGConductanceSubMatrix[each_conductance_submatrix.first] = tmpSparseMatrix;
	}
	mPG_EigenConductanceSubMatrix.clear();	

	// Construct G0_INV * Q (Dense Matrix)
	double *cG0INV_Q = new double[mPGEquivalentNodeSize * mCurrentSourceSize]();
	for(unsigned int each_row = 0 ; each_row < mPGEquivalentNodeSize ; ++each_row)
	{
		for(unsigned int each_col = 0 ; each_col < mCurrentSourceSize ; ++each_col)
		{
			double value = mEachRowPGInverseMatrix[each_row][mCurrentSourceIdx[each_col]];
			if(value != 0)	// positively dense
			{
				cG0INV_Q[each_row + each_col*mPGEquivalentNodeSize] = value;
			}
		}
	}

	// Construct INVGT * M_EtoT, find the grid index to calculate the temperature, ak
	EigenDenseMatrix AkMatrix = mDenseInverseThermalMatrix * Eigen::MatrixXd(mElectricalGridtoThermalGridIncidenceMatrix);
	
	// load AkMatrix to sparse vector data structure
	unsigned int thermalgridsize = mXCut * mYCut * mZlayer;
	unsigned int ArraySize = thermalgridsize * mPGEquivalentNodeSize;
	double *cAkMatrix = new double[ArraySize]();
	Eigen::Map<Eigen::MatrixXd>( cAkMatrix, AkMatrix.rows(), AkMatrix.cols() ) = AkMatrix;
	std::vector< std::vector<double> > DenseAkMatrix(thermalgridsize);	// Q * Ak

	// Factor
	double sFactor = -1 * mR_temp_coeff * mSupplyVoltage;	// -1 * alpha * Vdd 	

	for(unsigned int i=0 ; i<thermalgridsize ; ++i )
		DenseAkMatrix[i].resize(mCurrentSourceSize);
	
	for(unsigned int row = 0 ; row < thermalgridsize ; ++row)
	{
		for(unsigned int col = 0 ; col < mPGEquivalentNodeSize ; ++col)
		{
			unsigned int index = row + col * thermalgridsize;
			if(cAkMatrix[index] != 0)
			{
				DenseAkMatrix[row][mCurrentSourceIdxMapping[col]] = sFactor * cAkMatrix[index];
			}
		}
	}
	Dmatrixtime = (double)(clock()-start) / (CLOCKS_PER_SEC);
	std::cout << "Construct D matrix time = " << Dmatrixtime << std::endl;

	// start matlab engine

	Engine *ep;

	if (!(ep = engOpen(""))) {
		fprintf(stderr, "\nCan't start MATLAB engine\n");
		return ;
	}
	else {
		std::cout << "Started matlab engine !" << std::endl;
	}

	// engEvalString(ep, "LASTN = maxNumCompThreads(1);");

	// construct matlab QP A matrix
	//-----------------------------------------------------------------------------------------------------------
	// load mMatlabAmatrix & lower bound constraint to Amatrix  => 0 < Ui < Ig
	unsigned int Amatrixsize = 2 * mGlobalConstraintSize * mCurrentSourceSize;
	double *Amatrix = new double[Amatrixsize];
	Eigen::Map<Eigen::MatrixXd>( Amatrix, mMatlabAmatrix.rows(), mMatlabAmatrix.cols() ) = Eigen::MatrixXd(mMatlabAmatrix);
    
  	// load Amatrix to matAmatrix
	mxArray *matAmatrix = NULL;
	matAmatrix = mxCreateDoubleMatrix(2 * mGlobalConstraintSize, mCurrentSourceSize, mxREAL);
	memcpy((void *)mxGetPr(matAmatrix), (void *)Amatrix, sizeof(double) * Amatrixsize);
	engPutVariable(ep, "matAmatrix", matAmatrix);	// put variable to matlab workspace 
	//-----------------------------------------------------------------------------------------------------------
	
	// construct matlab QP b vector
	double *b = new double[2*mGlobalConstraintSize]();
	for(unsigned int i=0 ; i < mGlobalConstraintSize ; ++i)
	{
		b[i] = mGlobalCurrentConstraint[i];
	}
	mxArray *matb = NULL;
	matb = mxCreateDoubleMatrix(1, 2 * mGlobalConstraintSize, mxREAL);
	memcpy((void *)mxGetPr(matb), (void *)b, sizeof(double) * 2 * mGlobalConstraintSize);
	engPutVariable(ep, "matb", matb);

	// construct matlab QP local upper and lower variable constraints => 0 < i < Il
	//-----------------------------------------------------------------------------------------------------------
	double *clower = new double[mCurrentSourceSize]();
	double *cupper = new double[mCurrentSourceSize];
	for(unsigned int i=0 ; i < mCurrentSourceSize ; ++i)
	{
		cupper[i] = mLocalCurrentConstraint[i];
	}
	mxArray *matlower = NULL;
	mxArray *matupper = NULL;
	matlower = mxCreateDoubleMatrix(1, mCurrentSourceSize, mxREAL);
	matupper = mxCreateDoubleMatrix(1, mCurrentSourceSize, mxREAL);
	memcpy((void *)mxGetPr(matlower), (void *)clower, sizeof(double) * mCurrentSourceSize);
	memcpy((void *)mxGetPr(matupper), (void *)cupper, sizeof(double) * mCurrentSourceSize);
	engPutVariable(ep, "matlower", matlower);
	engPutVariable(ep, "matupper", matupper);
	//-----------------------------------------------------------------------------------------------------------
	
	//-----------------------------------------------------------------------------------------------------------
	// construct matlab linear terms 
	double *cf = new double[mCurrentSourceSize];
	mxArray *matf = NULL;
	matf = mxCreateDoubleMatrix(1, mCurrentSourceSize, mxREAL);

	// H matrix
	mxArray *matHmatrix = NULL;
	matHmatrix = mxCreateDoubleMatrix(mCurrentSourceSize, mCurrentSourceSize, mxREAL);
	double *cHbarmatrix = new double[mCurrentSourceSize*mCurrentSourceSize]();
	double *cHmatrix = new double[mCurrentSourceSize*mCurrentSourceSize]();
	

	

	// solution
	std::vector<double> VdropSolution, VdropSolution_Origin;

	
	engEvalString(ep, "qptime = 0;");	// for calculating matlab quadprog time 
	clock_t matrixStart, simulationStart; 
	double QPtime(0), MatrixTime(0), SimulationTime(0);

	EigenDenseVector DkRow(mCurrentSourceSize);

    //int debug_idx(0);

for(unsigned int each_node_idx = 0 ; each_node_idx < mPGEquivalentNodeSize ; ++each_node_idx)
{
	
if(mPGEquivalentNodes[each_node_idx]->is_connect_I)
{

	// construct matlab f vector
	matrixStart = clock();

	for(unsigned int i=0 ; i<mCurrentSourceSize ; ++i)
	{
		cf[i] = -1 * mEachRowPGInverseMatrix[each_node_idx][mCurrentSourceIdx[i]];
	}
	memcpy((void *)mxGetPr(matf), (void *)cf, sizeof(double) * mCurrentSourceSize);
	engPutVariable(ep, "matf", matf);	
	//std::cout << "row time = " << std::setprecision(15) << (double)(clock()-start) / (CLOCKS_PER_SEC) << std::endl;
	
	// Construct cHbarmatrix
	for(unsigned int i = 0 ; i < mCurrentSourceSize * mCurrentSourceSize ; ++i)
	{
		cHbarmatrix[i] = 0;
	}

	for(auto& each_grid : MyPGConductanceSubMatrix)		
	{
		MySparseMatrixClass::SparseVector MyPreRow; 	// G0INV * G0,k

		for(unsigned int each_col = 0 ; each_col < mPGEquivalentNodeSize ; ++each_col)
		{
			double adjointValue(0);
			for(unsigned int i = 0 ; i < each_grid.second.getSparseVector(each_col).getVector().size() ; ++i)
			{
				adjointValue += each_grid.second.getSparseVector(each_col).getVector()[i] * mEachRowPGInverseMatrix[each_node_idx][each_grid.second.getSparseVector(each_col).getLocation()[i]];
			}
			if(adjointValue != 0)
			{
				MyPreRow.PushBackElement(each_col, adjointValue);
			}
		}

		for(unsigned int each_col = 0 ; each_col < mCurrentSourceSize ; ++each_col)	  //   Dk = (G0INV * G0,k) * (G0INV * Q)
		{
			double adjointValue(0);
			for(unsigned int i = 0 ; i < MyPreRow.getVector().size() ; ++i)
			{
				unsigned int index = MyPreRow.getLocation()[i] + each_col * mPGEquivalentNodeSize;
				if( cG0INV_Q[index] != 0)
				{
					adjointValue += MyPreRow.getVector()[i] * cG0INV_Q[index];
				}
			}
			DkRow[each_col] = adjointValue;
		}
		
		for(unsigned int row = 0 ; row < mCurrentSourceSize ; ++row)  // ak * Dk
		{
			for(unsigned int col = 0 ; col < mCurrentSourceSize ; ++col)
			{
				cHbarmatrix[row + col * mCurrentSourceSize] += DenseAkMatrix[each_grid.first][row] * DkRow[col];
			}
		}

	}
	
	// load (cHbarmatrix + cHbarmatrix transpose) to cHmatrix
	for(unsigned int row = 0 ; row < mCurrentSourceSize ; ++row)
	{
		for(unsigned int col = 0 ; col < mCurrentSourceSize ; ++col)
		{
			unsigned int idx     =  row + col * mCurrentSourceSize;
			unsigned int tranidx =  col + row * mCurrentSourceSize;
			cHmatrix[idx] = (cHbarmatrix[idx] + cHbarmatrix[tranidx]);
		}
	}
	memcpy((void *)mxGetPr(matHmatrix), (void *)cHmatrix, sizeof(double) * mCurrentSourceSize * mCurrentSourceSize);
	engPutVariable(ep, "matHmatrix", matHmatrix);	
	MatrixTime += (double)(clock() - matrixStart) / (CLOCKS_PER_SEC);

	// Run Matlab Quadratic Program
	// engEvalString(ep, "tStart = tic;");
	engEvalString(ep, "tStart = cputime;");
	engEvalString(ep, "[x,fval,exitflag] = quadprog(matHmatrix,matf,matAmatrix,matb,[],[],matlower,matupper);");
	engEvalString(ep, "qptime = cputime - tStart;");	
	// engEvalString(ep, "qptime = toc(tStart);");	
	
	QPtime += *(double *)mxGetPr(engGetVariable(ep, "qptime"));

	VdropSolution_Origin.emplace_back(-1 * *(double *)mxGetPr(engGetVariable(ep, "fval")));

	
	// DC Simulation
	simulationStart = clock();

	double *xx = (double *)mxGetPr(engGetVariable(ep, "x"));		// current pattern solution of QP

	EigenDenseVector CurrentSourceVector(mPGEquivalentNodeSize), TemperatureVector(thermalgridsize);
	CurrentSourceVector.setZero();
	for(unsigned int i = 0 ; i < mCurrentSourceSize ; ++i)
	{
		CurrentSourceVector[mCurrentSourceIdx[i]] = xx[i];
	}
	TemperatureVector = AkMatrix * CurrentSourceVector;	// Obtain temperature vector

	// Stamping
	MATRIX tmpMatrix;
	tmpMatrix.resize(mPGEquivalentNodeSize);

	for (auto& stamp_node : mStampNodeList)
	{
		for (auto& nei_node : stamp_node.first->neighbor_NodeRes)    
		{
			if (nei_node.first->is_grounded == true)
			{
				for (auto& each_res : nei_node.second)
				{
					tmpMatrix[stamp_node.second][stamp_node.second] += 1 / (each_res->res_value * (1 + mR_temp_coeff * TemperatureVector[each_res->thermal_grid_index]));
				}
			}
			else if (nei_node.first->is_connect_V == true)
			{
				for (auto& each_res : nei_node.second)
				{				
					tmpMatrix[stamp_node.second][stamp_node.second] += 1 / (each_res->res_value * (1 + mR_temp_coeff * TemperatureVector[each_res->thermal_grid_index]));
				}
			}
			else if (nei_node.first->is_connect_V == false && nei_node.first->is_grounded == false)
			{
				for (auto& each_res : nei_node.second)
				{
					tmpMatrix[stamp_node.second][stamp_node.second] += 1 / (each_res->res_value * (1 + mR_temp_coeff * TemperatureVector[each_res->thermal_grid_index]));
					tmpMatrix[stamp_node.second][mStampNodeList[nei_node.first]] += -1 / (each_res->res_value * (1 + mR_temp_coeff * TemperatureVector[each_res->thermal_grid_index]));
				}
			}
		}
	}
	std::vector<triplet> tripletList;
	for(unsigned int row_idx = 0; row_idx < mPGEquivalentNodeSize; ++row_idx)
    {
        for(auto& row_data : tmpMatrix[row_idx])
        {           
            tripletList.push_back(triplet(row_idx, row_data.first, row_data.second));   
        }
    }
    Eigen::SparseMatrix<double> sys_mat(mPGEquivalentNodeSize, mPGEquivalentNodeSize);
	sys_mat.setFromTriplets(tripletList.begin(), tripletList.end());
	tripletList.clear();

	// Solve
	Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
	solver.analyzePattern(sys_mat);   
	solver.factorize(sys_mat);
	EigenDenseVector VdropVector = solver.solve(CurrentSourceVector);
	SimulationTime += (double)(clock() - simulationStart) / (CLOCKS_PER_SEC);
	VdropSolution.emplace_back(VdropVector[each_node_idx]);
}

}

	std::cout << "Complete Efficient QP Verification !" << std::endl;

	// Output solution 
	unsigned int violations(0);
	double sum(0);
	for(auto& vdrop : VdropSolution) { sum +=vdrop; }

	std::ofstream fout("../LPresult/" + mCktName + "_EfficientQPVdropMap_Simulation");
	
	int sol_idx(0);
	
	for(unsigned int i = 0 ; i < mPGEquivalentNodeSize ; ++i)
	{
		if(mPGEquivalentNodes[i]->is_connect_I)
		{
			fout << mPGEquivalentNodes[i]->x_loc_idx << " " << mPGEquivalentNodes[i]->y_loc_idx 
				<< " " << std::setprecision(15) << VdropSolution[sol_idx] << std::endl;
	
			if(VdropSolution[sol_idx] > mSupplyVoltage * 0.1)
			{
				++violations;
			}
			++sol_idx;
		}
	}
	fout.close();

	fout.open("../LPresult/" + mCktName + "_EfficientQP_TimeInfo_Simulation");
	fout << "Matrix Time = " << std::setprecision(15) << MatrixTime << " (S) " << std::endl;
	fout << "Matlab QP Time = " << std::setprecision(15) << QPtime << " (S) " << std::endl;
	fout << "Total time (Matrix + Matlab) = " << std::setprecision(15) << MatrixTime + QPtime << " (S) " << std::endl;
	fout << "Total verification time = " << mConstructDataTime + Dmatrixtime + MatrixTime + QPtime + SimulationTime << " (S) " << std::endl;
	fout << "Average voltage drop = " << std::setprecision(15) << sum / VdropSolution.size() << " (V) " << std::endl;
	fout << "Maximum voltage drop = " << std::setprecision(15) << *std::max_element(VdropSolution.begin(), VdropSolution.end()) << " (V) " << std::endl;
	fout << "Total violations = " << violations << std::endl;
	fout << "Simulation Time = " << SimulationTime << " (s) " << std::endl;
	fout.close();

	sum = 0;
	for(auto& vdrop : VdropSolution_Origin) { sum +=vdrop; }

	fout.open("../LPresult/" + mCktName + "_EfficientQPVdropMap");
	sol_idx = 0;
	violations = 0;
	for(unsigned int i = 0 ; i < mPGEquivalentNodeSize ; ++i)
	{
		if(mPGEquivalentNodes[i]->is_connect_I)
		{
			fout << mPGEquivalentNodes[i]->x_loc_idx << " " << mPGEquivalentNodes[i]->y_loc_idx 
				<< " " << std::setprecision(15) << VdropSolution_Origin[sol_idx] << std::endl;
	
			if(VdropSolution_Origin[sol_idx] > mSupplyVoltage * 0.1)
			{
				++violations;
			}
			++sol_idx;
		}
	}
	fout.close();

	fout.open("../LPresult/" + mCktName + "_EfficientQP_TimeInfo");
	fout << "Matrix Time = " << std::setprecision(15) << MatrixTime << " (S) " << std::endl;
	fout << "Matlab QP Time = " << std::setprecision(15) << QPtime << " (S) " << std::endl;
	fout << "Total time (Matrix + Matlab) = " << std::setprecision(15) << MatrixTime + QPtime << " (S) " << std::endl;
	fout << "Total verification time = " << mConstructDataTime + Dmatrixtime + MatrixTime + QPtime << " (S) " << std::endl;
	fout << "Average voltage drop = " << std::setprecision(15) << sum / VdropSolution_Origin.size() << " (V) " << std::endl;
	fout << "Maximum voltage drop = " << std::setprecision(15) << *std::max_element(VdropSolution_Origin.begin(), VdropSolution_Origin.end()) << " (V) " << std::endl;
	fout << "Total violations = " << violations << std::endl;
	fout.close();

	// Delete memory
	mxDestroyArray(matHmatrix);
	mxDestroyArray(matf);
	mxDestroyArray(matAmatrix);
	mxDestroyArray(matb);
	mxDestroyArray(matlower);
	mxDestroyArray(matupper);
	engClose(ep);

	delete [] cG0INV_Q;
	delete [] Amatrix;
	delete [] b;
	delete [] clower;
	delete [] cupper;
	delete [] cf;
	delete [] cHbarmatrix;
	delete [] cHmatrix;

}

// 2018.5.18
// 做兩次QP,第二次QP會在第一次QP找到的solution current pattern所造成的溫度之下展開 (not in thesis)
void PGSolver::PGVerificationQP_Efficient_Simulation_Verify()
{

	double Dmatrixtime(0);

	clock_t start;
	
	start = clock();

	// in order to perform sparse vector dot product
	std::map<unsigned int, MySparseMatrixClass::SparseMatrix> MyPGConductanceSubMatrix;	 
	for(auto& each_conductance_submatrix : mPG_EigenConductanceSubMatrix)
	{
		(each_conductance_submatrix.second).makeCompressed();	

		int    *aptrb    = (each_conductance_submatrix.second).outerIndexPtr();	  
		int    *asub     = (each_conductance_submatrix.second).innerIndexPtr();	 
		double *aval     = (each_conductance_submatrix.second).valuePtr();		
		unsigned int NonZerosIndexAdder(0);

		MySparseMatrixClass::SparseMatrix tmpSparseMatrix;
		for(unsigned int Col = 0 ; Col < mPGEquivalentNodeSize ; ++Col)
		{
			MySparseMatrixClass::SparseVector tmpSparseVector((aptrb[Col+1] - aptrb[Col]), asub, aval, &NonZerosIndexAdder);
			tmpSparseMatrix.PushBack(tmpSparseVector);
		}
		MyPGConductanceSubMatrix[each_conductance_submatrix.first] = tmpSparseMatrix;
	}
	// mPG_EigenConductanceSubMatrix.clear();	// clear origin Eigen matrix

	// Construct G0_INV * Q (Dense Matrix)
	double *cG0INV_Q = new double[mPGEquivalentNodeSize * mCurrentSourceSize]();
	for(unsigned int each_row = 0 ; each_row < mPGEquivalentNodeSize ; ++each_row)
	{
		for(unsigned int each_col = 0 ; each_col < mCurrentSourceSize ; ++each_col)
		{
			double value = mEachRowPGInverseMatrix[each_row][mCurrentSourceIdx[each_col]];
			if(value != 0)	// positively dense
			{
				cG0INV_Q[each_row + each_col*mPGEquivalentNodeSize] = value;
			}
		}
	}

	// Construct INVGT * M_EtoT, find the grid index to calculate the temperature, ak
	EigenDenseMatrix AkMatrix = mDenseInverseThermalMatrix * Eigen::MatrixXd(mElectricalGridtoThermalGridIncidenceMatrix);
	
	// load AkMatrix to sparse vector data structure
	unsigned int thermalgridsize = mXCut * mYCut * mZlayer;
	unsigned int ArraySize = thermalgridsize * mPGEquivalentNodeSize;
	double *cAkMatrix = new double[ArraySize]();
	Eigen::Map<Eigen::MatrixXd>( cAkMatrix, AkMatrix.rows(), AkMatrix.cols() ) = AkMatrix;
	std::vector< std::vector<double> > DenseAkMatrix(thermalgridsize);	// Q * Ak

	for(unsigned int i=0 ; i<thermalgridsize ; ++i )
		DenseAkMatrix[i].resize(mCurrentSourceSize);
	
	for(unsigned int row = 0 ; row < thermalgridsize ; ++row)
	{
		for(unsigned int col = 0 ; col < mPGEquivalentNodeSize ; ++col)
		{
			unsigned int index = row + col * thermalgridsize;
			if(cAkMatrix[index] != 0)
			{
				DenseAkMatrix[row][mCurrentSourceIdxMapping[col]] = cAkMatrix[index];
			}
		}
	}
	Dmatrixtime = (double)(clock()-start) / (CLOCKS_PER_SEC);
	std::cout << "Construct D matrix time = " << Dmatrixtime << std::endl;

	// start matlab engine

	Engine *ep;

	if (!(ep = engOpen(""))) {
		fprintf(stderr, "\nCan't start MATLAB engine\n");
		return ;
	}
	else {
		std::cout << "Started matlab engine !" << std::endl;
	}

	// engEvalString(ep, "LASTN = maxNumCompThreads(1);");

	// construct matlab QP A matrix
	//-----------------------------------------------------------------------------------------------------------
	// load mMatlabAmatrix & lower bound constraint to Amatrix  => 0 < Ui < Ig
	unsigned int Amatrixsize = 2 * mGlobalConstraintSize * mCurrentSourceSize;
	double *Amatrix = new double[Amatrixsize];
	Eigen::Map<Eigen::MatrixXd>( Amatrix, mMatlabAmatrix.rows(), mMatlabAmatrix.cols() ) = Eigen::MatrixXd(mMatlabAmatrix);
    
  	// load Amatrix to matAmatrix
	mxArray *matAmatrix = NULL;
	matAmatrix = mxCreateDoubleMatrix(2 * mGlobalConstraintSize, mCurrentSourceSize, mxREAL);
	memcpy((void *)mxGetPr(matAmatrix), (void *)Amatrix, sizeof(double) * Amatrixsize);
	engPutVariable(ep, "matAmatrix", matAmatrix);	// put variable to matlab workspace 
	//-----------------------------------------------------------------------------------------------------------
	
	// construct matlab QP b vector
	double *b = new double[2*mGlobalConstraintSize]();
	for(unsigned int i=0 ; i < mGlobalConstraintSize ; ++i)
	{
		b[i] = mGlobalCurrentConstraint[i];
	}
	mxArray *matb = NULL;
	matb = mxCreateDoubleMatrix(1, 2 * mGlobalConstraintSize, mxREAL);
	memcpy((void *)mxGetPr(matb), (void *)b, sizeof(double) * 2 * mGlobalConstraintSize);
	engPutVariable(ep, "matb", matb);

	// construct matlab QP local upper and lower variable constraints => 0 < i < Il
	//-----------------------------------------------------------------------------------------------------------
	double *clower = new double[mCurrentSourceSize]();
	double *cupper = new double[mCurrentSourceSize];
	for(unsigned int i=0 ; i < mCurrentSourceSize ; ++i)
	{
		cupper[i] = mLocalCurrentConstraint[i];
	}
	mxArray *matlower = NULL;
	mxArray *matupper = NULL;
	matlower = mxCreateDoubleMatrix(1, mCurrentSourceSize, mxREAL);
	matupper = mxCreateDoubleMatrix(1, mCurrentSourceSize, mxREAL);
	memcpy((void *)mxGetPr(matlower), (void *)clower, sizeof(double) * mCurrentSourceSize);
	memcpy((void *)mxGetPr(matupper), (void *)cupper, sizeof(double) * mCurrentSourceSize);
	engPutVariable(ep, "matlower", matlower);
	engPutVariable(ep, "matupper", matupper);
	//-----------------------------------------------------------------------------------------------------------
	
	//-----------------------------------------------------------------------------------------------------------
	// construct matlab linear terms 
	double *cf = new double[mCurrentSourceSize];
	mxArray *matf = NULL;
	matf = mxCreateDoubleMatrix(1, mCurrentSourceSize, mxREAL);

	// H matrix
	mxArray *matHmatrix = NULL;
	matHmatrix = mxCreateDoubleMatrix(mCurrentSourceSize, mCurrentSourceSize, mxREAL);
	double *cHbarmatrix = new double[mCurrentSourceSize*mCurrentSourceSize]();
	double *cHmatrix = new double[mCurrentSourceSize*mCurrentSourceSize]();
	

	// Factor
	double sFactor = -1 * mR_temp_coeff * mSupplyVoltage;	// -1 * alpha * Vdd 	

	// solution
	std::vector<double> VdropSolution;
	std::vector<double> CurrentErrorNorm;
	std::vector<double> CurrentRelativeErrorNorm;
	std::vector<double> CurrentErrorMaxElement;

	
	engEvalString(ep, "qptime = 0;");	// for calculating matlab quadprog time 
	clock_t matrixStart, simulationStart; 
	double QPtime(0), MatrixTime(0), SimulationTime(0);

	EigenDenseVector DkRow(mCurrentSourceSize);
	double *tmp_cG0INV_Q = new double[mPGEquivalentNodeSize * mCurrentSourceSize]();
    //int debug_idx(0);

    std::vector<unsigned int> veri_idx(17);

// PG1
/*
    veri_idx[0] = 79;
    veri_idx[1] = 159;
*/
// PG2
/*
    veri_idx[0] = 99;
    veri_idx[1] = 199;
*/
// PG3 (V2)
/*
veri_idx[0] = 139;
veri_idx[1] =144;
veri_idx[2] =149;
veri_idx[3] =289;
veri_idx[4] =294;
veri_idx[5] =439;
veri_idx[6] =444;
veri_idx[7] =589;
veri_idx[8] =594;
veri_idx[9] =744;
veri_idx[10] =894;
veri_idx[11] =1049;   
*/
// PG4
    
veri_idx[0] = 139;
veri_idx[1] =144;
veri_idx[2] =149;
veri_idx[3] =289;
veri_idx[4] =294;
veri_idx[5] =439;
veri_idx[6] =444;
veri_idx[7] =589;
veri_idx[8] =594;
veri_idx[9] =739;
veri_idx[10] =744;
veri_idx[11] =889;
veri_idx[12] =894;
veri_idx[13] =1044;
veri_idx[14] =1049;
veri_idx[15] =1199;
veri_idx[16] =1349;


//for(unsigned int each_node_idx = 0 ; each_node_idx < mPGEquivalentNodeSize ; ++each_node_idx)
//for(unsigned int each_node_idx = 5050 ; each_node_idx < 5050+1 ; ++each_node_idx)
for(auto &each_node_idx : veri_idx)
{
	
//if(each_node_idx == veri_idx[0] || each_node_idx == veri_idx[1])
//{

	// construct matlab f vector
	matrixStart = clock();

	for(unsigned int i=0 ; i<mCurrentSourceSize ; ++i)
	{
		cf[i] = -1 * mEachRowPGInverseMatrix[each_node_idx][mCurrentSourceIdx[i]];
	}
	memcpy((void *)mxGetPr(matf), (void *)cf, sizeof(double) * mCurrentSourceSize);
	engPutVariable(ep, "matf", matf);	
	//std::cout << "row time = " << std::setprecision(15) << (double)(clock()-start) / (CLOCKS_PER_SEC) << std::endl;
	
	// Construct cHbarmatrix
	for(unsigned int i = 0 ; i < mCurrentSourceSize * mCurrentSourceSize ; ++i)
	{
		cHbarmatrix[i] = 0;
	}

	for(auto& each_grid : MyPGConductanceSubMatrix)		
	{
		MySparseMatrixClass::SparseVector MyPreRow; 	// G0INV * G0,k

		for(unsigned int each_col = 0 ; each_col < mPGEquivalentNodeSize ; ++each_col)
		{
			double adjointValue(0);
			for(unsigned int i = 0 ; i < each_grid.second.getSparseVector(each_col).getVector().size() ; ++i)
			{
				adjointValue += each_grid.second.getSparseVector(each_col).getVector()[i] * mEachRowPGInverseMatrix[each_node_idx][each_grid.second.getSparseVector(each_col).getLocation()[i]];
			}
			if(adjointValue != 0)
			{
				MyPreRow.PushBackElement(each_col, adjointValue);
			}
		}

		for(unsigned int each_col = 0 ; each_col < mCurrentSourceSize ; ++each_col)	  //   Dk = (G0INV * G0,k) * (G0INV * Q)
		{
			double adjointValue(0);
			for(unsigned int i = 0 ; i < MyPreRow.getVector().size() ; ++i)
			{
				unsigned int index = MyPreRow.getLocation()[i] + each_col * mPGEquivalentNodeSize;
				if( cG0INV_Q[index] != 0)
				{
					adjointValue += MyPreRow.getVector()[i] * cG0INV_Q[index];
				}
			}
			DkRow[each_col] = adjointValue;
		}
		
		for(unsigned int row = 0 ; row < mCurrentSourceSize ; ++row)  // ak * Dk
		{
			for(unsigned int col = 0 ; col < mCurrentSourceSize ; ++col)
			{
				cHbarmatrix[row + col * mCurrentSourceSize] += DenseAkMatrix[each_grid.first][row] * DkRow[col];
			}
		}

	}
	
	// load sFactor * (cHbarmatrix + cHbarmatrix transpose) to cHmatrix
	for(unsigned int row = 0 ; row < mCurrentSourceSize ; ++row)
	{
		for(unsigned int col = 0 ; col < mCurrentSourceSize ; ++col)
		{
			unsigned int idx     =  row + col * mCurrentSourceSize;
			unsigned int tranidx =  col + row * mCurrentSourceSize;
			cHmatrix[idx] = sFactor * (cHbarmatrix[idx] + cHbarmatrix[tranidx]);
		}
	}
	memcpy((void *)mxGetPr(matHmatrix), (void *)cHmatrix, sizeof(double) * mCurrentSourceSize * mCurrentSourceSize);
	engPutVariable(ep, "matHmatrix", matHmatrix);	
	MatrixTime += (double)(clock() - matrixStart) / (CLOCKS_PER_SEC);

	// Run Matlab Quadratic Program
	// engEvalString(ep, "tStart = tic;");
	engEvalString(ep, "tStart = cputime;");
	engEvalString(ep, "[x,fval,exitflag] = quadprog(matHmatrix,matf,matAmatrix,matb,[],[],matlower,matupper);");
	engEvalString(ep, "qptime = cputime - tStart;");	
	// engEvalString(ep, "qptime = toc(tStart);");	

	std::cout <<  *(double *)mxGetPr(engGetVariable(ep, "exitflag")) << std::endl;
	
	QPtime += *(double *)mxGetPr(engGetVariable(ep, "qptime"));

	double *origin_current = (double *)mxGetPr(engGetVariable(ep, "x"));
	std::vector<double> origin_current_solution(mCurrentSourceSize);
	double origin_sum(0);
	for(unsigned int i = 0 ; i < mCurrentSourceSize ; ++i)
	{
		origin_current_solution[i] = origin_current[i];
		origin_sum += pow(origin_current_solution[i], 2);
	}

	// DC Simulation
	simulationStart = clock();

	double *xx = (double *)mxGetPr(engGetVariable(ep, "x"));		// current pattern solution of QP

	EigenDenseVector CurrentSourceVector(mPGEquivalentNodeSize), TemperatureVector(thermalgridsize);
	CurrentSourceVector.setZero();
	for(unsigned int i = 0 ; i < mCurrentSourceSize ; ++i)
	{
		CurrentSourceVector[mCurrentSourceIdx[i]] = xx[i];
	}
	TemperatureVector = AkMatrix * CurrentSourceVector;	// Obtain temperature vector

	// Stamping
	MATRIX tmpMatrix;
	tmpMatrix.resize(mPGEquivalentNodeSize);

	for (auto& stamp_node : mStampNodeList)
	{
		for (auto& nei_node : stamp_node.first->neighbor_NodeRes)    
		{
			if (nei_node.first->is_grounded == true)
			{
				for (auto& each_res : nei_node.second)
				{
					tmpMatrix[stamp_node.second][stamp_node.second] += 1 / (each_res->res_value * (1 + mR_temp_coeff * TemperatureVector[each_res->thermal_grid_index]));
				}
			}
			else if (nei_node.first->is_connect_V == true)
			{
				for (auto& each_res : nei_node.second)
				{				
					tmpMatrix[stamp_node.second][stamp_node.second] += 1 / (each_res->res_value * (1 + mR_temp_coeff * TemperatureVector[each_res->thermal_grid_index]));
				}
			}
			else if (nei_node.first->is_connect_V == false && nei_node.first->is_grounded == false)
			{
				for (auto& each_res : nei_node.second)
				{
					tmpMatrix[stamp_node.second][stamp_node.second] += 1 / (each_res->res_value * (1 + mR_temp_coeff * TemperatureVector[each_res->thermal_grid_index]));
					tmpMatrix[stamp_node.second][mStampNodeList[nei_node.first]] += -1 / (each_res->res_value * (1 + mR_temp_coeff * TemperatureVector[each_res->thermal_grid_index]));
				}
			}
		}
	}
	std::vector<triplet> tripletList;
	for(unsigned int row_idx = 0; row_idx < mPGEquivalentNodeSize; ++row_idx)
    {
        for(auto& row_data : tmpMatrix[row_idx])
        {           
            tripletList.push_back(triplet(row_idx, row_data.first, row_data.second));   
        }
    }
    Eigen::SparseMatrix<double> sys_mat(mPGEquivalentNodeSize, mPGEquivalentNodeSize);
	sys_mat.setFromTriplets(tripletList.begin(), tripletList.end());
	tripletList.clear();

	// Solve inverse matrix
	Eigen::VectorXd pattern(mPGEquivalentNodeSize);
	Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
	solver.analyzePattern(sys_mat);   
	solver.factorize(sys_mat);
	pattern.setZero();
	std::vector<EigenDenseVector> tmp_EachRowPGInverseMatrix;
	for(unsigned int i = 0 ; i < mPGEquivalentNodeSize ; ++i)
	{
		pattern[i] = 1;
		tmp_EachRowPGInverseMatrix.emplace_back(solver.solve(pattern));
		pattern[i] = 0;
	}
	// SubConductance Matrix
	for(unsigned int y = 0; y < mYCut ; ++y)
	{
		for(unsigned int x = 0; x < mXCut ; ++x)
		{
			MATRIX subConductanceMatrix(mPGEquivalentNodeSize);
			// Stamp each subconductance matrix
			for (auto& stamp_node : mStampNodeList)
			{
				for (auto& nei_node : stamp_node.first->neighbor_NodeRes)    
				{
					if (nei_node.first->is_grounded == true)
					{
						for (auto& each_res : nei_node.second)
						{
							if(each_res->is_in_submatrix)
							{
								subConductanceMatrix[stamp_node.second][stamp_node.second] += 1 / (each_res->res_value * ( 1 + mR_temp_coeff * TemperatureVector[each_res->thermal_grid_index])) ;
							}
						}
					}
					else if (nei_node.first->is_connect_V == true)
					{
						for (auto& each_res : nei_node.second)
						{				
							if(each_res->is_in_submatrix)
							{
								subConductanceMatrix[stamp_node.second][stamp_node.second] += 1 / (each_res->res_value * ( 1 + mR_temp_coeff * TemperatureVector[each_res->thermal_grid_index])) ;
							}
						}
					}
					else if (nei_node.first->is_connect_V == false && nei_node.first->is_grounded == false)
					{
						for (auto& each_res : nei_node.second)
						{
							if(each_res->is_in_submatrix)
							{
								subConductanceMatrix[stamp_node.second][stamp_node.second] += 1 / (each_res->res_value * ( 1 + mR_temp_coeff * TemperatureVector[each_res->thermal_grid_index])) ;
								subConductanceMatrix[stamp_node.second][mStampNodeList[nei_node.first]] += -1 / (each_res->res_value * ( 1 + mR_temp_coeff * TemperatureVector[each_res->thermal_grid_index])) ;
							}
						}
					}
				}
			}
			mPG_ConductanceSubMatrix[mTotalThermalGrids[mZlayer-1][y][x].mGridIndex] = subConductanceMatrix; // store
		}
	}
	std::unordered_map<unsigned int, EigenSparseMatrix> tmp_PG_EigenConductanceSubMatrix;
	for(auto& each_matrix : mPG_ConductanceSubMatrix)
	{
		EigenSparseMatrix Mat_PGSubMatrix(mPGEquivalentNodeSize, mPGEquivalentNodeSize);	
		std::vector<triplet> triplet_PGSubMatrix;

		for(unsigned int row_idx = 0; row_idx < each_matrix.second.size(); ++row_idx)
    	{
        	for(auto& row_data : each_matrix.second[row_idx])
        	{           
            	triplet_PGSubMatrix.push_back(triplet(row_idx, row_data.first, row_data.second));   
        	}
    	}
    	Mat_PGSubMatrix.setFromTriplets(triplet_PGSubMatrix.begin(), triplet_PGSubMatrix.end());
    	tmp_PG_EigenConductanceSubMatrix[each_matrix.first] = Mat_PGSubMatrix;
	}
	mPG_ConductanceSubMatrix.clear();

	std::map<unsigned int, MySparseMatrixClass::SparseMatrix> tmp_MyPGConductanceSubMatrix;	 
	for(auto& each_conductance_submatrix : tmp_PG_EigenConductanceSubMatrix)
	{
		each_conductance_submatrix.second.makeCompressed();	
		int    *aptrb    = each_conductance_submatrix.second.outerIndexPtr();	  
		int    *asub     = each_conductance_submatrix.second.innerIndexPtr();	 
		double *aval     = each_conductance_submatrix.second.valuePtr();		
		unsigned int NonZerosIndexAdder(0);

		MySparseMatrixClass::SparseMatrix tmpSparseMatrix;
		for(unsigned int Col = 0 ; Col < mPGEquivalentNodeSize ; ++Col)
		{
			MySparseMatrixClass::SparseVector tmpSparseVector((aptrb[Col+1] - aptrb[Col]), asub, aval, &NonZerosIndexAdder);
			tmpSparseMatrix.PushBack(tmpSparseVector);
		}
		tmp_MyPGConductanceSubMatrix[each_conductance_submatrix.first] = tmpSparseMatrix;
	}
	tmp_PG_EigenConductanceSubMatrix.clear();
	// G0INV_Q
	for(unsigned int each_row = 0 ; each_row < mPGEquivalentNodeSize ; ++each_row)
	{
		for(unsigned int each_col = 0 ; each_col < mCurrentSourceSize ; ++each_col)
		{
			double value = tmp_EachRowPGInverseMatrix[each_row][mCurrentSourceIdx[each_col]];
			if(value != 0)	// positively dense
			{
				tmp_cG0INV_Q[each_row + each_col*mPGEquivalentNodeSize] = value;
			}
		}
	}
	// Solve QP
	for(unsigned int i=0 ; i<mCurrentSourceSize ; ++i)
	{
		cf[i] = -1 * tmp_EachRowPGInverseMatrix[each_node_idx][mCurrentSourceIdx[i]];
	}
	memcpy((void *)mxGetPr(matf), (void *)cf, sizeof(double) * mCurrentSourceSize);
	engPutVariable(ep, "matf", matf);	

	for(unsigned int i = 0 ; i < mCurrentSourceSize * mCurrentSourceSize ; ++i)
	{
		cHbarmatrix[i] = 0;
	}

	for(auto& each_grid : tmp_MyPGConductanceSubMatrix)		
	{
		MySparseMatrixClass::SparseVector MyPreRow; 	// G0INV * G0,k

		for(unsigned int each_col = 0 ; each_col < mPGEquivalentNodeSize ; ++each_col)
		{
			double adjointValue(0);
			for(unsigned int i = 0 ; i < each_grid.second.getSparseVector(each_col).getVector().size() ; ++i)
			{
				adjointValue += each_grid.second.getSparseVector(each_col).getVector()[i] * tmp_EachRowPGInverseMatrix[each_node_idx][each_grid.second.getSparseVector(each_col).getLocation()[i]];
			}
			if(adjointValue != 0)
			{
				MyPreRow.PushBackElement(each_col, adjointValue);
			}
		}

		for(unsigned int each_col = 0 ; each_col < mCurrentSourceSize ; ++each_col)	  //   Dk = (G0INV * G0,k) * (G0INV * Q)
		{
			double adjointValue(0);
			for(unsigned int i = 0 ; i < MyPreRow.getVector().size() ; ++i)
			{
				unsigned int index = MyPreRow.getLocation()[i] + each_col * mPGEquivalentNodeSize;
				if( tmp_cG0INV_Q[index] != 0)
				{
					adjointValue += MyPreRow.getVector()[i] * tmp_cG0INV_Q[index];
				}
			}
			DkRow[each_col] = adjointValue;
		}
		
		for(unsigned int row = 0 ; row < mCurrentSourceSize ; ++row)  // ak * Dk
		{
			for(unsigned int col = 0 ; col < mCurrentSourceSize ; ++col)
			{
				cHbarmatrix[row + col * mCurrentSourceSize] += DenseAkMatrix[each_grid.first][row] * DkRow[col];
			}
		}

	}
	
	// load sFactor * (cHbarmatrix + cHbarmatrix transpose) to cHmatrix
	for(unsigned int row = 0 ; row < mCurrentSourceSize ; ++row)
	{
		for(unsigned int col = 0 ; col < mCurrentSourceSize ; ++col)
		{
			unsigned int idx     =  row + col * mCurrentSourceSize;
			unsigned int tranidx =  col + row * mCurrentSourceSize;
			cHmatrix[idx] = sFactor * (cHbarmatrix[idx] + cHbarmatrix[tranidx]);
		}
	}
	memcpy((void *)mxGetPr(matHmatrix), (void *)cHmatrix, sizeof(double) * mCurrentSourceSize * mCurrentSourceSize);
	engPutVariable(ep, "matHmatrix", matHmatrix);	
	//engEvalString(ep, "tStart = cputime;");
	engEvalString(ep, "[x,fval,exitflag] = quadprog(matHmatrix,matf,matAmatrix,matb,[],[],matlower,matupper);");
	//engEvalString(ep, "qptime = cputime - tStart;");	

	std::cout <<  *(double *)mxGetPr(engGetVariable(ep, "exitflag")) << std::endl;

	VdropSolution.emplace_back(-1 * *(double *)mxGetPr(engGetVariable(ep, "fval")));

	SimulationTime += (double)(clock() - simulationStart) / (CLOCKS_PER_SEC);

	// calculate current error
	double *solution_current = (double *)mxGetPr(engGetVariable(ep, "x"));
	double sum(0);
	std::vector<double> errors;
	for(unsigned int i = 0 ; i < mCurrentSourceSize ; ++i)
	{
		sum += pow(fabs(origin_current_solution[i] - solution_current[i]), 2);
		errors.emplace_back(fabs(origin_current_solution[i] - solution_current[i]));
	}
	CurrentErrorNorm.emplace_back(pow(sum, 0.5));
	CurrentErrorMaxElement.emplace_back(*std::max_element(errors.begin(), errors.end()));
	CurrentRelativeErrorNorm.emplace_back(pow(sum/origin_sum, 0.5) * 100);
//}

}

	std::cout << "Complete Efficient QP Verification !" << std::endl;

	// Output solution 
	std::ofstream fout("../LPresult/" + mCktName + "_EfficientQPVdropMap_Simulation_Verify");

	int sol_idx(0);
	for(auto& idx : veri_idx)
	{
		fout << mPGEquivalentNodes[idx]->x_loc_idx << " " << mPGEquivalentNodes[idx]->y_loc_idx 
			<< " " << std::setprecision(15) << VdropSolution[sol_idx] << std::endl;
		++sol_idx;
	}
	
	/*
	
	for(unsigned int i = 0 ; i < mPGEquivalentNodeSize ; ++i)
	{
		if(mPGEquivalentNodes[i]->is_connect_I)
		{
			fout << mPGEquivalentNodes[i]->x_loc_idx << " " << mPGEquivalentNodes[i]->y_loc_idx 
				<< " " << std::setprecision(15) << VdropSolution[sol_idx] << std::endl;
	
			if(VdropSolution[sol_idx] > mSupplyVoltage * 0.1)
			{
				++violations;
			}
			++sol_idx;
		}
	}
	*/

	fout.close();
	
	fout.open("../LPresult/" + mCktName + "_EfficientQP_CurrentError_Verify");
	for(unsigned int i = 0 ; i < CurrentErrorNorm.size() ; ++i)
	{
		fout << CurrentErrorNorm[i] << " " << CurrentRelativeErrorNorm[i] << " (%) " << CurrentErrorMaxElement[i] << std::endl;
	}
	fout.close();

	/*
	fout.open("../LPresult/" + mCktName + "_EfficientQP_TimeInfo_Simulation_Verify");
	fout << "Matrix Time = " << std::setprecision(15) << MatrixTime << " (S) " << std::endl;
	fout << "Matlab QP Time = " << std::setprecision(15) << QPtime << " (S) " << std::endl;
	fout << "Total time (Matrix + Matlab) = " << std::setprecision(15) << MatrixTime + QPtime + SimulationTime << " (S) " << std::endl;
	fout << "Total verification time = " << mConstructDataTime + Dmatrixtime + MatrixTime + QPtime << " (S) " << std::endl;
	fout << "Average voltage drop = " << std::setprecision(15) << sum / VdropSolution.size() << " (V) " << std::endl;
	fout << "Maximum voltage drop = " << std::setprecision(15) << *std::max_element(VdropSolution.begin(), VdropSolution.end()) << " (V) " << std::endl;
	fout << "Total violations = " << violations << std::endl;
	fout << "Simulation Time = " << SimulationTime << " (s) " << std::endl;
	fout.close();
	*/
	// Delete memory
	mxDestroyArray(matHmatrix);
	mxDestroyArray(matf);
	mxDestroyArray(matAmatrix);
	mxDestroyArray(matb);
	mxDestroyArray(matlower);
	mxDestroyArray(matupper);
	engClose(ep);

	delete [] cG0INV_Q;
	delete [] Amatrix;
	delete [] b;
	delete [] clower;
	delete [] cupper;
	delete [] cf;
	delete [] cHbarmatrix;
	delete [] cHmatrix;
	delete [] tmp_cG0INV_Q;

}

// 2018.5.24
void PGSolver::PGVerificationQP_Efficient_Simulation_Verify_Iterative()
{

	double Dmatrixtime(0);

	clock_t start;
	
	start = clock();

	// in order to perform sparse vector dot product
	std::map<unsigned int, MySparseMatrixClass::SparseMatrix> MyPGConductanceSubMatrix;	 
	for(auto& each_conductance_submatrix : mPG_EigenConductanceSubMatrix)
	{
		(each_conductance_submatrix.second).makeCompressed();	

		int    *aptrb    = (each_conductance_submatrix.second).outerIndexPtr();	  
		int    *asub     = (each_conductance_submatrix.second).innerIndexPtr();	 
		double *aval     = (each_conductance_submatrix.second).valuePtr();		
		unsigned int NonZerosIndexAdder(0);

		MySparseMatrixClass::SparseMatrix tmpSparseMatrix;
		for(unsigned int Col = 0 ; Col < mPGEquivalentNodeSize ; ++Col)
		{
			MySparseMatrixClass::SparseVector tmpSparseVector((aptrb[Col+1] - aptrb[Col]), asub, aval, &NonZerosIndexAdder);
			tmpSparseMatrix.PushBack(tmpSparseVector);
		}
		MyPGConductanceSubMatrix[each_conductance_submatrix.first] = tmpSparseMatrix;
	}
	// mPG_EigenConductanceSubMatrix.clear();	// clear origin Eigen matrix

	// Construct G0_INV * Q (Dense Matrix)
	double *cG0INV_Q = new double[mPGEquivalentNodeSize * mCurrentSourceSize]();
	for(unsigned int each_row = 0 ; each_row < mPGEquivalentNodeSize ; ++each_row)
	{
		for(unsigned int each_col = 0 ; each_col < mCurrentSourceSize ; ++each_col)
		{
			double value = mEachRowPGInverseMatrix[each_row][mCurrentSourceIdx[each_col]];
			if(value != 0)	// positively dense
			{
				cG0INV_Q[each_row + each_col*mPGEquivalentNodeSize] = value;
			}
		}
	}

	// Construct INVGT * M_EtoT, find the grid index to calculate the temperature, ak
	EigenDenseMatrix AkMatrix = mDenseInverseThermalMatrix * Eigen::MatrixXd(mElectricalGridtoThermalGridIncidenceMatrix);
	
	// load AkMatrix to sparse vector data structure
	unsigned int thermalgridsize = mXCut * mYCut * mZlayer;
	unsigned int ArraySize = thermalgridsize * mPGEquivalentNodeSize;
	double *cAkMatrix = new double[ArraySize]();
	Eigen::Map<Eigen::MatrixXd>( cAkMatrix, AkMatrix.rows(), AkMatrix.cols() ) = AkMatrix;
	std::vector< std::vector<double> > DenseAkMatrix(thermalgridsize);	// Q * Ak

	for(unsigned int i=0 ; i<thermalgridsize ; ++i )
		DenseAkMatrix[i].resize(mCurrentSourceSize);
	
	for(unsigned int row = 0 ; row < thermalgridsize ; ++row)
	{
		for(unsigned int col = 0 ; col < mPGEquivalentNodeSize ; ++col)
		{
			unsigned int index = row + col * thermalgridsize;
			if(cAkMatrix[index] != 0)
			{
				DenseAkMatrix[row][mCurrentSourceIdxMapping[col]] = cAkMatrix[index];
			}
		}
	}
	Dmatrixtime = (double)(clock()-start) / (CLOCKS_PER_SEC);
	std::cout << "Construct D matrix time = " << Dmatrixtime << std::endl;

	// start matlab engine

	Engine *ep;

	if (!(ep = engOpen(""))) {
		fprintf(stderr, "\nCan't start MATLAB engine\n");
		return ;
	}
	else {
		std::cout << "Started matlab engine !" << std::endl;
	}

	// engEvalString(ep, "LASTN = maxNumCompThreads(1);");

	// construct matlab QP A matrix
	//-----------------------------------------------------------------------------------------------------------
	// load mMatlabAmatrix & lower bound constraint to Amatrix  => 0 < Ui < Ig
	unsigned int Amatrixsize = 2 * mGlobalConstraintSize * mCurrentSourceSize;
	double *Amatrix = new double[Amatrixsize];
	Eigen::Map<Eigen::MatrixXd>( Amatrix, mMatlabAmatrix.rows(), mMatlabAmatrix.cols() ) = Eigen::MatrixXd(mMatlabAmatrix);
    
  	// load Amatrix to matAmatrix
	mxArray *matAmatrix = NULL;
	matAmatrix = mxCreateDoubleMatrix(2 * mGlobalConstraintSize, mCurrentSourceSize, mxREAL);
	memcpy((void *)mxGetPr(matAmatrix), (void *)Amatrix, sizeof(double) * Amatrixsize);
	engPutVariable(ep, "matAmatrix", matAmatrix);	// put variable to matlab workspace 
	//-----------------------------------------------------------------------------------------------------------
	
	// construct matlab QP b vector
	double *b = new double[2*mGlobalConstraintSize]();
	for(unsigned int i=0 ; i < mGlobalConstraintSize ; ++i)
	{
		b[i] = mGlobalCurrentConstraint[i];
	}
	mxArray *matb = NULL;
	matb = mxCreateDoubleMatrix(1, 2 * mGlobalConstraintSize, mxREAL);
	memcpy((void *)mxGetPr(matb), (void *)b, sizeof(double) * 2 * mGlobalConstraintSize);
	engPutVariable(ep, "matb", matb);

	// construct matlab QP local upper and lower variable constraints => 0 < i < Il
	//-----------------------------------------------------------------------------------------------------------
	double *clower = new double[mCurrentSourceSize]();
	double *cupper = new double[mCurrentSourceSize];
	for(unsigned int i=0 ; i < mCurrentSourceSize ; ++i)
	{
		cupper[i] = mLocalCurrentConstraint[i];
	}
	mxArray *matlower = NULL;
	mxArray *matupper = NULL;
	matlower = mxCreateDoubleMatrix(1, mCurrentSourceSize, mxREAL);
	matupper = mxCreateDoubleMatrix(1, mCurrentSourceSize, mxREAL);
	memcpy((void *)mxGetPr(matlower), (void *)clower, sizeof(double) * mCurrentSourceSize);
	memcpy((void *)mxGetPr(matupper), (void *)cupper, sizeof(double) * mCurrentSourceSize);
	engPutVariable(ep, "matlower", matlower);
	engPutVariable(ep, "matupper", matupper);
	//-----------------------------------------------------------------------------------------------------------
	
	//-----------------------------------------------------------------------------------------------------------
	// construct matlab linear terms 
	double *cf = new double[mCurrentSourceSize];
	mxArray *matf = NULL;
	matf = mxCreateDoubleMatrix(1, mCurrentSourceSize, mxREAL);

	// H matrix
	mxArray *matHmatrix = NULL;
	matHmatrix = mxCreateDoubleMatrix(mCurrentSourceSize, mCurrentSourceSize, mxREAL);
	double *cHbarmatrix = new double[mCurrentSourceSize*mCurrentSourceSize]();
	double *cHmatrix = new double[mCurrentSourceSize*mCurrentSourceSize]();
	

	// Factor
	double sFactor = -1 * mR_temp_coeff * mSupplyVoltage;	// -1 * alpha * Vdd 	

	// solution
	std::vector<double> VdropSolution;
	std::vector<double> CurrentErrorNorm;
	std::vector<double> CurrentErrorMaxElement;

	
	engEvalString(ep, "qptime = 0;");	// for calculating matlab quadprog time 
	clock_t matrixStart, simulationStart; 
	double QPtime(0), MatrixTime(0), SimulationTime(0);

	EigenDenseVector DkRow(mCurrentSourceSize);
	double *tmp_cG0INV_Q = new double[mPGEquivalentNodeSize * mCurrentSourceSize]();
    //int debug_idx(0);

    std::vector<unsigned int> veri_idx(2);

// PG1

    veri_idx[0] = 79;
    veri_idx[1] = 159;

// PG2
/*
    veri_idx[0] = 99;
    veri_idx[1] = 199;
*/
// PG3 (V2)
/*
veri_idx[0] = 139;
veri_idx[1] =144;
veri_idx[2] =149;
veri_idx[3] =289;
veri_idx[4] =294;
veri_idx[5] =439;
veri_idx[6] =444;
veri_idx[7] =589;
veri_idx[8] =594;
veri_idx[9] =744;
veri_idx[10] =894;
veri_idx[11] =1049;   
*/
// PG4
/*
veri_idx[0] = 139;
veri_idx[1] =144;
veri_idx[2] =149;
veri_idx[3] =289;
veri_idx[4] =294;
veri_idx[5] =439;
veri_idx[6] =444;
veri_idx[7] =589;
veri_idx[8] =594;
veri_idx[9] =739;
veri_idx[10] =744;
veri_idx[11] =889;
veri_idx[12] =894;
veri_idx[13] =1044;
veri_idx[14] =1049;
veri_idx[15] =1199;
veri_idx[16] =1349;
*/

//for(unsigned int each_node_idx = 0 ; each_node_idx < mPGEquivalentNodeSize ; ++each_node_idx)
//for(unsigned int each_node_idx = 5050 ; each_node_idx < 5050+1 ; ++each_node_idx)
for(auto &each_node_idx : veri_idx)
{
	
//if(each_node_idx == veri_idx[0] || each_node_idx == veri_idx[1])
//{

	// construct matlab f vector
	matrixStart = clock();

	for(unsigned int i=0 ; i<mCurrentSourceSize ; ++i)
	{
		cf[i] = -1 * mEachRowPGInverseMatrix[each_node_idx][mCurrentSourceIdx[i]];
	}
	memcpy((void *)mxGetPr(matf), (void *)cf, sizeof(double) * mCurrentSourceSize);
	engPutVariable(ep, "matf", matf);	
	//std::cout << "row time = " << std::setprecision(15) << (double)(clock()-start) / (CLOCKS_PER_SEC) << std::endl;
	
	// Construct cHbarmatrix
	for(unsigned int i = 0 ; i < mCurrentSourceSize * mCurrentSourceSize ; ++i)
	{
		cHbarmatrix[i] = 0;
	}

	for(auto& each_grid : MyPGConductanceSubMatrix)		
	{
		MySparseMatrixClass::SparseVector MyPreRow; 	// G0INV * G0,k

		for(unsigned int each_col = 0 ; each_col < mPGEquivalentNodeSize ; ++each_col)
		{
			double adjointValue(0);
			for(unsigned int i = 0 ; i < each_grid.second.getSparseVector(each_col).getVector().size() ; ++i)
			{
				adjointValue += each_grid.second.getSparseVector(each_col).getVector()[i] * mEachRowPGInverseMatrix[each_node_idx][each_grid.second.getSparseVector(each_col).getLocation()[i]];
			}
			if(adjointValue != 0)
			{
				MyPreRow.PushBackElement(each_col, adjointValue);
			}
		}

		for(unsigned int each_col = 0 ; each_col < mCurrentSourceSize ; ++each_col)	  //   Dk = (G0INV * G0,k) * (G0INV * Q)
		{
			double adjointValue(0);
			for(unsigned int i = 0 ; i < MyPreRow.getVector().size() ; ++i)
			{
				unsigned int index = MyPreRow.getLocation()[i] + each_col * mPGEquivalentNodeSize;
				if( cG0INV_Q[index] != 0)
				{
					adjointValue += MyPreRow.getVector()[i] * cG0INV_Q[index];
				}
			}
			DkRow[each_col] = adjointValue;
		}
		
		for(unsigned int row = 0 ; row < mCurrentSourceSize ; ++row)  // ak * Dk
		{
			for(unsigned int col = 0 ; col < mCurrentSourceSize ; ++col)
			{
				cHbarmatrix[row + col * mCurrentSourceSize] += DenseAkMatrix[each_grid.first][row] * DkRow[col];
			}
		}

	}
	
	// load sFactor * (cHbarmatrix + cHbarmatrix transpose) to cHmatrix
	for(unsigned int row = 0 ; row < mCurrentSourceSize ; ++row)
	{
		for(unsigned int col = 0 ; col < mCurrentSourceSize ; ++col)
		{
			unsigned int idx     =  row + col * mCurrentSourceSize;
			unsigned int tranidx =  col + row * mCurrentSourceSize;
			cHmatrix[idx] = sFactor * (cHbarmatrix[idx] + cHbarmatrix[tranidx]);
		}
	}
	memcpy((void *)mxGetPr(matHmatrix), (void *)cHmatrix, sizeof(double) * mCurrentSourceSize * mCurrentSourceSize);
	engPutVariable(ep, "matHmatrix", matHmatrix);	
	MatrixTime += (double)(clock() - matrixStart) / (CLOCKS_PER_SEC);

	// Run Matlab Quadratic Program
	// engEvalString(ep, "tStart = tic;");
	engEvalString(ep, "tStart = cputime;");
	engEvalString(ep, "[x,fval,exitflag] = quadprog(matHmatrix,matf,matAmatrix,matb,[],[],matlower,matupper);");
	engEvalString(ep, "qptime = cputime - tStart;");	
	// engEvalString(ep, "qptime = toc(tStart);");	

	std::cout <<  *(double *)mxGetPr(engGetVariable(ep, "exitflag")) << std::endl;
	
	QPtime += *(double *)mxGetPr(engGetVariable(ep, "qptime"));

	double *origin_current = (double *)mxGetPr(engGetVariable(ep, "x"));
	std::vector<double> origin_current_solution(mCurrentSourceSize);
	for(unsigned int i = 0 ; i < mCurrentSourceSize ; ++i)
	{
		origin_current_solution[i] = origin_current[i];
	}


	double sum(0);
	std::vector<double> errors;
	int MaxIteration = 15;
	int iteration = 0;

while( iteration < MaxIteration )
{
	sum = 0;
	errors.clear();

	// DC Simulation
	simulationStart = clock();

	double *xx = (double *)mxGetPr(engGetVariable(ep, "x"));		// current pattern solution of QP

	EigenDenseVector CurrentSourceVector(mPGEquivalentNodeSize), TemperatureVector(thermalgridsize);
	CurrentSourceVector.setZero();
	for(unsigned int i = 0 ; i < mCurrentSourceSize ; ++i)
	{
		CurrentSourceVector[mCurrentSourceIdx[i]] = xx[i];
	}
	TemperatureVector = AkMatrix * CurrentSourceVector;	// Obtain temperature vector

	// Stamping
	MATRIX tmpMatrix;
	tmpMatrix.resize(mPGEquivalentNodeSize);

	for (auto& stamp_node : mStampNodeList)
	{
		for (auto& nei_node : stamp_node.first->neighbor_NodeRes)    
		{
			if (nei_node.first->is_grounded == true)
			{
				for (auto& each_res : nei_node.second)
				{
					tmpMatrix[stamp_node.second][stamp_node.second] += 1 / (each_res->res_value * (1 + mR_temp_coeff * TemperatureVector[each_res->thermal_grid_index]));
				}
			}
			else if (nei_node.first->is_connect_V == true)
			{
				for (auto& each_res : nei_node.second)
				{				
					tmpMatrix[stamp_node.second][stamp_node.second] += 1 / (each_res->res_value * (1 + mR_temp_coeff * TemperatureVector[each_res->thermal_grid_index]));
				}
			}
			else if (nei_node.first->is_connect_V == false && nei_node.first->is_grounded == false)
			{
				for (auto& each_res : nei_node.second)
				{
					tmpMatrix[stamp_node.second][stamp_node.second] += 1 / (each_res->res_value * (1 + mR_temp_coeff * TemperatureVector[each_res->thermal_grid_index]));
					tmpMatrix[stamp_node.second][mStampNodeList[nei_node.first]] += -1 / (each_res->res_value * (1 + mR_temp_coeff * TemperatureVector[each_res->thermal_grid_index]));
				}
			}
		}
	}
	std::vector<triplet> tripletList;
	for(unsigned int row_idx = 0; row_idx < mPGEquivalentNodeSize; ++row_idx)
    {
        for(auto& row_data : tmpMatrix[row_idx])
        {           
            tripletList.push_back(triplet(row_idx, row_data.first, row_data.second));   
        }
    }
    Eigen::SparseMatrix<double> sys_mat(mPGEquivalentNodeSize, mPGEquivalentNodeSize);
	sys_mat.setFromTriplets(tripletList.begin(), tripletList.end());
	tripletList.clear();

	// Solve inverse matrix
	Eigen::VectorXd pattern(mPGEquivalentNodeSize);
	Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
	solver.analyzePattern(sys_mat);   
	solver.factorize(sys_mat);
	pattern.setZero();
	std::vector<EigenDenseVector> tmp_EachRowPGInverseMatrix;
	for(unsigned int i = 0 ; i < mPGEquivalentNodeSize ; ++i)
	{
		pattern[i] = 1;
		tmp_EachRowPGInverseMatrix.emplace_back(solver.solve(pattern));
		pattern[i] = 0;
	}
	// SubConductance Matrix
	for(unsigned int y = 0; y < mYCut ; ++y)
	{
		for(unsigned int x = 0; x < mXCut ; ++x)
		{
			MATRIX subConductanceMatrix(mPGEquivalentNodeSize);
			// Stamp each subconductance matrix
			for (auto& stamp_node : mStampNodeList)
			{
				for (auto& nei_node : stamp_node.first->neighbor_NodeRes)    
				{
					if (nei_node.first->is_grounded == true)
					{
						for (auto& each_res : nei_node.second)
						{
							if(each_res->is_in_submatrix)
							{
								subConductanceMatrix[stamp_node.second][stamp_node.second] += 1 / (each_res->res_value * ( 1 + mR_temp_coeff * TemperatureVector[each_res->thermal_grid_index])) ;
							}
						}
					}
					else if (nei_node.first->is_connect_V == true)
					{
						for (auto& each_res : nei_node.second)
						{				
							if(each_res->is_in_submatrix)
							{
								subConductanceMatrix[stamp_node.second][stamp_node.second] += 1 / (each_res->res_value * ( 1 + mR_temp_coeff * TemperatureVector[each_res->thermal_grid_index])) ;
							}
						}
					}
					else if (nei_node.first->is_connect_V == false && nei_node.first->is_grounded == false)
					{
						for (auto& each_res : nei_node.second)
						{
							if(each_res->is_in_submatrix)
							{
								subConductanceMatrix[stamp_node.second][stamp_node.second] += 1 / (each_res->res_value * ( 1 + mR_temp_coeff * TemperatureVector[each_res->thermal_grid_index])) ;
								subConductanceMatrix[stamp_node.second][mStampNodeList[nei_node.first]] += -1 / (each_res->res_value * ( 1 + mR_temp_coeff * TemperatureVector[each_res->thermal_grid_index])) ;
							}
						}
					}
				}
			}
			// mPG_ConductanceSubMatrix 有在PowerGrid.cpp裡面clear過
			mPG_ConductanceSubMatrix[mTotalThermalGrids[mZlayer-1][y][x].mGridIndex] = subConductanceMatrix; // store
		}
	}
	std::unordered_map<unsigned int, EigenSparseMatrix> tmp_PG_EigenConductanceSubMatrix;
	for(auto& each_matrix : mPG_ConductanceSubMatrix)
	{
		EigenSparseMatrix Mat_PGSubMatrix(mPGEquivalentNodeSize, mPGEquivalentNodeSize);	
		std::vector<triplet> triplet_PGSubMatrix;

		for(unsigned int row_idx = 0; row_idx < each_matrix.second.size(); ++row_idx)
    	{
        	for(auto& row_data : each_matrix.second[row_idx])
        	{           
            	triplet_PGSubMatrix.push_back(triplet(row_idx, row_data.first, row_data.second));   
        	}
    	}
    	Mat_PGSubMatrix.setFromTriplets(triplet_PGSubMatrix.begin(), triplet_PGSubMatrix.end());
    	tmp_PG_EigenConductanceSubMatrix[each_matrix.first] = Mat_PGSubMatrix;
	}
	mPG_ConductanceSubMatrix.clear();

	std::map<unsigned int, MySparseMatrixClass::SparseMatrix> tmp_MyPGConductanceSubMatrix;	 
	for(auto& each_conductance_submatrix : tmp_PG_EigenConductanceSubMatrix)
	{
		each_conductance_submatrix.second.makeCompressed();	
		int    *aptrb    = each_conductance_submatrix.second.outerIndexPtr();	  
		int    *asub     = each_conductance_submatrix.second.innerIndexPtr();	 
		double *aval     = each_conductance_submatrix.second.valuePtr();		
		unsigned int NonZerosIndexAdder(0);

		MySparseMatrixClass::SparseMatrix tmpSparseMatrix;
		for(unsigned int Col = 0 ; Col < mPGEquivalentNodeSize ; ++Col)
		{
			MySparseMatrixClass::SparseVector tmpSparseVector((aptrb[Col+1] - aptrb[Col]), asub, aval, &NonZerosIndexAdder);
			tmpSparseMatrix.PushBack(tmpSparseVector);
		}
		tmp_MyPGConductanceSubMatrix[each_conductance_submatrix.first] = tmpSparseMatrix;
	}
	tmp_PG_EigenConductanceSubMatrix.clear();
	// G0INV_Q
	for(unsigned int each_row = 0 ; each_row < mPGEquivalentNodeSize ; ++each_row)
	{
		for(unsigned int each_col = 0 ; each_col < mCurrentSourceSize ; ++each_col)
		{
			double value = tmp_EachRowPGInverseMatrix[each_row][mCurrentSourceIdx[each_col]];
			if(value != 0)	// positively dense
			{
				tmp_cG0INV_Q[each_row + each_col*mPGEquivalentNodeSize] = value;
			}
		}
	}
	// Solve QP
	for(unsigned int i=0 ; i<mCurrentSourceSize ; ++i)
	{
		cf[i] = -1 * tmp_EachRowPGInverseMatrix[each_node_idx][mCurrentSourceIdx[i]];
	}
	memcpy((void *)mxGetPr(matf), (void *)cf, sizeof(double) * mCurrentSourceSize);
	engPutVariable(ep, "matf", matf);	

	for(unsigned int i = 0 ; i < mCurrentSourceSize * mCurrentSourceSize ; ++i)
	{
		cHbarmatrix[i] = 0;
	}

	for(auto& each_grid : tmp_MyPGConductanceSubMatrix)		
	{
		MySparseMatrixClass::SparseVector MyPreRow; 	// G0INV * G0,k

		for(unsigned int each_col = 0 ; each_col < mPGEquivalentNodeSize ; ++each_col)
		{
			double adjointValue(0);
			for(unsigned int i = 0 ; i < each_grid.second.getSparseVector(each_col).getVector().size() ; ++i)
			{
				adjointValue += each_grid.second.getSparseVector(each_col).getVector()[i] * tmp_EachRowPGInverseMatrix[each_node_idx][each_grid.second.getSparseVector(each_col).getLocation()[i]];
			}
			if(adjointValue != 0)
			{
				MyPreRow.PushBackElement(each_col, adjointValue);
			}
		}

		for(unsigned int each_col = 0 ; each_col < mCurrentSourceSize ; ++each_col)	  //   Dk = (G0INV * G0,k) * (G0INV * Q)
		{
			double adjointValue(0);
			for(unsigned int i = 0 ; i < MyPreRow.getVector().size() ; ++i)
			{
				unsigned int index = MyPreRow.getLocation()[i] + each_col * mPGEquivalentNodeSize;
				if( tmp_cG0INV_Q[index] != 0)
				{
					adjointValue += MyPreRow.getVector()[i] * tmp_cG0INV_Q[index];
				}
			}
			DkRow[each_col] = adjointValue;
		}
		
		for(unsigned int row = 0 ; row < mCurrentSourceSize ; ++row)  // ak * Dk
		{
			for(unsigned int col = 0 ; col < mCurrentSourceSize ; ++col)
			{
				cHbarmatrix[row + col * mCurrentSourceSize] += DenseAkMatrix[each_grid.first][row] * DkRow[col];
			}
		}

	}
	
	// load sFactor * (cHbarmatrix + cHbarmatrix transpose) to cHmatrix
	for(unsigned int row = 0 ; row < mCurrentSourceSize ; ++row)
	{
		for(unsigned int col = 0 ; col < mCurrentSourceSize ; ++col)
		{
			unsigned int idx     =  row + col * mCurrentSourceSize;
			unsigned int tranidx =  col + row * mCurrentSourceSize;
			cHmatrix[idx] = sFactor * (cHbarmatrix[idx] + cHbarmatrix[tranidx]);
		}
	}
	memcpy((void *)mxGetPr(matHmatrix), (void *)cHmatrix, sizeof(double) * mCurrentSourceSize * mCurrentSourceSize);
	engPutVariable(ep, "matHmatrix", matHmatrix);	
	engEvalString(ep, "[x,fval,exitflag] = quadprog(matHmatrix,matf,matAmatrix,matb,[],[],matlower,matupper);");
	

	std::cout <<  *(double *)mxGetPr(engGetVariable(ep, "exitflag")) << std::endl;
	SimulationTime += (double)(clock() - simulationStart) / (CLOCKS_PER_SEC);

	// calculate current error
	double *solution_current = (double *)mxGetPr(engGetVariable(ep, "x"));
	for(unsigned int i = 0 ; i < mCurrentSourceSize ; ++i)
	{
		sum += pow(fabs(origin_current_solution[i] - solution_current[i]), 2);
		errors.emplace_back(fabs(origin_current_solution[i] - solution_current[i]));
	}

	if(pow(sum, 0.5) < 1e-3)
	{
		std::cout << "Converge at iteration " << iteration << " ! " << std::endl;
		break;
	}
	++iteration;

}

	VdropSolution.emplace_back(-1 * *(double *)mxGetPr(engGetVariable(ep, "fval")));
	CurrentErrorNorm.emplace_back(pow(sum, 0.5));
	CurrentErrorMaxElement.emplace_back(*std::max_element(errors.begin(), errors.end()));
//}

}

	std::cout << "Complete Efficient QP Verification !" << std::endl;

	// Output solution 
	std::ofstream fout("../LPresult/" + mCktName + "_EfficientQPVdropMap_Simulation_Verify_Iterative");

	int sol_idx(0);
	for(auto& idx : veri_idx)
	{
		fout << mPGEquivalentNodes[idx]->x_loc_idx << " " << mPGEquivalentNodes[idx]->y_loc_idx 
			<< " " << std::setprecision(15) << VdropSolution[sol_idx] << std::endl;
		++sol_idx;
	}
	
	/*
	
	for(unsigned int i = 0 ; i < mPGEquivalentNodeSize ; ++i)
	{
		if(mPGEquivalentNodes[i]->is_connect_I)
		{
			fout << mPGEquivalentNodes[i]->x_loc_idx << " " << mPGEquivalentNodes[i]->y_loc_idx 
				<< " " << std::setprecision(15) << VdropSolution[sol_idx] << std::endl;
	
			if(VdropSolution[sol_idx] > mSupplyVoltage * 0.1)
			{
				++violations;
			}
			++sol_idx;
		}
	}
	*/

	fout.close();
	
	fout.open("../LPresult/" + mCktName + "_EfficientQP_CurrentError_Verify_Iterative");
	for(unsigned int i = 0 ; i < CurrentErrorNorm.size() ; ++i)
	{
		fout << CurrentErrorNorm[i] << " " << CurrentErrorMaxElement[i] << std::endl;
	}
	fout.close();

	/*
	fout.open("../LPresult/" + mCktName + "_EfficientQP_TimeInfo_Simulation_Verify");
	fout << "Matrix Time = " << std::setprecision(15) << MatrixTime << " (S) " << std::endl;
	fout << "Matlab QP Time = " << std::setprecision(15) << QPtime << " (S) " << std::endl;
	fout << "Total time (Matrix + Matlab) = " << std::setprecision(15) << MatrixTime + QPtime + SimulationTime << " (S) " << std::endl;
	fout << "Total verification time = " << mConstructDataTime + Dmatrixtime + MatrixTime + QPtime << " (S) " << std::endl;
	fout << "Average voltage drop = " << std::setprecision(15) << sum / VdropSolution.size() << " (V) " << std::endl;
	fout << "Maximum voltage drop = " << std::setprecision(15) << *std::max_element(VdropSolution.begin(), VdropSolution.end()) << " (V) " << std::endl;
	fout << "Total violations = " << violations << std::endl;
	fout << "Simulation Time = " << SimulationTime << " (s) " << std::endl;
	fout.close();
	*/
	// Delete memory
	mxDestroyArray(matHmatrix);
	mxDestroyArray(matf);
	mxDestroyArray(matAmatrix);
	mxDestroyArray(matb);
	mxDestroyArray(matlower);
	mxDestroyArray(matupper);
	engClose(ep);

	delete [] cG0INV_Q;
	delete [] Amatrix;
	delete [] b;
	delete [] clower;
	delete [] cupper;
	delete [] cf;
	delete [] cHbarmatrix;
	delete [] cHmatrix;
	delete [] tmp_cG0INV_Q;

}


// Expansion Point Method without IR drop recalculation
// using Mosek for LP
void PGSolver::PGVerificationExpandPoint_Efficient()
{
	double TotalTime(0);

	clock_t start;
	
	start = clock();

	std::map<unsigned int, MySparseMatrixClass::SparseMatrix> MyPGConductanceSubMatrix;	 

	for(auto& each_conductance_submatrix : mPG_EigenConductanceSubMatrix)
	{
		(each_conductance_submatrix.second).makeCompressed();	

		int    *aptrb    = (each_conductance_submatrix.second).outerIndexPtr();	  
		int    *asub     = (each_conductance_submatrix.second).innerIndexPtr();	 
		double *aval     = (each_conductance_submatrix.second).valuePtr();		
		unsigned int NonZerosIndexAdder(0);

		MySparseMatrixClass::SparseMatrix tmpSparseMatrix;
		for(unsigned int Col = 0 ; Col < mPGEquivalentNodeSize ; ++Col)
		{
			MySparseMatrixClass::SparseVector tmpSparseVector((aptrb[Col+1] - aptrb[Col]), asub, aval, &NonZerosIndexAdder);
			tmpSparseMatrix.PushBack(tmpSparseVector);
		}
		MyPGConductanceSubMatrix[each_conductance_submatrix.first] = tmpSparseMatrix;
	}
	mPG_EigenConductanceSubMatrix.clear();


	double *cG0INV_Q = new double[mPGEquivalentNodeSize * mCurrentSourceSize]();   // G0_inv * Q

	for(unsigned int each_row = 0 ; each_row < mPGEquivalentNodeSize ; ++each_row)
	{
		for(unsigned int each_col = 0 ; each_col < mCurrentSourceSize ; ++each_col)
		{
			double value = mEachRowPGInverseMatrix[each_row][mCurrentSourceIdx[each_col]];
			if(value != 0)
			{
				cG0INV_Q[each_row + each_col*mPGEquivalentNodeSize] = value;
			}
		}
	}

	// Construct INVGT * M_EtoT, find the grid index to calculate the temperature, ak
	EigenDenseMatrix AkMatrix = mDenseInverseThermalMatrix * Eigen::MatrixXd(mElectricalGridtoThermalGridIncidenceMatrix);
	
	// load AkMatrix to sparse vector data structure
	unsigned int thermalgridsize = mXCut * mYCut * mZlayer;
	
	unsigned int ArraySize = thermalgridsize * mPGEquivalentNodeSize;
	
	double *cAkMatrix = new double[ArraySize]();
	
	Eigen::Map<Eigen::MatrixXd>( cAkMatrix, AkMatrix.rows(), AkMatrix.cols() ) = AkMatrix;
	
	double thermalResistanceFactor = mSupplyVoltage * mR_temp_coeff; 	

	std::vector< std::vector<double> > DenseAkMatrix(thermalgridsize);  // alpha * Vdd * Gt_inv  * M_EtoT * Q
	
	for(unsigned int i = 0 ; i < thermalgridsize ; ++i )
	{	
		DenseAkMatrix[i].resize(mCurrentSourceSize);
	}
	
	for(unsigned int row = 0 ; row < thermalgridsize ; ++row)
	{
		for(unsigned int col = 0 ; col < mPGEquivalentNodeSize ; ++col)
		{
			unsigned int index = row + col * thermalgridsize;
			if(cAkMatrix[index] != 0)
			{
				DenseAkMatrix[row][mCurrentSourceIdxMapping[col]] = thermalResistanceFactor * cAkMatrix[index];
			}
		}
	}
	

	int NUMCON = mGlobalConstraintSize;  	
	int NUMVAR = mCurrentSourceSize;	   
	
	MSKboundkeye *bkc = new MSKboundkeye[NUMCON];
	for(int i=0 ; i<NUMCON ; ++i) { bkc[i] = MSK_BK_RA; }
	double *blc = new double[NUMCON]();  
	double *buc = new double[NUMCON]; 
	for(int i=0 ; i<NUMCON ; ++i) { buc[i] = mGlobalCurrentConstraint[i]; }

	mGlobalCurrentConstraintIncidenceMatrix.makeCompressed();
	
	int    *aptrb    = mGlobalCurrentConstraintIncidenceMatrix.outerIndexPtr();
	int    *asub     = mGlobalCurrentConstraintIncidenceMatrix.innerIndexPtr();
	double *aval     = mGlobalCurrentConstraintIncidenceMatrix.valuePtr();
	
	
	MSKboundkeye *bkx = new MSKboundkeye[NUMVAR];
	for(int i=0 ; i<NUMVAR ; ++i) { bkx[i] = MSK_BK_RA; }
	double *blx = new double[NUMVAR](); 
	double *bux = new double[NUMVAR];
	for(int i=0 ; i<NUMVAR ; ++i) { bux[i] = mLocalCurrentConstraint[i]; }
	
	
	Eigen::VectorXd LPC(NUMVAR); 	
	//double *cHjmatrix = new double[mCurrentSourceSize * mCurrentSourceSize]();
	EigenDenseVector DkRow(mCurrentSourceSize);
	

	double *xx = new double[NUMVAR];  
	std::vector<double> VdropSolution;


	MSKenv_t env = NULL;
	MSKtask_t task = NULL;
	MSKrescodee r;
	MSKrescodee trmcode;

	r = MSK_makeenv(&env, NULL);
	r = MSK_maketask(env, NUMCON, NUMVAR, &task);
	r = MSK_appendcons(task, NUMCON); 
	r = MSK_appendvars(task, NUMVAR);

	for(int j=0; j<NUMVAR && r == MSK_RES_OK; ++j)
	{	
		r = MSK_putvarbound(task, j, bkx[j], blx[j], bux[j]); 
		r = MSK_putacol(task, j, aptrb[j+1]-aptrb[j], asub+aptrb[j], aval+aptrb[j]);
	}
	
	for(int i=0; i<NUMCON && r==MSK_RES_OK; ++i)
		r = MSK_putconbound(task, i, bkc[i], blc[i], buc[i]); 

	
	r = MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MAXIMIZE);  

	std::cout << "Optimizing ..." << std::endl; 
	

	// unsigned int specificIdx(5050);
	std::unordered_map<unsigned int, std::vector<double> > ExpandPoint;
	std::vector<unsigned int> ThermalGridIndexContainer;
	std::unordered_map<unsigned int, std::vector<unsigned int> > NodeContainer;
	std::vector<unsigned int> VdropIndex;

for(auto & submatrix : MyPGConductanceSubMatrix)
{
	
	for(int i=0 ; i<NUMVAR ; ++i)
	{
		r = MSK_putcj(task, i, DenseAkMatrix[submatrix.first][i]);	
	}
	
	r = MSK_optimizetrm(task, &trmcode);
	MSK_getxx(task, MSK_SOL_ITR, xx);

	// store expand point
	std::vector<double> exp(NUMVAR);
	for(int i=0 ; i<NUMVAR ; ++i)
	{
		exp[i] = xx[i];
	}
	ExpandPoint[submatrix.first] = exp;


	// Node Partitioning, stored in NodeContainer
	for(auto& node : mPGEquivalentNodes)
	{
		if(node->is_connect_I && node->thermal_grid_index == submatrix.first)
		{
			NodeContainer[submatrix.first].emplace_back((node->x_loc_idx - 1) * mY_loc_max + (node->y_loc_idx - 1));
		}
	}

	ThermalGridIndexContainer.emplace_back(submatrix.first);

	// Temperature map
/*
	if(mPGEquivalentNodes[specificIdx]->thermal_grid_index == submatrix.first)
		OnePointThermalAnalysis(specificIdx, xx, AkMatrix, 1);
*/
}

	TotalTime += (double)(clock()-start) / (CLOCKS_PER_SEC);

	std::cout << "Storing expand point complete !" << std::endl;

	double MatrixTime(0);
	double MosekTime(0);
	clock_t matrixstart;

	EigenDenseVector Hj_trans_row(mCurrentSourceSize);
	

// Run each node for LP
for(auto& each_container : NodeContainer)
{

	std::vector<double> exp = ExpandPoint[each_container.first];
	std::vector<double> TemperatureAdjointValue(thermalgridsize);

	for(auto& thermalGridIndex : ThermalGridIndexContainer)
	{
		double adjointValue(0);
		for(unsigned int row = 0 ; row < mCurrentSourceSize ; ++row) 
		{
			adjointValue += DenseAkMatrix[thermalGridIndex][row] * exp[row];
		}
		TemperatureAdjointValue[thermalGridIndex] = adjointValue;
	}

for(auto& each_node_idx : each_container.second)
{

	start = clock();

	matrixstart = clock();

	for(int i = 0 ; i < NUMVAR ; ++i)
	{
		unsigned int index = mPGEquivalentNodeSize * i + each_node_idx;	// rowsize * column index + row index
		LPC[i] = cG0INV_Q[index];
	}
	/*
	for(unsigned int i = 0 ; i < mCurrentSourceSize * mCurrentSourceSize ; ++i)  
	{
		cHjmatrix[i] = 0;
	}
	*/
	Hj_trans_row.setZero();
	
	for(auto& each_grid : MyPGConductanceSubMatrix)		
	{
		MySparseMatrixClass::SparseVector MyPreRow; 	// row each_node_idx of 0INV * G0,k

		for(unsigned int each_col = 0 ; each_col < mPGEquivalentNodeSize ; ++each_col)
		{
			double adjointValue(0);
			for(unsigned int i = 0 ; i < each_grid.second.getSparseVector(each_col).getVector().size() ; ++i)
			{
				adjointValue += each_grid.second.getSparseVector(each_col).getVector()[i] * mEachRowPGInverseMatrix[each_node_idx][each_grid.second.getSparseVector(each_col).getLocation()[i]];
			}
			if(adjointValue != 0)
			{
				MyPreRow.PushBackElement(each_col, adjointValue);
			}
		}

		for(unsigned int each_col = 0 ; each_col < mCurrentSourceSize ; ++each_col)	  // row each_node_idx of Dk = (G0INV * G0,k) * (G0INV * Q)
		{
			double adjointValue(0);
			for(unsigned int i = 0 ; i < MyPreRow.getVector().size() ; ++i)
			{
				unsigned int index = MyPreRow.getLocation()[i] + each_col * mPGEquivalentNodeSize;
				if( cG0INV_Q[index] != 0)
				{
					adjointValue += MyPreRow.getVector()[i] * cG0INV_Q[index];
				}
			}
			DkRow[each_col] = adjointValue;
		}
		/*
		for(unsigned int row = 0 ; row < mCurrentSourceSize ; ++row)  // ak * Dk
		{
			for(unsigned int col = 0 ; col < mCurrentSourceSize ; ++col)
			{
				cHjmatrix[row + col * mCurrentSourceSize] += DenseAkMatrix[each_grid.first][row] * DkRow[col];
			}
		}
		*/
		double TransAdjointValue(0);
		for(unsigned int row = 0 ; row < mCurrentSourceSize ; ++row) 
		{
			TransAdjointValue += DkRow[row] * exp[row];
		}

		for(unsigned int row = 0 ; row < mCurrentSourceSize ; ++row) 
		{
			double tmp = TransAdjointValue * DenseAkMatrix[each_grid.first][row];
			Hj_trans_row[row] += tmp;
			LPC[row] += TemperatureAdjointValue[each_grid.first] * DkRow[row] + tmp;

		}
	}

	double VoltageShift(0);
	for(unsigned int i = 0 ; i < mCurrentSourceSize ; ++i)
	{
		VoltageShift += exp[i] * Hj_trans_row[i];
	}

	/*
	std::vector<double> exp = ExpandPoint[mPGEquivalentNodes[each_node_idx]->thermal_grid_index];
	for(unsigned int row = 0 ; row < mCurrentSourceSize ; ++row)
	{
		double adjointValue(0);
		double adjointValueShift(0);
		for(unsigned int col = 0 ; col < mCurrentSourceSize ; ++col)
		{
			unsigned int idx     =  col + row * mCurrentSourceSize;
			unsigned int tranidx =  row + col * mCurrentSourceSize; 
			double SumHjtran = exp[col] * cHjmatrix[tranidx];
			double SumHj_Hjtran = exp[col] * cHjmatrix[idx] + SumHjtran;
			adjointValue += SumHj_Hjtran;
			adjointValueShift += SumHjtran;
		}
		LPC[row] += adjointValue;
		VoltageShift += adjointValueShift * exp[row];
	}
	*/
	MatrixTime += (double)(clock()-matrixstart) / (CLOCKS_PER_SEC);
	//-----------------------------------------------------------------------------------------------
	// Mosek
	
	clock_t mosek_start;
	mosek_start = clock();

	for(unsigned int row = 0 ; row < mCurrentSourceSize ; ++row)
	{
		r = MSK_putcj(task, row, LPC[row]);
	}

	double IRdrop;
	r = MSK_optimizetrm(task, &trmcode);
	MSK_getprimalobj(task, MSK_SOL_ITR, &IRdrop);

	MosekTime += (double)(clock()-mosek_start) / (CLOCKS_PER_SEC);
	TotalTime += (double)(clock()-start) / (CLOCKS_PER_SEC);

	double Vdrop(0);
	Vdrop = IRdrop - VoltageShift;
	VdropSolution.emplace_back(Vdrop);
	VdropIndex.emplace_back(each_node_idx);

	// OnePointThermalAnalysis(specificIdx, xx, AkMatrix, 2);
}

}

	// std::cout << "Matrix Time = " << MatrixTime << std::endl;

	std::cout << "Complete Efficient Expand Point Verification !" << std::endl;

// Output solution 
	
	unsigned int violations(0);
	
	double sum(0);
	for(auto& vdrop : VdropSolution) { sum +=vdrop; }

	for(unsigned int i = 0 ; i < VdropIndex.size() ; ++i)
	{
		for(unsigned int j = i + 1 ; j < VdropIndex.size() ; ++j)
		{
			if(VdropIndex[i] > VdropIndex[j])
			{
				unsigned int idx_tmp = VdropIndex[i];
				VdropIndex[i] = VdropIndex[j];
				VdropIndex[j] = idx_tmp;

				double sol_tmp = VdropSolution[i];
				VdropSolution[i] = VdropSolution[j];
				VdropSolution[j] = sol_tmp;
			}
		}
	}

	std::ofstream fout("../LPresult/" + mCktName + "_EfficientExpandLPVdropMap");
	int sol_idx(0);
	for(unsigned int i = 0 ; i < mPGEquivalentNodeSize ; ++i)
	{
		if(mPGEquivalentNodes[i]->is_connect_I)
		{
			fout << mPGEquivalentNodes[i]->x_loc_idx << " " << mPGEquivalentNodes[i]->y_loc_idx 
				 << " " << std::setprecision(15) << VdropSolution[sol_idx] << std::endl;
		
			if(VdropSolution[sol_idx] > mSupplyVoltage * 0.1)
			{
				++violations;
			}
			++sol_idx;
		}
	}
	fout.close();

	fout.open("../LPresult/" + mCktName + "_EfficientExpandLP_TimeInfo");
	fout << "Total verification time = " << mConstructDataTime + TotalTime << " (s) " << std::endl;
	fout << "Average voltage drop = " << std::setprecision(15) << sum / VdropSolution.size() << " (V) " << std::endl;
	fout << "Maximum voltage drop = " << std::setprecision(15) << *std::max_element(VdropSolution.begin(), VdropSolution.end()) << " (V) " << std::endl;
	fout << "Total violations = " << violations << std::endl;
	fout << "Matrix Time = " << MatrixTime << " (s) " << std::endl;
	fout << "Mosek Time  = " << MosekTime << " (s) " << std::endl;
	fout.close();
	
// delete memory

	MSK_deletetask(&task);
	MSK_deleteenv(&env);

	delete [] cAkMatrix;
	delete [] cG0INV_Q;
	//delete [] cHjmatrix;
	delete [] bkc;
	delete [] blc;	
	delete [] buc;
	delete [] bkx;
	delete [] blx;	
	delete [] bux;
	delete [] xx;


}


// Expansion Point Method with IR drop recalculation
void PGSolver::PGVerificationExpandPoint_Efficient_Simulation()
{
	double TotalTime(0);

	clock_t start;
	
	start = clock();

	std::map<unsigned int, MySparseMatrixClass::SparseMatrix> MyPGConductanceSubMatrix;	 

	for(auto& each_conductance_submatrix : mPG_EigenConductanceSubMatrix)
	{
		(each_conductance_submatrix.second).makeCompressed();	

		int    *aptrb    = (each_conductance_submatrix.second).outerIndexPtr();	  
		int    *asub     = (each_conductance_submatrix.second).innerIndexPtr();	 
		double *aval     = (each_conductance_submatrix.second).valuePtr();		
		unsigned int NonZerosIndexAdder(0);

		MySparseMatrixClass::SparseMatrix tmpSparseMatrix;
		for(unsigned int Col = 0 ; Col < mPGEquivalentNodeSize ; ++Col)
		{
			MySparseMatrixClass::SparseVector tmpSparseVector((aptrb[Col+1] - aptrb[Col]), asub, aval, &NonZerosIndexAdder);
			tmpSparseMatrix.PushBack(tmpSparseVector);
		}
		MyPGConductanceSubMatrix[each_conductance_submatrix.first] = tmpSparseMatrix;
	}
	mPG_EigenConductanceSubMatrix.clear();


	double *cG0INV_Q = new double[mPGEquivalentNodeSize * mCurrentSourceSize]();   // G0_inv * Q

	for(unsigned int each_row = 0 ; each_row < mPGEquivalentNodeSize ; ++each_row)
	{
		for(unsigned int each_col = 0 ; each_col < mCurrentSourceSize ; ++each_col)
		{
			double value = mEachRowPGInverseMatrix[each_row][mCurrentSourceIdx[each_col]];
			if(value != 0)
			{
				cG0INV_Q[each_row + each_col*mPGEquivalentNodeSize] = value;
			}
		}
	}

	// Construct Gt_inv * M_EtoT, find the grid index to calculate the temperature, ak
	EigenDenseMatrix AkMatrix = mDenseInverseThermalMatrix * Eigen::MatrixXd(mElectricalGridtoThermalGridIncidenceMatrix);
	
	// load AkMatrix to sparse vector data structure
	unsigned int thermalgridsize = mXCut * mYCut * mZlayer;
	
	unsigned int ArraySize = thermalgridsize * mPGEquivalentNodeSize;
	
	double *cAkMatrix = new double[ArraySize]();
	
	Eigen::Map<Eigen::MatrixXd>( cAkMatrix, AkMatrix.rows(), AkMatrix.cols() ) = AkMatrix;
	
	double thermalResistanceFactor = mSupplyVoltage * mR_temp_coeff; 	

	std::vector< std::vector<double> > DenseAkMatrix(thermalgridsize);		// alpha * Vdd * Gt_inv  * M_EtoT * Q
	
	for(unsigned int i = 0 ; i < thermalgridsize ; ++i )
	{	
		DenseAkMatrix[i].resize(mCurrentSourceSize);
	}
	
	for(unsigned int row = 0 ; row < thermalgridsize ; ++row)
	{
		for(unsigned int col = 0 ; col < mPGEquivalentNodeSize ; ++col)
		{
			unsigned int index = row + col * thermalgridsize;
			if(cAkMatrix[index] != 0)
			{
				DenseAkMatrix[row][mCurrentSourceIdxMapping[col]] = thermalResistanceFactor * cAkMatrix[index];
			}
		}
	}
	

	int NUMCON = mGlobalConstraintSize;  	
	int NUMVAR = mCurrentSourceSize;	   
	
	MSKboundkeye *bkc = new MSKboundkeye[NUMCON];
	for(int i=0 ; i<NUMCON ; ++i) { bkc[i] = MSK_BK_RA; }
	double *blc = new double[NUMCON]();  
	double *buc = new double[NUMCON]; 
	for(int i=0 ; i<NUMCON ; ++i) { buc[i] = mGlobalCurrentConstraint[i]; }

	mGlobalCurrentConstraintIncidenceMatrix.makeCompressed();
	
	int    *aptrb    = mGlobalCurrentConstraintIncidenceMatrix.outerIndexPtr();
	int    *asub     = mGlobalCurrentConstraintIncidenceMatrix.innerIndexPtr();
	double *aval     = mGlobalCurrentConstraintIncidenceMatrix.valuePtr();
	
	
	MSKboundkeye *bkx = new MSKboundkeye[NUMVAR];
	for(int i=0 ; i<NUMVAR ; ++i) { bkx[i] = MSK_BK_RA; }
	double *blx = new double[NUMVAR](); 
	double *bux = new double[NUMVAR];
	for(int i=0 ; i<NUMVAR ; ++i) { bux[i] = mLocalCurrentConstraint[i]; }
	
	
	Eigen::VectorXd LPC(NUMVAR); 	
	// double *cHjmatrix = new double[mCurrentSourceSize * mCurrentSourceSize]();
	EigenDenseVector DkRow(mCurrentSourceSize);
	

	double *xx = new double[NUMVAR];  
	std::vector<double> VdropSolution, VdropSolution_Origin;


	MSKenv_t env = NULL;
	MSKtask_t task = NULL;
	MSKrescodee r;
	MSKrescodee trmcode;

	r = MSK_makeenv(&env, NULL);
	r = MSK_maketask(env, NUMCON, NUMVAR, &task);
	r = MSK_appendcons(task, NUMCON); 
	r = MSK_appendvars(task, NUMVAR);

	for(int j=0; j<NUMVAR && r == MSK_RES_OK; ++j)
	{	
		r = MSK_putvarbound(task, j, bkx[j], blx[j], bux[j]); 
		r = MSK_putacol(task, j, aptrb[j+1]-aptrb[j], asub+aptrb[j], aval+aptrb[j]);
	}
	
	for(int i=0; i<NUMCON && r==MSK_RES_OK; ++i)
		r = MSK_putconbound(task, i, bkc[i], blc[i], buc[i]); 

	
	r = MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MAXIMIZE);  

	std::cout << "Optimizing ..." << std::endl; 
	

	// unsigned int specificIdx(5050);

	std::unordered_map<unsigned int, std::vector<double> > ExpandPoint;		// Container for expansion point current patterns
	std::vector<unsigned int> ThermalGridIndexContainer;  // recording the thermal grid index in interconnect layer
	std::unordered_map<unsigned int, std::vector<unsigned int> > NodeContainer;	// Storing the nodes in each thermal grid, key=thermal grid index, value=grouop of power grid nodes
	std::vector<unsigned int> VdropIndex;	// recording the node index

for(auto & submatrix : MyPGConductanceSubMatrix)		
{
	
	for(int i=0 ; i<NUMVAR ; ++i)
	{
		r = MSK_putcj(task, i, DenseAkMatrix[submatrix.first][i]);	
	}
	
	r = MSK_optimizetrm(task, &trmcode);
	MSK_getxx(task, MSK_SOL_ITR, xx);

	// store expand point
	std::vector<double> exp(NUMVAR);
	for(int i=0 ; i<NUMVAR ; ++i)
	{
		exp[i] = xx[i];		// storing expansion point current pattern 
	}
	ExpandPoint[submatrix.first] = exp;	


	// Node Partitioning, stored in NodeContainer
	for(auto& node : mPGEquivalentNodes)
	{
		if(node->is_connect_I && node->thermal_grid_index == submatrix.first)
		{
			NodeContainer[submatrix.first].emplace_back((node->x_loc_idx - 1) * mY_loc_max + (node->y_loc_idx - 1));
		}
	}

	ThermalGridIndexContainer.emplace_back(submatrix.first);

/*
	if(mPGEquivalentNodes[specificIdx]->thermal_grid_index == submatrix.first)
		OnePointThermalAnalysis(specificIdx, xx, AkMatrix, 1);
*/
}

	TotalTime += (double)(clock()-start) / (CLOCKS_PER_SEC);

	std::cout << "Storing expand point complete !" << std::endl;

	double MatrixTime(0);
	double MosekTime(0);
	clock_t matrixstart;
	clock_t simulationStart;
	double SimulationTime(0);

// Run each node for LP
// for(unsigned int each_node_idx = 0 ; each_node_idx < mPGEquivalentNodeSize ; ++each_node_idx)
// for(unsigned int each_node_idx = specificIdx ; each_node_idx < specificIdx+1 ; ++each_node_idx)
for(auto& each_container : NodeContainer)
{

	std::vector<double> exp = ExpandPoint[each_container.first];
	std::vector<double> TemperatureAdjointValue(thermalgridsize);

	for(auto& thermalGridIndex : ThermalGridIndexContainer)
	{
		double adjointValue(0);
		for(unsigned int row = 0 ; row < mCurrentSourceSize ; ++row) 
		{
			adjointValue += DenseAkMatrix[thermalGridIndex][row] * exp[row];
		}
		TemperatureAdjointValue[thermalGridIndex] = adjointValue;
	}

for(auto& each_node_idx : each_container.second)
{
//if(mPGEquivalentNodes[each_node_idx]->is_connect_I)
//{
	start = clock();

	matrixstart = clock();

	for(int i = 0 ; i < NUMVAR ; ++i)
	{
		unsigned int index = mPGEquivalentNodeSize * i + each_node_idx;	// rowsize * column index + row index
		LPC[i] = cG0INV_Q[index];
	}

	/*
	for(unsigned int i = 0 ; i < mCurrentSourceSize * mCurrentSourceSize ; ++i)  
	{
		cHjmatrix[i] = 0;
	}
	*/

	for(auto& each_grid : MyPGConductanceSubMatrix)		
	{
		MySparseMatrixClass::SparseVector MyPreRow; 	// row each_node_idx of G0INV * G0,k

		for(unsigned int each_col = 0 ; each_col < mPGEquivalentNodeSize ; ++each_col)
		{
			double adjointValue(0);
			for(unsigned int i = 0 ; i < each_grid.second.getSparseVector(each_col).getVector().size() ; ++i)
			{
				adjointValue += each_grid.second.getSparseVector(each_col).getVector()[i] * mEachRowPGInverseMatrix[each_node_idx][each_grid.second.getSparseVector(each_col).getLocation()[i]];
			}
			if(adjointValue != 0)
			{
				MyPreRow.PushBackElement(each_col, adjointValue);
			}
		}

		for(unsigned int each_col = 0 ; each_col < mCurrentSourceSize ; ++each_col)	  // row each_node_idx of Dk = (G0INV * G0,k) * (G0INV * Q)
		{
			double adjointValue(0);
			for(unsigned int i = 0 ; i < MyPreRow.getVector().size() ; ++i)
			{
				unsigned int index = MyPreRow.getLocation()[i] + each_col * mPGEquivalentNodeSize;
				if( cG0INV_Q[index] != 0)
				{
					adjointValue += MyPreRow.getVector()[i] * cG0INV_Q[index];
				}
			}
			DkRow[each_col] = adjointValue;
		}
		/*
		for(unsigned int row = 0 ; row < mCurrentSourceSize ; ++row)  // ak * Dk
		{
			for(unsigned int col = 0 ; col < mCurrentSourceSize ; ++col)
			{
				cHjmatrix[row + col * mCurrentSourceSize] += DenseAkMatrix[each_grid.first][row] * DkRow[col];
			}
		}
		*/
		double TransAdjointValue(0);
		for(unsigned int row = 0 ; row < mCurrentSourceSize ; ++row) 
		{
			TransAdjointValue += DkRow[row] * exp[row];
		}

		for(unsigned int row = 0 ; row < mCurrentSourceSize ; ++row) 
		{
			LPC[row] += TemperatureAdjointValue[each_grid.first] * DkRow[row] + TransAdjointValue * DenseAkMatrix[each_grid.first][row];
		}
	}

	// std::vector<double> exp = ExpandPoint[mPGEquivalentNodes[each_node_idx]->thermal_grid_index];

	//double VoltageShift(0);
	/*
	for(unsigned int row = 0 ; row < mCurrentSourceSize ; ++row)
	{
		double adjointValue(0);
		//double adjointValueShift(0);
		for(unsigned int col = 0 ; col < mCurrentSourceSize ; ++col)
		{
			unsigned int idx     =  col + row * mCurrentSourceSize;
			unsigned int tranidx =  row + col * mCurrentSourceSize; 
			// double SumHjtran = exp[col] * cHjmatrix[tranidx];
			double SumHj_Hjtran = exp[col] * (cHjmatrix[idx] + cHjmatrix[tranidx]);
			adjointValue += SumHj_Hjtran;
			//adjointValueShift += SumHjtran;
		}
		LPC[row] += adjointValue;
		//VoltageShift += adjointValueShift * exp[row];
	}
	*/
	MatrixTime += (double)(clock()-matrixstart) / (CLOCKS_PER_SEC);
	//-----------------------------------------------------------------------------------------------
	// Mosek
	
	clock_t mosek_start;
	mosek_start = clock();
	for(unsigned int row = 0 ; row < mCurrentSourceSize ; ++row)
	{
		r = MSK_putcj(task, row, LPC[row]);
	}
	r = MSK_optimizetrm(task, &trmcode);
	MSK_getxx(task, MSK_SOL_ITR, xx);
	MosekTime += (double)(clock()-mosek_start) / (CLOCKS_PER_SEC);
	TotalTime += (double)(clock()-start) / (CLOCKS_PER_SEC);
	
	//double IRdrop;
	//MSK_getprimalobj(task, MSK_SOL_ITR, &IRdrop);
	//VdropSolution_Origin.emplace_back(IRdrop - VoltageShift);

	simulationStart = clock();
	//pattern[each_node_idx] = 1;
	Eigen::SparseMatrix<double> sys_mat;
	EigenDenseVector CurrentSourceVector(mPGEquivalentNodeSize);
	sys_mat = SimulationForConductanceMatrix<double>(AkMatrix, xx, thermalgridsize, mStampNodeList, CurrentSourceVector);
	// Solve
	Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> tmp_solver;
	tmp_solver.analyzePattern(sys_mat);   
	tmp_solver.factorize(sys_mat);
	EigenDenseVector VdropVector = tmp_solver.solve(CurrentSourceVector);
	SimulationTime += (double)(clock() - simulationStart) / (CLOCKS_PER_SEC);
	VdropSolution.emplace_back(VdropVector[each_node_idx]);
	// OnePointThermalAnalysis(specificIdx, xx, AkMatrix, 2);
	VdropIndex.emplace_back(each_node_idx);
}

}

	// std::cout << "Matrix Time = " << MatrixTime << std::endl;

	std::cout << "Complete Efficient Expand Point Verification !" << std::endl;

// Output solution 
	
	unsigned int violations(0);
	
	double sum(0);
	for(auto& vdrop : VdropSolution) { sum +=vdrop; }

	for(unsigned int i = 0 ; i < VdropIndex.size() ; ++i)
	{
		for(unsigned int j = i + 1 ; j < VdropIndex.size() ; ++j)
		{
			if(VdropIndex[i] > VdropIndex[j])
			{
				unsigned int idx_tmp = VdropIndex[i];
				VdropIndex[i] = VdropIndex[j];
				VdropIndex[j] = idx_tmp;

				double sol_tmp = VdropSolution[i];
				VdropSolution[i] = VdropSolution[j];
				VdropSolution[j] = sol_tmp;
			}
		}
	}


	std::ofstream fout("../LPresult/" + mCktName + "_EfficientExpandLPVdropMap_Simulation");
	int sol_idx(0);
	for(unsigned int i = 0 ; i < mPGEquivalentNodeSize ; ++i)
	{
		if(mPGEquivalentNodes[i]->is_connect_I)
		{
			fout << mPGEquivalentNodes[i]->x_loc_idx << " " << mPGEquivalentNodes[i]->y_loc_idx 
				 << " " << std::setprecision(15) << VdropSolution[sol_idx] << std::endl;
		
			if(VdropSolution[sol_idx] > mSupplyVoltage * 0.1)
			{
				++violations;
			}
			++sol_idx;
		}
	}
	fout.close();

	fout.open("../LPresult/" + mCktName + "_EfficientExpandLP_TimeInfo_Simulation");
	fout << "Total verification time = " << mConstructDataTime + TotalTime + SimulationTime << " (s) " << std::endl;
	fout << "Average voltage drop = " << std::setprecision(15) << sum / VdropSolution.size() << " (V) " << std::endl;
	fout << "Maximum voltage drop = " << std::setprecision(15) << *std::max_element(VdropSolution.begin(), VdropSolution.end()) << " (V) " << std::endl;
	fout << "Total violations = " << violations << std::endl;
	fout << "Matrix Time = " << MatrixTime << " (s) " << std::endl;
	fout << "Mosek Time  = " << MosekTime << " (s) " << std::endl;
	fout << "Simulation Time = " << SimulationTime << " (s) " << std::endl;
	fout.close();

	/*
	
	sum = 0;
	for(auto& vdrop : VdropSolution_Origin) { sum +=vdrop; }

	fout.open("../LPresult/" + mCktName + "_EfficientExpandLPVdropMap");
	sol_idx = 0;
	violations = 0;
	for(unsigned int i = 0 ; i < mPGEquivalentNodeSize ; ++i)
	{
		if(mPGEquivalentNodes[i]->is_connect_I)
		{
			fout << mPGEquivalentNodes[i]->x_loc_idx << " " << mPGEquivalentNodes[i]->y_loc_idx 
				 << " " << std::setprecision(15) << VdropSolution_Origin[sol_idx] << std::endl;
		
			if(VdropSolution_Origin[sol_idx] > mSupplyVoltage * 0.1)
			{
				++violations;
			}
			++sol_idx;
		}
	}
	fout.close();


	fout.open("../LPresult/" + mCktName + "_EfficientExpandLP_TimeInfo");
	fout << "Total verification time = " << mConstructDataTime + TotalTime << " (s) " << std::endl;
	fout << "Average voltage drop = " << std::setprecision(15) << sum / VdropSolution_Origin.size() << " (V) " << std::endl;
	fout << "Maximum voltage drop = " << std::setprecision(15) << *std::max_element(VdropSolution_Origin.begin(), VdropSolution_Origin.end()) << " (V) " << std::endl;
	fout << "Total violations = " << violations << std::endl;
	fout << "Matrix Time = " << MatrixTime << " (s) " << std::endl;
	fout << "Mosek Time  = " << MosekTime << " (s) " << std::endl;
	fout.close();

	*/
	
// delete memory

	MSK_deletetask(&task);
	MSK_deleteenv(&env);

	delete [] cAkMatrix;
	delete [] cG0INV_Q;
	// delete [] cHjmatrix;
	delete [] bkc;
	delete [] blc;	
	delete [] buc;
	delete [] bkx;
	delete [] blx;	
	delete [] bux;
	delete [] xx;

}

// 2018.5.9
// A function for calculating system matrix
template<typename RetType>
Eigen::SparseMatrix<RetType> PGSolver::
SimulationForConductanceMatrix
(const Eigen::MatrixXd & AkMatrix, 
 const double *xx,	
 unsigned int thermalgridsize,
 const STAMPLIST & stampList,
 EigenDenseVector & CurrentSourceVector)
{
	EigenDenseVector TemperatureVector(thermalgridsize);
	CurrentSourceVector.setZero();
	for(unsigned int i = 0 ; i < mCurrentSourceSize ; ++i)
	{
		CurrentSourceVector[mCurrentSourceIdx[i]] = xx[i];
	}
	TemperatureVector = AkMatrix * CurrentSourceVector;

	// Stamping
	MATRIX tmpMatrix;
	tmpMatrix.resize(mPGEquivalentNodeSize);
	for (auto& stamp_node : mStampNodeList)
	{
		for (auto& nei_node : stamp_node.first->neighbor_NodeRes)    
		{
			if (nei_node.first->is_grounded == true)
			{
				for (auto& each_res : nei_node.second)
				{
					tmpMatrix[stamp_node.second][stamp_node.second] += 1 / (each_res->res_value * (1 + mR_temp_coeff * TemperatureVector[each_res->thermal_grid_index]));
				}
			}
			else if (nei_node.first->is_connect_V == true)
			{
				for (auto& each_res : nei_node.second)
				{				
					tmpMatrix[stamp_node.second][stamp_node.second] += 1 / (each_res->res_value * (1 + mR_temp_coeff * TemperatureVector[each_res->thermal_grid_index]));
				}
			}
			else if (nei_node.first->is_connect_V == false && nei_node.first->is_grounded == false)
			{
				for (auto& each_res : nei_node.second)
				{
					tmpMatrix[stamp_node.second][stamp_node.second] += 1 / (each_res->res_value * (1 + mR_temp_coeff * TemperatureVector[each_res->thermal_grid_index]));
					tmpMatrix[stamp_node.second][mStampNodeList[nei_node.first]] += -1 / (each_res->res_value * (1 + mR_temp_coeff * TemperatureVector[each_res->thermal_grid_index]));
				}
			}
		}
	}
	std::vector<triplet> tripletList;
	for(unsigned int row_idx = 0; row_idx < mPGEquivalentNodeSize; ++row_idx)
    {
        for(auto& row_data : tmpMatrix[row_idx])
        {           
            tripletList.push_back(triplet(row_idx, row_data.first, row_data.second));   
        }
    }
    Eigen::SparseMatrix<RetType> sys_mat(mPGEquivalentNodeSize, mPGEquivalentNodeSize);
	sys_mat.setFromTriplets(tripletList.begin(), tripletList.end());
	tripletList.clear();
	return sys_mat;
}

// 2018.4.21
// expand point method 2
// 1. Direct
// 2. Iterative
// 3. Direct expand point
void PGSolver::PGVerificationExpandPoint_Efficient_Method2(int Method)
{
	clock_t verificationStart; 
	double verificationTime(0);
	
	verificationStart = clock();
	
	EigenDenseMatrix AkMatrix = mDenseInverseThermalMatrix * Eigen::MatrixXd(mElectricalGridtoThermalGridIncidenceMatrix);
	
	unsigned int thermalgridsize = mXCut * mYCut * mZlayer;
	
	unsigned int ArraySize = thermalgridsize * mPGEquivalentNodeSize;
	
	double *cAkMatrix = new double[ArraySize]();
	
	Eigen::Map<Eigen::MatrixXd>( cAkMatrix, AkMatrix.rows(), AkMatrix.cols() ) = AkMatrix;
	
	double thermalResistanceFactor = mSupplyVoltage * mR_temp_coeff; 	

	std::vector< std::vector<double> > DenseAkMatrix(thermalgridsize);
	
	for(unsigned int i = 0 ; i < thermalgridsize ; ++i )
	{	
		DenseAkMatrix[i].resize(mCurrentSourceSize);
	}
	
	for(unsigned int row = 0 ; row < thermalgridsize ; ++row)
	{
		for(unsigned int col = 0 ; col < mPGEquivalentNodeSize ; ++col)
		{
			unsigned int index = row + col * thermalgridsize;
			if(cAkMatrix[index] != 0)
			{
				DenseAkMatrix[row][mCurrentSourceIdxMapping[col]] = thermalResistanceFactor * cAkMatrix[index];
			}
		}
	}
	
	unsigned int NUMCON = mGlobalConstraintSize;  	
	unsigned int NUMVAR = mCurrentSourceSize;	   
	
	MSKboundkeye *bkc = new MSKboundkeye[NUMCON];
	for(unsigned int i=0 ; i<NUMCON ; ++i) { bkc[i] = MSK_BK_RA; }
	double *blc = new double[NUMCON]();  
	double *buc = new double[NUMCON]; 
	for(unsigned int i=0 ; i<NUMCON ; ++i) { buc[i] = mGlobalCurrentConstraint[i]; }

	mGlobalCurrentConstraintIncidenceMatrix.makeCompressed();
	
	int    *aptrb    = mGlobalCurrentConstraintIncidenceMatrix.outerIndexPtr();
	int    *asub     = mGlobalCurrentConstraintIncidenceMatrix.innerIndexPtr();
	double *aval     = mGlobalCurrentConstraintIncidenceMatrix.valuePtr();
	
	MSKboundkeye *bkx = new MSKboundkeye[NUMVAR];
	for(unsigned int i=0 ; i<NUMVAR ; ++i) { bkx[i] = MSK_BK_RA; }
	double *blx = new double[NUMVAR](); 
	double *bux = new double[NUMVAR];
	for(unsigned int i=0 ; i<NUMVAR ; ++i) { bux[i] = mLocalCurrentConstraint[i]; }
	
	double *xx = new double[NUMVAR];  // for retriving solution (current pattern)
	
	MSKenv_t env = NULL;
	MSKtask_t task = NULL;
	MSKrescodee r;
	MSKrescodee trmcode;

	r = MSK_makeenv(&env, NULL);
	r = MSK_maketask(env, NUMCON, NUMVAR, &task);
	r = MSK_appendcons(task, NUMCON); 
	r = MSK_appendvars(task, NUMVAR);

	for(unsigned int j=0; j<NUMVAR && r == MSK_RES_OK; ++j)
	{	
		r = MSK_putvarbound(task, j, bkx[j], blx[j], bux[j]); 
		r = MSK_putacol(task, j, aptrb[j+1]-aptrb[j], asub+aptrb[j], aval+aptrb[j]);
	}
	
	for(unsigned int i=0; i<NUMCON && r==MSK_RES_OK; ++i)
		r = MSK_putconbound(task, i, bkc[i], blc[i], buc[i]); 

	
	r = MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MAXIMIZE);  

	std::cout << "Finding expand point ..." << std::endl; 
	
	// unsigned int specificIdx(5050);
	std::unordered_map<unsigned int, EigenSparseMatrix > SystemMatrix;
	std::unordered_map<unsigned int, std::vector<unsigned int> > NodeContainer;
	std::unordered_map<unsigned int, EigenDenseVector> ExpandPointCurrentPattern;

	// for stamping
	unsigned int mapping_idx(0);
	for (auto& node : mPGEquivalentNodes)
	{
		mStampNodeList.insert(std::pair<Node*, unsigned int>(node, mapping_idx));
		mapping_idx++;
	}

for(auto & thermalGridIndex : mThermalGridIndexConstainer)
{
	for(unsigned int i = 0 ; i < NUMVAR ; ++i)
	{
		r = MSK_putcj(task, i, DenseAkMatrix[thermalGridIndex][i]);	
	}
	
	r = MSK_optimizetrm(task, &trmcode);
	MSK_getxx(task, MSK_SOL_ITR, xx);

	EigenDenseVector CurrentSourceVector(mPGEquivalentNodeSize);
	SystemMatrix[thermalGridIndex] = 
		SimulationForConductanceMatrix<double>(AkMatrix, xx, thermalgridsize, mStampNodeList, CurrentSourceVector);

	if(Method == 3) {	// for Direct expand point method 2 
		ExpandPointCurrentPattern[thermalGridIndex] = CurrentSourceVector;
	}

	// Node Partitioning, stored in NodeContainer
	for(auto& node : mPGEquivalentNodes)
	{
		if(node->is_connect_I && node->thermal_grid_index == thermalGridIndex)
		{
			NodeContainer[thermalGridIndex].emplace_back((node->x_loc_idx - 1) * mY_loc_max + (node->y_loc_idx - 1));
		}
	}

	/*
	if(mPGEquivalentNodes[specificIdx]->thermal_grid_index == thermalGridIndex)
		OnePointThermalAnalysis(specificIdx, xx, AkMatrix, 1);
	*/
}

	// temporary verification datas

	Eigen::VectorXd pattern(mPGEquivalentNodeSize);	// determine which row of inverse matrix to solve
	pattern.setZero();

	Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;	// solve Row of inverse matrix
	EigenSparseMatrix sysMatrix;	// approximated conductance matrix
	EigenDenseVector Row;
	double IRdrop;	

	std::vector<double> VdropSolution;
	std::vector<unsigned int> VdropIndex;
	clock_t simulationStart;
	double SimulationTime(0);

	std::cout << "Optimizing ..." << std::endl; 

// Method = 1 : Direct
if(Method == 1)
{

for(auto & each_container : NodeContainer)	// nodes in each thermal grid
{
	sysMatrix = SystemMatrix[each_container.first];
	solver.analyzePattern(sysMatrix);   
	solver.factorize(sysMatrix);

	for(auto & each_node : each_container.second)	// verify each node in each thermal grid
	{
		pattern[each_node] = 1;
		Row = solver.solve(pattern);
		for(unsigned int i = 0 ; i < NUMVAR ; ++i)
		{
			r = MSK_putcj(task, i, Row[mCurrentSourceIdx[i]]);	// cj^T
		}
		r = MSK_optimizetrm(task, &trmcode);
		MSK_getxx(task, MSK_SOL_ITR, xx);
		
		simulationStart = clock();
		Eigen::SparseMatrix<double> sys_mat;
		EigenDenseVector CurrentSourceVector(mPGEquivalentNodeSize);
		sys_mat = SimulationForConductanceMatrix<double>(AkMatrix, xx, thermalgridsize, mStampNodeList, CurrentSourceVector);
		// Solve
		Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> tmp_solver;
		tmp_solver.analyzePattern(sys_mat);   
		tmp_solver.factorize(sys_mat);
		Row = tmp_solver.solve(pattern);
		double simVdrop = Row.adjoint() * CurrentSourceVector;
		SimulationTime += (double)(clock() - simulationStart) / (CLOCKS_PER_SEC);
		VdropSolution.emplace_back(simVdrop);
		VdropIndex.emplace_back(each_node);
		pattern[each_node] = 0;
		// OnePointThermalAnalysis(specificIdx, xx, AkMatrix, 2);
	}
}

}
// Method = 2 : Iterative
else if(Method == 2)
{
double *next_sol = new double[NUMVAR];
EigenDenseVector CurrentSourceVector(mPGEquivalentNodeSize);
for(auto & each_container : NodeContainer)	// nodes in each thermal grid
{
	sysMatrix = SystemMatrix[each_container.first];
	solver.analyzePattern(sysMatrix);   
	solver.factorize(sysMatrix);

	for(auto & each_node : each_container.second)	// verify each node in each thermal grid
	{
		pattern[each_node] = 1;
		Row = solver.solve(pattern);
		
		unsigned int MaxIteration(20), Iter(0);
		double Tolerance(1e-6), Max_Tolerance(1e-6);
		for(unsigned int i = 0 ; i < NUMVAR ; ++i)
		{
			r = MSK_putcj(task, i, Row[mCurrentSourceIdx[i]]);	// cj^T
		}
		r = MSK_optimizetrm(task, &trmcode);
		MSK_getprimalobj(task, MSK_SOL_ITR, &IRdrop);	// Obtain first solution
		MSK_getxx(task, MSK_SOL_ITR, xx);
		while(Iter <= MaxIteration)
		{
			for(unsigned int i = 0 ; i < NUMVAR ; ++i)
			{
				r = MSK_putcj(task, i, Row[mCurrentSourceIdx[i]]);	// cj^T
			}
			r = MSK_optimizetrm(task, &trmcode);
			MSK_getprimalobj(task, MSK_SOL_ITR, &IRdrop);	// Obtain first solution
			MSK_getxx(task, MSK_SOL_ITR, next_sol);

			double err_sum(0);
			std::vector<double> each_err;
			for(unsigned int i = 0 ; i < NUMVAR ; ++i)
			{
				err_sum += pow(next_sol[i] - xx[i], 2);
				each_err.emplace_back(fabs(next_sol[i] - xx[i]));
			}
			err_sum = pow(err_sum, 0.5);
			if(err_sum <= Tolerance && 
				*std::max_element(each_err.begin(), each_err.end()) <= Max_Tolerance)
			{
				break;
			}
			Eigen::SparseMatrix<double> sys_mat = 
				SimulationForConductanceMatrix<double>(AkMatrix, xx, thermalgridsize, mStampNodeList, CurrentSourceVector);
			solver.analyzePattern(sys_mat);   
			solver.factorize(sys_mat);
			Row = solver.solve(pattern); 
			for(unsigned int i = 0 ; i < NUMVAR ; ++i)
			{
				xx[i] = next_sol[i];
			}	
			++Iter;
		}
		VdropSolution.emplace_back(IRdrop);
		VdropIndex.emplace_back(each_node);
		pattern[each_node] = 0;
	}
}
delete [] next_sol;
}

// Direct Expansion Point Method with IR drop recalculation (in thesis)
else if (Method == 3) 
{

for(auto & each_container : NodeContainer)	// nodes in each thermal grid
{
	sysMatrix = SystemMatrix[each_container.first];
	solver.analyzePattern(sysMatrix);   
	solver.factorize(sysMatrix);

	for(auto & each_node : each_container.second)	// verify each node in each thermal grid
	{
		pattern[each_node] = 1;
		Row = solver.solve(pattern);
		IRdrop = Row.adjoint() * ExpandPointCurrentPattern[each_container.first];
		VdropSolution.emplace_back(IRdrop);
		VdropIndex.emplace_back(each_node);
		pattern[each_node] = 0;
		// OnePointThermalAnalysis(specificIdx, xx, AkMatrix, 2);
	}
}

}

	verificationTime = (double) ( clock() - verificationStart ) / CLOCKS_PER_SEC;

	std::cout << "Complete Efficient Expand Point Method 2 Verification !" << std::endl;

// Output solution 
	
	unsigned int violations(0);
	
	double sum(0);
	for(auto& vdrop : VdropSolution) { sum +=vdrop; }

	for(unsigned int i = 0 ; i < VdropIndex.size() ; ++i)
	{
		for(unsigned int j = i + 1 ; j < VdropIndex.size() ; ++j)
		{
			if(VdropIndex[i] > VdropIndex[j])
			{
				unsigned int idx_tmp = VdropIndex[i];
				VdropIndex[i] = VdropIndex[j];
				VdropIndex[j] = idx_tmp;

				double sol_tmp = VdropSolution[i];
				VdropSolution[i] = VdropSolution[j];
				VdropSolution[j] = sol_tmp;
			}
		}
	}

	std::string T_method, V_method;
	if(Method == 1) {
		V_method = "_EfficientExpandLPVdropMap_Method2_Direct";
		T_method = "_EfficientExpandLP_TimeInfo_Method2_Direct";
	}
	else if (Method == 2)
	{
		V_method = "_EfficientExpandLPVdropMap_Method2_Iterative";
		T_method = "_EfficientExpandLP_TimeInfo_Method2_Iterative";
	}
	else if (Method == 3)
	{
		V_method = "_EfficientExpandLPVdropMap_DirectExpandPointMethod_Simulation";
		T_method = "_EfficientExpandLP_TimeInfo_DirectExpandPointMethod_Simulation";
	}
	// output solution from small node index to large node index
	std::ofstream fout("../LPresult/" + mCktName + V_method);
	
	for(unsigned int i = 0 ; i < VdropSolution.size() ; ++i)
	{
		fout << mPGEquivalentNodes[VdropIndex[i]]->x_loc_idx << " " << mPGEquivalentNodes[VdropIndex[i]]->y_loc_idx 
			 << " " << std::setprecision(15) << VdropSolution[i] << std::endl;
		
		if(VdropSolution[i] > mSupplyVoltage * 0.1)
		{
			++violations;
		}
	}
	fout.close();

	fout.open("../LPresult/" + mCktName + T_method);
	fout << "Total verification time = " << mConstructDataTime + verificationTime << " (s) " << std::endl;
	fout << "Average voltage drop = " << std::setprecision(15) << sum / VdropSolution.size() << " (V) " << std::endl;
	fout << "Maximum voltage drop = " << std::setprecision(15) << *std::max_element(VdropSolution.begin(), VdropSolution.end()) << " (V) " << std::endl;
	fout << "Total violations = " << violations << std::endl;
	fout << "Simulation Time = " << SimulationTime << " (s) " << std::endl;
	fout.close();
	
// delete memory

	MSK_deletetask(&task);
	MSK_deleteenv(&env);

	delete [] cAkMatrix;
	delete [] bkc;
	delete [] blc;	
	delete [] buc;
	delete [] bkx;
	delete [] blx;	
	delete [] bux;
	delete [] xx;

}


// Direct Expansion Point Method without IR drop recalculation (not in thesis)
void PGSolver::PGVerificationExpandPoint_EfficientExpand()
{

	
	std::map<unsigned int, MySparseMatrixClass::SparseMatrix> MyPGConductanceSubMatrix;	 

	for(auto& each_conductance_submatrix : mPG_EigenConductanceSubMatrix)
	{
		(each_conductance_submatrix.second).makeCompressed();	

		int    *aptrb    = (each_conductance_submatrix.second).outerIndexPtr();	  
		int    *asub     = (each_conductance_submatrix.second).innerIndexPtr();	 
		double *aval     = (each_conductance_submatrix.second).valuePtr();		
		unsigned int NonZerosIndexAdder(0);

		MySparseMatrixClass::SparseMatrix tmpSparseMatrix;
		for(unsigned int Col = 0 ; Col < mPGEquivalentNodeSize ; ++Col)
		{
			MySparseMatrixClass::SparseVector tmpSparseVector((aptrb[Col+1] - aptrb[Col]), asub, aval, &NonZerosIndexAdder);
			tmpSparseMatrix.PushBack(tmpSparseVector);
		}
		MyPGConductanceSubMatrix[each_conductance_submatrix.first] = tmpSparseMatrix;
	}
	mPG_EigenConductanceSubMatrix.clear();


	double *cG0INV_Q = new double[mPGEquivalentNodeSize * mCurrentSourceSize]();

	for(unsigned int each_row = 0 ; each_row < mPGEquivalentNodeSize ; ++each_row)
	{
		for(unsigned int each_col = 0 ; each_col < mCurrentSourceSize ; ++each_col)
		{
			double value = mEachRowPGInverseMatrix[each_row][mCurrentSourceIdx[each_col]];
			if(value != 0)
			{
				cG0INV_Q[each_row + each_col*mPGEquivalentNodeSize] = value;
			}
		}
	}

	// Construct INVGT * M_EtoT, find the grid index to calculate the temperature, ak
	EigenDenseMatrix AkMatrix = mDenseInverseThermalMatrix * Eigen::MatrixXd(mElectricalGridtoThermalGridIncidenceMatrix);
	
	// load AkMatrix to sparse vector data structure
	unsigned int thermalgridsize = mXCut * mYCut * mZlayer;
	
	unsigned int ArraySize = thermalgridsize * mPGEquivalentNodeSize;
	
	double *cAkMatrix = new double[ArraySize]();
	
	Eigen::Map<Eigen::MatrixXd>( cAkMatrix, AkMatrix.rows(), AkMatrix.cols() ) = AkMatrix;
	
	std::vector< std::vector<double> > DenseAkMatrix(thermalgridsize);

	double thermalResistanceFactor = mSupplyVoltage * mR_temp_coeff;
	
	for(unsigned int i = 0 ; i < thermalgridsize ; ++i )
	{	
		DenseAkMatrix[i].resize(mCurrentSourceSize);
	}
	
	for(unsigned int row = 0 ; row < thermalgridsize ; ++row)
	{
		for(unsigned int col = 0 ; col < mPGEquivalentNodeSize ; ++col)
		{
			unsigned int index = row + col * thermalgridsize;
			if(cAkMatrix[index] != 0)
			{
				DenseAkMatrix[row][mCurrentSourceIdxMapping[col]] = thermalResistanceFactor * cAkMatrix[index];
			}
		}
	}
	

	int NUMCON = mGlobalConstraintSize;  	
	int NUMVAR = mCurrentSourceSize;	   
	
	MSKboundkeye *bkc = new MSKboundkeye[NUMCON];
	for(int i=0 ; i<NUMCON ; ++i) { bkc[i] = MSK_BK_RA; }
	double *blc = new double[NUMCON]();  
	double *buc = new double[NUMCON]; 
	for(int i=0 ; i<NUMCON ; ++i) { buc[i] = mGlobalCurrentConstraint[i]; }

	mGlobalCurrentConstraintIncidenceMatrix.makeCompressed();
	
	int    *aptrb    = mGlobalCurrentConstraintIncidenceMatrix.outerIndexPtr();
	int    *asub     = mGlobalCurrentConstraintIncidenceMatrix.innerIndexPtr();
	double *aval     = mGlobalCurrentConstraintIncidenceMatrix.valuePtr();
	

	MSKboundkeye *bkx = new MSKboundkeye[NUMVAR];
	for(int i=0 ; i<NUMVAR ; ++i) { bkx[i] = MSK_BK_RA; }
	double *blx = new double[NUMVAR](); 
	double *bux = new double[NUMVAR];
	for(int i=0 ; i<NUMVAR ; ++i) { bux[i] = mLocalCurrentConstraint[i]; }
	
	
	Eigen::VectorXd LPC(NUMVAR); 	
	double *cHjmatrix = new double[mCurrentSourceSize * mCurrentSourceSize]();
	EigenDenseVector DkRow(mCurrentSourceSize);
	 	

	double *xx = new double[NUMVAR];  
	std::vector<double> VdropSolution;


	MSKenv_t env = NULL;
	MSKtask_t task = NULL;
	MSKrescodee r;
	MSKrescodee trmcode;

	r = MSK_makeenv(&env, NULL);
	r = MSK_maketask(env, NUMCON, NUMVAR, &task);
	r = MSK_appendcons(task, NUMCON); 
	r = MSK_appendvars(task, NUMVAR);

	for(int j=0; j<NUMVAR && r == MSK_RES_OK; ++j)
	{	
		r = MSK_putvarbound(task, j, bkx[j], blx[j], bux[j]); 
		r = MSK_putacol(task, j, aptrb[j+1]-aptrb[j], asub+aptrb[j], aval+aptrb[j]);
	}
	
	for(int i=0; i<NUMCON && r==MSK_RES_OK; ++i)
		r = MSK_putconbound(task, i, bkc[i], blc[i], buc[i]); 

	r = MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MAXIMIZE);  

	std::unordered_map<unsigned int, std::vector<double> > ExpandPoint;


for(auto & submatrix : MyPGConductanceSubMatrix)
{
	
	for(int i=0 ; i<NUMVAR ; ++i)
	{
		r = MSK_putcj(task, i, DenseAkMatrix[submatrix.first][i]);	
	}
	
	r = MSK_optimizetrm(task, &trmcode);
	MSK_getxx(task, MSK_SOL_ITR, xx);

	// store expand point
	std::vector<double> exp(NUMVAR);
	for(int i=0 ; i<NUMVAR ; ++i)
	{
		exp[i] = xx[i];
	}
	ExpandPoint[submatrix.first] = exp;

/*
	if(mPGEquivalentNodes[specificIdx]->thermal_grid_index == submatrix.first)
		OnePointThermalAnalysis(specificIdx, xx, AkMatrix, 1);
*/

}


// Run each node for LP
for(unsigned int each_node_idx = 0 ; each_node_idx < mPGEquivalentNodeSize ; ++each_node_idx)
{

if(mPGEquivalentNodes[each_node_idx]->is_connect_I)
{	

	for(int i = 0 ; i < NUMVAR ; ++i)
	{
		unsigned int index = mPGEquivalentNodeSize * i + each_node_idx;	// rowsize * column index + row index
		LPC[i] = cG0INV_Q[index];
	}

	for(unsigned int i = 0 ; i < mCurrentSourceSize * mCurrentSourceSize ; ++i)  
	{
		cHjmatrix[i] = 0;
	}

	for(auto& each_grid : MyPGConductanceSubMatrix)		
	{
		MySparseMatrixClass::SparseVector MyPreRow; 	// G0INV * G0,k

		for(unsigned int each_col = 0 ; each_col < mPGEquivalentNodeSize ; ++each_col)
		{
			double adjointValue(0);
			for(unsigned int i = 0 ; i < each_grid.second.getSparseVector(each_col).getVector().size() ; ++i)
			{
				adjointValue += each_grid.second.getSparseVector(each_col).getVector()[i] * mEachRowPGInverseMatrix[each_node_idx][each_grid.second.getSparseVector(each_col).getLocation()[i]];
			}
			if(adjointValue != 0)
			{
				MyPreRow.PushBackElement(each_col, adjointValue);
			}
		}

		for(unsigned int each_col = 0 ; each_col < mCurrentSourceSize ; ++each_col)	  // Dk = (G0INV * G0,k) * (G0INV * Q)
		{
			double adjointValue(0);
			for(unsigned int i = 0 ; i < MyPreRow.getVector().size() ; ++i)
			{
				unsigned int index = MyPreRow.getLocation()[i] + each_col * mPGEquivalentNodeSize;
				if( cG0INV_Q[index] != 0)
				{
					adjointValue += MyPreRow.getVector()[i] * cG0INV_Q[index];
				}
			}
			DkRow[each_col] = adjointValue;
		}
		
		for(unsigned int row = 0 ; row < mCurrentSourceSize ; ++row)  // ak * Dk
		{
			for(unsigned int col = 0 ; col < mCurrentSourceSize ; ++col)
			{
				cHjmatrix[row + col * mCurrentSourceSize] += DenseAkMatrix[each_grid.first][row] * DkRow[col];
			}
		}

	}

	for(unsigned int i = 0 ; i < mCurrentSourceSize * mCurrentSourceSize ; ++i)
	{
		cHjmatrix[i] = cHjmatrix[i] * thermalResistanceFactor;
	}

	std::vector<double> exp = ExpandPoint[mPGEquivalentNodes[each_node_idx]->thermal_grid_index];

	for(unsigned int row = 0 ; row < mCurrentSourceSize ; ++row)
	{
		for(unsigned int col = 0 ; col < mCurrentSourceSize ; ++col)
		{
			unsigned int idx =  col + row * mCurrentSourceSize;
			LPC[row] += exp[col] * cHjmatrix[idx];
		}
	}

	double Vdrop(0);
	for(int i = 0 ; i < NUMVAR ; ++i) { Vdrop += LPC[i] * exp[i]; } 
	VdropSolution.emplace_back(Vdrop);
	
	// OnePointThermalAnalysis(specificIdx, xx, AkMatrix, 2);
}
	
}

	std::cout << "Complete Efficient Expand Point Verification !" << std::endl;

// Output solution 
	
	unsigned int violations(0);
	double sum(0);
	for(auto& vdrop : VdropSolution) { sum +=vdrop; }

	std::ofstream fout("../LPresult/" + mCktName + "_EfficientExpandLPVdropMapExpand");
	int sol_idx(0);
	for(unsigned int i = 0 ; i < mPGEquivalentNodeSize ; ++i)
	{
		if(mPGEquivalentNodes[i]->is_connect_I)
		{
			fout << mPGEquivalentNodes[i]->x_loc_idx << " " << mPGEquivalentNodes[i]->y_loc_idx 
				 << " " << std::setprecision(15) << VdropSolution[sol_idx] << std::endl;
		
			if(VdropSolution[sol_idx] > mSupplyVoltage * 0.1)
			{
				++violations;
			}
			++sol_idx;
		}
	}
	fout.close();

// delete memory

	MSK_deletetask(&task);
	MSK_deleteenv(&env);

	delete [] cAkMatrix;
	delete [] cG0INV_Q;
	delete [] cHjmatrix;
	delete [] bkc;
	delete [] blc;	
	delete [] buc;
	delete [] bkx;
	delete [] blx;	
	delete [] bux;
	delete [] xx;

}



//	Original LP Verification 
void PGSolver::PGVerificationLP()
{
	
	clock_t start = clock();

	EigenSparseMatrix EigenMatrix(mPGEquivalentNodeSize, mPGEquivalentNodeSize);
	
	std::vector<triplet> tripletList;
	for(unsigned int row_idx = 0; row_idx < mPGEquivalentNodeSize; ++row_idx)
    {
        for(auto& row_data : mPG_OriginConductanceMatrix[row_idx])
        {           
            tripletList.push_back(triplet(row_idx, row_data.first, row_data.second));   
        }
    }
	EigenMatrix.setFromTriplets(tripletList.begin(), tripletList.end());
	tripletList.clear();

	Eigen::VectorXd pattern(mPGEquivalentNodeSize), Row;
	pattern.setZero();

	Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
	solver.analyzePattern(EigenMatrix);   
	solver.factorize(EigenMatrix);
	

	int NUMCON = mGlobalConstraintSize; 
	int NUMVAR = mPGEquivalentNodeSize;			  
	
	MSKboundkeye *bkc = new MSKboundkeye[NUMCON];
	for(int i=0 ; i<NUMCON ; ++i) { bkc[i] = MSK_BK_RA; }

	double *blc = new double[NUMCON]();  
	double *buc = new double[NUMCON]; 
	for(int i=0 ; i<NUMCON ; ++i) { buc[i] = mGlobalCurrentConstraint[i]; }
	 
	mGlobalCurrentConstraintIncidenceMatrix.makeCompressed();

	int    *aptrb    = mGlobalCurrentConstraintIncidenceMatrix.outerIndexPtr();
	int    *asub     = mGlobalCurrentConstraintIncidenceMatrix.innerIndexPtr();
	double *aval     = mGlobalCurrentConstraintIncidenceMatrix.valuePtr();
	
	MSKboundkeye *bkx = new MSKboundkeye[NUMVAR];
	for(int i=0 ; i<NUMVAR ; ++i) { bkx[i] = MSK_BK_RA; }
	double *blx = new double[NUMVAR](); 
	double *bux = new double[NUMVAR];
	for(int i=0 ; i<NUMVAR ; ++i) { bux[i] = mLocalCurrentConstraint[i]; }
	
	
	Eigen::VectorXd LPC; 
	
	double *xx = new double[NUMVAR];    
	std::vector<double> VdropSolution;

	MSKenv_t env = NULL;
	MSKtask_t task = NULL;
	MSKrescodee r;
	MSKrescodee trmcode;	

	r = MSK_makeenv(&env, NULL);
	r = MSK_maketask(env, NUMCON, NUMVAR, &task);
	r = MSK_appendcons(task, NUMCON); 
	r = MSK_appendvars(task, NUMVAR);
	r = MSK_putcfix(task, 0.0);
	r = MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MAXIMIZE); 
	
	for(int j=0 ; j<NUMVAR && r == MSK_RES_OK ; ++j)
	{	
		r = MSK_putvarbound(task, j, bkx[j], blx[j], bux[j]); 
		r = MSK_putacol(task, j, aptrb[j+1]-aptrb[j], asub+aptrb[j], aval+aptrb[j]);
	}

	for(int i=0 ; i<NUMCON && r == MSK_RES_OK ; ++i)
		r = MSK_putconbound(task, i, bkc[i], blc[i], buc[i]); 


for(unsigned int each_node_idx = 0 ; each_node_idx < mPGEquivalentNodeSize ; ++each_node_idx)
{
	
	pattern[each_node_idx] = 1;
	Row = solver.solve(pattern);
	pattern[each_node_idx] = 0;

	
	for(int j=0 ; j<NUMVAR ; ++j)
	{
		r = MSK_putcj(task, j, Row[j]);
	}

	r = MSK_optimizetrm(task, &trmcode);	
	MSK_getxx(task, MSK_SOL_ITR, xx);
	
	double Vdrop(0);
	for(int i = 0 ; i < NUMVAR ; ++i) { Vdrop += Row[i] * xx[i]; } 
	VdropSolution.emplace_back(Vdrop);
	
}
	double LPTime = (double)(clock()-start) / (CLOCKS_PER_SEC);

	std::cout << "Complete Thermal-less Verification !" << std::endl;
	
	unsigned int violations(0);

	double sum(0);
	for(auto& vdrop : VdropSolution) { sum +=vdrop; }

	std::ofstream fout("../LPresult/" + mCktName + "_LPVdropMap");
	for(unsigned int i = 0 ; i < mPGEquivalentNodeSize ; ++i)
	{
		fout << mPGEquivalentNodes[i]->x_loc_idx << " " << mPGEquivalentNodes[i]->y_loc_idx << " " << std::setprecision(15) << VdropSolution[i] << std::endl;
		
		if(VdropSolution[i] > mSupplyVoltage * 0.1)
		{
			++violations;
		}

	}
	fout.close();

	fout.open("../LPresult/" + mCktName + "_LP_TimeInfo");
	fout << "Total verification time = " << mConstructDataTime + LPTime << " (S) " << std::endl;
	fout << "Average voltage drop = " << std::setprecision(15) << sum / VdropSolution.size() << " (V) " << std::endl;
	fout << "Maximum voltage drop = " << std::setprecision(15) << *std::max_element(VdropSolution.begin(), VdropSolution.end()) << " (V) " << std::endl;
	fout << "Total violations = " << violations << std::endl;
	fout.close();

	MSK_deletetask(&task);
	MSK_deleteenv(&env);

	delete [] bkc;
	delete [] blc;	
	delete [] buc;
	delete [] bkx;
	delete [] blx;	
	delete [] bux;
	delete [] xx;
	
}


// Efficient Thermal-less LP
void PGSolver::PGVerificationLP_Efficient()
{

	clock_t start;

	EigenSparseMatrix EigenMatrix(mPGEquivalentNodeSize, mPGEquivalentNodeSize);
	
	std::vector<triplet> tripletList;
	for(unsigned int row_idx = 0; row_idx < mPGEquivalentNodeSize; ++row_idx)
    {
        for(auto& row_data : mPG_OriginConductanceMatrix[row_idx])
        {           
            tripletList.push_back(triplet(row_idx, row_data.first, row_data.second));   
        }
    }
	EigenMatrix.setFromTriplets(tripletList.begin(), tripletList.end());
	tripletList.clear();

	Eigen::VectorXd pattern(mPGEquivalentNodeSize), Row;
	pattern.setZero();

	Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
	solver.analyzePattern(EigenMatrix);   
	solver.factorize(EigenMatrix);
	
	
	int NUMCON = mGlobalConstraintSize; 
	int NUMVAR = mCurrentSourceSize;			  
	
	MSKboundkeye *bkc = new MSKboundkeye[NUMCON];
	for(int i=0 ; i<NUMCON ; ++i) { bkc[i] = MSK_BK_RA; }

	double *blc = new double[NUMCON]();  
	double *buc = new double[NUMCON]; 
	for(int i=0 ; i<NUMCON ; ++i) { buc[i] = mGlobalCurrentConstraint[i]; }
	 
	mGlobalCurrentConstraintIncidenceMatrix.makeCompressed();
	
	int    *aptrb    = mGlobalCurrentConstraintIncidenceMatrix.outerIndexPtr();
	int    *asub     = mGlobalCurrentConstraintIncidenceMatrix.innerIndexPtr();
	double *aval     = mGlobalCurrentConstraintIncidenceMatrix.valuePtr();
	
	MSKboundkeye *bkx = new MSKboundkeye[NUMVAR];
	for(int i=0 ; i<NUMVAR ; ++i) { bkx[i] = MSK_BK_RA; }
	double *blx = new double[NUMVAR](); 
	double *bux = new double[NUMVAR];
	for(int i=0 ; i<NUMVAR ; ++i) { bux[i] = mLocalCurrentConstraint[i]; }
	
	
	Eigen::VectorXd LPC; 
	
	double *xx = new double[NUMVAR];    
	std::vector<double> VdropSolution;

	MSKenv_t env = NULL;
	MSKtask_t task = NULL;
	MSKrescodee r;
	MSKrescodee trmcode;	

	r = MSK_makeenv(&env, NULL);
	r = MSK_maketask(env, NUMCON, NUMVAR, &task);
	r = MSK_appendcons(task, NUMCON); 
	r = MSK_appendvars(task, NUMVAR);
	r = MSK_putcfix(task, 0.0);
	r = MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MAXIMIZE); 
	
	for(int j=0 ; j<NUMVAR && r == MSK_RES_OK ; ++j)
	{	
		r = MSK_putvarbound(task, j, bkx[j], blx[j], bux[j]); 
		r = MSK_putacol(task, j, aptrb[j+1]-aptrb[j], asub+aptrb[j], aval+aptrb[j]);
	}

	for(int i=0 ; i<NUMCON && r == MSK_RES_OK ; ++i)
		r = MSK_putconbound(task, i, bkc[i], blc[i], buc[i]); 

double TotalTime(0);

for(unsigned int each_node_idx = 0 ; each_node_idx < mPGEquivalentNodeSize ; ++each_node_idx)
{
	start = clock();

if(mPGEquivalentNodes[each_node_idx]->is_connect_I)
{	
	pattern[each_node_idx] = 1;
	Row = solver.solve(pattern);
	pattern[each_node_idx] = 0;

	for(int j=0 ; j<NUMVAR ; ++j)
	{
		r = MSK_putcj(task, j, Row[mCurrentSourceIdx[j]]);
	}

	r = MSK_optimizetrm(task, &trmcode);	
	double IRdrop;
	MSK_getprimalobj(task, MSK_SOL_ITR, &IRdrop);
	TotalTime += (double)(clock()-start) / (CLOCKS_PER_SEC);
	VdropSolution.emplace_back(IRdrop);
}

}

	std::cout << "Complete Efficient Thermal-less Verification !" << std::endl;
	
	unsigned int violations(0);
	double sum(0);
	for(auto& vdrop : VdropSolution) { sum +=vdrop; }

	std::ofstream fout("../LPresult/" + mCktName + "_EfficientLPVdropMap");
	int sol_idx(0);
	for(unsigned int i = 0 ; i < mPGEquivalentNodeSize ; ++i)
	{
		if(mPGEquivalentNodes[i]->is_connect_I)
		{
			fout << mPGEquivalentNodes[i]->x_loc_idx << " " << mPGEquivalentNodes[i]->y_loc_idx 
				 << " " << std::setprecision(15) << VdropSolution[sol_idx] << std::endl;
		
			if(VdropSolution[sol_idx] > mSupplyVoltage * 0.1)
			{
				++violations;
			}
			++sol_idx;
		}
	}
	fout.close();

	fout.open("../LPresult/" + mCktName + "_EfficientLP_TimeInfo");
	fout << "Total verification time = " << mConstructDataTime + TotalTime << " (S) " << std::endl;
	fout << "Average voltage drop = " << std::setprecision(15) << sum / VdropSolution.size() << " (V) " << std::endl;
	fout << "Maximum voltage drop = " << std::setprecision(15) << *std::max_element(VdropSolution.begin(), VdropSolution.end()) << " (V) " << std::endl;
	fout << "Total violations = " << violations << std::endl;
	fout.close();

	MSK_deletetask(&task);
	MSK_deleteenv(&env);

	delete [] bkc;
	delete [] blc;	
	delete [] buc;
	delete [] bkx;
	delete [] blx;	
	delete [] bux;
	delete [] xx;
	
}


void PGSolver::MonteCarlo()
{
	

	EigenSparseMatrix EigenPGConductanceMatrix(mPGEquivalentNodeSize, mPGEquivalentNodeSize);

	Eigen::VectorXd CurrentPattern(mPGEquivalentNodeSize);
	
	Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
	
	std::vector<triplet> tripletList;
	for(unsigned int row_idx = 0; row_idx < mPGEquivalentNodeSize; ++row_idx)
    {
        for(auto& row_data : mPG_OriginConductanceMatrix[row_idx])
        {           
            tripletList.push_back(triplet(row_idx, row_data.first, row_data.second));   
        }
    }
	EigenPGConductanceMatrix.setFromTriplets(tripletList.begin(), tripletList.end());
	tripletList.clear();

	// solver.analyzePattern(EigenPGConductanceMatrix);   
	// solver.factorize(EigenPGConductanceMatrix);

	std::vector< std::vector<double> > CurrentSetRecorder;
	CurrentSetRecorder.resize(mCurrentSourceSize);


	for(unsigned int times = 0 ; times < 1000 ; ++times)
	{

		std::vector<double> GlobalCurrentConstraintRecorder(mGlobalCurrentConstraint);
		
		// construct current set in random
		for(unsigned int each_source = 0 ; each_source < mCurrentSourceSize ; ++each_source)
		{
			double CurrentValue;
			bool exit = false;
			
			while(!exit) {
				
				// CurrentValue = randMToN(0, mLocalCurrentConstraint[each_source]);
				unsigned int counter(0);
				for(unsigned int i = 0 ; i < CurrentSetRecorder[each_source].size() ; ++i)
				{
					if(CurrentSetRecorder[each_source][i] == CurrentValue) {
						break;
					}
					else {
						++counter;
					}
				}

				if(counter == CurrentSetRecorder[each_source].size()-1) {
					exit = true;
				}

			}
			CurrentSetRecorder[each_source].emplace_back(CurrentValue);

		}
		// verify each node
		for(unsigned int each_node_idx = 0 ; each_node_idx < mPGEquivalentNodeSize ; ++each_node_idx)
		{


		}

	}
	
}



// EXP (iterative) with IR drop recalculation using Matlab
void PGSolver::PGVerificationExpandPoint_Efficient_Simulation_Matlab()
{

	clock_t matrixStart, simulationStart; 
	double LPTime(0), MatrixTime(0), SimulationTime(0);

	matrixStart = clock();

	std::map<unsigned int, MySparseMatrixClass::SparseMatrix> MyPGConductanceSubMatrix;	 

	for(auto& each_conductance_submatrix : mPG_EigenConductanceSubMatrix)
	{
		(each_conductance_submatrix.second).makeCompressed();	

		int    *aptrb    = (each_conductance_submatrix.second).outerIndexPtr();	  
		int    *asub     = (each_conductance_submatrix.second).innerIndexPtr();	 
		double *aval     = (each_conductance_submatrix.second).valuePtr();		
		unsigned int NonZerosIndexAdder(0);

		MySparseMatrixClass::SparseMatrix tmpSparseMatrix;
		for(unsigned int Col = 0 ; Col < mPGEquivalentNodeSize ; ++Col)
		{
			MySparseMatrixClass::SparseVector tmpSparseVector((aptrb[Col+1] - aptrb[Col]), asub, aval, &NonZerosIndexAdder);
			tmpSparseMatrix.PushBack(tmpSparseVector);
		}
		MyPGConductanceSubMatrix[each_conductance_submatrix.first] = tmpSparseMatrix;
	}
	mPG_EigenConductanceSubMatrix.clear();


	double *cG0INV_Q = new double[mPGEquivalentNodeSize * mCurrentSourceSize]();

	for(unsigned int each_row = 0 ; each_row < mPGEquivalentNodeSize ; ++each_row)
	{
		for(unsigned int each_col = 0 ; each_col < mCurrentSourceSize ; ++each_col)
		{
			double value = mEachRowPGInverseMatrix[each_row][mCurrentSourceIdx[each_col]];
			if(value != 0)
			{
				cG0INV_Q[each_row + each_col*mPGEquivalentNodeSize] = value;
			}
		}
	}

	// Construct INVGT * M_EtoT, find the grid index to calculate the temperature, ak
	EigenDenseMatrix AkMatrix = mDenseInverseThermalMatrix * Eigen::MatrixXd(mElectricalGridtoThermalGridIncidenceMatrix);
	
	// load AkMatrix to sparse vector data structure
	unsigned int thermalgridsize = mXCut * mYCut * mZlayer;
	
	unsigned int ArraySize = thermalgridsize * mPGEquivalentNodeSize;
	
	double *cAkMatrix = new double[ArraySize]();
	
	Eigen::Map<Eigen::MatrixXd>( cAkMatrix, AkMatrix.rows(), AkMatrix.cols() ) = AkMatrix;
	

	std::vector< std::vector<double> > DenseAkMatrix(thermalgridsize);	// Q * Ak
	
	for(unsigned int i=0 ; i<thermalgridsize ; ++i )
		DenseAkMatrix[i].resize(mCurrentSourceSize);

	// Factor
	double sFactor = -1 * mR_temp_coeff * mSupplyVoltage;	
	
	for(unsigned int row = 0 ; row < thermalgridsize ; ++row)
	{
		for(unsigned int col = 0 ; col < mPGEquivalentNodeSize ; ++col)
		{
			unsigned int index = row + col * thermalgridsize;
			if(cAkMatrix[index] != 0)
			{
				DenseAkMatrix[row][mCurrentSourceIdxMapping[col]] = sFactor * cAkMatrix[index];
			}
		}
	}
	
	
	// start matlab engine

	Engine *ep;

	if (!(ep = engOpen(""))) {
		fprintf(stderr, "\nCan't start MATLAB engine\n");
		return ;
	}
	else {
		std::cout << "Started matlab engine !" << std::endl;
	}

	// engEvalString(ep, "LASTN = maxNumCompThreads(1);");

	// construct matlab A matrix
	//-----------------------------------------------------------------------------------------------------------
	// load mMatlabAmatrix & lower bound constraint to Amatrix  => 0 < Ui < Ig
	unsigned int Amatrixsize = 2 * mGlobalConstraintSize * mCurrentSourceSize;
	double *Amatrix = new double[Amatrixsize];
	Eigen::Map<Eigen::MatrixXd>( Amatrix, mMatlabAmatrix.rows(), mMatlabAmatrix.cols() ) = Eigen::MatrixXd(mMatlabAmatrix);
    
  	// load Amatrix to matAmatrix
	mxArray *matAmatrix = NULL;
	matAmatrix = mxCreateDoubleMatrix(2 * mGlobalConstraintSize, mCurrentSourceSize, mxREAL);
	memcpy((void *)mxGetPr(matAmatrix), (void *)Amatrix, sizeof(double) * Amatrixsize);
	engPutVariable(ep, "matAmatrix", matAmatrix);	// put variable to matlab workspace 
	//-----------------------------------------------------------------------------------------------------------
	
	// construct matlab b vector
	double *b = new double[2*mGlobalConstraintSize]();
	for(unsigned int i=0 ; i < mGlobalConstraintSize ; ++i)
	{
		b[i] = mGlobalCurrentConstraint[i];
	}
	mxArray *matb = NULL;
	matb = mxCreateDoubleMatrix(1, 2 * mGlobalConstraintSize, mxREAL);
	memcpy((void *)mxGetPr(matb), (void *)b, sizeof(double) * 2 * mGlobalConstraintSize);
	engPutVariable(ep, "matb", matb);

	// construct matlab local upper and lower variable constraints => 0 < i < Il
	//-----------------------------------------------------------------------------------------------------------
	double *clower = new double[mCurrentSourceSize]();
	double *cupper = new double[mCurrentSourceSize];
	for(unsigned int i=0 ; i < mCurrentSourceSize ; ++i)
	{
		cupper[i] = mLocalCurrentConstraint[i];
	}
	mxArray *matlower = NULL;
	mxArray *matupper = NULL;
	matlower = mxCreateDoubleMatrix(1, mCurrentSourceSize, mxREAL);
	matupper = mxCreateDoubleMatrix(1, mCurrentSourceSize, mxREAL);
	memcpy((void *)mxGetPr(matlower), (void *)clower, sizeof(double) * mCurrentSourceSize);
	memcpy((void *)mxGetPr(matupper), (void *)cupper, sizeof(double) * mCurrentSourceSize);
	engPutVariable(ep, "matlower", matlower);
	engPutVariable(ep, "matupper", matupper);
	//-----------------------------------------------------------------------------------------------------------
	
	//-----------------------------------------------------------------------------------------------------------
	// construct matlab linear terms 
	double *cf = new double[mCurrentSourceSize];
	mxArray *matf = NULL;
	matf = mxCreateDoubleMatrix(1, mCurrentSourceSize, mxREAL);

	// solution
	double *xx;	// for retrieving current pattern
	std::vector<double> VdropSolution;

	
	engEvalString(ep, "lptime = 0;");	// for calculating matlab time 
	EigenDenseVector DkRow(mCurrentSourceSize);


// Finding expansion point ------------------------------------------------------------------------
	std::unordered_map<unsigned int, std::vector<double> > ExpandPoint;
	std::vector<unsigned int> ThermalGridIndexContainer;
	std::unordered_map<unsigned int, std::vector<unsigned int> > NodeContainer;
	std::vector<unsigned int> VdropIndex;

for(auto & submatrix : MyPGConductanceSubMatrix)
{
	
	for(unsigned int i=0 ; i<mCurrentSourceSize ; ++i)
	{
		cf[i] = DenseAkMatrix[submatrix.first][i];	
	}

	memcpy((void *)mxGetPr(matf), (void *)cf, sizeof(double) * mCurrentSourceSize);
	engPutVariable(ep, "matf", matf);	
	engEvalString(ep, "tStart = cputime;");
	engEvalString(ep, "xx = linprog(matf,matAmatrix,matb,[],[],matlower,matupper);");
	engEvalString(ep, "lptime = cputime - tStart;");	
	LPTime += *(double *)mxGetPr(engGetVariable(ep, "lptime"));

	
	// store expand point
	xx = (double *)mxGetPr(engGetVariable(ep, "xx"));
	std::vector<double> exp(mCurrentSourceSize);
	for(unsigned int i=0 ; i < mCurrentSourceSize ; ++i)
	{
		exp[i] = xx[i];
		//std::cout << exp[i] << std::endl;
	}
	ExpandPoint[submatrix.first] = exp;


	// Node Partitioning, stored in NodeContainer
	for(auto& node : mPGEquivalentNodes)
	{
		if(node->is_connect_I && node->thermal_grid_index == submatrix.first)
		{
			NodeContainer[submatrix.first].emplace_back((node->x_loc_idx - 1) * mY_loc_max + (node->y_loc_idx - 1));
		}
	}
	ThermalGridIndexContainer.emplace_back(submatrix.first);

	// Temperature map
/*
	if(mPGEquivalentNodes[specificIdx]->thermal_grid_index == submatrix.first)
		OnePointThermalAnalysis(specificIdx, xx, AkMatrix, 1);
*/
}
//  ---------------------------------------------------------------------------------------------------------------------
	
	MatrixTime += (double)(clock()-matrixStart) / (CLOCKS_PER_SEC);

	std::cout << "Optimizing ..." << std::endl;

	std::vector<unsigned int> ConvergeContainer;	// record the number of iterations
	double *inv_row = new double[mCurrentSourceSize]; // for storing the inverse conductance row at initial temperature
	std::unordered_map<unsigned int, double*> DkRowsContainer;	// for storing each Dk rows in each thermal grid
	for(auto& thermalGridIndex : ThermalGridIndexContainer)
	{
		double *tmpDkRow = new double[mCurrentSourceSize];
		// storing pointer to each DkRow
		DkRowsContainer[thermalGridIndex] = tmpDkRow;	// storing the pointer 
	}
	
for(auto& each_container : NodeContainer)
{
	matrixStart = clock();

	std::vector<double> exp = ExpandPoint[each_container.first];  // select initial expansion point
	std::vector<double> TemperatureAdjointValue(thermalgridsize); // for storing initial temperature induced by expansion point

	for(auto& thermalGridIndex : ThermalGridIndexContainer)
	{
		double adjointValue(0);
		for(unsigned int row = 0 ; row < mCurrentSourceSize ; ++row) 
		{
			adjointValue += DenseAkMatrix[thermalGridIndex][row] * exp[row];
		}
		TemperatureAdjointValue[thermalGridIndex] = adjointValue;
	}

	MatrixTime += (double)(clock()-matrixStart) / (CLOCKS_PER_SEC);

for(auto& each_node_idx : each_container.second)
{
//if(mPGEquivalentNodes[each_node_idx]->is_connect_I)
//{
	
	matrixStart = clock();

	// initial conditions
	unsigned int MaxIteration = 10;
	double Tolerance = 1e-5;
	//std::vector<double> InitialCurrentPattern;
	std::vector<double> InitialTemperatureAdjointValue(thermalgridsize);  // temperature terms for iterations
	double InitialVdrop; // objective solutions for iterations
	double *tmpDkRow;    // temporary pointer for retriving each DkRow

	for(unsigned int i = 0 ; i < mCurrentSourceSize ; ++i)
	{
		unsigned int index = mPGEquivalentNodeSize * i + each_node_idx;	// rowsize * column index + row index
		cf[i] = -1 * cG0INV_Q[index];
		inv_row[i] = -1 * cG0INV_Q[index];
	}

	for(auto& each_grid : MyPGConductanceSubMatrix)		
	{
		MySparseMatrixClass::SparseVector MyPreRow; 	// row each_node_idx of G0INV * G0,k

		tmpDkRow = DkRowsContainer[each_grid.first];	// retriving the precomputed DkRow in DkRowsContainer

		for(unsigned int each_col = 0 ; each_col < mPGEquivalentNodeSize ; ++each_col)
		{
			double adjointValue(0);
			for(unsigned int i = 0 ; i < each_grid.second.getSparseVector(each_col).getVector().size() ; ++i)
			{
				adjointValue += each_grid.second.getSparseVector(each_col).getVector()[i] * mEachRowPGInverseMatrix[each_node_idx][each_grid.second.getSparseVector(each_col).getLocation()[i]];
			}
			if(adjointValue != 0)
			{
				MyPreRow.PushBackElement(each_col, adjointValue);
			}
		}

		for(unsigned int each_col = 0 ; each_col < mCurrentSourceSize ; ++each_col)	  // row each_node_idx of Dk = (G0INV * G0,k) * (G0INV * Q)
		{
			double adjointValue(0);
			for(unsigned int i = 0 ; i < MyPreRow.getVector().size() ; ++i)
			{
				unsigned int index = MyPreRow.getLocation()[i] + each_col * mPGEquivalentNodeSize;
				if( cG0INV_Q[index] != 0)
				{
					adjointValue += MyPreRow.getVector()[i] * cG0INV_Q[index];
				}
			}
			// DkRow[each_col] = adjointValue;
			tmpDkRow[each_col] = adjointValue;
		}

		double TransAdjointValue(0);
		for(unsigned int row = 0 ; row < mCurrentSourceSize ; ++row) 
		{
			TransAdjointValue += tmpDkRow[row] * exp[row];
		}

		for(unsigned int row = 0 ; row < mCurrentSourceSize ; ++row) 
		{
			cf[row] += TemperatureAdjointValue[each_grid.first] * tmpDkRow[row] + TransAdjointValue * DenseAkMatrix[each_grid.first][row];
		}
		
	}

	memcpy((void *)mxGetPr(matf), (void *)cf, sizeof(double) * mCurrentSourceSize);
	engPutVariable(ep, "matf", matf);	
	engEvalString(ep, "tStart = cputime;");
	engEvalString(ep, "[xx fval] = linprog(matf,matAmatrix,matb,[],[],matlower,matupper);");
	engEvalString(ep, "lptime = cputime - tStart;");	
	LPTime += *(double *)mxGetPr(engGetVariable(ep, "lptime"));

	InitialVdrop = -1 * *(double*)mxGetPr(engGetVariable(ep, "fval")); // first objective solution
	xx = (double*)mxGetPr(engGetVariable(ep, "xx"));	// first solution current pattern

	for(auto& thermalGridIndex : ThermalGridIndexContainer)
	{
		double adjointValue(0);
		for(unsigned int row = 0 ; row < mCurrentSourceSize ; ++row) 
		{
			adjointValue += DenseAkMatrix[thermalGridIndex][row] * xx[row];
		}
		InitialTemperatureAdjointValue[thermalGridIndex] = adjointValue;	// for the first iteration
	}

	// Iterations, solutions from the first iteration until converge
	for(unsigned int iter = 1 ; iter <= MaxIteration ; ++iter)
	{

		for(unsigned int i = 0 ; i < mCurrentSourceSize ; ++i)
		{
			cf[i] = 0;
		}

		for(auto& each_grid : DkRowsContainer)
		{
			double *curDkRow = each_grid.second;
			double TransAdjointValue(0);
			for(unsigned int row = 0 ; row < mCurrentSourceSize ; ++row) 
			{
				TransAdjointValue += curDkRow[row] * xx[row];
			}

			for(unsigned int row = 0 ; row < mCurrentSourceSize ; ++row) 
			{
				cf[row] += InitialTemperatureAdjointValue[each_grid.first] * curDkRow[row] + TransAdjointValue * DenseAkMatrix[each_grid.first][row];
			}
		}

		for(unsigned int i = 0 ; i < mCurrentSourceSize ; ++i)
		{
			cf[i] += inv_row[i];
		}

		// Matlab
		memcpy((void *)mxGetPr(matf), (void *)cf, sizeof(double) * mCurrentSourceSize);
		engPutVariable(ep, "matf", matf);	
		engEvalString(ep, "tStart = cputime;");
		engEvalString(ep, "[xx fval] = linprog(matf,matAmatrix,matb,[],[],matlower,matupper);");
		engEvalString(ep, "lptime = cputime - tStart;");	
		LPTime += *(double *)mxGetPr(engGetVariable(ep, "lptime"));

		if( fabs(InitialVdrop - (-1 * *(double *)mxGetPr(engGetVariable(ep, "fval")))) <= Tolerance )	// if objective solution error < tolerance, then converge
		{
			//std::cout << "Converge at iteration " << iter << std::endl;
			xx = (double *)mxGetPr(engGetVariable(ep, "xx"));  // get solution
			ConvergeContainer.emplace_back(iter);	// record the iterations 
			break;
		}
		else {
			InitialVdrop = -1 * *(double *)mxGetPr(engGetVariable(ep, "fval")); // update objective solution
			xx = (double *)mxGetPr(engGetVariable(ep, "xx"));  // update solution current pattern for next iteration
			for(auto& thermalGridIndex : ThermalGridIndexContainer)
			{
				double adjointValue(0);
				for(unsigned int row = 0 ; row < mCurrentSourceSize ; ++row) 
				{
					adjointValue += DenseAkMatrix[thermalGridIndex][row] * xx[row];
				}
				InitialTemperatureAdjointValue[thermalGridIndex] = adjointValue;  // update temperature terms for next iteration
			}
		}

		if(iter == MaxIteration)
		{
			std::cout << "exceed max iteration" << std::endl;
		}

	}

	MatrixTime += (double)(clock()-matrixStart) / (CLOCKS_PER_SEC);
	
	simulationStart = clock();
	//pattern[each_node_idx] = 1;
	Eigen::SparseMatrix<double> sys_mat;
	EigenDenseVector CurrentSourceVector(mPGEquivalentNodeSize);
	sys_mat = SimulationForConductanceMatrix<double>(AkMatrix, xx, thermalgridsize, mStampNodeList, CurrentSourceVector);
	// Solve
	Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> tmp_solver;
	tmp_solver.analyzePattern(sys_mat);   
	tmp_solver.factorize(sys_mat);
	EigenDenseVector VdropVector = tmp_solver.solve(CurrentSourceVector);
	SimulationTime += (double)(clock() - simulationStart) / (CLOCKS_PER_SEC);
	VdropSolution.emplace_back(VdropVector[each_node_idx]);
	// OnePointThermalAnalysis(specificIdx, xx, AkMatrix, 2);
	VdropIndex.emplace_back(each_node_idx);

}

}

	std::cout << "Complete Efficient Expand Point Verification !" << std::endl;

// Output solution 
	
	unsigned int violations(0);
	
	double sum(0);
	for(auto& vdrop : VdropSolution) { sum +=vdrop; }

	for(unsigned int i = 0 ; i < VdropIndex.size() ; ++i)
	{
		for(unsigned int j = i + 1 ; j < VdropIndex.size() ; ++j)
		{
			if(VdropIndex[i] > VdropIndex[j])
			{
				unsigned int idx_tmp = VdropIndex[i];
				VdropIndex[i] = VdropIndex[j];
				VdropIndex[j] = idx_tmp;

				double sol_tmp = VdropSolution[i];
				VdropSolution[i] = VdropSolution[j];
				VdropSolution[j] = sol_tmp;
			}
		}
	}

	std::ofstream fout("../LPresult/" + mCktName + "_EfficientExpandLPVdropMap_Simulation_Matlab");
	int sol_idx(0);
	for(unsigned int i = 0 ; i < mPGEquivalentNodeSize ; ++i)
	{
		if(mPGEquivalentNodes[i]->is_connect_I)
		{
			fout << mPGEquivalentNodes[i]->x_loc_idx << " " << mPGEquivalentNodes[i]->y_loc_idx 
				 << " " << std::setprecision(15) << VdropSolution[sol_idx] << std::endl;
		
			if(VdropSolution[sol_idx] > mSupplyVoltage * 0.1)
			{
				++violations;
			}
			++sol_idx;
		}
	}
	fout.close();

	fout.open("../LPresult/" + mCktName + "_EfficientExpandLP_TimeInfo_Simulation_Matlab");
	fout << "Total verification time = " << mConstructDataTime + LPTime + MatrixTime + SimulationTime << " (s) " << std::endl;
	fout << "Average voltage drop = " << std::setprecision(15) << sum / VdropSolution.size() << " (V) " << std::endl;
	fout << "Maximum voltage drop = " << std::setprecision(15) << *std::max_element(VdropSolution.begin(), VdropSolution.end()) << " (V) " << std::endl;
	fout << "Total violations = " << violations << std::endl;
	fout << "Matrix Time = " << MatrixTime << " (s) " << std::endl;
	fout << "Matlab Time  = " << LPTime << " (s) " << std::endl;
	fout << "Simulation Time = " << SimulationTime << " (s) " << std::endl;
	fout.close();


	fout.open(mCktName + "_Converge");
	for(auto& i : ConvergeContainer)
	{
		fout << i << std::endl;
	}
	fout.close();

	mxDestroyArray(matf);
	mxDestroyArray(matAmatrix);
	mxDestroyArray(matb);
	mxDestroyArray(matlower);
	mxDestroyArray(matupper);
	engClose(ep);
	std::cout << "Closed Matlab engine" << std::endl;

	for(auto& each_grid : DkRowsContainer)
	{
		delete each_grid.second;
	}

	delete [] cG0INV_Q;
	delete [] cAkMatrix;
	delete [] Amatrix;
	delete [] b;
	delete [] clower;
	delete [] cupper;
	delete [] cf;
	delete [] inv_row;
}



// EXP (iterative) without IR drop recalculation using matlab
void PGSolver::PGVerificationExpandPoint_Efficient_Matlab()
{

	clock_t matrixStart; 
	double LPTime(0), MatrixTime(0);

	matrixStart = clock();

	std::map<unsigned int, MySparseMatrixClass::SparseMatrix> MyPGConductanceSubMatrix;	 

	for(auto& each_conductance_submatrix : mPG_EigenConductanceSubMatrix)
	{
		(each_conductance_submatrix.second).makeCompressed();	

		int    *aptrb    = (each_conductance_submatrix.second).outerIndexPtr();	  
		int    *asub     = (each_conductance_submatrix.second).innerIndexPtr();	 
		double *aval     = (each_conductance_submatrix.second).valuePtr();		
		unsigned int NonZerosIndexAdder(0);

		MySparseMatrixClass::SparseMatrix tmpSparseMatrix;
		for(unsigned int Col = 0 ; Col < mPGEquivalentNodeSize ; ++Col)
		{
			MySparseMatrixClass::SparseVector tmpSparseVector((aptrb[Col+1] - aptrb[Col]), asub, aval, &NonZerosIndexAdder);
			tmpSparseMatrix.PushBack(tmpSparseVector);
		}
		MyPGConductanceSubMatrix[each_conductance_submatrix.first] = tmpSparseMatrix;
	}
	mPG_EigenConductanceSubMatrix.clear();


	double *cG0INV_Q = new double[mPGEquivalentNodeSize * mCurrentSourceSize]();  // G0_inv * Q

	for(unsigned int each_row = 0 ; each_row < mPGEquivalentNodeSize ; ++each_row)
	{
		for(unsigned int each_col = 0 ; each_col < mCurrentSourceSize ; ++each_col)
		{
			double value = mEachRowPGInverseMatrix[each_row][mCurrentSourceIdx[each_col]];
			if(value != 0)
			{
				cG0INV_Q[each_row + each_col*mPGEquivalentNodeSize] = value;
			}
		}
	}

	// Construct INVGT * M_EtoT, find the grid index to calculate the temperature, ak
	EigenDenseMatrix AkMatrix = mDenseInverseThermalMatrix * Eigen::MatrixXd(mElectricalGridtoThermalGridIncidenceMatrix);
	
	// load AkMatrix to sparse vector data structure
	unsigned int thermalgridsize = mXCut * mYCut * mZlayer;
	
	unsigned int ArraySize = thermalgridsize * mPGEquivalentNodeSize;
	
	double *cAkMatrix = new double[ArraySize]();
	
	Eigen::Map<Eigen::MatrixXd>( cAkMatrix, AkMatrix.rows(), AkMatrix.cols() ) = AkMatrix;
	

	std::vector< std::vector<double> > DenseAkMatrix(thermalgridsize);	 // alpha * Vdd * Gt_inv  * M_EtoT * Q
	
	for(unsigned int i=0 ; i<thermalgridsize ; ++i )
		DenseAkMatrix[i].resize(mCurrentSourceSize);

	// Factor
	double sFactor = -1 * mR_temp_coeff * mSupplyVoltage;	
	
	for(unsigned int row = 0 ; row < thermalgridsize ; ++row)
	{
		for(unsigned int col = 0 ; col < mPGEquivalentNodeSize ; ++col)
		{
			unsigned int index = row + col * thermalgridsize;
			if(cAkMatrix[index] != 0)
			{
				DenseAkMatrix[row][mCurrentSourceIdxMapping[col]] = sFactor * cAkMatrix[index];		// alpha * Vdd * Gt_inv  * M_EtoT * Q
			}
		}
	}
	
	
	// start matlab engine

	Engine *ep;

	if (!(ep = engOpen(""))) {
		fprintf(stderr, "\nCan't start MATLAB engine\n");
		return ;
	}
	else {
		std::cout << "Started matlab engine !" << std::endl;
	}

	// engEvalString(ep, "LASTN = maxNumCompThreads(1);");

	// construct matlab A matrix
	//-----------------------------------------------------------------------------------------------------------
	// load mMatlabAmatrix & lower bound constraint to Amatrix  => 0 < Ui < Ig
	unsigned int Amatrixsize = 2 * mGlobalConstraintSize * mCurrentSourceSize;
	double *Amatrix = new double[Amatrixsize];
	Eigen::Map<Eigen::MatrixXd>( Amatrix, mMatlabAmatrix.rows(), mMatlabAmatrix.cols() ) = Eigen::MatrixXd(mMatlabAmatrix);
    
  	// load Amatrix to matAmatrix
	mxArray *matAmatrix = NULL;
	matAmatrix = mxCreateDoubleMatrix(2 * mGlobalConstraintSize, mCurrentSourceSize, mxREAL);
	memcpy((void *)mxGetPr(matAmatrix), (void *)Amatrix, sizeof(double) * Amatrixsize);
	engPutVariable(ep, "matAmatrix", matAmatrix);	// put variable to matlab workspace 
	//-----------------------------------------------------------------------------------------------------------
	
	// construct matlab b vector
	double *b = new double[2*mGlobalConstraintSize]();
	for(unsigned int i=0 ; i < mGlobalConstraintSize ; ++i)
	{
		b[i] = mGlobalCurrentConstraint[i];
	}
	mxArray *matb = NULL;
	matb = mxCreateDoubleMatrix(1, 2 * mGlobalConstraintSize, mxREAL);
	memcpy((void *)mxGetPr(matb), (void *)b, sizeof(double) * 2 * mGlobalConstraintSize);
	engPutVariable(ep, "matb", matb);

	// construct matlab local upper and lower variable constraints => 0 < i < Il
	//-----------------------------------------------------------------------------------------------------------
	double *clower = new double[mCurrentSourceSize]();
	double *cupper = new double[mCurrentSourceSize];
	for(unsigned int i=0 ; i < mCurrentSourceSize ; ++i)
	{
		cupper[i] = mLocalCurrentConstraint[i];
	}
	mxArray *matlower = NULL;
	mxArray *matupper = NULL;
	matlower = mxCreateDoubleMatrix(1, mCurrentSourceSize, mxREAL);
	matupper = mxCreateDoubleMatrix(1, mCurrentSourceSize, mxREAL);
	memcpy((void *)mxGetPr(matlower), (void *)clower, sizeof(double) * mCurrentSourceSize);
	memcpy((void *)mxGetPr(matupper), (void *)cupper, sizeof(double) * mCurrentSourceSize);
	engPutVariable(ep, "matlower", matlower);
	engPutVariable(ep, "matupper", matupper);
	//-----------------------------------------------------------------------------------------------------------
	
	//-----------------------------------------------------------------------------------------------------------
	// construct matlab linear terms 
	double *cf = new double[mCurrentSourceSize];
	mxArray *matf = NULL;
	matf = mxCreateDoubleMatrix(1, mCurrentSourceSize, mxREAL);

	// solution
	double *xx;	// for retrieving current pattern
	std::vector<double> VdropSolution;

	
	engEvalString(ep, "lptime = 0;");	// for calculating matlab time 
	EigenDenseVector DkRow(mCurrentSourceSize);


// Finding expansion point ------------------------------------------------------------------------
	std::unordered_map<unsigned int, std::vector<double> > ExpandPoint;
	std::vector<unsigned int> ThermalGridIndexContainer;
	std::unordered_map<unsigned int, std::vector<unsigned int> > NodeContainer;
	std::vector<unsigned int> VdropIndex;

for(auto & submatrix : MyPGConductanceSubMatrix)
{
	
	for(unsigned int i=0 ; i<mCurrentSourceSize ; ++i)
	{
		cf[i] = DenseAkMatrix[submatrix.first][i];	
	}

	memcpy((void *)mxGetPr(matf), (void *)cf, sizeof(double) * mCurrentSourceSize);
	engPutVariable(ep, "matf", matf);	
	engEvalString(ep, "tStart = cputime;");
	engEvalString(ep, "xx = linprog(matf,matAmatrix,matb,[],[],matlower,matupper);");
	engEvalString(ep, "lptime = cputime - tStart;");	
	LPTime += *(double *)mxGetPr(engGetVariable(ep, "lptime"));

	
	// store expansion point current pattern
	xx = (double *)mxGetPr(engGetVariable(ep, "xx"));
	std::vector<double> exp(mCurrentSourceSize);
	for(unsigned int i=0 ; i < mCurrentSourceSize ; ++i)
	{
		exp[i] = xx[i];
		//std::cout << exp[i] << std::endl;
	}
	ExpandPoint[submatrix.first] = exp;


	// Node Partitioning, stored in NodeContainer
	for(auto& node : mPGEquivalentNodes)
	{
		if(node->is_connect_I && node->thermal_grid_index == submatrix.first)
		{
			NodeContainer[submatrix.first].emplace_back((node->x_loc_idx - 1) * mY_loc_max + (node->y_loc_idx - 1));
		}
	}
	ThermalGridIndexContainer.emplace_back(submatrix.first);

	// Temperature map
/*
	if(mPGEquivalentNodes[specificIdx]->thermal_grid_index == submatrix.first)
		OnePointThermalAnalysis(specificIdx, xx, AkMatrix, 1);
*/
}
//  ---------------------------------------------------------------------------------------------------------------------
	
	MatrixTime += (double)(clock()-matrixStart) / (CLOCKS_PER_SEC);

	std::cout << "Optimizing ..." << std::endl;
	
	EigenDenseVector Hj_trans_row(mCurrentSourceSize);
	std::vector<unsigned int> ConvergeContainer;	// record the number of iterations
	double *inv_row = new double[mCurrentSourceSize]; // for storing the inverse conductance row at initial temperature
	std::unordered_map<unsigned int, double*> DkRowsContainer;	// for storing each Dk rows in each thermal grid
	for(auto& thermalGridIndex : ThermalGridIndexContainer)
	{
		double *tmpDkRow = new double[mCurrentSourceSize];
		// storing pointer to each DkRow
		DkRowsContainer[thermalGridIndex] = tmpDkRow;	// storing the pointer 
	}
	
for(auto& each_container : NodeContainer)
{
	matrixStart = clock();

	std::vector<double> exp = ExpandPoint[each_container.first];  // select initial expansion point
	std::vector<double> TemperatureAdjointValue(thermalgridsize); // for storing initial temperature terms induced by expansion point

	for(auto& thermalGridIndex : ThermalGridIndexContainer)
	{
		double adjointValue(0);
		for(unsigned int row = 0 ; row < mCurrentSourceSize ; ++row) 
		{
			adjointValue += DenseAkMatrix[thermalGridIndex][row] * exp[row];
		}
		TemperatureAdjointValue[thermalGridIndex] = adjointValue;
	}

	MatrixTime += (double)(clock()-matrixStart) / (CLOCKS_PER_SEC);

for(auto& each_node_idx : each_container.second)
{
//if(mPGEquivalentNodes[each_node_idx]->is_connect_I)
//{
	
	matrixStart = clock();

	// initial conditions
	unsigned int MaxIteration = 10;  // maximum iterations
	double Tolerance = 1e-5;	// tolerance of the objective solution of LP (1e-5 V)
	//std::vector<double> InitialCurrentPattern;
	std::vector<double> InitialTemperatureAdjointValue(thermalgridsize);  // temperature terms for iterations
	double InitialVdrop; // objective solutions for iterations
	double *tmpDkRow;    // temporary pointer for retriving each DkRow
	Hj_trans_row.setZero();

	for(unsigned int i = 0 ; i < mCurrentSourceSize ; ++i)
	{
		unsigned int index = mPGEquivalentNodeSize * i + each_node_idx;	// rowsize * column index + row index
		cf[i] = -1 * cG0INV_Q[index];
		inv_row[i] = -1 * cG0INV_Q[index];
	}

	for(auto& each_grid : MyPGConductanceSubMatrix)		
	{
		MySparseMatrixClass::SparseVector MyPreRow; 	// row each_node_idx of G0INV * G0,k

		tmpDkRow = DkRowsContainer[each_grid.first];	// retriving the precomputed DkRow in DkRowsContainer

		for(unsigned int each_col = 0 ; each_col < mPGEquivalentNodeSize ; ++each_col)
		{
			double adjointValue(0);
			for(unsigned int i = 0 ; i < each_grid.second.getSparseVector(each_col).getVector().size() ; ++i)
			{
				adjointValue += each_grid.second.getSparseVector(each_col).getVector()[i] * mEachRowPGInverseMatrix[each_node_idx][each_grid.second.getSparseVector(each_col).getLocation()[i]];
			}
			if(adjointValue != 0)
			{
				MyPreRow.PushBackElement(each_col, adjointValue);
			}
		}

		for(unsigned int each_col = 0 ; each_col < mCurrentSourceSize ; ++each_col)	  // row each_node_idx of Dk = (G0INV * G0,k) * (G0INV * Q)
		{
			double adjointValue(0);
			for(unsigned int i = 0 ; i < MyPreRow.getVector().size() ; ++i)
			{
				unsigned int index = MyPreRow.getLocation()[i] + each_col * mPGEquivalentNodeSize;
				if( cG0INV_Q[index] != 0)
				{
					adjointValue += MyPreRow.getVector()[i] * cG0INV_Q[index];
				}
			}
			// DkRow[each_col] = adjointValue;
			tmpDkRow[each_col] = adjointValue;
		}

		double TransAdjointValue(0);
		for(unsigned int row = 0 ; row < mCurrentSourceSize ; ++row) 
		{
			TransAdjointValue += tmpDkRow[row] * exp[row];
		}

		for(unsigned int row = 0 ; row < mCurrentSourceSize ; ++row) 
		{
			double tmp = TransAdjointValue * DenseAkMatrix[each_grid.first][row];
			Hj_trans_row[row] += tmp;
			cf[row] += TemperatureAdjointValue[each_grid.first] * tmpDkRow[row] + tmp;
		}
		
	}

	memcpy((void *)mxGetPr(matf), (void *)cf, sizeof(double) * mCurrentSourceSize);
	engPutVariable(ep, "matf", matf);	
	engEvalString(ep, "tStart = cputime;");
	engEvalString(ep, "[xx fval] = linprog(matf,matAmatrix,matb,[],[],matlower,matupper);");
	engEvalString(ep, "lptime = cputime - tStart;");	
	LPTime += *(double *)mxGetPr(engGetVariable(ep, "lptime"));

	InitialVdrop = -1 * *(double*)mxGetPr(engGetVariable(ep, "fval")); // first objective solution
	xx = (double*)mxGetPr(engGetVariable(ep, "xx"));	// first solution current pattern

	for(auto& thermalGridIndex : ThermalGridIndexContainer)
	{
		double adjointValue(0);
		for(unsigned int row = 0 ; row < mCurrentSourceSize ; ++row) 
		{
			adjointValue += DenseAkMatrix[thermalGridIndex][row] * xx[row];
		}
		InitialTemperatureAdjointValue[thermalGridIndex] = adjointValue;	// for the first iteration
	}

	// Iterations, solutions from the first iteration until converge
	for(unsigned int iter = 1 ; iter <= MaxIteration ; ++iter)
	{
		Hj_trans_row.setZero();
		for(unsigned int i = 0 ; i < mCurrentSourceSize ; ++i)
		{
			cf[i] = 0;
		}

		for(auto& each_grid : DkRowsContainer)
		{
			double *curDkRow = each_grid.second;
			double TransAdjointValue(0);
			for(unsigned int row = 0 ; row < mCurrentSourceSize ; ++row) 
			{
				TransAdjointValue += curDkRow[row] * xx[row];
			}

			for(unsigned int row = 0 ; row < mCurrentSourceSize ; ++row) 
			{
				double tmp = TransAdjointValue * DenseAkMatrix[each_grid.first][row];
				Hj_trans_row[row] += tmp;
				cf[row] += InitialTemperatureAdjointValue[each_grid.first] * curDkRow[row] + tmp;
			}
		}

		for(unsigned int i = 0 ; i < mCurrentSourceSize ; ++i)
		{
			cf[i] += inv_row[i];
		}

		// Matlab
		memcpy((void *)mxGetPr(matf), (void *)cf, sizeof(double) * mCurrentSourceSize);
		engPutVariable(ep, "matf", matf);	
		engEvalString(ep, "tStart = cputime;");
		engEvalString(ep, "[xx fval] = linprog(matf,matAmatrix,matb,[],[],matlower,matupper);");
		engEvalString(ep, "lptime = cputime - tStart;");	
		LPTime += *(double *)mxGetPr(engGetVariable(ep, "lptime"));
		double curSolution = -1 * *(double *)mxGetPr(engGetVariable(ep, "fval"));
		if( fabs(InitialVdrop - curSolution) <= Tolerance )	// if objective solution error < tolerance, then converge
		{
			//std::cout << "Converge at iteration " << iter << std::endl;
			xx = (double *)mxGetPr(engGetVariable(ep, "xx"));  // get solution
			double VoltageShift(0);
			for(unsigned int i = 0 ; i < mCurrentSourceSize ; ++i)
			{
				VoltageShift += Hj_trans_row[i] * xx[i];
			}
			VdropSolution.emplace_back(curSolution + VoltageShift);
			ConvergeContainer.emplace_back(iter);	// record the iterations 
			break;
		}
		else {
			InitialVdrop = curSolution; // update objective solution
			xx = (double *)mxGetPr(engGetVariable(ep, "xx"));  // update solution current pattern for next iteration
			for(auto& thermalGridIndex : ThermalGridIndexContainer)
			{
				double adjointValue(0);
				for(unsigned int row = 0 ; row < mCurrentSourceSize ; ++row) 
				{
					adjointValue += DenseAkMatrix[thermalGridIndex][row] * xx[row];
				}
				InitialTemperatureAdjointValue[thermalGridIndex] = adjointValue;  // update temperature terms for next iteration
			}
		}

		if(iter == MaxIteration)
		{
			std::cout << "exceed max iteration" << std::endl;
		}

	}

	MatrixTime += (double)(clock()-matrixStart) / (CLOCKS_PER_SEC);
	// OnePointThermalAnalysis(specificIdx, xx, AkMatrix, 2);
	VdropIndex.emplace_back(each_node_idx);
}

}

	std::cout << "Complete Efficient Expand Point Verification !" << std::endl;

// Output solution 
	
	unsigned int violations(0);
	
	double sum(0);
	for(auto& vdrop : VdropSolution) { sum +=vdrop; }

	for(unsigned int i = 0 ; i < VdropIndex.size() ; ++i)
	{
		for(unsigned int j = i + 1 ; j < VdropIndex.size() ; ++j)
		{
			if(VdropIndex[i] > VdropIndex[j])
			{
				unsigned int idx_tmp = VdropIndex[i];
				VdropIndex[i] = VdropIndex[j];
				VdropIndex[j] = idx_tmp;

				double sol_tmp = VdropSolution[i];
				VdropSolution[i] = VdropSolution[j];
				VdropSolution[j] = sol_tmp;
			}
		}
	}

	std::ofstream fout("../LPresult/" + mCktName + "_EfficientExpandLPVdropMap_Matlab");
	int sol_idx(0);
	for(unsigned int i = 0 ; i < mPGEquivalentNodeSize ; ++i)
	{
		if(mPGEquivalentNodes[i]->is_connect_I)
		{
			fout << mPGEquivalentNodes[i]->x_loc_idx << " " << mPGEquivalentNodes[i]->y_loc_idx 
				 << " " << std::setprecision(15) << VdropSolution[sol_idx] << std::endl;
		
			if(VdropSolution[sol_idx] > mSupplyVoltage * 0.1)
			{
				++violations;
			}
			++sol_idx;
		}
	}
	fout.close();

	fout.open("../LPresult/" + mCktName + "_EfficientExpandLP_TimeInfo_Matlab");
	fout << "Total verification time = " << mConstructDataTime + LPTime + MatrixTime << " (s) " << std::endl;
	fout << "Average voltage drop = " << std::setprecision(15) << sum / VdropSolution.size() << " (V) " << std::endl;
	fout << "Maximum voltage drop = " << std::setprecision(15) << *std::max_element(VdropSolution.begin(), VdropSolution.end()) << " (V) " << std::endl;
	fout << "Total violations = " << violations << std::endl;
	fout << "Matrix Time = " << MatrixTime << " (s) " << std::endl;
	fout << "Matlab Time  = " << LPTime << " (s) " << std::endl;
	fout.close();

	fout.open(mCktName + "_Converge_WithoutSim");
	for(auto& i : ConvergeContainer)
	{
		fout << i << std::endl;
	}
	fout.close();

	mxDestroyArray(matf);
	mxDestroyArray(matAmatrix);
	mxDestroyArray(matb);
	mxDestroyArray(matlower);
	mxDestroyArray(matupper);
	engClose(ep);
	std::cout << "Closed Matlab engine" << std::endl;

	for(auto& each_grid : DkRowsContainer)
	{
		delete each_grid.second;
	}

	delete [] cG0INV_Q;
	delete [] cAkMatrix;
	delete [] Amatrix;
	delete [] b;
	delete [] clower;
	delete [] cupper;
	delete [] cf;
	
}


} // end namespace


	