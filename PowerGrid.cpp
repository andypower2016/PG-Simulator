// PoweGrid.cpp 
// 1. Power grid DC simulation (thermal simulation)
// 2. Constructing power grid verification datas
#include "PowerGrid.hpp"

namespace PowerGridPlane {

PowerGrid::PowerGrid()
 {}


void PowerGrid::ConstructPowerGridPlane()
{

	mTotalNodes["0"]->node_voltage = 0;
	mTotalNodes["0"]->is_grounded = true;
	mTotalNodes["0"]->is_stamp    = false;
	
	for(auto& volt : mtotal_supply_Vsource)
	{
		volt.second->node_voltage = volt.first->volt_value;
		volt.second->is_stamp = false;	// no need to stamp the nodes connect to vdd
	}

	if(mtotal_supply_Vsource.size() != 0)	// record Vdd
	{
		auto volt = mtotal_supply_Vsource.begin();
		mSupplyVoltage = (*volt).first->volt_value;	
	}

	mPGEquivalentNodeSize = (mTotalNodes.size() - 1 - mtotal_supply_Vsource.size());
	mPGEquivalentNodes.resize(mPGEquivalentNodeSize);
	for(auto& PNode : mTotalNodes)
	{
		if(PNode.second->node_name != "0" && !PNode.second->is_connect_V)
		{	
			unsigned int idx((PNode.second->x_loc_idx - 1) * mY_loc_max + (PNode.second->y_loc_idx - 1));
			mPGEquivalentNodes[idx] = PNode.second;	// store the nodes that needs to stamp in mPGEquivalentNodes container

			if(PNode.second->is_connect_I)
			{
				mtotal_current += PNode.second->connected_Isource["pos"][0]->curr_value;	// roughly estimate the total current 
			}
		}	
	}

	std::vector<double> PeakCurrentValue;
	for(auto& i : mtotal_Isource)
	{
		PeakCurrentValue.emplace_back(i->curr_value);
	}

	std::cout << "Voltage Sources = " << mtotal_supply_Vsource.size() << std::endl;
	std::cout << "Current Sources = " << mtotal_Isource.size() << std::endl;
	std::cout << "Max Peak Current Value = " << *std::max_element(PeakCurrentValue.begin(), PeakCurrentValue.end()) << std::endl;
	std::cout << "Min Peak Current Value = " << *std::min_element(PeakCurrentValue.begin(), PeakCurrentValue.end()) << std::endl;

}


// For power grid DC simulation
// Construct Power Map and thermal grids
void PowerGrid::ConstructThermalPlane()
{
	m2DPowerMap = new double*[(int)mXCut];
	for(int i=0;i<(int)mXCut;++i) m2DPowerMap[i] = new double[(int)mYCut]();
    
	for(auto& curr : mtotal_Isource)
	{	
			unsigned int xidx(0), yidx(0);
			for(unsigned int i=1 ; i <= mXCut ; ++i)
			{
				if( (curr->x_loc_idx-mX_loc_min) * mwire_seg_Xlen <= (double)(i*mPowerPlaneXLength/mXCut) )
				{
					xidx = i-1;
					break;
				}
			}
			for(unsigned int i=1 ; i <= mYCut ; ++i)
			{
				if( (curr->y_loc_idx-mY_loc_min) * mwire_seg_Ylen <= (double)(i*mPowerPlaneYLength/mYCut) )
				{
					yidx = i-1;
					break;
				}
			}
			double power = mTotalNodes[curr->pos_node]->node_voltage*curr->curr_value;
			curr->power = power;
			m2DPowerMap[xidx][yidx] += power;
			mTotalPower += power;
	}	
	
	// Construct thermal grid
	unsigned int idx(0);
	mZlayer = mZ_metalCut + 1 + mZ_subCut;   // 1 is for active layer 
	
	for(double z = 0 ; z < mZlayer ; ++z)
	{
		std::vector<std::vector<ThermalGrid> > XYGrids;
		for(double y = 0 ; y < mYCut ; ++y)
		{	
			std::vector<ThermalGrid> XGrids;
			for(double x = 0 ; x < mXCut ; ++x)
			{
				if( z < mZ_subCut )		// Silicon
				{
					if(z == 0)	// at boundary
					{
						ThermalGrid Grid
						( mX_loc_min+(x*mPowerPlaneXLength/mXCut), 
						  mY_loc_min+(y*mPowerPlaneYLength/mYCut), 
					  	  0 + z * mZ_silicon_len/mZ_subCut, 
					  	  mX_loc_min+( (x+1)*mPowerPlaneXLength/mXCut ), 
					  	  mY_loc_min+( (y+1)*mPowerPlaneYLength/mYCut ), 
					  	  0 + (z+1) * mZ_silicon_len/mZ_subCut,
					  	  130, 8000);

						Grid.mGridIndex = idx;
						XGrids.emplace_back(Grid);	
					}
					else
					{
						ThermalGrid Grid
						( mX_loc_min+(x*mPowerPlaneXLength/mXCut), 
						  mY_loc_min+(y*mPowerPlaneYLength/mYCut), 
					  	  0 + z * mZ_silicon_len/mZ_subCut, 
					  	  mX_loc_min+( (x+1)*mPowerPlaneXLength/mXCut ), 
					  	  mY_loc_min+( (y+1)*mPowerPlaneYLength/mYCut ), 
					  	  0 + (z+1) * mZ_silicon_len/mZ_subCut,
					  	  130, 0);

						Grid.mGridIndex = idx;
						XGrids.emplace_back(Grid);	
					}
				}
				else if( z == mZ_subCut )	// cut thermal grid for active layer 
				{
					ThermalGrid Grid
					( mX_loc_min+(x*mPowerPlaneXLength/mXCut), 
					  mY_loc_min+(y*mPowerPlaneYLength/mYCut), 
					  mZ_silicon_len, 
					  mX_loc_min+( (x+1)*mPowerPlaneXLength/mXCut ), 
					  mY_loc_min+( (y+1)*mPowerPlaneYLength/mYCut ), 
					  mZ_silicon_len + mZ_powerWidth,
					  130, 0);

					Grid.mGridIndex = idx;  
					XGrids.emplace_back(Grid);	

				}
				else if( z > mZ_subCut )   //  Metal
				{
					if(z == mZlayer-1)
					{
						ThermalGrid Grid
						( mX_loc_min+(x*mPowerPlaneXLength/mXCut), 
					  	  mY_loc_min+(y*mPowerPlaneYLength/mYCut), 
					  	  mZ_silicon_len + mZ_powerWidth + (z - mZ_subCut - 1) * mZ_metal_len/mZ_metalCut, 
					  	  mX_loc_min+( (x+1)*mPowerPlaneXLength/mXCut ), 
					  	  mY_loc_min+( (y+1)*mPowerPlaneYLength/mYCut ), 
					  	  mZ_silicon_len + mZ_powerWidth + (z - mZ_subCut) * mZ_metal_len/mZ_metalCut,
					  	  101, 2000);
					  
						Grid.mGridIndex = idx;  
						XGrids.emplace_back(Grid);	
					}
					else
					{
						ThermalGrid Grid
						( mX_loc_min+(x*mPowerPlaneXLength/mXCut), 
					  	  mY_loc_min+(y*mPowerPlaneYLength/mYCut), 
					  	  mZ_silicon_len + mZ_powerWidth + (z - mZ_subCut - 1) * mZ_metal_len/mZ_metalCut, 
					  	  mX_loc_min+( (x+1)*mPowerPlaneXLength/mXCut ), 
					  	  mY_loc_min+( (y+1)*mPowerPlaneYLength/mYCut ), 
					  	  mZ_silicon_len + mZ_powerWidth + (z - mZ_subCut) * mZ_metal_len/mZ_metalCut,
					  	  101, 0);
					  
						Grid.mGridIndex = idx;  
						XGrids.emplace_back(Grid);	
					}
				}
				idx++;
			}
			XYGrids.emplace_back(XGrids);	
		}
		mTotalThermalGrids.emplace_back(XYGrids);
	}
	
	for(unsigned int y = 0; y < mYCut ; ++y)
	{
		for(unsigned int x = 0; x < mXCut ; ++x)
		{
			mTotalThermalGrids[mZ_subCut][y][x].mPower = m2DPowerMap[x][y];
		}
	}
#if TEST		
	std::cout << "Construct Thermal Plane Complete !" << std::endl;
#endif
}	  


// Ouputs
void PowerGrid::OutputVdrop()
{
	std::fstream fout( "../result/" + mCktName + ".drop", std::ios::out);
	for (auto& node : mVdropList)
	{
		fout << node.first << " " << node.second << std::endl;
	}
	fout.close();
}


void PowerGrid::Output_DC_Solution()
{
	std::fstream fout( "../result/" + mCktName + ".sol", std::ios::out);
	
	for (auto& node : mTotalNodes)
	{
		if(node.first != "0")
			fout << node.first << " " << node.second->node_voltage << std::endl;
		else
			fout << "G" << " " << node.second->node_voltage << std::endl;
	}
	std::cout << "Output solution complete !" << std::endl;
	fout.close();
}


void PowerGrid::Output_Element_Solution()
{

		std::fstream fout("../result/" + mCktName + "_" + "List", std::ios::out);

		fout << "#--------------------------- Voltage Source ---------------------------" << std::endl;
		for (auto& volt_source : mtotal_supply_Vsource)
		{
			fout << volt_source.first->volt_name 
			<< std::setw(20) 
			<< "voltage value = " << volt_source.first->volt_value << "V"
			<< std::setw(20) 
			<< "current = " << volt_source.first->curr_value << "A"
			<< std::setw(20) 
			<< "power = " << volt_source.first->power_value << "W" 
			<< std::endl;
		}
		fout << "#----------------------------------------------------------------------" << std::endl << std::endl;

		fout << "#--------------------------- Total Resistors ---------------------------" << std::endl;
		for (auto& each_res : mtotal_unshorted_Res)
		{
			fout << each_res->res_name 
			<< std::setw(20) 
			<< "resistor value = " << each_res->res_value << "ohm" 
			<< std::setw(20) 
			<< "Vdrop = " << each_res->Vdrop << "V" 
			<< std::setw(20) 
			<< "current = " << each_res->curr_value << "A" 
			<< std::setw(20) 
			<< "power = " << each_res->power_value << "W" << std::endl;
		}
		fout << "#----------------------------------------------------------------------" << std::endl;
#if TEST
		std::cout << "Output element list file complete !" << std::endl;
#endif
}

void PowerGrid::OutputThermalSolution()
{
	std::fstream fout;
	/*
	fout.open("../result/STMap", std::ios::out);

	for(unsigned int y = 0; y < mYCut ; ++y)
		for(unsigned int x = 0; x < mXCut ; ++x)
		{
			//if(mThermalSolution[mTotalThermalGrids[0][y][x].mGridIndex] != 0)
			fout << x
		         << " " 
			     << y
			     << " " 
			     << 23 + mThermalSolution[mTotalThermalGrids[0][y][x].mGridIndex] << std::endl;	
		}

	fout.close();
	*/

	fout.open("../result/" + mCktName + "_MTMap", std::ios::out);

	for(unsigned int y = 0; y < mYCut ; ++y)
		for(unsigned int x = 0; x < mXCut ; ++x)
		{
			//if(mThermalSolution[mTotalThermalGrids[1][y][x].mGridIndex] != 0)
			fout << x
		         << " " 
			     << y
			     << " " 
			     << 23 + mThermalSolution[mTotalThermalGrids[mZlayer-1][y][x].mGridIndex] << std::endl;	
		}

	fout.close();
	

	// for icepack
	/*
	fout.open("../result/" + mCktName + "_" + "Temperature", std::ios::out);

	for(int z = 0 ; z < mZlayer ; ++z)
		for(int y = 0; y < mYCut ; ++y)
			for(int x = 0; x < mXCut ; ++x)
			{
				fout << std::setprecision(10) << mTotalThermalGrids[z][y][x].getXmidpoint() << " " 
				     << std::setprecision(10) << mTotalThermalGrids[z][y][x].getYmidpoint() << " " 
				     << std::setprecision(10) << mTotalThermalGrids[z][y][x].getZmidpoint() << " "
				     << std::setprecision(10) << 23 + mThermalSolution[mTotalThermalGrids[z][y][x].mGridIndex] << std::endl;
			}

	fout.close();
	*/
}

void PowerGrid::OutuptDistributionMaps()
{
	std::fstream fout;

	fout.open("../result/" + mCktName + "_" + "VoltageMap", std::ios::out);
	for(auto& node : mPGEquivalentNodes)
	{
		fout << node->x_loc_idx << " " 
			 << node->y_loc_idx << " " 
		     << node->node_voltage << std::endl;
	}
	fout.close();
	

	/*
	fout.open("../result/" file_name + "_" + "CurrentMap", std::ios::out);
	for(auto& curr : mtotal_Isource)
	{	
		fout << curr->x_loc_idx << " " 
			 << curr->y_loc_idx << " " 
			 << curr->curr_value << std::endl;	
	}
	
	for(auto& node : mTotal_Equivalent_PG_Nodes)
	{
			fout << node.second->x_loc_idx << " " 
			 	 << node.second->y_loc_idx << " " 
			 	 << 0 << std::endl;	
	}
	fout.close();
	*/

	
	fout.open("../result/" + mCktName + "_" + "PowerMap", std::ios::out);
	
	for(unsigned int y = 0; y < mYCut ; ++y)
		for(unsigned int x = 0; x < mXCut ; ++x)
		{
			fout << x
		         << " " 
			     << y
			     << " " 
			     << m2DPowerMap[x][y]<< std::endl;
		}
	fout.close();
	

}


void PowerGrid::OutputDesignInfo()
{
	std::fstream fout;
	
	fout.open("../PGGen_result/" + mCktName + "_" + "DesignInfo", std::ios::out);

	fout << "#----------------------------------------------------------------------------------------------#" << std::endl << std::endl;
	fout << "# PG geometry : " << std::endl << std::endl;
	fout << "X length           = " << mPowerPlaneXLength << "m" << std::endl;
	fout << "Y length           = " << mPowerPlaneYLength << "m" << std::endl;
	//fout << "Z metal length     = " << mZ_metal_len  << "m" << std::endl;
	//fout << "Z substrate length = " << mZ_silicon_len  << "m" << std::endl;
	fout << std::endl;
	fout << "#----------------------------------------------------------------------------------------------#" << std::endl << std::endl;

	for(auto& node : mPGEquivalentNodes)
	{
		if(node->node_voltage < mSupplyVoltage * 0.9)
		{
			++mViolations;
		}
	}

	std::cout << "total violations = " << mViolations << std::endl;
	
	fout << "# PG specs : " << std::endl << std::endl;
	fout << "Total Nodes             = " << mTotalNodes.size() << std::endl;
	fout << "Equivalent Nodes        = " << mPGEquivalentNodes.size() << std::endl;
	fout << "Total Voltage Source    = " << mtotal_supply_Vsource.size() << std::endl;
	fout << "Supply Voltage          = " << mSupplyVoltage << " V" << std::endl;
	fout << "Total Current Source    = " << mtotal_Isource.size() << std::endl;
	fout << "Total Current           = " << mtotal_current << " A" << std::endl;
	mTotalPower = 0;
	for(auto& i : mtotal_Isource) { mTotalPower += mSupplyVoltage * i->curr_value; }
	fout << "Total Power Consumption = " << mTotalPower << " W" <<  std::endl;
	//fout << "Total IR drop violation = " << mViolations << std::endl;
	fout << std::endl;
	fout << "#----------------------------------------------------------------------------------------------#" << std::endl;
	
	fout.close();

	std::cout << "Output Design Info Complete !" << std::endl;
}


void PowerGrid::OutputThermalGrids()
{
	std::ofstream fout("../result/" + mCktName + "_" + "ThermalGrids", std::ios::out);
	for(int z = 0 ; z < mZlayer ; ++z)
		for(int y = 0; y < mYCut ; ++y)
			for(int x = 0; x < mXCut ; ++x)
			{
				fout << "(" << mTotalThermalGrids[z][y][x].dx_loc_idx << ", " 
				            << mTotalThermalGrids[z][y][x].dy_loc_idx << ", " 
				            << mTotalThermalGrids[z][y][x].dz_loc_idx << ")" << " "
				     << "(" << mTotalThermalGrids[z][y][x].ux_loc_idx << ", " 
				            << mTotalThermalGrids[z][y][x].uy_loc_idx << ", " 
				            << mTotalThermalGrids[z][y][x].uz_loc_idx << ")" << " "				           
				            << "power = " << mTotalThermalGrids[z][y][x].mPower << " "
				            << "T = " << 27 + mThermalSolution[mTotalThermalGrids[z][y][x].mGridIndex] << " "
				            << "Mid point = " 
				            << " (" << (mTotalThermalGrids[z][y][x].dx_loc_idx + mTotalThermalGrids[z][y][x].ux_loc_idx)/2  << ", " 
				            << (mTotalThermalGrids[z][y][x].dy_loc_idx + mTotalThermalGrids[z][y][x].uy_loc_idx)/2 << ", " 
				            << (mTotalThermalGrids[z][y][x].dz_loc_idx + mTotalThermalGrids[z][y][x].uz_loc_idx)/2 << ")"
				            << std::endl;
			

			}

    fout.close();
}


// for icepack 
void PowerGrid::OutputPointReport()
{
	std::ofstream fout("../result/" + mCktName + "_" 
						+ std::to_string((int)mXCut) + "_" 
						+ std::to_string((int)mYCut) + "_" 
						+ std::to_string((int)mZlayer) 
						+ "_pointreport.rpt", std::ios::out);

	fout << "Report point data"     << std::endl;
	fout << mZlayer * mYCut * mXCut << std::endl;
	for(int z = 0 ; z < mZlayer ; ++z)
		for(int y = 0; y < mYCut ; ++y)
			for(int x = 0; x < mXCut ; ++x)
			{
				fout << "point {" << std::setprecision(10) << mTotalThermalGrids[z][y][x].getXmidpoint() << " " 
				                  << std::setprecision(10) << mTotalThermalGrids[z][y][x].getYmidpoint() << " "
				                  << std::setprecision(10) << mTotalThermalGrids[z][y][x].getZmidpoint() 
				                  << " } value Temperature active 1" << std::endl;
			}
	fout.close();
}



// 2017.7.19 ~ 
// For verification 

// 2017.8.1 ~
// VERIFICATION_MODE
// 0 LP
// 1 QP
// 2 Expand Point

// METHOD_MODE
// 0 Origin 
// 1 Efficient

void PowerGrid::ConstructCurrentConstraint(int VERIFICATION_MODE, int METHOD_MODE)
{

// Parse constraint profile 

	std::ifstream fin("../input/constraint.txt");
	std::string name, x, y;
	
	if(fin.is_open())
	{
		fin >> name >> x >> y ;
		mXBlockDim = stod(x);
		mYBlockDim = stod(y);
		fin >> name >> x ;
		mGlobalCCSF = stod(x);
		fin.close();
	}
	else {
		std::cerr << "Open constraint.txt failed !" << std::endl;
		fin.close();
	}

/*
#if TEST
	std::cout << "Block x dim = " << mXBlockDim  << std::endl;
	std::cout << "Block y dim = " << mYBlockDim  << std::endl;
	std::cout << "Scale factor= " << mGlobalCCSF << std::endl;
#endif
*/


// Set local current constraint
if(METHOD_MODE == 0)		// Origin method
{
	mLocalCurrentConstraint.resize(mPGEquivalentNodeSize, 0);
	for(auto& PNode : mPGEquivalentNodes)
	{
		if(PNode->is_connect_I)
		{
			unsigned int idx((PNode->x_loc_idx - 1) * mY_loc_max + (PNode->y_loc_idx - 1));
			mLocalCurrentConstraint[idx] = PNode->connected_Isource["pos"][0]->curr_value;
		}
	}
}
else if (METHOD_MODE == 1)	// Efficient Method
{
	mCurrentSourceSize = mtotal_Isource.size();
	mLocalCurrentConstraint.resize(mCurrentSourceSize, 0);
	std::vector<double> idx_source(mPGEquivalentNodeSize);
	for(auto& PNode : mPGEquivalentNodes)
	{
		if(PNode->is_connect_I)
		{
			unsigned int idx((PNode->x_loc_idx - 1) * mY_loc_max + (PNode->y_loc_idx - 1));
			idx_source[idx] = PNode->connected_Isource["pos"][0]->curr_value;
			mCurrentSourceIdx.emplace_back(idx);	
		}
	}
	std::sort(mCurrentSourceIdx.begin(), mCurrentSourceIdx.end());

	for(unsigned int node_idx = 0 ; node_idx < mCurrentSourceSize ; ++node_idx)
	{
		mLocalCurrentConstraint[node_idx] = idx_source[mCurrentSourceIdx[node_idx]];
		mCurrentSourceIdxMapping[mCurrentSourceIdx[node_idx]] = node_idx;
	}
}

/*
#if TEST
	std::cout << "Construct Local Current Constraint Complete ! " << std::endl;
#endif	
*/

if((VERIFICATION_MODE == 0 || VERIFICATION_MODE == 2)  && METHOD_MODE == 0)	// Origin LP & Expand Point
{
// Partition Power Grid and set global constraint

	mGlobalConstraintSize = mXBlockDim * mYBlockDim;
	EigenSparseMatrix Mat_GlobalCurrentConstraintIncidenceMatrix(mGlobalConstraintSize, mPGEquivalentNodeSize);	
	std::vector<triplet> triplet_GlobalCurrentConstraint;
	mGlobalCurrentConstraint.resize(mGlobalConstraintSize, 0);

	double Xidxlen(mX_loc_max-mX_loc_min);
	double Yidxlen(mY_loc_max-mY_loc_min);
	double Xlength(Xidxlen/mXBlockDim);
	double Ylength(Yidxlen/mYBlockDim);
	unsigned int blockCount(0);
	std::vector<bool> visit(mPGEquivalentNodeSize, false);
	for(unsigned int i = 1 ; i <= (unsigned int)mXBlockDim ; ++i)
	{
		for(unsigned int j = 1 ; j <= (unsigned int)mYBlockDim ; ++j)
		{
			double blockCurrentSum(0);
			for(auto& node : mPGEquivalentNodes)
			{

				unsigned int idx((node->x_loc_idx - 1) * mY_loc_max + (node->y_loc_idx - 1));

				if( node->x_loc_idx <= 1 + i * Xlength && node->x_loc_idx >= 1 + (i-1) * Xlength &&
					node->y_loc_idx <= 1 + j * Ylength && node->y_loc_idx >= 1 + (j-1) * Ylength &&
					node->is_connect_I && 
					!visit[idx] )
				{
					visit[idx] = true;
					blockCurrentSum += node->connected_Isource["pos"][0]->curr_value;
			
					triplet_GlobalCurrentConstraint.emplace_back(triplet(blockCount, idx, 1));
				}
			}
			
			mGlobalCurrentConstraint[blockCount] = blockCurrentSum * mGlobalCCSF;
		
			++blockCount;
		}
	}
	Mat_GlobalCurrentConstraintIncidenceMatrix.setFromTriplets(triplet_GlobalCurrentConstraint.begin(), triplet_GlobalCurrentConstraint.end());
	mGlobalCurrentConstraintIncidenceMatrix = Mat_GlobalCurrentConstraintIncidenceMatrix;

}
else if((VERIFICATION_MODE == 0 || VERIFICATION_MODE == 2)  && METHOD_MODE == 1) // Efficient LP & Expand Point
{
	std::vector<int> SourceCount;



	mGlobalConstraintSize = mXBlockDim * mYBlockDim;
	EigenSparseMatrix Mat_GlobalCurrentConstraintIncidenceMatrix(mGlobalConstraintSize, mCurrentSourceSize);	
	std::vector<triplet> triplet_GlobalCurrentConstraint;
	mGlobalCurrentConstraint.resize(mGlobalConstraintSize, 0);

	double Xidxlen(mX_loc_max-mX_loc_min);
	double Yidxlen(mY_loc_max-mY_loc_min);
	double Xlength(Xidxlen/mXBlockDim);
	double Ylength(Yidxlen/mYBlockDim);
	unsigned int blockCount(0);
	std::vector<bool> visit(mPGEquivalentNodeSize, false);
	for(unsigned int i = 1 ; i <= (unsigned int)mXBlockDim ; ++i)
	{
		for(unsigned int j = 1 ; j <= (unsigned int)mYBlockDim ; ++j)
		{
			int source_cnt(0);
			double blockCurrentSum(0);
			for(auto& node : mPGEquivalentNodes)
			{

				unsigned int idx((node->x_loc_idx - 1) * mY_loc_max + (node->y_loc_idx - 1));

				if( node->x_loc_idx <= 1 + i * Xlength && node->x_loc_idx >= 1 + (i-1) * Xlength &&
					node->y_loc_idx <= 1 + j * Ylength && node->y_loc_idx >= 1 + (j-1) * Ylength &&
					node->is_connect_I && 
					!visit[idx] )
				{
					visit[idx] = true;
					blockCurrentSum += node->connected_Isource["pos"][0]->curr_value;
				
					triplet_GlobalCurrentConstraint.emplace_back(triplet(blockCount, mCurrentSourceIdxMapping[idx], 1));
				
					++source_cnt;
				}
			}
			mGlobalCurrentConstraint[blockCount] = blockCurrentSum * mGlobalCCSF;
			++blockCount;

			SourceCount.emplace_back(source_cnt);
		}
	}
	Mat_GlobalCurrentConstraintIncidenceMatrix.setFromTriplets(triplet_GlobalCurrentConstraint.begin(), triplet_GlobalCurrentConstraint.end());
	mGlobalCurrentConstraintIncidenceMatrix = Mat_GlobalCurrentConstraintIncidenceMatrix;

/*
	for(auto& src : SourceCount)
		std::cout <<src << std::endl;
*/
	/*
	std::ofstream fout(mCktName + "_ConstraintBLocks");
	for(unsigned int i = 0 ; i < 4 ; ++i)
	{
		fout << SourceCount[i] << " " << mGlobalCurrentConstraint[i] << std::endl;
	}
	fout.close();
	*/

}
else if( VERIFICATION_MODE == 1  && METHOD_MODE == 0 )				// Origin QP
{

	mGlobalConstraintSize = mXBlockDim*mYBlockDim;
	
	EigenSparseMatrix MatlabAMatrix(2 * mGlobalConstraintSize, mPGEquivalentNodeSize);	
	std::vector<triplet> triplet_MatlabAmatrix;
	mGlobalCurrentConstraint.resize(mGlobalConstraintSize, 0);

	double Xidxlen(mX_loc_max-mX_loc_min);
	double Yidxlen(mY_loc_max-mY_loc_min);
	double Xlength(Xidxlen/mXBlockDim);
	double Ylength(Yidxlen/mYBlockDim);
	unsigned int blockCount(0);
	std::vector<bool> visit(mPGEquivalentNodeSize, false);
	for(unsigned int i = 1 ; i <= (unsigned int)mXBlockDim ; ++i)
	{
		for(unsigned int j = 1 ; j <= (unsigned int)mYBlockDim ; ++j)
		{
			double blockCurrentSum(0);
			for(auto& node : mPGEquivalentNodes)
			{
				unsigned int idx((node->x_loc_idx - 1) * mY_loc_max + (node->y_loc_idx - 1));

				if( node->x_loc_idx <= 1 + i * Xlength && node->x_loc_idx >= 1 + (i-1) * Xlength &&
					node->y_loc_idx <= 1 + j * Ylength && node->y_loc_idx >= 1 + (j-1) * Ylength &&
					node->is_connect_I && 
					!visit[idx] )
				{

					visit[idx] = true;
					blockCurrentSum += node->connected_Isource["pos"][0]->curr_value;
			
					// for matlab 
					triplet_MatlabAmatrix.emplace_back(triplet(blockCount, idx, 1));
					triplet_MatlabAmatrix.emplace_back(triplet(blockCount + mGlobalConstraintSize, idx, -1));
				}
			}
			mGlobalCurrentConstraint[blockCount] = blockCurrentSum * mGlobalCCSF;
			++blockCount;
		}
	}

	MatlabAMatrix.setFromTriplets(triplet_MatlabAmatrix.begin(), triplet_MatlabAmatrix.end());
	mMatlabAmatrix = MatlabAMatrix;	


}
else if( VERIFICATION_MODE == 1  && METHOD_MODE == 1 )				// Efficient QP
{


	mGlobalConstraintSize = mXBlockDim * mYBlockDim;
	EigenSparseMatrix MatlabAMatrix(2 * mGlobalConstraintSize, mCurrentSourceSize);	
	std::vector<triplet> triplet_MatlabAmatrix;
	mGlobalCurrentConstraint.resize(mGlobalConstraintSize, 0);

	double Xidxlen(mX_loc_max-mX_loc_min);
	double Yidxlen(mY_loc_max-mY_loc_min);
	double Xlength(Xidxlen/mXBlockDim);
	double Ylength(Yidxlen/mYBlockDim);
	unsigned int blockCount(0);
	std::vector<bool> visit(mPGEquivalentNodeSize, false);

	
	for(unsigned int i = 1 ; i <= (unsigned int)mXBlockDim ; ++i)
	{
		for(unsigned int j = 1 ; j <= (unsigned int)mYBlockDim ; ++j)
		{
			double blockCurrentSum(0);
			for(auto& node : mPGEquivalentNodes)
			{

				unsigned int idx((node->x_loc_idx - 1) * mY_loc_max + (node->y_loc_idx - 1));

				if( node->x_loc_idx <= 1 + i * Xlength && node->x_loc_idx >= 1 + (i-1) * Xlength &&
					node->y_loc_idx <= 1 + j * Ylength && node->y_loc_idx >= 1 + (j-1) * Ylength &&
					node->is_connect_I && 
					!visit[idx] )
				{
					
					visit[idx] = true;
					blockCurrentSum += node->connected_Isource["pos"][0]->curr_value;
					
					// for matlab 
					triplet_MatlabAmatrix.emplace_back(triplet(blockCount, mCurrentSourceIdxMapping[idx], 1));
					triplet_MatlabAmatrix.emplace_back(triplet(blockCount + mGlobalConstraintSize, mCurrentSourceIdxMapping[idx], -1));
				}
			}

			mGlobalCurrentConstraint[blockCount] = blockCurrentSum * mGlobalCCSF;
		
			++blockCount;
		}
	}
	
	MatlabAMatrix.setFromTriplets(triplet_MatlabAmatrix.begin(), triplet_MatlabAmatrix.end());
	mMatlabAmatrix = MatlabAMatrix;

}

/*
#if TEST
	std::ofstream fout("eigen_incidence.txt");
	fout << mGlobalCurrentConstraintIncidenceMatrix;
	fout.close();
	fout.open("local.txt");
	int idx(0);
	for(auto& node : mPGEquivalentNodes)
	{
		fout << (node->x_loc_idx - 1) * mY_loc_max + (node->y_loc_idx - 1) << " " << mLocalCurrentConstraint[(node->x_loc_idx - 1) * mY_loc_max + (node->y_loc_idx - 1)] << std::endl;
		idx++;
	}
	fout.close();
#endif	
*/

}


// 2017.7.24 ~ 
void PowerGrid::ConstructThermalGrids()
{

	// Construct thermal grid
	unsigned int idx(0);
	mZlayer = mZ_metalCut + 1 + mZ_subCut;   // 1 is for active layer 
	
	for(double z = 0 ; z < mZlayer ; ++z)
	{
		std::vector<std::vector<ThermalGrid> > XYGrids;
		for(double y = 0 ; y < mYCut ; ++y)
		{	
			std::vector<ThermalGrid> XGrids;
			for(double x = 0 ; x < mXCut ; ++x)
			{
				if( z < mZ_subCut )		// Silicon
				{
					if(z == 0)	// at boundary, primary heat flow path
					{
						ThermalGrid Grid
						( mX_loc_min+(x*mPowerPlaneXLength/mXCut), 
						  mY_loc_min+(y*mPowerPlaneYLength/mYCut), 
					  	  0 + z * mZ_silicon_len/mZ_subCut, 
					  	  mX_loc_min+( (x+1)*mPowerPlaneXLength/mXCut ), 
					  	  mY_loc_min+( (y+1)*mPowerPlaneYLength/mYCut ), 
					  	  0 + (z+1) * mZ_silicon_len/mZ_subCut,
					  	  130, 8000);

						Grid.mGridIndex = idx;
						XGrids.emplace_back(Grid);	
					}
					else  // not at boundary 
					{
						ThermalGrid Grid
						( mX_loc_min+(x*mPowerPlaneXLength/mXCut), 
						  mY_loc_min+(y*mPowerPlaneYLength/mYCut), 
					  	  0 + z * mZ_silicon_len/mZ_subCut, 
					  	  mX_loc_min+( (x+1)*mPowerPlaneXLength/mXCut ), 
					  	  mY_loc_min+( (y+1)*mPowerPlaneYLength/mYCut ), 
					  	  0 + (z+1) * mZ_silicon_len/mZ_subCut,
					  	  130, 0);

						Grid.mGridIndex = idx;
						XGrids.emplace_back(Grid);	
					}
				}
				else if( z == mZ_subCut )	// cut thermal grid for active layer 
				{
					ThermalGrid Grid
					( mX_loc_min+(x*mPowerPlaneXLength/mXCut), 
					  mY_loc_min+(y*mPowerPlaneYLength/mYCut), 
					  mZ_silicon_len, 
					  mX_loc_min+( (x+1)*mPowerPlaneXLength/mXCut ), 
					  mY_loc_min+( (y+1)*mPowerPlaneYLength/mYCut ), 
					  mZ_silicon_len + mZ_powerWidth,
					  130, 0);

					Grid.mGridIndex = idx;  
					XGrids.emplace_back(Grid);	

				}
				else if( z > mZ_subCut )   //  Metal
				{
					if(z == mZlayer-1)	// at boundary, secondary heat flow path
					{
						ThermalGrid Grid
						( mX_loc_min+(x*mPowerPlaneXLength/mXCut), 
					  	  mY_loc_min+(y*mPowerPlaneYLength/mYCut), 
					  	  mZ_silicon_len + mZ_powerWidth + (z - mZ_subCut - 1) * mZ_metal_len/mZ_metalCut, 
					  	  mX_loc_min+( (x+1)*mPowerPlaneXLength/mXCut ), 
					  	  mY_loc_min+( (y+1)*mPowerPlaneYLength/mYCut ), 
					  	  mZ_silicon_len + mZ_powerWidth + (z - mZ_subCut) * mZ_metal_len/mZ_metalCut,
					  	  101, 2000);
					  
						Grid.mGridIndex = idx;  
						XGrids.emplace_back(Grid);	
					}
					else  // not at boundary 
					{
						ThermalGrid Grid
						( mX_loc_min+(x*mPowerPlaneXLength/mXCut), 
					  	  mY_loc_min+(y*mPowerPlaneYLength/mYCut), 
					  	  mZ_silicon_len + mZ_powerWidth + (z - mZ_subCut - 1) * mZ_metal_len/mZ_metalCut, 
					  	  mX_loc_min+( (x+1)*mPowerPlaneXLength/mXCut ), 
					  	  mY_loc_min+( (y+1)*mPowerPlaneYLength/mYCut ), 
					  	  mZ_silicon_len + mZ_powerWidth + (z - mZ_subCut) * mZ_metal_len/mZ_metalCut,
					  	  101, 0);
					  
						Grid.mGridIndex = idx;  
						XGrids.emplace_back(Grid);	
					}
				}
				idx++;
			}
			XYGrids.emplace_back(XGrids);	
		}
		mTotalThermalGrids.emplace_back(XYGrids);
	}
/*
#if TEST
	std::cout << "Construct Thermal Grids Complete !" << std::endl;
#endif	
*/

}



// 2017.7.24 ~ 
void PowerGrid::StampPGConductanceMatrix()
{

// Stamp origin conductance matrix G0

	mPG_OriginConductanceMatrix.resize(mPGEquivalentNodeSize);

	
	unsigned int mapping_idx(0);
	for (auto& node : mPGEquivalentNodes)
	{
		mStampNodeList.insert(std::pair<Node*, unsigned int>(node, mapping_idx));
		mapping_idx++;
	}

	for (auto& stamp_node : mStampNodeList)
	{
		for (auto& nei_node : stamp_node.first->neighbor_NodeRes)    
		{
			if (nei_node.first->is_grounded == true)
			{
				for (auto& each_res : nei_node.second)
				{
					mPG_OriginConductanceMatrix[stamp_node.second][stamp_node.second] += 1 / each_res->res_value;
				}
			}
			else if (nei_node.first->is_connect_V == true)
			{
				for (auto& each_res : nei_node.second)
				{				
					mPG_OriginConductanceMatrix[stamp_node.second][stamp_node.second] += 1 / each_res->res_value;
				}
			}
			else if (nei_node.first->is_connect_V == false && nei_node.first->is_grounded == false)
			{
				for (auto& each_res : nei_node.second)
				{
					mPG_OriginConductanceMatrix[stamp_node.second][stamp_node.second] += 1 / each_res->res_value;
					mPG_OriginConductanceMatrix[stamp_node.second][mStampNodeList[nei_node.first]] += -1 / each_res->res_value;
				}
			}
		}
	}
	
/*
#if TEST
	std::cout << "Stamping origin conductance matrix complete !" << std::endl;
#endif
*/
}


// for QP and Expand point method 1
void PowerGrid::StampThermalAwareMatrix()
{

// Stamp thermal matrix
	
	mPG_ThermalMatrix.resize(mXCut*mYCut*mZlayer);

	for(int z = 0 ; z < mZlayer ; ++z)
		for(int y = 0; y < mYCut ; ++y)
			for(int x = 0; x < mXCut ; ++x)
	{
				
				if(x-1 < 0)   // at x plane min bound
				{
					mPG_ThermalMatrix[mTotalThermalGrids[z][y][x].mGridIndex][mTotalThermalGrids[z][y][x].mGridIndex]
					+= 1/(mTotalThermalGrids[z][y][x].x_res + mTotalThermalGrids[z][y][x+1].x_res);
					
					mPG_ThermalMatrix[mTotalThermalGrids[z][y][x].mGridIndex][mTotalThermalGrids[z][y][x+1].mGridIndex]
					+= -1/(mTotalThermalGrids[z][y][x].x_res + mTotalThermalGrids[z][y][x+1].x_res);
				}
				else if(x-1 >= 0 && x+1 < mXCut) // at x plane between
				{
					mPG_ThermalMatrix[mTotalThermalGrids[z][y][x].mGridIndex][mTotalThermalGrids[z][y][x].mGridIndex]
					+= 1/(mTotalThermalGrids[z][y][x].x_res + mTotalThermalGrids[z][y][x-1].x_res)+
					   1/(mTotalThermalGrids[z][y][x].x_res + mTotalThermalGrids[z][y][x+1].x_res);
					
					mPG_ThermalMatrix[mTotalThermalGrids[z][y][x].mGridIndex][mTotalThermalGrids[z][y][x+1].mGridIndex]
					+= -1/(mTotalThermalGrids[z][y][x].x_res + mTotalThermalGrids[z][y][x+1].x_res);
					
					mPG_ThermalMatrix[mTotalThermalGrids[z][y][x].mGridIndex][mTotalThermalGrids[z][y][x-1].mGridIndex]
					+= -1/(mTotalThermalGrids[z][y][x].x_res + mTotalThermalGrids[z][y][x-1].x_res);
				}
				else if(x+1 >= mXCut) // at x plane max bound
				{
					mPG_ThermalMatrix[mTotalThermalGrids[z][y][x].mGridIndex][mTotalThermalGrids[z][y][x].mGridIndex]
					+= 1/(mTotalThermalGrids[z][y][x].x_res + mTotalThermalGrids[z][y][x-1].x_res);
					
					mPG_ThermalMatrix[mTotalThermalGrids[z][y][x].mGridIndex][mTotalThermalGrids[z][y][x-1].mGridIndex]
					+= -1/(mTotalThermalGrids[z][y][x].x_res + mTotalThermalGrids[z][y][x-1].x_res);	
				}


				
				if(y-1 < 0)   
				{
					mPG_ThermalMatrix[mTotalThermalGrids[z][y][x].mGridIndex][mTotalThermalGrids[z][y][x].mGridIndex]
					+= 1/(mTotalThermalGrids[z][y][x].y_res + mTotalThermalGrids[z][y+1][x].y_res);
					
					mPG_ThermalMatrix[mTotalThermalGrids[z][y][x].mGridIndex][mTotalThermalGrids[z][y+1][x].mGridIndex]
					+= -1/(mTotalThermalGrids[z][y][x].y_res + mTotalThermalGrids[z][y+1][x].y_res);
				}
				else if(y-1 >= 0 && y+1 < mYCut) 
				{
					mPG_ThermalMatrix[mTotalThermalGrids[z][y][x].mGridIndex][mTotalThermalGrids[z][y][x].mGridIndex]
					+= 1/(mTotalThermalGrids[z][y][x].y_res + mTotalThermalGrids[z][y-1][x].y_res)+
					   1/(mTotalThermalGrids[z][y][x].y_res + mTotalThermalGrids[z][y+1][x].y_res);
					
					mPG_ThermalMatrix[mTotalThermalGrids[z][y][x].mGridIndex][mTotalThermalGrids[z][y+1][x].mGridIndex]
					+= -1/(mTotalThermalGrids[z][y][x].y_res + mTotalThermalGrids[z][y+1][x].y_res);
					
					mPG_ThermalMatrix[mTotalThermalGrids[z][y][x].mGridIndex][mTotalThermalGrids[z][y-1][x].mGridIndex]
					+= -1/(mTotalThermalGrids[z][y][x].y_res + mTotalThermalGrids[z][y-1][x].y_res);	
				}
				else if(y+1 >= mYCut) 
				{
					mPG_ThermalMatrix[mTotalThermalGrids[z][y][x].mGridIndex][mTotalThermalGrids[z][y][x].mGridIndex]
					+= 1/(mTotalThermalGrids[z][y][x].y_res + mTotalThermalGrids[z][y-1][x].y_res);
					
					mPG_ThermalMatrix[mTotalThermalGrids[z][y][x].mGridIndex][mTotalThermalGrids[z][y-1][x].mGridIndex]
					+= -1/(mTotalThermalGrids[z][y][x].y_res + mTotalThermalGrids[z][y-1][x].y_res);		
				}
					
				if(z-1 < 0)   
				{
					mPG_ThermalMatrix[mTotalThermalGrids[z][y][x].mGridIndex][mTotalThermalGrids[z][y][x].mGridIndex]
					+= 1/(mTotalThermalGrids[z][y][x].z_res + mTotalThermalGrids[z][y][x].z_bound_res)+
					   1/(mTotalThermalGrids[z][y][x].z_res + mTotalThermalGrids[z+1][y][x].z_res);
					
					mPG_ThermalMatrix[mTotalThermalGrids[z][y][x].mGridIndex][mTotalThermalGrids[z+1][y][x].mGridIndex]
					+= -1/(mTotalThermalGrids[z][y][x].z_res + mTotalThermalGrids[z+1][y][x].z_res);
				}
				else if(z-1 >= 0 && z+1 < mZlayer) 
				{
					mPG_ThermalMatrix[mTotalThermalGrids[z][y][x].mGridIndex][mTotalThermalGrids[z][y][x].mGridIndex]
					+= 1/(mTotalThermalGrids[z][y][x].z_res + mTotalThermalGrids[z-1][y][x].z_res)+
					   1/(mTotalThermalGrids[z][y][x].z_res + mTotalThermalGrids[z+1][y][x].z_res);
					
					mPG_ThermalMatrix[mTotalThermalGrids[z][y][x].mGridIndex][mTotalThermalGrids[z+1][y][x].mGridIndex]
					+= -1/(mTotalThermalGrids[z][y][x].z_res + mTotalThermalGrids[z+1][y][x].z_res);
					
					mPG_ThermalMatrix[mTotalThermalGrids[z][y][x].mGridIndex][mTotalThermalGrids[z-1][y][x].mGridIndex]
					+= -1/(mTotalThermalGrids[z][y][x].z_res + mTotalThermalGrids[z-1][y][x].z_res);	
				}
				else if(z+1 >= mZlayer) 
				{
					mPG_ThermalMatrix[mTotalThermalGrids[z][y][x].mGridIndex][mTotalThermalGrids[z][y][x].mGridIndex]
					+= 1/(mTotalThermalGrids[z][y][x].z_res + mTotalThermalGrids[z][y][x].z_bound_res)+
					   1/(mTotalThermalGrids[z][y][x].z_res + mTotalThermalGrids[z-1][y][x].z_res);
					
					mPG_ThermalMatrix[mTotalThermalGrids[z][y][x].mGridIndex][mTotalThermalGrids[z-1][y][x].mGridIndex]
					+= -1/(mTotalThermalGrids[z][y][x].z_res + mTotalThermalGrids[z-1][y][x].z_res);	
				}
	}


/*
#if TEST
	std::cout << "Stamping thermal grids complete !" << std::endl;
#endif
*/

// Stamp sub conductance matrix	due to each thermal grids in metal layer

	std::vector<bool> res_visit(mtotal_unshorted_Res.size(), false);
	std::vector<bool> node_visit(mPGEquivalentNodeSize, false);
	unsigned int count(0);
	
	for(unsigned int y = 0; y < mYCut ; ++y)
	{
		for(unsigned int x = 0; x < mXCut ; ++x)
		{
			MATRIX subConductanceMatrix(mPGEquivalentNodeSize);

			// mark each resistors in each thermal grid
			unsigned int idx(0);
			for(auto & each_res : mtotal_unshorted_Res)
			{
				each_res->is_in_submatrix = false;	// turn of every resistors
				double x_loc_idx = 1 + ( (each_res->shorted_node.first->x_loc_idx-1)*mwire_seg_Xlen + (each_res->shorted_node.second->x_loc_idx-1)*mwire_seg_Xlen ) / 2;
				double y_loc_idx = 1 + ( (each_res->shorted_node.first->y_loc_idx-1)*mwire_seg_Ylen + (each_res->shorted_node.second->y_loc_idx-1)*mwire_seg_Ylen ) / 2;

				if(  x_loc_idx >= mTotalThermalGrids[mZlayer-1][y][x].dx_loc_idx && x_loc_idx <= mTotalThermalGrids[mZlayer-1][y][x].ux_loc_idx && 
					 y_loc_idx >= mTotalThermalGrids[mZlayer-1][y][x].dy_loc_idx && y_loc_idx <= mTotalThermalGrids[mZlayer-1][y][x].uy_loc_idx &&
					 !res_visit[idx])
				{
					res_visit[idx] = true;
					each_res->is_in_submatrix = true; // turn on each resistors in each thermal grid
					each_res->thermal_grid_index = mTotalThermalGrids[mZlayer-1][y][x].mGridIndex;
				}
				else if(  (x == mXCut-1 ) && (y != mYCut-1) && ( x_loc_idx >= 1 + (mX_loc_max-1)*mwire_seg_Xlen ) && 	// PG nodes at X boundary
					      (y_loc_idx <= mTotalThermalGrids[mZlayer-1][y][x].uy_loc_idx ) && !res_visit[idx]) 
				{
					res_visit[idx] = true;
					each_res->is_in_submatrix = true; 
					each_res->thermal_grid_index = mTotalThermalGrids[mZlayer-1][y][x].mGridIndex;
				}
				else if(  (x != mXCut-1) && (y == mYCut-1) && ( x_loc_idx <= mTotalThermalGrids[mZlayer-1][y][x].ux_loc_idx ) && 	// PG nodes at Y boundary
				     	  (y_loc_idx >= 1 + (mY_loc_max-1)*mwire_seg_Ylen ) && !res_visit[idx])
				{
					res_visit[idx] = true;
					each_res->is_in_submatrix = true; 
					each_res->thermal_grid_index = mTotalThermalGrids[mZlayer-1][y][x].mGridIndex;
				}
				else if(  (x == mXCut-1) && (y == mYCut-1) && ( x_loc_idx >= 1 + (mX_loc_max-1)*mwire_seg_Xlen ) && 	// resistor at the last thermal grid (corner)
				     	  (y_loc_idx >= 1 + (mY_loc_max-1)*mwire_seg_Ylen ) && !res_visit[idx])
				{
					res_visit[idx] = true;
					each_res->is_in_submatrix = true; 
					each_res->thermal_grid_index = mTotalThermalGrids[mZlayer-1][y][x].mGridIndex;
				}
				else if(  (x == mXCut-1) && (y == mYCut-1) && ( x_loc_idx <= 1 + (mX_loc_max-1)*mwire_seg_Xlen ) && 	// resistors at the last thermal grid 
				     	  (y_loc_idx <= 1 + (mY_loc_max-1)*mwire_seg_Ylen ) && !res_visit[idx])
				{
					res_visit[idx] = true;
					each_res->is_in_submatrix = true; 
					each_res->thermal_grid_index = mTotalThermalGrids[mZlayer-1][y][x].mGridIndex;
				}
				
				++idx;
			}

			// 2017.8.5
			
			// for delta G matrix
			idx = 0;

			for(auto & each_node : mPGEquivalentNodes)
			{
				if(!node_visit[idx]){

					double x_loc_idx = 1 + (each_node->x_loc_idx-1)*mwire_seg_Xlen;
					double y_loc_idx = 1 + (each_node->y_loc_idx-1)*mwire_seg_Ylen;

					if(  x_loc_idx >= mTotalThermalGrids[mZlayer-1][y][x].dx_loc_idx && x_loc_idx <= mTotalThermalGrids[mZlayer-1][y][x].ux_loc_idx && 
					     y_loc_idx >= mTotalThermalGrids[mZlayer-1][y][x].dy_loc_idx && y_loc_idx <= mTotalThermalGrids[mZlayer-1][y][x].uy_loc_idx )
					{
						each_node->thermal_grid_index = mTotalThermalGrids[mZlayer-1][y][x].mGridIndex;	// mark grid index to each PG node 
						node_visit[idx] = true;	
						++count;
						
					}
					else if(  (x == mXCut-1 ) && (y != mYCut-1) && ( x_loc_idx >= 1 + (mX_loc_max-1)*mwire_seg_Xlen ) && 	// PG nodes at X boundary
					     	  (y_loc_idx <= mTotalThermalGrids[mZlayer-1][y][x].uy_loc_idx ) ) 
					{
						each_node->thermal_grid_index = mTotalThermalGrids[mZlayer-1][y][x].mGridIndex;		
						node_visit[idx] = true;	
						++count;
						
					}
					else if(  (x != mXCut-1) && (y == mYCut-1) && ( x_loc_idx <= mTotalThermalGrids[mZlayer-1][y][x].ux_loc_idx ) && 	// PG nodes at Y boundary
					     	  (y_loc_idx >= 1 + (mY_loc_max-1)*mwire_seg_Ylen ) )
					{
						each_node->thermal_grid_index = mTotalThermalGrids[mZlayer-1][y][x].mGridIndex;	
						node_visit[idx] = true;	
						++count;
						
					}
					else if(  (x == mXCut-1) && (y == mYCut-1) && ( x_loc_idx >= 1 + (mX_loc_max-1)*mwire_seg_Xlen ) && 	 // node at the last thermal grid
					     	  (y_loc_idx >= 1 + (mY_loc_max-1)*mwire_seg_Ylen ) )
					{
						each_node->thermal_grid_index = mTotalThermalGrids[mZlayer-1][y][x].mGridIndex;
						node_visit[idx] = true;	
						++count;
						
					}
					else if(  (x == mXCut-1) && (y == mYCut-1) && ( x_loc_idx <= 1 + (mX_loc_max-1)*mwire_seg_Xlen ) && 	// nodes at the last thermal grid	
					     	  (y_loc_idx <= 1 + (mY_loc_max-1)*mwire_seg_Ylen ) )
					{
						each_node->thermal_grid_index = mTotalThermalGrids[mZlayer-1][y][x].mGridIndex;
						node_visit[idx] = true;	
						++count;
						
					}


				}
				++idx;
			}
			

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
								subConductanceMatrix[stamp_node.second][stamp_node.second] += 1 / each_res->res_value;
							}
						}
					}
					else if (nei_node.first->is_connect_V == true)
					{
						for (auto& each_res : nei_node.second)
						{				
							if(each_res->is_in_submatrix)
							{
								subConductanceMatrix[stamp_node.second][stamp_node.second] += 1 / each_res->res_value;
							}
						}
					}
					else if (nei_node.first->is_connect_V == false && nei_node.first->is_grounded == false)
					{
						for (auto& each_res : nei_node.second)
						{
							if(each_res->is_in_submatrix)
							{
								subConductanceMatrix[stamp_node.second][stamp_node.second] += 1 / each_res->res_value;
								subConductanceMatrix[stamp_node.second][mStampNodeList[nei_node.first]] += -1 / each_res->res_value;
							}
						}
					}
				}
			}
			mPG_ConductanceSubMatrix[mTotalThermalGrids[mZlayer-1][y][x].mGridIndex] = subConductanceMatrix; // store
			// mThermalGridIndexConstainer.emplace_back(mTotalThermalGrids[mZlayer-1][y][x].mGridIndex);
		}
	}

	
	
#if DEBUG	
	
	for(unsigned int i=0 ;i < res_visit.size();++i)
	{
		if(!res_visit[i]) std::cerr << "resistor error !" << std::endl;
	}
	
	std::cout << "Counted nodes = " << count << std::endl;
	
	if(count != mPGEquivalentNodeSize) 
	{
		std::cerr << "Counted Nodes error !" << std::endl;
		for(unsigned int i=0 ; i < node_visit.size() ; ++i)
		{
			if(!node_visit[i])
			{
				std::cout << mPGEquivalentNodes[i]->x_loc_idx << " " << mPGEquivalentNodes[i]->y_loc_idx << std::endl;
			}
		}	
	}
	/*
	for(auto & each_node : mPGEquivalentNodes)
	{
		if(each_node->x_loc_idx == 1 && each_node->y_loc_idx==100)
		{
			std::cout << "thermal grid idx = " << each_node->thermal_grid_index << std::endl;
		}
	}
	*/
#endif
	
	// 2017.8.3 ~ 
	// Load sub-conductance matrix to Eigen
	
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
    	mPG_EigenConductanceSubMatrix[each_matrix.first] = Mat_PGSubMatrix;
	}
	mPG_ConductanceSubMatrix.clear();

/*
#if TEST
	
	std::cout << "Stamping sub conductance matrix complete !" << std::endl;
	std::ofstream fout("origin_conductance_matrix");
	for(unsigned int row_idx = 0; row_idx < mPG_OriginConductanceMatrix.size(); ++row_idx)
    {
        for(auto& row_data : mPG_OriginConductanceMatrix[row_idx])
        {
        	fout << "(" << row_idx << ", " << row_data.first << ")" << " = " << row_data.second << std::endl;
        }
    }
	fout.close();
	fout.open("verify_conductance_matrix");
	EigenSparseMatrix Mat(mPGEquivalentNodeSize, mPGEquivalentNodeSize);
	Mat.setZero();
	for(auto& each_matrix : mPG_EigenConductanceSubMatrix)
	{
		Mat+=each_matrix.second;
	}
	for(int i=0;i<mPGEquivalentNodeSize;++i)
	{
		for(int j=0;j<mPGEquivalentNodeSize;++j)
		{
			if(Mat.coeff(i, j) != 0) 
				fout << "(" << i << ", " << j << ")" << " = " << Mat.coeff(i, j) << std::endl;
		}
	}
	fout.close();
	
#endif
*/
	

//Construct electrical grid to thermal grid incidence matrix 

	for(unsigned int i=0 ;i < node_visit.size() ; ++i)
	{
		node_visit[i] = false;
	}
	
	EigenSparseMatrix Mat_EtoT(mXCut*mYCut*mZlayer, mPGEquivalentNodeSize);	
	std::vector<triplet> triplet_EtoT;
		
	for(unsigned int y = 0 ; y < mYCut ; ++y)
	{
		for(unsigned int x = 0 ; x < mXCut ; ++x)
		{		
			unsigned int node_idx(0);
			for(auto& node : mPGEquivalentNodes)
			{
				if(node->is_connect_I)
				{	
					
					double x_loc_idx = 1 + (node->x_loc_idx-1)*mwire_seg_Xlen;
					double y_loc_idx = 1 + (node->y_loc_idx-1)*mwire_seg_Ylen;	
					unsigned int pgnodeidx((node->x_loc_idx - 1) * mY_loc_max + (node->y_loc_idx - 1));

			 		if(  x_loc_idx >= mTotalThermalGrids[mZlayer-1][y][x].dx_loc_idx && x_loc_idx <= mTotalThermalGrids[mZlayer-1][y][x].ux_loc_idx && 
			 		     y_loc_idx >= mTotalThermalGrids[mZlayer-1][y][x].dy_loc_idx && y_loc_idx <= mTotalThermalGrids[mZlayer-1][y][x].uy_loc_idx 
			 		     && !node_visit[node_idx])
				    {
				    	node_visit[node_idx] = true;
				    	// load data to matrix 
				    	triplet_EtoT.emplace_back(triplet(mTotalThermalGrids[mZ_subCut][y][x].mGridIndex, pgnodeidx, 1));	// grid index of active layer
				    }
				    else if(  (x == mXCut-1 ) && (y != mYCut-1) && ( x_loc_idx >= 1 + (mX_loc_max-1)*mwire_seg_Xlen ) && 	// PG nodes at X boundary
					     	  (y_loc_idx <= mTotalThermalGrids[mZlayer-1][y][x].uy_loc_idx ) 
					     	   && !node_visit[node_idx]) 
					{
						node_visit[node_idx] = true;
						triplet_EtoT.emplace_back(triplet(mTotalThermalGrids[mZ_subCut][y][x].mGridIndex, pgnodeidx, 1));	// grid index of active layer
				    	
					}
					else if(  (x != mXCut-1) && (y == mYCut-1) && ( x_loc_idx <= mTotalThermalGrids[mZlayer-1][y][x].ux_loc_idx ) && 	// PG nodes at Y boundary
					     	  (y_loc_idx >= 1 + (mY_loc_max-1)*mwire_seg_Ylen ) 
					     	   && !node_visit[node_idx]) 
					{
						node_visit[node_idx] = true;
						triplet_EtoT.emplace_back(triplet(mTotalThermalGrids[mZ_subCut][y][x].mGridIndex, pgnodeidx, 1));	// grid index of active layer
					}
					else if(  (x == mXCut-1) && (y == mYCut-1) && ( x_loc_idx >= 1 + (mX_loc_max-1)*mwire_seg_Xlen ) && 	 // node at the last thermal grid
					     	  (y_loc_idx >= 1 + (mY_loc_max-1)*mwire_seg_Ylen ) 
					     	   && !node_visit[node_idx]) 
					{
						node_visit[node_idx] = true;
						triplet_EtoT.emplace_back(triplet(mTotalThermalGrids[mZ_subCut][y][x].mGridIndex, pgnodeidx, 1));	// grid index of active layer
				    	
					}
					else if(  (x == mXCut-1) && (y == mYCut-1) && ( x_loc_idx <= 1 + (mX_loc_max-1)*mwire_seg_Xlen ) && 	// nodes at the last thermal grid	
					     	  (y_loc_idx <= 1 + (mY_loc_max-1)*mwire_seg_Ylen ) 
					     	   && !node_visit[node_idx]) 
					{
						node_visit[node_idx] = true;
						triplet_EtoT.emplace_back(triplet(mTotalThermalGrids[mZ_subCut][y][x].mGridIndex, pgnodeidx, 1));	// grid index of active layer
				    	
					}
					

			    }
			    else if(!node_visit[node_idx])	// nodes that doesn't connect I
			    {
			    	node_visit[node_idx] = true;
			    }
			    ++node_idx;
			}
		}

	}
	Mat_EtoT.setFromTriplets(triplet_EtoT.begin(), triplet_EtoT.end());
	mElectricalGridtoThermalGridIncidenceMatrix = Mat_EtoT;

	
#if DEBUG
	/*
	fout.open("eigen_EtoT_incidence.txt");
	fout << mElectricalGridtoThermalGridIncidenceMatrix;
	fout.close();
	*/
	for(unsigned int i=0 ;i < node_visit.size();++i)
	{
		if(!node_visit[i]) { std::cerr << "Incidence error" << std::endl; }
	}
	std::cout << "Incidence Matrix Complete" << std::endl;
#endif


}



// 2018.4.26
// for expand point method 2
void PowerGrid::StampThermalAwareMatrix_Method2()
{

// Stamp thermal matrix
	
	mPG_ThermalMatrix.resize(mXCut*mYCut*mZlayer);

	for(int z = 0 ; z < mZlayer ; ++z)
		for(int y = 0; y < mYCut ; ++y)
			for(int x = 0; x < mXCut ; ++x)
	{
				
				if(x-1 < 0)   // at x plane min bound
				{
					mPG_ThermalMatrix[mTotalThermalGrids[z][y][x].mGridIndex][mTotalThermalGrids[z][y][x].mGridIndex]
					+= 1/(mTotalThermalGrids[z][y][x].x_res + mTotalThermalGrids[z][y][x+1].x_res);
					
					mPG_ThermalMatrix[mTotalThermalGrids[z][y][x].mGridIndex][mTotalThermalGrids[z][y][x+1].mGridIndex]
					+= -1/(mTotalThermalGrids[z][y][x].x_res + mTotalThermalGrids[z][y][x+1].x_res);
				}
				else if(x-1 >= 0 && x+1 < mXCut) // at x plane between
				{
					mPG_ThermalMatrix[mTotalThermalGrids[z][y][x].mGridIndex][mTotalThermalGrids[z][y][x].mGridIndex]
					+= 1/(mTotalThermalGrids[z][y][x].x_res + mTotalThermalGrids[z][y][x-1].x_res)+
					   1/(mTotalThermalGrids[z][y][x].x_res + mTotalThermalGrids[z][y][x+1].x_res);
					
					mPG_ThermalMatrix[mTotalThermalGrids[z][y][x].mGridIndex][mTotalThermalGrids[z][y][x+1].mGridIndex]
					+= -1/(mTotalThermalGrids[z][y][x].x_res + mTotalThermalGrids[z][y][x+1].x_res);
					
					mPG_ThermalMatrix[mTotalThermalGrids[z][y][x].mGridIndex][mTotalThermalGrids[z][y][x-1].mGridIndex]
					+= -1/(mTotalThermalGrids[z][y][x].x_res + mTotalThermalGrids[z][y][x-1].x_res);
				}
				else if(x+1 >= mXCut) // at x plane max bound
				{
					mPG_ThermalMatrix[mTotalThermalGrids[z][y][x].mGridIndex][mTotalThermalGrids[z][y][x].mGridIndex]
					+= 1/(mTotalThermalGrids[z][y][x].x_res + mTotalThermalGrids[z][y][x-1].x_res);
					
					mPG_ThermalMatrix[mTotalThermalGrids[z][y][x].mGridIndex][mTotalThermalGrids[z][y][x-1].mGridIndex]
					+= -1/(mTotalThermalGrids[z][y][x].x_res + mTotalThermalGrids[z][y][x-1].x_res);	
				}

				if(y-1 < 0)   
				{
					mPG_ThermalMatrix[mTotalThermalGrids[z][y][x].mGridIndex][mTotalThermalGrids[z][y][x].mGridIndex]
					+= 1/(mTotalThermalGrids[z][y][x].y_res + mTotalThermalGrids[z][y+1][x].y_res);
					
					mPG_ThermalMatrix[mTotalThermalGrids[z][y][x].mGridIndex][mTotalThermalGrids[z][y+1][x].mGridIndex]
					+= -1/(mTotalThermalGrids[z][y][x].y_res + mTotalThermalGrids[z][y+1][x].y_res);
				}
				else if(y-1 >= 0 && y+1 < mYCut) 
				{
					mPG_ThermalMatrix[mTotalThermalGrids[z][y][x].mGridIndex][mTotalThermalGrids[z][y][x].mGridIndex]
					+= 1/(mTotalThermalGrids[z][y][x].y_res + mTotalThermalGrids[z][y-1][x].y_res)+
					   1/(mTotalThermalGrids[z][y][x].y_res + mTotalThermalGrids[z][y+1][x].y_res);
					
					mPG_ThermalMatrix[mTotalThermalGrids[z][y][x].mGridIndex][mTotalThermalGrids[z][y+1][x].mGridIndex]
					+= -1/(mTotalThermalGrids[z][y][x].y_res + mTotalThermalGrids[z][y+1][x].y_res);
					
					mPG_ThermalMatrix[mTotalThermalGrids[z][y][x].mGridIndex][mTotalThermalGrids[z][y-1][x].mGridIndex]
					+= -1/(mTotalThermalGrids[z][y][x].y_res + mTotalThermalGrids[z][y-1][x].y_res);	
				}
				else if(y+1 >= mYCut) 
				{
					mPG_ThermalMatrix[mTotalThermalGrids[z][y][x].mGridIndex][mTotalThermalGrids[z][y][x].mGridIndex]
					+= 1/(mTotalThermalGrids[z][y][x].y_res + mTotalThermalGrids[z][y-1][x].y_res);
					
					mPG_ThermalMatrix[mTotalThermalGrids[z][y][x].mGridIndex][mTotalThermalGrids[z][y-1][x].mGridIndex]
					+= -1/(mTotalThermalGrids[z][y][x].y_res + mTotalThermalGrids[z][y-1][x].y_res);		
				}
					
				if(z-1 < 0)   
				{
					mPG_ThermalMatrix[mTotalThermalGrids[z][y][x].mGridIndex][mTotalThermalGrids[z][y][x].mGridIndex]
					+= 1/(mTotalThermalGrids[z][y][x].z_res + mTotalThermalGrids[z][y][x].z_bound_res)+
					   1/(mTotalThermalGrids[z][y][x].z_res + mTotalThermalGrids[z+1][y][x].z_res);
					
					mPG_ThermalMatrix[mTotalThermalGrids[z][y][x].mGridIndex][mTotalThermalGrids[z+1][y][x].mGridIndex]
					+= -1/(mTotalThermalGrids[z][y][x].z_res + mTotalThermalGrids[z+1][y][x].z_res);
				}
				else if(z-1 >= 0 && z+1 < mZlayer) 
				{
					mPG_ThermalMatrix[mTotalThermalGrids[z][y][x].mGridIndex][mTotalThermalGrids[z][y][x].mGridIndex]
					+= 1/(mTotalThermalGrids[z][y][x].z_res + mTotalThermalGrids[z-1][y][x].z_res)+
					   1/(mTotalThermalGrids[z][y][x].z_res + mTotalThermalGrids[z+1][y][x].z_res);
					
					mPG_ThermalMatrix[mTotalThermalGrids[z][y][x].mGridIndex][mTotalThermalGrids[z+1][y][x].mGridIndex]
					+= -1/(mTotalThermalGrids[z][y][x].z_res + mTotalThermalGrids[z+1][y][x].z_res);
					
					mPG_ThermalMatrix[mTotalThermalGrids[z][y][x].mGridIndex][mTotalThermalGrids[z-1][y][x].mGridIndex]
					+= -1/(mTotalThermalGrids[z][y][x].z_res + mTotalThermalGrids[z-1][y][x].z_res);	
				}
				else if(z+1 >= mZlayer) 
				{
					mPG_ThermalMatrix[mTotalThermalGrids[z][y][x].mGridIndex][mTotalThermalGrids[z][y][x].mGridIndex]
					+= 1/(mTotalThermalGrids[z][y][x].z_res + mTotalThermalGrids[z][y][x].z_bound_res)+
					   1/(mTotalThermalGrids[z][y][x].z_res + mTotalThermalGrids[z-1][y][x].z_res);
					
					mPG_ThermalMatrix[mTotalThermalGrids[z][y][x].mGridIndex][mTotalThermalGrids[z-1][y][x].mGridIndex]
					+= -1/(mTotalThermalGrids[z][y][x].z_res + mTotalThermalGrids[z-1][y][x].z_res);	
				}
	}


/*
#if TEST
	std::cout << "Stamping thermal grids complete !" << std::endl;
#endif
*/

// Stamp sub conductance matrix	due to each thermal grids in metal layer

	std::vector<bool> res_visit(mtotal_unshorted_Res.size(), false);
	std::vector<bool> node_visit(mPGEquivalentNodeSize, false);
	unsigned int count(0);
	
	for(unsigned int y = 0; y < mYCut ; ++y)
	{
		for(unsigned int x = 0; x < mXCut ; ++x)
		{
			// mark each resistors in each thermal grid
			unsigned int idx(0);
			for(auto & each_res : mtotal_unshorted_Res)
			{
				each_res->is_in_submatrix = false;	// turn of every resistors
				double x_loc_idx = 1 + ( (each_res->shorted_node.first->x_loc_idx-1)*mwire_seg_Xlen + (each_res->shorted_node.second->x_loc_idx-1)*mwire_seg_Xlen ) / 2;
				double y_loc_idx = 1 + ( (each_res->shorted_node.first->y_loc_idx-1)*mwire_seg_Ylen + (each_res->shorted_node.second->y_loc_idx-1)*mwire_seg_Ylen ) / 2;

				if(  x_loc_idx >= mTotalThermalGrids[mZlayer-1][y][x].dx_loc_idx && x_loc_idx <= mTotalThermalGrids[mZlayer-1][y][x].ux_loc_idx && 
					 y_loc_idx >= mTotalThermalGrids[mZlayer-1][y][x].dy_loc_idx && y_loc_idx <= mTotalThermalGrids[mZlayer-1][y][x].uy_loc_idx &&
					 !res_visit[idx])
				{
					res_visit[idx] = true;
					each_res->is_in_submatrix = true; // turn on each resistors in each thermal grid
					each_res->thermal_grid_index = mTotalThermalGrids[mZlayer-1][y][x].mGridIndex;
				}
				else if(  (x == mXCut-1 ) && (y != mYCut-1) && ( x_loc_idx >= 1 + (mX_loc_max-1)*mwire_seg_Xlen ) && 	// PG nodes at X boundary
					      (y_loc_idx <= mTotalThermalGrids[mZlayer-1][y][x].uy_loc_idx ) && !res_visit[idx]) 
				{
					res_visit[idx] = true;
					each_res->is_in_submatrix = true; 
					each_res->thermal_grid_index = mTotalThermalGrids[mZlayer-1][y][x].mGridIndex;
				}
				else if(  (x != mXCut-1) && (y == mYCut-1) && ( x_loc_idx <= mTotalThermalGrids[mZlayer-1][y][x].ux_loc_idx ) && 	// PG nodes at Y boundary
				     	  (y_loc_idx >= 1 + (mY_loc_max-1)*mwire_seg_Ylen ) && !res_visit[idx])
				{
					res_visit[idx] = true;
					each_res->is_in_submatrix = true; 
					each_res->thermal_grid_index = mTotalThermalGrids[mZlayer-1][y][x].mGridIndex;
				}
				else if(  (x == mXCut-1) && (y == mYCut-1) && ( x_loc_idx >= 1 + (mX_loc_max-1)*mwire_seg_Xlen ) && 	// resistor at the last thermal grid (corner)
				     	  (y_loc_idx >= 1 + (mY_loc_max-1)*mwire_seg_Ylen ) && !res_visit[idx])
				{
					res_visit[idx] = true;
					each_res->is_in_submatrix = true; 
					each_res->thermal_grid_index = mTotalThermalGrids[mZlayer-1][y][x].mGridIndex;
				}
				else if(  (x == mXCut-1) && (y == mYCut-1) && ( x_loc_idx <= 1 + (mX_loc_max-1)*mwire_seg_Xlen ) && 	// resistors at the last thermal grid 
				     	  (y_loc_idx <= 1 + (mY_loc_max-1)*mwire_seg_Ylen ) && !res_visit[idx])
				{
					res_visit[idx] = true;
					each_res->is_in_submatrix = true; 
					each_res->thermal_grid_index = mTotalThermalGrids[mZlayer-1][y][x].mGridIndex;
				}
				
				++idx;
			}

			
			idx = 0;

			for(auto & each_node : mPGEquivalentNodes)
			{
				if(!node_visit[idx]){

					double x_loc_idx = 1 + (each_node->x_loc_idx-1)*mwire_seg_Xlen;
					double y_loc_idx = 1 + (each_node->y_loc_idx-1)*mwire_seg_Ylen;

					if(  x_loc_idx >= mTotalThermalGrids[mZlayer-1][y][x].dx_loc_idx && x_loc_idx <= mTotalThermalGrids[mZlayer-1][y][x].ux_loc_idx && 
					     y_loc_idx >= mTotalThermalGrids[mZlayer-1][y][x].dy_loc_idx && y_loc_idx <= mTotalThermalGrids[mZlayer-1][y][x].uy_loc_idx )
					{
						each_node->thermal_grid_index = mTotalThermalGrids[mZlayer-1][y][x].mGridIndex;	// mark grid index to each PG node 
						node_visit[idx] = true;	
						++count;
						
					}
					else if(  (x == mXCut-1 ) && (y != mYCut-1) && ( x_loc_idx >= 1 + (mX_loc_max-1)*mwire_seg_Xlen ) && 	// PG nodes at X boundary
					     	  (y_loc_idx <= mTotalThermalGrids[mZlayer-1][y][x].uy_loc_idx ) ) 
					{
						each_node->thermal_grid_index = mTotalThermalGrids[mZlayer-1][y][x].mGridIndex;		
						node_visit[idx] = true;	
						++count;
						
					}
					else if(  (x != mXCut-1) && (y == mYCut-1) && ( x_loc_idx <= mTotalThermalGrids[mZlayer-1][y][x].ux_loc_idx ) && 	// PG nodes at Y boundary
					     	  (y_loc_idx >= 1 + (mY_loc_max-1)*mwire_seg_Ylen ) )
					{
						each_node->thermal_grid_index = mTotalThermalGrids[mZlayer-1][y][x].mGridIndex;	
						node_visit[idx] = true;	
						++count;
						
					}
					else if(  (x == mXCut-1) && (y == mYCut-1) && ( x_loc_idx >= 1 + (mX_loc_max-1)*mwire_seg_Xlen ) && 	 // node at the last thermal grid
					     	  (y_loc_idx >= 1 + (mY_loc_max-1)*mwire_seg_Ylen ) )
					{
						each_node->thermal_grid_index = mTotalThermalGrids[mZlayer-1][y][x].mGridIndex;
						node_visit[idx] = true;	
						++count;
						
					}
					else if(  (x == mXCut-1) && (y == mYCut-1) && ( x_loc_idx <= 1 + (mX_loc_max-1)*mwire_seg_Xlen ) && 	// nodes at the last thermal grid	
					     	  (y_loc_idx <= 1 + (mY_loc_max-1)*mwire_seg_Ylen ) )
					{
						each_node->thermal_grid_index = mTotalThermalGrids[mZlayer-1][y][x].mGridIndex;
						node_visit[idx] = true;	
						++count;
						
					}

				}
				++idx;
			}
			mThermalGridIndexConstainer.emplace_back(mTotalThermalGrids[mZlayer-1][y][x].mGridIndex);
		}
	}

	
	
#if DEBUG	
	
	for(unsigned int i=0 ;i < res_visit.size();++i)
	{
		if(!res_visit[i]) std::cerr << "resistor error !" << std::endl;
	}
	
	std::cout << "Counted nodes = " << count << std::endl;
	
	if(count != mPGEquivalentNodeSize) 
	{
		std::cerr << "Counted Nodes error !" << std::endl;
		for(unsigned int i=0 ; i < node_visit.size() ; ++i)
		{
			if(!node_visit[i])
			{
				std::cout << mPGEquivalentNodes[i]->x_loc_idx << " " << mPGEquivalentNodes[i]->y_loc_idx << std::endl;
			}
		}	
	}
	/*
	for(auto & each_node : mPGEquivalentNodes)
	{
		if(each_node->x_loc_idx == 1 && each_node->y_loc_idx==100)
		{
			std::cout << "thermal grid idx = " << each_node->thermal_grid_index << std::endl;
		}
	}
	*/

#endif
	

//Construct electrical grid to thermal grid incidence matrix 

	for(unsigned int i=0 ;i < node_visit.size() ; ++i)
	{
		node_visit[i] = false;
	}
	
	EigenSparseMatrix Mat_EtoT(mXCut*mYCut*mZlayer, mPGEquivalentNodeSize);	
	std::vector<triplet> triplet_EtoT;
		
	for(unsigned int y = 0 ; y < mYCut ; ++y)
	{
		for(unsigned int x = 0 ; x < mXCut ; ++x)
		{		
			unsigned int node_idx(0);
			for(auto& node : mPGEquivalentNodes)
			{
				if(node->is_connect_I)
				{	
					
					double x_loc_idx = 1 + (node->x_loc_idx-1)*mwire_seg_Xlen;
					double y_loc_idx = 1 + (node->y_loc_idx-1)*mwire_seg_Ylen;	
					unsigned int pgnodeidx((node->x_loc_idx - 1) * mY_loc_max + (node->y_loc_idx - 1));

			 		if(  x_loc_idx >= mTotalThermalGrids[mZlayer-1][y][x].dx_loc_idx && x_loc_idx <= mTotalThermalGrids[mZlayer-1][y][x].ux_loc_idx && 
			 		     y_loc_idx >= mTotalThermalGrids[mZlayer-1][y][x].dy_loc_idx && y_loc_idx <= mTotalThermalGrids[mZlayer-1][y][x].uy_loc_idx 
			 		     && !node_visit[node_idx])
				    {
				    	node_visit[node_idx] = true;
				    	// load data to matrix 
				    	triplet_EtoT.emplace_back(triplet(mTotalThermalGrids[mZ_subCut][y][x].mGridIndex, pgnodeidx, 1));	// grid index of active layer
				    }
				    else if(  (x == mXCut-1 ) && (y != mYCut-1) && ( x_loc_idx >= 1 + (mX_loc_max-1)*mwire_seg_Xlen ) && 	// PG nodes at X boundary
					     	  (y_loc_idx <= mTotalThermalGrids[mZlayer-1][y][x].uy_loc_idx ) 
					     	   && !node_visit[node_idx]) 
					{
						node_visit[node_idx] = true;
						triplet_EtoT.emplace_back(triplet(mTotalThermalGrids[mZ_subCut][y][x].mGridIndex, pgnodeidx, 1));	// grid index of active layer
				    	
					}
					else if(  (x != mXCut-1) && (y == mYCut-1) && ( x_loc_idx <= mTotalThermalGrids[mZlayer-1][y][x].ux_loc_idx ) && 	// PG nodes at Y boundary
					     	  (y_loc_idx >= 1 + (mY_loc_max-1)*mwire_seg_Ylen ) 
					     	   && !node_visit[node_idx]) 
					{
						node_visit[node_idx] = true;
						triplet_EtoT.emplace_back(triplet(mTotalThermalGrids[mZ_subCut][y][x].mGridIndex, pgnodeidx, 1));	// grid index of active layer
					}
					else if(  (x == mXCut-1) && (y == mYCut-1) && ( x_loc_idx >= 1 + (mX_loc_max-1)*mwire_seg_Xlen ) && 	 // node at the last thermal grid
					     	  (y_loc_idx >= 1 + (mY_loc_max-1)*mwire_seg_Ylen ) 
					     	   && !node_visit[node_idx]) 
					{
						node_visit[node_idx] = true;
						triplet_EtoT.emplace_back(triplet(mTotalThermalGrids[mZ_subCut][y][x].mGridIndex, pgnodeidx, 1));	// grid index of active layer
				    	
					}
					else if(  (x == mXCut-1) && (y == mYCut-1) && ( x_loc_idx <= 1 + (mX_loc_max-1)*mwire_seg_Xlen ) && 	// nodes at the last thermal grid	
					     	  (y_loc_idx <= 1 + (mY_loc_max-1)*mwire_seg_Ylen ) 
					     	   && !node_visit[node_idx]) 
					{
						node_visit[node_idx] = true;
						triplet_EtoT.emplace_back(triplet(mTotalThermalGrids[mZ_subCut][y][x].mGridIndex, pgnodeidx, 1));	// grid index of active layer
				    	
					}
					

			    }
			    else if(!node_visit[node_idx])	// nodes that doesn't connect I
			    {
			    	node_visit[node_idx] = true;
			    }
			    ++node_idx;
			}
		}

	}
	Mat_EtoT.setFromTriplets(triplet_EtoT.begin(), triplet_EtoT.end());
	mElectricalGridtoThermalGridIncidenceMatrix = Mat_EtoT;

	
#if DEBUG
	/*
	fout.open("eigen_EtoT_incidence.txt");
	fout << mElectricalGridtoThermalGridIncidenceMatrix;
	fout.close();
	*/
	for(unsigned int i=0 ;i < node_visit.size();++i)
	{
		if(!node_visit[i]) { std::cerr << "Incidence error" << std::endl; }
	}
	std::cout << "Incidence Matrix Complete" << std::endl;
#endif


}



void PowerGrid::ComputeInversePGConductanceMatrix()
{

	EigenSparseMatrix EigenMatrix(mPGEquivalentNodeSize, mPGEquivalentNodeSize);
					 		
	// Load matrix data to eigen 
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
	for(unsigned int i = 0 ; i < mPGEquivalentNodeSize ; ++i)
	{
		pattern[i] = 1;
		mEachRowPGInverseMatrix.emplace_back(solver.solve(pattern));
		pattern[i] = 0;
	}

}


void PowerGrid::ComputeInverseThermalConductanceMatrix()
{

	unsigned int thermalsize = mPG_ThermalMatrix.size();
	EigenSparseMatrix EigenMatrix(thermalsize, thermalsize);
						  
	// Load matrix data to eigen 
	std::vector<triplet> tripletList;
	for(unsigned int row_idx = 0; row_idx < thermalsize; ++row_idx)
    {
        for(auto& row_data : mPG_ThermalMatrix[row_idx])
        {           
            tripletList.push_back(triplet(row_idx, row_data.first, row_data.second));   
        }
    }
	EigenMatrix.setFromTriplets(tripletList.begin(), tripletList.end());
	tripletList.clear();

	mDenseInverseThermalMatrix = Eigen::MatrixXd(EigenMatrix).inverse();

}



// 2017.7.31

// Thermal-less
void PowerGrid::ConstructThermalLessVerificationData(int METHOD)
{
	
	ConstructCurrentConstraint(0, METHOD);
	
	StampPGConductanceMatrix();		

}


// QP
void PowerGrid::ConstructThermalAwareVerificationData_QP(int METHOD)
{
	
	ConstructCurrentConstraint(1, METHOD);		
	
	ConstructThermalGrids();		
	
	StampPGConductanceMatrix();		

	StampThermalAwareMatrix(); 

	ComputeInversePGConductanceMatrix();

	ComputeInverseThermalConductanceMatrix();

}


// Expand Point
void PowerGrid::ConstructThermalAwareVerificationData_ExpandPoint(int METHOD)
{
	
	ConstructCurrentConstraint(2, METHOD);
	
	ConstructThermalGrids();		
	
	StampPGConductanceMatrix();		/* no use if using expand point method 2 */  

	StampThermalAwareMatrix(); 

	ComputeInversePGConductanceMatrix();		/* no use if using expand point method 2 */

	ComputeInverseThermalConductanceMatrix();

}

// Expand Point Method 2
void PowerGrid::ConstructThermalAwareVerificationData_ExpandPoint_Method2(int METHOD)
{
	
	ConstructCurrentConstraint(2, METHOD);
	
	ConstructThermalGrids();		
	
	StampThermalAwareMatrix_Method2(); 

	ComputeInverseThermalConductanceMatrix();

}


void PowerGrid::ConstructMonteCarloData()
{

	ConstructCurrentConstraint(0, 1);

	ConstructThermalGrids();	
	
	StampPGConductanceMatrix();	

	StampThermalAwareMatrix(); 

}

// added 2018.8.18 ~
// Expand Point using matlab
void PowerGrid::ConstructThermalAwareVerificationData_ExpandPoint_Matlab(int METHOD)
{
	
	ConstructCurrentConstraint(1, METHOD);		//  construct matrix for matlab  
	
	ConstructThermalGrids();		
	
	StampPGConductanceMatrix();		/* no use if using expand point method 2 */  

	StampThermalAwareMatrix(); 

	ComputeInversePGConductanceMatrix();		/* no use if using expand point method 2 */

	ComputeInverseThermalConductanceMatrix();

}


}	

