#include "Parser.hpp"


Parser::Parser()
{}

Parser::Parser(const std::string& file_path) 
{} 

void Parser::ParseLine(const std::string& line, const std::string& seperators, std::vector<std::string>& line_element)
{
	line_element.clear();
	std::string token;
	for(size_t i=0;i<line.length();i++){
		if( seperators.find(line[i]) != std::string::npos ){ 
			if(token.size()!=0){
				line_element.push_back(token);
				token.clear();
			}	
		}
		else token.push_back(line[i]);
	}
	if(token.size()!=0){
		line_element.push_back(token);
		token.clear();
	}	
}

void Parser::ParseInput(const std::string &file_path, const std::string &cktname)
{ 
	
	
// Parse .sp file 

	std::ifstream in(file_path, std::ios::in);

if(!in.is_open()){
	
	std::cerr << "no such file under current directory !" << std::endl;	 
	in.close();
	exit(1);

}
else{	
	
	std::cout << "Open " << cktname << " success !" << std::endl;

	mCktName = cktname;
	
	std::vector<std::string> line_elements;	 
	std::string line;
	while( std::getline(in, line) )
	{
	
	ParseLine(line,"	 \n", line_elements);   
	
	if(line_elements.size()!=0)
		switch( char( toupper(line_elements[0].at(0)) ) ) 
		{
		 	   

case 'V':
		 	   	       																				
				if( stod(line_elements[3]) !=0 )
				{		
						volt = new PowerGridModel::VoltageSource
						(line_elements[0], line_elements[1], line_elements[2], stod(line_elements[3]));
						
						if(mTotalNodes[line_elements[1]] == NULL)
						{
							node = new PowerGridModel::Node(line_elements[1]);
							mTotalNodes[line_elements[1]] = node;
						}
						
						if(mTotalNodes[line_elements[2]] == NULL)
						{
							node = new PowerGridModel::Node(line_elements[2]);
							mTotalNodes[line_elements[2]] = node;						
						}
											
						mTotalNodes[line_elements[1]]->Add_Vsource("pos", volt);
						mTotalNodes[line_elements[2]]->Add_Vsource("neg", volt);
						mtotal_supply_Vsource[volt] = mTotalNodes[line_elements[1]];
												
				}									
				break;
		 	   			 	   	
case 'I':

				curr = new PowerGridModel::CurrentSource
					(line_elements[0], line_elements[1], line_elements[2], stod(line_elements[3])); 
				
				if (mTotalNodes[line_elements[1]] == NULL)
				{
					node = new PowerGridModel::Node(line_elements[1]);
					mTotalNodes[line_elements[1]] = node;
				}
				

				if (mTotalNodes[line_elements[2]] == NULL)
				{
					node = new PowerGridModel::Node(line_elements[2]);
					mTotalNodes[line_elements[2]] = node;
				}
				
				mTotalNodes[line_elements[1]]->Add_Isource("pos", curr);
				mTotalNodes[line_elements[1]]->is_connect_I = true;
				mTotalNodes[line_elements[2]]->Add_Isource("neg", curr);											
				mtotal_Isource.emplace_back(curr);

		 	   	break;	 	   	     
case 'R':
		 	   	
				res = new PowerGridModel::Resistor
					(line_elements[0], line_elements[1], line_elements[2], stod(line_elements[3]));         	
					
				if (mTotalNodes[line_elements[1]] == NULL)
				{
					node = new PowerGridModel::Node(line_elements[1]);
					mTotalNodes[line_elements[1]] = node;
				}
				

				if (mTotalNodes[line_elements[2]] == NULL)
				{
					node = new PowerGridModel::Node(line_elements[2]);
					mTotalNodes[line_elements[2]] = node;
				}
				
				mTotalNodes[line_elements[1]]->
					Add_Neighbor_Node(mTotalNodes[line_elements[2]], res);
				
				mTotalNodes[line_elements[2]]->
					Add_Neighbor_Node(mTotalNodes[line_elements[1]], res);
				
				
				res->Add_Node
					(mTotalNodes[line_elements[1]], mTotalNodes[line_elements[2]]);

				mtotal_unshorted_Res.emplace_back(res);
				
		 	   	break;		 	   	     
case '*' : 		
				if( line_elements[0] == "*ROL_COL" )
				{
					mX_loc_max = stod(line_elements[1]);
					mY_loc_max = stod(line_elements[2]);
				}	
				else if(line_elements[0] == "*WIRE_LENGTH")
				{
					mwire_seg_Xlen = stod(line_elements[1]);
					mwire_seg_Ylen = stod(line_elements[2]);
					//std::cout << "Wire length = " << mwire_seg_len << std::endl;
				}
				else if(line_elements[0] == "*PG" && line_elements[1] == "x" && line_elements[2] == "length")
				{
					mPowerPlaneXLength = stod(line_elements[4]);
					//std::cout << "X length = " << mPowerPlaneXLength << std::endl;
				}
				else if(line_elements[0] == "*PG" && line_elements[1] == "y" && line_elements[2] == "length")
				{
					mPowerPlaneYLength = stod(line_elements[4]);
					//std::cout << "Y length = " << mPowerPlaneYLength << std::endl;	
				}
				break;
case '.':				
				break;						
default :		 		   	        
		 	    break;
			}				 
	}

	std::cout << "Parse input file complete ! " << std::endl;
	in.close();
}

}

