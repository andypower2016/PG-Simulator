#include "PowerGridModel.hpp"

namespace PowerGridModel {

Resistor::Resistor() 
 {}

Resistor::Resistor(std::string n, std::string pos, std::string neg, double value)
 : res_name(n), pos_node(pos), neg_node(neg), res_value(value) 
 {
 	is_in_submatrix = false;
 }

void Resistor::Add_Node(Node* n1, Node* n2)
{
	shorted_node = std::make_pair(n1, n2);
}

VoltageSource::VoltageSource() 
 {}
	
VoltageSource::VoltageSource(std::string n, std::string pos, std::string neg, double v)
 : volt_name(n), pos_node(pos), neg_node(neg), volt_value(v) 
 {}

 
CurrentSource::CurrentSource(){}
	
CurrentSource::CurrentSource(std::string n, std::string pos, std::string neg, double c)
 : curr_name(n), pos_node(pos), neg_node(neg), curr_value(c)
 {
	if(neg == "0") ParseNodeInfo(pos, "_");  
	else 		   ParseNodeInfo(neg, "_");  
 }

void CurrentSource::ParseNodeInfo
(const std::string &line, const std::string &seperators)
{
	std::vector<std::string> line_element;
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
	
	if(line_element.size() == 2){	//  x_y => (x, y)
		x_loc_idx = std::stoul(line_element[0]);
		y_loc_idx = std::stoul(line_element[1]);
	}
}

Node::Node()
 {}
	
Node::Node(std::string n) 
{
 	node_name    = n;
 	is_connect_V = false;
 	is_grounded  = false;
	is_connect_I = false;
	is_stamp     = true;
	ParseNodeInfo(n, "n_");	
}
	

void Node::ParseNodeInfo 
(const std::string &line, const std::string &seperators)
{
	std::vector<std::string> line_element;
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
	
	if(line_element.size() == 2){	//  x_y => (x, y)
		x_loc_idx = std::stoul(line_element[0]);
		y_loc_idx = std::stoul(line_element[1]);
	}
	
}


void Node::Add_Neighbor_Node(Node* nei, Resistor* res)		
{
		if(nei!=NULL && res!=NULL) neighbor_NodeRes[nei].emplace_back(res);
		else std::cerr<<"Null Neighbor Node, or resistor !"<<std::endl; 
}
	
void Node::Add_Vsource(std::string n, VoltageSource* v)
{
		if(v!=NULL) {
			connected_Vsource[n].emplace_back(v);
			is_connect_V = true;
		}
		else std::cerr<<"Null voltage source, or resistor !"<<std::endl; 
}
	
void Node::Add_Isource(std::string n, CurrentSource* c)
{
		if(c!=NULL) connected_Isource[n].emplace_back(c);
		else std::cerr<<"Null current source, or resistor !"<<std::endl; 
}
	


PowerGridData::PowerGridData() 
 {}

}

