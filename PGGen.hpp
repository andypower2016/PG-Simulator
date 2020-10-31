//* PG circuit generator
//* 2017.7.15 by 楊卓諭
//*
//* * * * Program specifications * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *                                                                    
//*  
//*  1. input file format : 
//*
//*  cktPara.txt is the circuit parameter .txt file constructed as the input format below.
//*  
//*  line1 : [CKT FILE NAME] [ckt_name] 
//*  line2 : [DIMENSION] [Row Dimension] [Column Dimension] 
//*  line3 : [VMOD and CMOD] [v location mod] [i location mod] 
//*  line4 : [Vdd] [voltage] 
//*  line5 : [Wire Reisstance] [ohms] 
//*  line6 : [Current value] [load current] 
//*  line7 : [PG length] [x length(m)] [y length(m)]
//*
//*  Note : Vdd connects a internal resistance to the power grid
//*         The node name direct connects to Vdd is specified as n_x1_y1
//*  
//*  
//* 
//*  2. PG circuit example : 
//*
//*  ex. 3x3 power grid
//*             
//*         (1, 1) --- (1, 2) --- (1, 3)
//*           |          |          |
//*         (2, 1) --- (2, 2) --- (2, 3)
//*           |          |          |
//*         (3, 1) --- (3, 2) --- (3, 3) 
//*             
//* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
//* last update : 2017.8.8 


#ifndef PGGEN_H
#define PGGEN_H

#include "Parser.hpp"

using namespace std;

class PGGen : public Parser
{

public : 

PGGen() {}

void GenerateCircuit()
{
	
	ifstream fin("../PGGen_input/cktPara.txt");

    string outputName = "";

    /*  cktPara.txt informations */

    int RowDimension(0); int ColDimension(0); 

    int Vmod(0);  int Cmod(0);

    double VDD(0);

    double minimumRes(0);

    double CURR_LOAD(0);

    double PGXlength(0), PGYlength(0);

    /*  Parse cktPara.txt  */
    
    int lineNum(0);
    while( std::getline(fin, line) )
    {
    
        ParseLine(line,"    \n\r", line_elements);  
           
        if(line_elements.size()!=0 && line_elements[0].at(0) != '#')
        {
            switch(lineNum)
            {
                case 0 : 
                    outputName = line_elements[1];
                    break;                 
                case 1 :                 
                    RowDimension = stoi(line_elements[1]);
                    ColDimension = stoi(line_elements[2]);
                    break;                          
                case 2 :                
                    Vmod = stoi(line_elements[1]);
                    Cmod = stoi(line_elements[2]); 
                    break;  
                case 3 :
                    VDD = stod(line_elements[1]);
                    break;
                case 4 :  
                    minimumRes = stod(line_elements[1]);
                    break;
                case 5 :  
                    CURR_LOAD = stod(line_elements[1]);
                    break;
                case 6 : 
                    PGXlength = stod(line_elements[1]);
                    PGYlength = stod(line_elements[2]);
                    break;
                default : 
                    break;
            }
            ++lineNum;
        }       
    }
    fin.close();

  
    // Generate PG Circuit 

    ofstream fout("../PGGen_result/" + outputName);
    ofstream VDfout("../PGGen_result/V_D");
    ofstream CDfout("../PGGen_result/C_D");

    fout << "*CKT_FILE_NAME    " << outputName << " " << endl;
    fout << "*ROL_COL          " << RowDimension << " " << ColDimension << endl;
    fout << "*VMOD_CMOD        " << Vmod << " " << Cmod << endl;
    fout << "*VDD_VOLT         " << VDD << endl;
    fout << "*MIN_WIRE_RES     " << minimumRes << endl;       
    fout << "*CURRENT_VALUE    " << CURR_LOAD << endl;  
    fout << "*PG x length =    " << PGXlength << " m " << endl;
    fout << "*PG y length =    " << PGYlength << " m " << endl;
    fout << "*WIRE_LENGTH      " << PGXlength/(RowDimension-1) << " " << PGYlength/(ColDimension-1) << endl;  // x, y wire segment length(m)  
    fout << endl;

    cout << "Minimum wire resistance = " << minimumRes << " ohm" << endl;

    int r_idx(0);
    for (int i = 1 ; i <= RowDimension ; i++)   // Generate row resistors
    {
        for (int j = 1 ; j < ColDimension ; j++)
        {
            fout<<"R" << r_idx << "_row" << " " 
            << i << "_" << j << " " 
            << i << "_" << j + 1 << " " 
            << randMToN(minimumRes, 2 * minimumRes) << endl;
            ++r_idx;
        }
    }


    for (int i = 1 ; i <= ColDimension ; i++)   // Generate column resistors 
    {
        for (int j = 1 ; j < RowDimension ; j++)
        {
            fout<<"R" << r_idx << "_col" << " " 
            << j << "_" << i << " " 
            << j + 1 << "_" << i << " " 
            << randMToN(minimumRes, 2 * minimumRes) << endl;
            ++r_idx;
        }
    }

  
    // Generate Vdds and current loads
   
    int loadCnt(0);
    double totalPower(0), totalCurrent(0);

    for (int i=1 ; i < RowDimension * ColDimension ; i++)   
    {
        
        if ( (i % Vmod) == 0) {

            int X(0), Y(0);
            bool exit(false);
            for (int x=1 ; x <= RowDimension ; ++x){            
                for (int y=1 ; y <= ColDimension ; ++y)
                {
                    if( (x-1)*(ColDimension) + y == i )
                    {
                        X = x;
                        Y = y;
                        exit = true;
                        break;
                    }
                }
                if(exit) break;
            }
            
            fout << "R" << i << " " 
                 << X << "_" << Y << " " 
                 << "n_" 
                 << X << "_" << Y << " " 
                 << "0.001" << endl;    


            fout << "V" << i << " " 
                 << "n_"
                 << X << "_" << Y << " " 
                 << "0" << " " 
                 << VDD << endl; 

            VDfout << X << " " << Y << endl;

        }
        else if ( (i % Cmod) == 0 ){

            int X(0), Y(0);
            bool exit(false);
            for (int x=1 ; x <= RowDimension ; ++x){            
                for (int y=1 ; y <= ColDimension ; ++y)
                {
                    if( (x-1)*(ColDimension) + y == i )
                    {
                        X = x;
                        Y = y;
                        exit = true;
                        break;
                    }
                }
                if(exit) break;
            }

            double current( randMToN(CURR_LOAD * 0.5, CURR_LOAD * 2) );
           
            fout << "I" << i << " " 
                 << X << "_" << Y << " " 
                 << "0" << " " 
                 <<  current << endl; 

            ++loadCnt;
            totalCurrent += current;
            totalPower += VDD * current;
            
            CDfout << X << " " << Y << endl;
        }
    }

    fout << ".op" << endl;
    fout << ".end" << endl;

    fout.close();
    VDfout.close();
    CDfout.close();

    std::cout << "total load    = " << loadCnt << std::endl;
    std::cout << "total current = " << totalCurrent << std::endl;
    std::cout << "total power   = " << totalPower << std::endl;

    ParseInput("../PGGen_result/" + outputName, outputName);
    
    ConstructPowerGridPlane();
    SolvePowerGridPlane(mTotalNodes);
	OutputDesignInfo();

}
		

double randMToN(double M, double N)
{
    return M + (rand() / ( RAND_MAX / (N-M) ) ) ;  
}

private : 
    
    string line;
    
    vector<string> line_elements;

};

#endif