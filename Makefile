CXX = g++

CXXFLAGS = -O3 -std=c++11 -Wall -g -pthread -DTARGET_API_VERSION=700 -DUSE_MEX_CMD -D_GNU_SOURCE -DMATLAB_MEX_FILE -fexceptions -fPIC -fno-omit-frame-pointer

VPATH = src

Mat_INCLUDE = -I"/usr/local/Matlab/R2017a/extern/include" -I"/usr/local/Matlab/R2017a/simulink/include"  

Mat_LIB_DIR = -Wl,--no-undefined "/usr/local/Matlab/R2017a/sys/os/glnxa64/libstdc++.so.6" -Wl,-rpath-link,/usr/local/Matlab/R2017a/bin/glnxa64 -L"/usr/local/Matlab/R2017a/bin/glnxa64" 

Mat_LIBS    = -lmx -lmex -lmat -lm -leng -lstdc++ 

MOSEK_LIB_PATH = /home/edausers/jacob24788/mosek/bin

MOSEK_INCLUDE = -I"/home/edausers/jacob24788/mosek/h"

MOSEK_LIBS = $(MOSEK_LIB_PATH)/libiomp5.so $(MOSEK_LIB_PATH)/libmosek64.so $(MOSEK_LIB_PATH)/libcilkrts.so.5 $(MOSEK_LIB_PATH)/libmosek64.so.8.0

Eigen_INCLUDE = -I"/home/edausers/jacob24788/Eigen"

EXE = ./bin/run	

OBJ = ./obj

MAIN = ./main

SRC = ./src

all:$(OBJ)/Parser.o $(OBJ)/PowerGridModel.o $(OBJ)/PowerGrid.o $(OBJ)/PGSolver.o $(OBJ)/Run.o $(OBJ)/ThermalModel.o
	$(CXX) $(CXXFLAGS) $(Eigen_INCLUDE) $(MOSEK_INCLUDE) $(Mat_INCLUDE) $(MAIN)/main.cpp $^ $(MOSEK_LIBS) $(Mat_LIB_DIR) $(Mat_LIBS) -o $(EXE)
	rm $(OBJ)/*.o
	@echo -e "Make target [$@] is complete !";

$(OBJ)/Run.o: Run.cpp
	$(CXX) $(CXXFLAGS) -c $(Eigen_INCLUDE) $(MOSEK_INCLUDE) $(Mat_INCLUDE) $< -o $@
	@echo -e "Make target [$@] is complete !";
	
$(OBJ)/Parser.o: Parser.cpp 
	$(CXX) $(CXXFLAGS) -c $(Eigen_INCLUDE) $(MOSEK_INCLUDE) $(Mat_INCLUDE) $< -o $@
	@echo -e "Make target [$@] is complete !";
	
$(OBJ)/PGSolver.o: PGSolver.cpp
	$(CXX) $(CXXFLAGS) -c $(Eigen_INCLUDE) $(MOSEK_INCLUDE) $(Mat_INCLUDE) $< -o $@
	@echo -e "Make target [$@] is complete !";

$(OBJ)/PowerGrid.o: PowerGrid.cpp
	$(CXX) $(CXXFLAGS) -c $(Eigen_INCLUDE) $< -o $@
	@echo -e "Make target [$@] is complete !";	
	
$(OBJ)/PowerGridModel.o: PowerGridModel.cpp
	$(CXX) $(CXXFLAGS) -c $(Eigen_INCLUDE) $< -o $@
	@echo -e "Make target [$@] is complete !";
	
$(OBJ)/ThermalModel.o: ThermalModel.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
	@echo -e "Make target [$@] is complete !";
	
StaticLib:	
	@echo -e "Make target [$@] is complete !";
	
clean:
	rm -f $(EXE) $(OBJ)/*.o
	@echo -e "Make target [$@] is complete !";