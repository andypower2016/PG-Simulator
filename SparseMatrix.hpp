#ifndef SPARSE_Matrix
#define SPARSE_Matrix

#include <vector>

namespace MySparseMatrixClass {

class SparseVector;

class SparseMatrix {

public : 

	SparseMatrix() {}

	SparseMatrix(const SparseMatrix & InputMatrix) 
	{
		mEachVector = InputMatrix.getMatrix();
	}

	~SparseMatrix() {}

	void PushBack(SparseVector & InputVector)
	{
		mEachVector.emplace_back(InputVector);
	}

	const SparseVector & getSparseVector(unsigned int idx) const { return mEachVector[idx]; }

	const std::vector<SparseVector> & getMatrix() const { return mEachVector; }
	
private : 

	std::vector<SparseVector> mEachVector;

};


class SparseVector {

public : 

	SparseVector() {}
	
	SparseVector(int NonZeros, int *Positions, double *NonzeroElements, unsigned int *IndexAdder) 
	{
		if(NonZeros != 0)
		{
			for(int i = 0 ; i < NonZeros ; ++i)
			{
				mLocation.emplace_back(Positions[*IndexAdder]);
				mVector.emplace_back(NonzeroElements[*IndexAdder]);
				(*IndexAdder)++;
			}
		}
	}

	~SparseVector() {}

	void PushBackElement(int Position, double Value)
	{
		mLocation.emplace_back(Position);
		mVector.emplace_back(Value);
	}

	const std::vector<double> & getVector() const { return mVector; }

	const std::vector<unsigned int> & getLocation() const { return mLocation; }

	
private : 

	std::vector<double> mVector;

	std::vector<unsigned int> mLocation;


};

}

#endif














