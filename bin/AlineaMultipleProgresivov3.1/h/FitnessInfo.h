#pragma once
#include "../h/TFitness.h"

class FitnessInfo :
	public TFitness
{
public:
	FitnessInfo(void);
	~FitnessInfo(void);
	FitnessInfo(vector<Secuencia*> sD,int mLD);
	double homoColumns(vector<EC> matrixAux);
	double homoPerColumn(COLUMN col);
	double homoBlocks(map<int, EC> alignMatrix,int point, int blockSize);
	double calculate(map<int, EC> matrix, vector<EC> matrixAux);	
	double penalty(vector<EC> matrixAux);
};
