#pragma once
#include "../h/TFitness.h"

class FitnessSumaPares :
	public TFitness
{
public:
	FitnessSumaPares(void);
	~FitnessSumaPares(void);
	FitnessSumaPares(vector<Secuencia*> sD,int mLD);
	double homoColumns(vector<EC> matrixAux);
	double homoPerColumn(COLUMN col);
	double homoBlocks(map<int, EC> alignMatrix,int point, int blockSize);
	double penalty(vector<EC> matrixAux);
	double calculate(map<int, EC> matrix, vector<EC> matrixAux);	
};
