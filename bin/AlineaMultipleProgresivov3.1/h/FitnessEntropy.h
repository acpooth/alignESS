#pragma once
#include "../h/TFitness.h"

class FitnessEntropy :
	public TFitness
{
private:
	double maxValue;
	double penGap;
	double penCol;
	double facHomo;
public:
	FitnessEntropy(void);
	~FitnessEntropy(void);
	FitnessEntropy(vector<Secuencia*> sD,int mLD);
	FitnessEntropy(vector<Secuencia*> sD,int mLD, double gap, double col, double homo);
	double calculate(map<int, EC> matrix, vector<EC> matrixAux);	
	double homoColumns(vector<EC> matrixAux);
	double homoPerColumn(COLUMN col);
	double penalty(vector<EC> matrixAux);
	double homoBlocks(map<int, EC> alignMatrix,int point, int blockSize);
};
