#pragma once
#include <cstdlib>
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include "../h/Secuencia.h"
#include "../h/TFitness.h"

class Individuo
{
private:
	int numSeq;
	int maxLength;
	vector<int> seqLengths;
	double fitness;	
	vector<EC> printMatrix;
	vector<EC> alineados;
	vector<string> printSecuencias;
	int numSeqSinAlinear;
public:
	map<int, EC> alignMatrix;
	Individuo(void);
	~Individuo(void);
	//Individuo(const Individuo& copy);
	Individuo(vector<int> sL, int mL);
	Individuo(vector<EC> a, int sinAlinear, vector<int> sL, int mL);
	double calculateFitness(TFitness* f,vector<Secuencia*> seqs);
	void setMaxLength(int m);
	double getFitness();
	int getMaxLength();
	double getRelFitness(double tF);	
	void initializeMatrix();
	void initializeMatrixAlineados();
	void print(vector<Secuencia*> seqs);
	void printAux();
	double calculateFitnessPerBlock(TFitness* f,vector<Secuencia*> seqs,int p,int bz);
	void convertirMatrix(vector<EC> pM);
	vector<string> getPrintSecuencias(vector<Secuencia *> seqs);
};
