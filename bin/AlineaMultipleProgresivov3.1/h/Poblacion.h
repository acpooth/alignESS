#pragma once
#include <cstdlib>
#include <iostream>
#include <vector>
#include "../h/Secuencia.h"
#include "../h/Individuo.h"
#include "../h/TFitness.h"

class Poblacion
{
private:
	double totalFitness;
	int sizePop;
	Individuo best;	
	vector<Secuencia*> sequences;
	vector<int> seqLengths;
	int maxLength;
	int numSeqs;
	TFitness* fitness;
	int numSeqSinAlinear;
public:
	vector<Individuo> individuals;
	Poblacion(void);
	~Poblacion(void);
	Poblacion(int sP, vector<Secuencia*> s,vector<int> sL, int mL, int fT);
	Poblacion(int sP, vector<Secuencia*> s,vector<int> sL, int mL, int fT, double gap, double col, double homo);
	Poblacion(vector<EC> convertidos, int sP, vector<Secuencia*> s,vector<int> sL, int mL, int fT, int nSSA);
	Poblacion(vector<EC> convertidos, int sP, vector<Secuencia*> s,vector<int> sL, int mL, int fT, int nSSA,double g, double c,double h);
	void createInitialPopulation();
	void createInitialPopulationAlineados(vector<EC> alineadas);
	vector<int> getSeqLengths();
	vector<Secuencia*> getSequences();
	int getNumSeqs();
	int getSizePop();
	int getMaxLength();
	void setMaxLength(int m);
	double getTotalFitness();	
	Individuo getBestIndividual();
	void calculateFitness();
	double getAverageFitness();	
	//////////////
	double fitnessPerBlock(int point,int blockSize);
	TFitness* getTipoFitness();
};
