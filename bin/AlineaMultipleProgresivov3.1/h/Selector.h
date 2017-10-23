#pragma once
#include "../h/Poblacion.h"

class Selector
{
private:
	double fitTotal;
	int totalPuntos;
	vector<double> fitRelativo;
	vector<double> probAcumulada;
	map<int,double> fitnessPuntos;
	map<int,double> probAcumuladaPuntos;        
public:
	Selector(void);
	~Selector(void);
	int length;
	int blockSize;
	vector<int> points;
	vector<double> fitness;
	Selector(int blockSize, int length);
	int getTotalPuntos();
	double getFitTotal();
	void actualizaFitness(Poblacion& pob, int blockSize);		
	int selectPunto(double r);		
	void initializePoints(int blockSize, int length);
};
