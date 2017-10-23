#pragma once
#include <cstdlib>
#include <iostream>
#include <vector>
#include <valarray>
#include "../h/Alineamiento.h"
#include "../h/Secuencia.h"
#include "../h/Individuo.h"
#include "../h/Poblacion.h"
#include "../h/Selector.h"

inline double getRandom()
{
	return (double)rand () / (double) RAND_MAX;
}
class AG
{
private:
	int maxGen;
	int sizePop;
	double crossRate;
	double mutRate;		
	Poblacion population;
	vector<Individuo> individuals;
	int newMaxLength;
	int numSeqSinAlinear;
	Selector select;
public:
	AG(void);
	~AG(void);
	AG(int mG, double cR, double mR, Poblacion p);
	AG(int mG, double cR, double mR, Poblacion p, int numSSA);
	double getCrossRate();
	double getMutRate();
	int getMaxGen();
	int getSizePop();
	vector<Individuo> getPopulation();
	Individuo actualizarMaxLengthBest(Individuo b);
	//Selección
	vector<Individuo> selection(vector<Individuo> oldPop,int tmsize,int elitist,Individuo best);
	Individuo tournament(vector<Individuo> oldPop,int tmsize);
	vector<Individuo> sus(vector<Individuo> oldPop, int size,Individuo best);
	//Mutacion..
	vector<Individuo> mutationGaps(vector<Individuo> oldPop);
	vector<Individuo> mutationGaps2(vector<Individuo> oldPop);
	vector<Individuo> ajustarMaxLength(vector<Individuo> oldPop, int max);
	void mutationPermutate(vector<Individuo> oldPop);
	vector<Individuo> mutationCircular(vector<Individuo> oldPop);
	//Recombinación
	vector<Individuo> crossoverColumnas2(vector<Individuo> oldPop);
	vector<int> crossoverColumnas2Child(valarray<int> circular,valarray<int> circular2,int maxUnos,int maxLength,int maxLength2,int random);
	vector<int> crossoverColumnas3Child(valarray<int> circular,valarray<int> circular2,int maxUnos,int maxLength,int random);
	vector<Individuo> crossoverExchange(vector<Individuo> oldPop);
	vector<Individuo> crossoverColumns(vector<Individuo> oldPop);
	void checkSequences(vector<int>& c1,int numEnz,int maxLength);
	vector<Individuo> crossoverColumnsSelectPoint(vector<Individuo> oldPop, Selector select);
	//Ciclo principal
	Individuo start();
        Individuo startMax();
        Individuo start2();//Usa los operadores por bloques
        Individuo start2Max();//Usa los operadores por bloques
};
