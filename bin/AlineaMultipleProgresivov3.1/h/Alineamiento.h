#pragma once
#include <cstdlib>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include "../h/Secuencia.h"
#include "../h/Individuo.h"

class Alineamiento
{
private:
	int numSeq;
	int numExp;
	double pc;
	double pm;
	int numGen;
	int sizePob;
	int seqId;//consecutivo de todas las secuencias	
	vector<Secuencia*> sequences;
	vector<int> seqLengths;
	vector<Secuencia*> sequencesParciales;
	vector<int> seqLengthsParciales;
	int maxLength;
	Individuo best;
	int iSeqMasLarga;
	double penGap;
	double penCol;
	double facHomo;
public:
	map<string, double> matrixCalis;
	Alineamiento(void);
	~Alineamiento(void);
	Alineamiento(string nameFile, string fileCalis,double c,double m, int sP, int nG, int fT, int cT, double gap, double col, double homo, int nExp);
	//Getters and setters
	void setGAParams(double c,double m, int sP, int nG);
	void setNumSeq(int n);
	void setSeqLengths();
	int getNumSeq();
	double getPc();
	double getPm();
	int getSizePob();
	int getNumGen();
	vector<int> getSeqLengths();
	//Auxiliares
	void readFile(string n);
	void readMatrixCalis(string n);
	void createAlineamientoProgresivo(int fT, int cT);
	vector<EC> convertirAlineadosBin(vector<string> alineados);
	Secuencia* createSequence(string nMap, string seqID, string lineSeq,string nameCluster);
	Individuo createGA(int fT, int cT);
	Individuo createGAAlineados(vector<EC> convertidos, int fitnessType, int crossType, int numSeqSinAlinear);
	Individuo getBest();
	vector<Secuencia*> getSequences();
	void ordenarSecuencias();
};
