#pragma once
#include "../h/Secuencia.h"
#include <cstdlib>
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <sstream>
#include <boost/regex.hpp>

typedef vector<EC> COLUMN;
//Utilerias....
inline string trim_right(const string &source , const string& t = "0")
{
	string str = source;
	return str.erase( str.find_last_not_of(t) + 1);
}

inline string trim_left( const string& source, const string& t = "0")
{
	string str = source;
	return str.erase(0 , source.find_first_not_of(t) );
}

inline string trim(const string& source, const string& t = "0")
{
	string str = source;
	return trim_left( trim_right( str , t) , t );
} 

inline string int_ToString(int n) 
{
	std::ostringstream result;
        result << n;
        return result.str();

}
//////////////////////////////

class TFitness
{
protected:
	vector<Secuencia*> sequences;
	int maxLength;
	int numSeqs;
	double openingGap;//Para penalizar gaps en SumaPares
	double extGap;	//Para penalizar gaps en SumaPares
	vector<int> factors;//Para ponderar en FitnessEntropia
public:
	TFitness(void);
	~TFitness(void);
	TFitness(vector<Secuencia*> s,int mL);
	void setMaxLength(int m);
	double incrementoColumnas(vector<EC> matrixAux);
	double calculateBlocks(map<int, EC> matrix,int point, int blockSize);
	virtual double calculate(map<int, EC> matrix, vector<EC> matrixAux) = 0;	
	virtual double homoColumns(vector<EC> matrixAux) = 0;
	virtual double homoPerColumn(COLUMN col) = 0;
	virtual double penalty(vector<EC> matrixAux) = 0;
	virtual double homoBlocks(map<int, EC> alignMatrix,int point, int blockSize) = 0;
};
