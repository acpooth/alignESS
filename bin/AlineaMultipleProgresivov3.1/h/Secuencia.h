#pragma once

#include <iostream>
#include <vector>
#include <map>
#include <string>

using namespace std;
typedef vector<int> EC;
typedef vector<EC> SEQUENCE;

class Secuencia
{
private: 
	int length;
	SEQUENCE enzymes;
	string mapId;
	string seqId;
	string name;
	string sequence;
public:
	Secuencia(void);
	Secuencia(string n,string all,string map,string seqI);
	Secuencia(string all);//Sin nombre del mapa, se genera uno default
	~Secuencia(void);
	//Getters and setters
	void setLength(int l);
	void setEnzymes(vector<EC> e);
	void setMapId(string m);
	void setSeqId(string s);
	void setName(string n);
	void setSequence(string seq);
	int getLength();
	vector<EC> getEnzymes();
	EC getEnzyme(int pos);
	string getEnzymeStr(int pos);
	string getMapId();
	string getSeqId();
	string getName();
	string getSequence();
	void print();
};

