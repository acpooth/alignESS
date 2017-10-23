#include <sstream>
#include <cstring>
#include <cstdlib>
#include "../h/Secuencia.h"

using namespace std;

///////////////////****Utilerias...***////////////
vector<string> &split(const string &s, char delim, vector<string> &elems) 
{
	stringstream ss(s);
	string item;
	while(getline(ss, item, delim)) {
		elems.push_back(item);
	}
	return elems;
}


vector<string> split(const string &s, char delim) 
{
	vector<string> elems;
	return split(s, delim, elems);
}

string intToString(int n )
{
	std::ostringstream result;
	result << n;
	return result.str();
}
////////////////////////////////////////////////


Secuencia::Secuencia(void)
{
}

Secuencia::~Secuencia(void)
{
}


/*
	Constructor
	s = Cadena completa como esta en el archivo de texto: 2.1.3:4.3.2: ...
	n = Nombre del cluster >Nombre
	s = secuencia de ec numbers separados por :
	map = nombre del mapa
	seqI = id de la secuencia, en 3 digitos 001 - 452
	map = nombre del mapa: >Mapa1
*/
Secuencia::Secuencia(string n, string s,string map,string seqI)
{
	setName(n);
	setSequence(s);
	setMapId(map);
	setSeqId(seqI);
	vector<string> ecs = split(s,':');//Me devuelve el vector con los ec numbers
	//////////////////////
	//Remover los gaps de la secuencia..
	vector<string>::iterator it;
	for(it = ecs.begin(); it != ecs.end(); it++)
	{
		if( (*it).compare("-.-.-") == 0)
			ecs.erase(it);
	}
	/////////////////////
	setLength(ecs.size());
	for (int j = 0; j < length; j++)
	{			
		EC ec;
		vector<string> num = split(ecs[j],'.'); //Me devuelve los 3 numeros
		ec.push_back(atoi(num[0].c_str())); ec.push_back(atoi(num[1].c_str())); ec.push_back(atoi(num[2].c_str()));
		enzymes.push_back(ec);		
	}
}

/*
	Constructor
	s = Cadena completa como esta en el archivo de texto: 2.1.3:4.3.2: ...
	Se usa un nombre por default para el mapa
*/
Secuencia::Secuencia(string s)
{	
}

/*
	l = número de enzimas
*/
void Secuencia::setLength(int l)
{
	length = l;
}

/*
	e = la lista de enzimas ya divididas, cada elemento de la lista es un vector de 3 elementos
*/
void Secuencia::setEnzymes(vector<EC> e)
{
	enzymes = e;
}

/*
	m = nombre del mapa
*/
void Secuencia::setMapId(string m)
{
	mapId = m;
}

/*
	s = id de la secuencia
*/
void Secuencia::setSeqId(string s)
{
	seqId = s;
}

/*
	nombre del cluster
*/
void Secuencia::setName(string n)
{
	name = n;
}

/*
	Secuencia de ec numbersi, sin dividir
*/
void Secuencia::setSequence(string s)
{
	sequence = s;
}

/*
	numero de enzimas que hay, osea, bloques de 3 enteros
*/
int Secuencia::getLength()
{
	return length;
}

/*
	Aqui ya estan divididas en bloques de 3
*/
SEQUENCE Secuencia::getEnzymes()
{
	return enzymes;
}

/*
	Aqui se devuelve una enzima en particular
*/
EC Secuencia::getEnzyme(int pos)
{
	return enzymes[pos];
}

/*
	Se convierte el arreglo de 3 enteros a cadena
*/
string Secuencia::getEnzymeStr(int pos)
{
	string s;
	s = intToString((enzymes[pos])[0]) + "." +
		intToString((enzymes[pos])[1]) + "." +
		intToString((enzymes[pos])[2]);
	return s;
}

/*
	El id del mapa
*/
string Secuencia::getMapId()
{
	return mapId;
}

/*
	Identificador de la secuencia
*/
string Secuencia::getSeqId()
{
	return seqId;
}

/*
	Nombre del clusters
*/
string Secuencia::getName()
{
	return name;
}

/*
	Cadena que contiene la secuencia de ec numbers
*/
string Secuencia::getSequence()
{
	return sequence;
}

/*
	Se imprime toda la secuencia
*/
void Secuencia::print()
{
	cout << getSequence() << endl;
}
