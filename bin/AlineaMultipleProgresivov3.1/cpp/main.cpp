#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <ctime>
#include <cmath>
#include <map>
#include <sys/stat.h> 
#include "../h/Alineamiento.h"
#include "../h/Individuo.h"

using namespace std;

bool validateInput(int argc,char* argv[])
{
	bool ban = false;
	if (argc < 11)
	{
		cout << "	Uso:"<<endl;
		cout << "		Alinea archivo archCalificaciones tamPob numGen tasaCruza tasaMuta penalizacionGap factorHomogeneidad penalizacionIncCol numExp" <<endl;
		cout << "	donde:"<<endl;
		cout << "		archivo: Nombre del archivo que tiene las secuencias para alinear"<<endl;
		cout << "		archCalif: Nombre del archivo que tiene las evaluacion por pares"<<endl;
		cout << "		tamaPob: número de individuos en la población"<<endl;
		cout << "		NumGen:  número de generaciones máximo del algoritmo" << endl;
		cout << "		tasaCruza: probabilidad con la que se realiza la cruza"<<endl;
		cout << "		tasaMuta:  probabilidad con la que se realiza la mutacion" << endl;
		cout << "		penGap:  penalizacion por apertura de gap (.05)" << endl;
		cout << "		factorHomo:  factor de ponderacion para la homogeneidad (.9)" << endl;
		cout << "		penIncCol:  penalización por incremento de columnas (.05)" << endl;
		cout << "		numExp:  numero de repeticiones del AG" << endl;
		
		
	}
	else{
		if (argc == 11)
		return true;
	}
	return ban;
}

int main(int argc, char* argv[])
{
	//Iniciamos la semilla aleatoria
	time_t t = time(0);	
	srand((unsigned)t);
	if (!validateInput(argc,argv))
		return 0;
	string file = argv[1]; //archivo secuencias
	string fileCalis = argv[2]; //archivo evaluaciones por pares
	double crossRate =  atof(argv[5]);//Tamaño poblacion
	double mutRate = atof(argv[6]);//Numero de generaciones
	int sizePop =  atoi(argv[3]);//Numero de generaciones
	int numGen =  atoi(argv[4]);//Numero de generaciones
	double penGap = atof(argv[7]);
	double homo = atof(argv[8]);
	double incremento = atof(argv[9]);
	int numExp = atoi(argv[10]);
	//Alineamiento* a = new Alineamiento(file, crossRate, mutRate, sizePop, numGen, fit, crossoverType);
	int fit = 1; // se deja el default Entropia
	int crossoverType = 1; // se deja el default el normal, no por bloques
        Individuo b;
        Alineamiento* a;
	a = new Alineamiento(file, fileCalis,crossRate, mutRate, sizePop, numGen, fit, crossoverType, penGap, incremento, homo, numExp);
        return 0;

}

