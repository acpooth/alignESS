#include <cstdlib>
#include <iostream>
#include <vector>
#include "../h/Poblacion.h"
#include "../h/Secuencia.h"
#include "../h/Individuo.h"
#include "../h/FitnessInfo.h"
#include "../h/FitnessEntropy.h"
#include "../h/FitnessSumaPares.h"

/*
	Clase que crea las poblaciones para el AG, y controla las operaciones sobre los individuos
	de la misma.	
*/

Poblacion::Poblacion(void)
{
}

Poblacion::~Poblacion(void)
{
}

/*
	sp = Tamaño de la poblacion
	s = vector con las secuencias que se van a alinear (objetos Secuencia)
	sL = vector con las longitudes de las secuencias
	mL = Cual es la longitud máxima inicial
	fT = tipo de Fitness (1,2,3)
	Este es el constructor normal.
*/
Poblacion::Poblacion(int sP, vector<Secuencia*> s,vector<int> sL, int mL, int fT, double g, double c, double h)
{
	sequences = s;
	sizePop = sP;
	seqLengths = sL;
	numSeqSinAlinear = 1;
	numSeqs = sL.size();
	maxLength = mL;
	createInitialPopulation();
	switch(fT)
	{
		case 1:
			fitness = new FitnessEntropy(sequences, maxLength,g,c,h);
			break;
		case 2:
			fitness = new FitnessInfo(sequences, maxLength);
			break;
		case 3:
			fitness = new FitnessSumaPares(sequences, maxLength);
			break;
		default:
			fitness = new FitnessEntropy(sequences, maxLength,g,c,h);			
	}
	
}
/*
	Agregada 210911 : Solo le agregue los parametros de ponderación para FitnessEntropy
	sp = Tamaño de la poblacion
	s = vector con las secuencias que se van a alinear (objetos Secuencia)
	sL = vector con las longitudes de las secuencias
	mL = Cual es la longitud máxima inicial
	fT = tipo de Fitness (1,2,3)
	convertidos = vector con las secuencias que ya estan alineadas (en 0's y 1's)
	nSSA = el número de secuencias que ya se alinearon
*/
Poblacion::Poblacion(vector<EC> convertidos, int sP, vector<Secuencia*> s,vector<int> sL, int mL, int fT, int nSSA,double g, double c, double h)
{
	numSeqSinAlinear = nSSA;
	sequences = s;
	sizePop = sP;
	seqLengths = sL;
	numSeqs = sL.size();
	maxLength = mL;
	createInitialPopulationAlineados(convertidos);
	switch(fT)
	{
		case 1:
			fitness = new FitnessEntropy(sequences, maxLength,g,c,h);
			break;
		case 2:
			fitness = new FitnessInfo(sequences, maxLength);
			break;
		case 3:
			fitness = new FitnessSumaPares(sequences, maxLength);
			break;
		default:
			fitness = new FitnessEntropy(sequences, maxLength,g,c,h);			
	}
	
}
/*
	sp = Tamaño de la poblacion
	s = vector con las secuencias que se van a alinear (objetos Secuencia)
	sL = vector con las longitudes de las secuencias
	mL = Cual es la longitud máxima inicial
	fT = tipo de Fitness (1,2,3)
	convertidos = vector con las secuencias que ya estan alineadas (en 0's y 1's)
	nSSA = el número de secuencias que ya se alinearon
*/
Poblacion::Poblacion(vector<EC> convertidos, int sP, vector<Secuencia*> s,vector<int> sL, int mL, int fT, int nSSA)
{
	numSeqSinAlinear = nSSA;
	sequences = s;
	sizePop = sP;
	seqLengths = sL;
	numSeqs = sL.size();
	maxLength = mL;
	createInitialPopulationAlineados(convertidos);
	switch(fT)
	{
		case 1:
			fitness = new FitnessEntropy(sequences, maxLength);
			break;
		case 2:
			fitness = new FitnessInfo(sequences, maxLength);
			break;
		case 3:
			fitness = new FitnessSumaPares(sequences, maxLength);
			break;
		default:
			fitness = new FitnessEntropy(sequences, maxLength);			
	}
	
}

/*
	sp = Tamaño de la poblacion
	s = vector con las secuencias que se van a alinear (objetos Secuencia)
	sL = vector con las longitudes de las secuencias
	mL = Cual es la longitud máxima inicial
	fT = tipo de Fitness (1,2,3)
	Este es el constructor normal.
*/
Poblacion::Poblacion(int sP, vector<Secuencia*> s,vector<int> sL, int mL, int fT)
{
	sequences = s;
	sizePop = sP;
	seqLengths = sL;
	numSeqSinAlinear = 1;
	numSeqs = sL.size();
	maxLength = mL;
	createInitialPopulation();
	switch(fT)
	{
		case 1:
			fitness = new FitnessEntropy(sequences, maxLength);
			break;
		case 2:
			fitness = new FitnessInfo(sequences, maxLength);
			break;
		case 3:
			fitness = new FitnessSumaPares(sequences, maxLength);
			break;
		default:
			fitness = new FitnessEntropy(sequences, maxLength);			
	}
	
}

/*
	Se crea la población inicial, aleatoria
*/
void Poblacion::createInitialPopulation()
{	
	for(int i = 0; i < sizePop; i++)
	{
		Individuo ind(seqLengths, maxLength);
		individuals.push_back(ind);
	}
}

/*
	Se crea la población inicial, aleatoria, para el AG iterativo
*/
void Poblacion::createInitialPopulationAlineados(vector<EC> alineadas)
{	
	
	for(int i = 0; i < sizePop; i++)
	{
		Individuo ind(alineadas, numSeqSinAlinear, seqLengths, maxLength);
		individuals.push_back(ind);
	}
}

/*
	Devuel el numero de individuos en la población
*/
int Poblacion::getSizePop()
{
	return sizePop;
}

/*
	Devuelve el vector de longitudes de las secuencias
*/
vector<int> Poblacion::getSeqLengths()
{
	return seqLengths;
}

/*
	Devuelve el vector con los objetos Secuencia
*/
vector<Secuencia*> Poblacion::getSequences()
{
	return sequences;
}

/*
	Devuelve el numero de secuencias
*/
int Poblacion::getNumSeqs()
{
	return numSeqs;
}

/*
	Devuelve la longitud máxima de la secuencia, es el tamaño de la matriz. Pero no siempre es confiable. Depende de los 
	operadores que se esten usando.
*/
int Poblacion::getMaxLength()
{
	return maxLength;
}

/*
	Esta es una función auxiliar para definir la máxima longitud de la secuencia (permitida), sirve cuando se usan
	los operadores que modifican las longitudes.
*/
void Poblacion::setMaxLength(int m)
{
	maxLength = m;
	fitness->setMaxLength(m);
}

/*
	Regresa el fitness total de la población, Esta variable se actualiza cuando se llama a calcularFitness
*/
double Poblacion::getTotalFitness()
{
	return totalFitness;
}

/*
	Regresa el mejor individuo, dependiendo del tipo de Fitness puede ser el que tenga menor o mayor fitness
	en toda la población.
*/
Individuo Poblacion::getBestIndividual()
{
 vector<Individuo>::iterator i;
 if (dynamic_cast<FitnessEntropy*>(fitness) == NULL)
 {
	double MAX = -999999999999999999999.0;
 	int j = 0;

	 for (i = individuals.begin(); i != individuals.end(); i++)
     	{
                if ((*i).getFitness() > MAX)
                {
                        MAX = (*i).getFitness();
                        best = (*i);
                }
                j++;
     	}

 }
else{
	double MIN = 999999999999999999999.0;
	 int j = 0;
 
	 for (i = individuals.begin(); i != individuals.end(); i++)
	     {
		if ((*i).getFitness() < MIN)
		{
			MIN = (*i).getFitness();
			best = (*i);
		}
		j++;
	     }	

}
  return best;		
}

/*
	Se devuelve un apuntador al objeto fitness
*/
TFitness* Poblacion::getTipoFitness()
{
	return fitness;
}

/*
	Esta función manda llamar al calculo del fitness de cada individuo y actualiza la variable totalFitness
*/
void Poblacion::calculateFitness()
{
	totalFitness = 0.0;
	for(int i = 0; i < sizePop; i++)
	{
		totalFitness += individuals[i].calculateFitness(fitness,sequences);
	}
}

/*
	Aqui se devuelve el fitness promedio de la población.
*/
double Poblacion::getAverageFitness()
{
	return (double)totalFitness / (double)sizePop;
}

/*
	Es la función que se llama cuando se van a usar bloques para el fitness
*/
double Poblacion::fitnessPerBlock(int point,int blockSize)
{
	double totalFitness = 0.0;
	for(int i = 0; i < sizePop; i++)
	{
		totalFitness += individuals[i].calculateFitnessPerBlock(fitness,sequences,point,blockSize);
		//if (individuals[i].calculateFitnessPerBlock(fitness,sequences,point,blockSize) > .8)
		//	return 1;

	}
	return totalFitness;
}
