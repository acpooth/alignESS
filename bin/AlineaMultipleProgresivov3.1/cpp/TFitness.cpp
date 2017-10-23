#include "../h/TFitness.h"

/*
	Clase base para las funciones fitness
*/


TFitness::TFitness(void)
{
}

TFitness::~TFitness(void)
{
}


/*
	Constructor
*/
TFitness::TFitness(vector<Secuencia*> s,int mL)
{
	sequences = s;
	maxLength = mL;
	numSeqs = s.size();
}

/*
	Auxiliar para la longitud máxima, cuando los individuos cambian por los operadores
*/
void TFitness::setMaxLength(int m)
{
	maxLength = m;
}


/*
	Para calcular el fitness por bloques
	Le falta revisión
*/
double TFitness::calculateBlocks(map<int, EC> matrix,int point, int blockSize)
{
	
	double fitness = homoBlocks(matrix,point, blockSize);
					//- penalty(matrixWithoutGaps);
	return fitness;
}


/*
	Penalización por incremento de columnas..
	Es un valor positivo,  de 0 a 1, mientras mas cercano a 1 es peor
*/
double TFitness::incrementoColumnas(vector<EC> matrix)
{
        double p = 0.0;
        int sizePrintMatrix = matrix.size();
        int MAX = -999;
        int aux = 0;
        for(int iSeq = 0; iSeq < numSeqs; iSeq++)//Columnas
        {
                aux = (sequences[iSeq])->getLength();
                if ( aux > MAX)
                {
                        MAX = aux;
                }
        }

        p = (double)MAX / (double)sizePrintMatrix;
        p -= 1;
        p = p * -1;
        return p;
}
