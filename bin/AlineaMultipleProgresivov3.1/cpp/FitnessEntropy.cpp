#include "../h/FitnessEntropy.h"


/*
	Clase que usa la entropia en las columnas para calcular la homogeneidad. Mientras mas cercana a 0 es mejor
	Se regresa un promedio de la entropia en las columnas.
	Por columna penalizamos la proporción de gaps.
*/


FitnessEntropy::FitnessEntropy(void)
{
}

FitnessEntropy::~FitnessEntropy(void)
{
}

/*	
	Estos son los factores que se usan para ponderar los 3 niveles de los Ec numbers 
*/
FitnessEntropy::FitnessEntropy(vector<Secuencia*> s,int mL) : TFitness(s,mL)
{
	factors.push_back(15);
	factors.push_back(10);
	factors.push_back(1);
}

/*	
	Estos son los factores que se usan para ponderar los 3 niveles de los Ec numbers 
*/
FitnessEntropy::FitnessEntropy(vector<Secuencia*> s,int mL, double gap, double col, double homo) : TFitness(s,mL)
{
	factors.push_back(15);
	factors.push_back(10);
	factors.push_back(1);
	penGap = gap;
	penCol = col;
	facHomo = homo;
}

/*
	Función que pondera las distintas partes que se califican..
	Homogeneidad, gaps e incremento de columnas
*/
double FitnessEntropy::calculate(map<int, EC> matrix, vector<EC> matrixWithoutGaps)
{
	
	double homo = homoColumns(matrixWithoutGaps);
	double penal =  penalty(matrixWithoutGaps);
	double incremento = incrementoColumnas(matrixWithoutGaps);
	incremento = (incremento);
	double fitness = 0.0;
	//Esta es para la 1.1 
	fitness = homo*facHomo + penal * penGap + incremento*penCol;
	//cout << "HOMO = " << homo << " FACHOMO = " << facHomo << " PENAL = " << penal << " penGap = " << penGap << " INC = " << penCol << " F = " << fitness << endl;
	//fitness = homo;
	return fitness;
}

/*
	La calificación es la suma de la entropia de todas las columnas
	Una menor calificación indica mas homogeneidad en la columna..
	Si toda la columna es homogenea la calificación es 0.
*/
double FitnessEntropy::homoColumns(vector<EC> matrix)
{
	vector<int> counts(numSeqs,0);
	vector<int> gap(3,0);
	double entropy = 0.0;
	int sizePrintMatrix = matrix.size();
	for(int iCol = 0; iCol < sizePrintMatrix; iCol++)//Columnas
	{
		COLUMN aux;
		EC line;
		
		for(int iSeq = 0; iSeq < numSeqs; iSeq++)//Filas - Secuencias
		{
			int element = matrix[iCol][iSeq];
			if(element == 1)//Enzima
			{				
				line = sequences[iSeq]->getEnzyme(counts[iSeq]);
				counts[iSeq]++;
			}
			else//Gap
			{
				line = gap;				
			}
			aux.push_back(line);
		}
		//Aqui ya se tiene la columna traducida..
		//La enviamos para que se calcule la entropía..		
		entropy += homoPerColumn(aux);
	}
	//Promedio..
	entropy = entropy / (double)sizePrintMatrix;
	return entropy;
}

/*
	La entropia de cada columna se calcula como la entropia para cada nivel del ec number
	multiplicada por un factor E. Col = EN1 * 15 + EN2 * 10 + EN3 * 5
*/
double FitnessEntropy::homoPerColumn(COLUMN col)
{
	double entropy = 0.0;
	int total = numSeqs;		
	int gaps = 0;
	double penal = 0.0;
	for(int iCol = 0; iCol < 3; iCol++)
	{
		map<int,int> counts;
		map<int,int>::iterator it;
		double entropyAux = 0.0;
		//CAlcular l amaxima entropia..
                //double maxima = (((double)1)*log((double)1/ (double)numSeqs)) * numSeqs * -1;

		for(int iRow = 0; iRow < numSeqs; iRow++)
		{
			int symbol = col[iRow][iCol];
			counts[symbol] = counts[symbol] + 1;			
		}		
		//Eliminamos los gaps...no los tomamos en cuenta, o hay que ver
		//como penalizarlos para que no afecten..
		//240210 - Se toman en cuenta solo en el denominador..
		if (counts.find(0) != counts.end())
		{
			gaps = counts[0];
			//counts.erase(counts.find(0));
		}
		for(it = counts.begin(); it != counts.end(); it++)
		{
			//	entropyAux +=  ((double)it->second)*log((double)it->second / (double)total);
			entropyAux +=  ((double)it->second / (double)total)*log((double)it->second / (double)total);
			//entropyAux +=  ((double)it->second / (double)(total-gaps))*log((double)it->second / (double)(total-gaps));
		}
		entropyAux *= -1;
		if (gaps > 0)
		{	//Penalizar los gaps..
			double porcen = (double) gaps / (double) (total-1);
			penal = porcen * 1.0;//maxima;

		}
		//Hay que normalizarla antes...
        //int contador = (int)(counts.size());
        //double denominador = (double)contador;
        //if (entropyAux != 0.0 && denominador != 1)
	        entropyAux = entropyAux / log((double)total);
                //cout << "Entropy norm " << entropyAux << endl;
                /////////////////////////////
		//Ponderar segun la columna...
		entropyAux *=factors[iCol];
		//entropyAux += penal; 
		entropy += entropyAux; /// 2.0;// + penal;
		//cout << "Entropy col " << entropyAux << endl;
	}
	//Primero calculamos la entropia de la columna completa
	entropy = entropy / (double)(factors[0] + factors[1] + factors[2]);
	//y ahora penalizamos el porcetaje de gaps en la columna
	//Aqui hacemos el calculo para la penalización de gaps...El mayor numero de gaps que puede tener una columna
        //es total - 1 , eso equivale a un 1, hacemos una regla de tres para ver a cuanto equivale el numero de gaps
        //en una escala del 0 al 1
        double penaGaps = (double)gaps / (double)(total - 1); // Esta es la pena en un valor de 0 a 1, 1 siendo el peor y 0 el mejor

	//entropy += penal;
	//entropy = (entropy + penaGaps) / 2.0;
	entropy = (entropy*.6) + (.4*penaGaps);
	return entropy;
}

/*
	Función para calcular la homogeneidad por bloques
*/
double FitnessEntropy::homoBlocks(map<int, EC> matrix,int point, int blockSize)
{
	///////////////////
	//Factor para block
	double factor =  numSeqs*log( 1.0 / (double)numSeqs)*15
			+ numSeqs*log( 1.0 / (double)numSeqs)*10
			+ numSeqs*log( 1.0 / (double)numSeqs)*5 ;
	factor *= -1;
	//////////////////
	vector<int> counts(numSeqs,0);
	vector<int> gap(3,0);
	double entropy = 0.0;	
	for(int iCol = point; iCol < point+blockSize; iCol++)//Columnas
	{
		COLUMN aux;
		EC line;		
		for(int iSeq = 1; iSeq <= numSeqs; iSeq++)//Filas - Secuencias
		{
			int element = matrix[iSeq][iCol];
			if(element == 1)//Enzima
			{				
				line = sequences[iSeq - 1]->getEnzyme(counts[iSeq - 1]);
				counts[iSeq - 1]++;
			}
			else//Gap
			{
				line = gap;				
			}
			aux.push_back(line);
		}
		//Aqui ya se tiene la columna traducida..
		//La enviamos para que se calcule la entropía..				
		//entropy += (homoPerColumn(aux) - factor) *-1;
		entropy += (homoPerColumn(aux));
	}
	//Hay que cambiarla para maximizar por la ruleta..
	entropy = entropy / (double)blockSize;
	//entropy = (entropy - 1 ) * -1;
	return entropy;
}

/*
	La penalización por los gaps horizontales, se usa la fórmula de Edgar
*/
double FitnessEntropy::penalty(vector<EC> matrix)
{
	double CG = 0.0;
	int sizePrintMatrix = matrix.size();
	int totalGapsInd = 0;
	int totalBloquesGaps = 0;
	//Hay que eliminar antes las columnas que solo tengan gaps..
	for(int iSeq = 0; iSeq < numSeqs; iSeq++)//Columnas
	{
		string seq = "";
		for(int iCol = 0; iCol < sizePrintMatrix; iCol++)//Filas - Secuencias
		{
			seq += int_ToString(matrix[iCol][iSeq]);			
		}
		//Aqui ya tengo los 0's y 1's concatenados..
		//Hay que hacerles el trim
		string aux = trim(seq);
		//string aux = seq;
		//Se eliminan los gaps del inicio y el final..
		vector<string> gaps;
		boost::regex re("1+");
		boost::sregex_token_iterator i(aux.begin(), aux.end(), re, -1);
		boost::sregex_token_iterator j;
		unsigned count = 0;
		*i++;//El primer elemento esta en blanco..
		int suma = 0;
		while(i != j)
		{
			string a = *i++;
			gaps.push_back(a);
			suma += a.length();
			totalGapsInd += a.length();
			totalBloquesGaps++;
			count++;
		}
		
	}
	//Esta penalización tiene 1 si es un solo bloque y es menor si hay varios bloques,
	//tenemos que invertirla..
	double aux = 0.0;
	if (totalGapsInd > 0)
	{
		CG = (double)totalGapsInd / (double)totalBloquesGaps;
		CG = CG / (double)totalGapsInd;
		aux = (CG - 1) * -1;
	}
	return aux;
}

