#include "../h/FitnessInfo.h"

/*
	Aqui se usa la formula de la tesis de Edgar
	También se maximiza..
*/
FitnessInfo::FitnessInfo(void)
{
}

FitnessInfo::~FitnessInfo(void)
{
}

/*
	Se inicializan los factores para ponderar los niveles de los ec numbers
*/
FitnessInfo::FitnessInfo(vector<Secuencia*> s,int mL) : TFitness(s,mL)
{
	factors.push_back(1);
	factors.push_back(1);
	factors.push_back(1);
}

double FitnessInfo::calculate(map<int, EC> matrix, vector<EC> matrixWithoutGaps)
{

        double homo = homoColumns(matrixWithoutGaps);
        double penal =  penalty(matrixWithoutGaps);
        double incremento = incrementoColumnas(matrixWithoutGaps);
        incremento = (incremento);
        double fitness = 0.0;
        fitness = homo * (penal);
        fitness = homo - (homo * penal);
        fitness = fitness - (fitness* (1.0 + incremento));
        return fitness;
}

/*
	Para calcular la homogeneidad por columna
*/
double FitnessInfo::homoPerColumn(COLUMN col)
{
	double info = 0.0;
	
	for(int iCol = 0; iCol < 3; iCol++)
	{
		map<int,int> counts;
		map<int,int>::iterator it;
		double infoAux = 0.0;

		for(int iRow = 0; iRow < numSeqs; iRow++)
		{
			int symbol = col[iRow][iCol];
			counts[symbol] = counts[symbol] + 1;			
		}		
		//int total = numSeqs;		
		//Eliminamos los gaps...no los tomamos en cuenta, o hay que ver
		//como penalizarlos para que no afecten..
		int gaps = 0;
		if (counts.find(0) != counts.end())
		{
			gaps = counts[0];
			counts.erase(counts.find(0));
		}
			
		//////////////////////////////
		double sumBelow = 0.0;
		for(it = counts.begin(); it != counts.end(); it++)
		{
			sumBelow += ((double)it->second);
			infoAux +=  ((double)it->second)*((double)it->second);
		}
		sumBelow += gaps;
		infoAux = ((double)infoAux / (double)(sumBelow * sumBelow));//multiplicamos por -1 pq vamos a minimizar
		//Ponderar segun la columna...
		//YA no se pondera el nivel del ec number..infoAux *=factors[iCol];
		info += infoAux;// + (iCol * 5);
	}	
	return info;
}

/*
	Es la que va enviando las columnas para que se calcule el fitness
*/
double FitnessInfo::homoColumns(vector<EC> matrix)
{
	vector<int> counts(numSeqs,0);
	vector<int> gap(3,0);
	double info = 0.0;
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
		//La enviamos para que se calcule la informacion..				
		info += homoPerColumn(aux);
	}
	return info / (double)sizePrintMatrix; //Número de columnas..
}

/*
	Se usa para calcular el fitness por bloques
	PENDIENTE!!
*/
double FitnessInfo::homoBlocks(map<int, EC> matrix,int point, int blockSize)
{
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
		entropy += homoPerColumn(aux) * -1; //Para regresar al positivo...mayor homogeneidad mayor fitness..
	}
	return entropy;
}

/*
	La penalización por los gaps horizontales, se usa la fórmula de Edgar
*/
double FitnessInfo::penalty(vector<EC> matrix)
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

