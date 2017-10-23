#include "../h/FitnessSumaPares.h"

/*
	Clase que usa la función Suma de Pares para calificar la homogeneidad de las columnas
*/

FitnessSumaPares::FitnessSumaPares(void)
{
}

FitnessSumaPares::~FitnessSumaPares(void)
{
}

/*
	Estos son los factores que se usan en la suma de pares, y en la penalización de los gaps
*/
FitnessSumaPares::FitnessSumaPares(vector<Secuencia*> s,int mL) : TFitness(s,mL)
{
	factors.push_back(2);//match
	factors.push_back(-1);//mismatch
	factors.push_back(-2);//gap (simbolo,gap)
	openingGap = -2;
	extGap = -.125;

}

/*
	Esta es la función principal, es la que pondera los factores involucrados
*/
double FitnessSumaPares::calculate(map<int, EC> matrix, vector<EC> matrixWithoutGaps)
{
	
	double homo = homoColumns(matrixWithoutGaps);
	double penal =  penalty(matrixWithoutGaps);
	double incremento = incrementoColumnas(matrixWithoutGaps);
	double fitness = 0.0;
	fitness = homo+penal-incremento;
	return fitness;
}

/*
	La calificación es la suma de la entropia de todas las columnas
	Una menor calificación indica mas homogeneidad en la columna..
	Si toda la columna es homogenea la calificación es 0.
*/
double FitnessSumaPares::homoColumns(vector<EC> matrix)
{
	vector<int> counts(numSeqs,0);
	vector<int> gap(3,0);
	double suma = 0.0;
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
		suma += homoPerColumn(aux);
	}
	//Promedio..
	suma = suma;// / (double)sizePrintMatrix;
	return suma;
}

/*
	Se calcular los matches y mismatches por columna, el (-,X) es un gap, (-,-) no se califica
	(x,y) mismatch, (x,x) match
	Se hace primero para los 3 niveles, y asi para cada columna, no se divide, ni se saca promedios.
*/
double FitnessSumaPares::homoPerColumn(COLUMN col)
{
	double sumaPares = 0.0;
		
	for(int iCol = 0; iCol < 3; iCol++)
	{
		map<int,int> counts;
		map<int,int>::iterator it;
		double sumaParesAux = 0.0;

		for(int iRow = 0; iRow < numSeqs - 1; iRow++)
		{
			for(int iRow2 = iRow + 1; iRow2 < numSeqs; iRow2++)
			{
				int symbol1 = col[iRow][iCol];
				int symbol2 = col[iRow2][iCol];
				if ((symbol1 == 0 && symbol2 == 0) )
                                {       //gap
                                        sumaParesAux += 0;
                                }
				else
				{
					if ((symbol1 == 0 || symbol2 == 0) )
					{	//gap
						sumaParesAux += factors[2];
					}
					else
					{
						if ( symbol1 != symbol2 )
						{	//mismatch
							sumaParesAux += factors[1];
						}
						else
						{	//match
							sumaParesAux += factors[0];					
						}
					}
				}
				
			}	
		}
		
		//Ponderar segun la columna...
		//sumaParesAux *=factors[iCol];		
		sumaPares += sumaParesAux; 
		//cout << "Entropy col " << entropyAux << endl;
	}
	//Primero calculamos la entropia de la columna completa
	sumaPares = sumaPares;// / 3.0;
	return sumaPares;
}

/*
	PENDIENTE, para calcular por bloques
*/
double FitnessSumaPares::homoBlocks(map<int, EC> matrix,int point, int blockSize)
{
	return 0.0;
}


/*
	Penalización de los gaps, por inserción e incremento
*/
double FitnessSumaPares::penalty(vector<EC> matrix)
{
	double p = 0.0;
	int sizePrintMatrix = matrix.size();
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
		//Se eliminan los gaps del inicio y el final..
		vector<string> gaps;
		boost::regex re("1+");
		boost::sregex_token_iterator i(aux.begin(), aux.end(), re, -1);
		boost::sregex_token_iterator j;
		unsigned count = 0;
		*i++;//El primer elemento esta en blanco..
		while(i != j)
		{
			string a = *i++;
			gaps.push_back(a);
			count++;
		}
		vector<string>::iterator iGaps;		
		for(iGaps = gaps.begin(); iGaps != gaps.end(); iGaps++)
		{
			p += openingGap + ((*iGaps).length() - 1) * (extGap);		
		}
	}
	return p;
}
