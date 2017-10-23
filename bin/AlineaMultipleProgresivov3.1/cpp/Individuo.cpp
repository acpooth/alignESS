#include "../h/Individuo.h"
#include "../h/Secuencia.h"

/*
	Clase que representa a los individuos, cada individuo representa un posible alineamiento, usando una matriz. Las matrices
	tienen en cada fila, una secuencia, todas las secuencias de la matriz, deben de ser del mismo tamaño.
	Puede ser que haya individuos dentro de una población que tengan matrices de distintos tamaños, según lso operadores que se 
	esten usando en el AG
	Se usan 0's para representar gaps, y 1's para representar enzimas. Debemos de cuidar que el numero de 1's que aparece en cada
	fila sea igual al número de enzimas que hay en la secuencia.
	Esta clase tiene una modificación para poder manejar los AG's iterativos
	donde los individuos tienen cierta parte de la matriz fija.
*/

Individuo::Individuo(void)
{
}

Individuo::~Individuo(void)
{
}

/*
	Constructor, aqui se inicializa la matriz aleatoriamente.
*/
Individuo::Individuo(vector<int> sL, int mL)
{
	seqLengths = sL;
	maxLength = mL;
	numSeq = sL.size();
	initializeMatrix();
}

/*
	Este constructor es para los individuos que usan los AG iterativos
	Hay que mantener cierta parte de la matriz fija, y el resto es aleatorio.
*/
Individuo::Individuo(vector<EC> a, int sinAlinear, vector<int> sL, int mL)
{
	alineados = a;
	numSeqSinAlinear = sinAlinear;
	seqLengths = sL;
	maxLength = mL;
	numSeq = sL.size();
	initializeMatrixAlineados();
}

/*
	Esta función es auxiliar para aumentar la longitud de las secuencias, cuando se usan ciertos operadores en el AG
*/
void Individuo::setMaxLength(int m)
{
	maxLength = m;
}

/*
	Devuelve la longitud de las secuencias de la matriz en este individuo
*/
int Individuo::getMaxLength()
{
	return maxLength;
}

/*
	Inicializa la matriz aleatoriamente, sirve para crear la población inicial
*/
void Individuo::initializeMatrix()
{
	vector<int> aux;
	
	//Vamos a crear un vector para cada secuencia. Primero generamos un vector de posiciones 
	for(int i = 0; i < numSeq; i++)
	{
		int l = seqLengths[i];//Número de enzimas
		for(int j = 0; j < l; j++) // Primero agregamos 1's (Enzimas)
		{
			aux.push_back(1);				
		}
		for(int j = l; j < maxLength; j++) // Despues acompletamos con 0's (Gaps)
		{
			aux.push_back(0);				
		}
		//Ahora generamos el arreglo aleatorio..
		for (int k = maxLength - 2; k  >= 0 ; k--)
          {
               int random = rand() % maxLength;
               int element = aux[random];
               aux[random] = aux[k];
               aux[k] = element;
          }		
		alignMatrix[i + 1] = aux;
		aux.clear();	
	}
}

/*
	Inicializa la matriz aleatoriamente, sirve para crear la población inicial
	Esta función tiene cierta parte de la matriz fija, y el resto es el que
	se crea aleatorio.
	El vector alineados, contiene las secuencias que se van a mantener fijas en la matriz.	
*/
void Individuo::initializeMatrixAlineados()
{
	vector<int> aux;
	
	//Vamos a crear un vector para cada secuencia. Primero generamos un vector de posiciones 
	//Primero agregamos las secuencias  ya alineadas..
	int iAlineados = 0;
	for(iAlineados = 0; iAlineados < (int)alineados.size(); iAlineados++)
	{
		while((int)alineados[iAlineados].size() < maxLength )
			alineados[iAlineados].push_back(0);
		alignMatrix[iAlineados + 1] = alineados[iAlineados];		
	}
	
	int l = seqLengths[numSeqSinAlinear - 1];//Número de enzimas
	for(int j = 0; j < l; j++) // Primero agregamos 1's (Enzimas)
	{
		aux.push_back(1);				
	}
	for(int j = l; j < maxLength; j++) // Despues acompletamos con 0's (Gaps)
	{
		aux.push_back(0);				
	}
	//Ahora generamos el arreglo aleatorio..
	for (int k = maxLength - 2; k  >= 0 ; k--)
      {
           int random = rand() % maxLength;
           int element = aux[random];
           aux[random] = aux[k];
           aux[k] = element;
      }		
	alignMatrix[iAlineados + 1] = aux;
	
}

/*
	Sirve para llamar a la función adecuada que calcula el fitness del individuo, se envia, la matriz de unos y ceros..
	y la matriz que ya esta convertida a ec numbers.
	La función printAux, crea esa última matriz. 
	La función convertirMatrix, es opcional, lo que hace es eliminar las columnas que solamente tiene gaps, para que 
	los operadores se apliquen sobre la matriz ya depurada de estas columnas.
*/
double Individuo::calculateFitness(TFitness* f,vector<Secuencia*> seqs)
{
	printAux();
	fitness = f->calculate(alignMatrix, printMatrix);
	convertirMatrix(printMatrix);
	return fitness;
}

/*
	Crea la matriz de 0's y 1's sin columnas que tengan puros ceros.
	Se copia la matriz printMatrix, que ya no tiene las columnas con gaps, a la matriz que sirve para aplicar los operadores
	alignMatrix.
*/
void Individuo::convertirMatrix(vector<EC> pM)
{
	int sizeMatrix = pM.size();
	alignMatrix.clear();
	vector<int> aux;
	for(int iCol = 0; iCol < numSeq; iCol++)
	{
		aux.clear();
		for(int iFila = 0; iFila < sizeMatrix; iFila++)
		{
			aux.push_back(pM[iFila][iCol]);
		}
		alignMatrix[iCol + 1] = aux;
	}	
}
/*
	Regresa el fitness del individuo
*/
double Individuo::getFitness()
{
	return fitness;
}

/*
	Fitness relativo del individuo
*/
double Individuo::getRelFitness(double totalFitness)
{
	return (double)fitness / (double) totalFitness;
}

/*
	Esta función, elimina las columnas que solo tienen gaps, de la matriz alignMatrix, y 
	esta nueva matriz se guarda en printMatrix.
*/
void Individuo::printAux()
{ //En la matriz para imprimir las  columnas de la alignMatrix se almacenan en las filas
	printMatrix.clear();
	maxLength = alignMatrix[1].size();
	for(int iCol = 0; iCol < maxLength; iCol++)//Columnas
	{
		int gaps = 0;
		vector<int> column;
		for(int iSeq = 1; iSeq <= numSeq; iSeq++)
		{
			gaps += alignMatrix[iSeq][iCol];
			column.push_back(alignMatrix[iSeq][iCol]);
		}
		if (gaps != 0)//Se eliminan las columnas que solo tienen gaps..
		{
			printMatrix.push_back(column);
		}
	}
}

/*
	Imprime a pantalla el individuo convertido a enzimas
*/
void Individuo::print(vector<Secuencia*> seqs)
{	
	string s = "";
	printSecuencias.clear();
	vector<int> counts(numSeq,0);
	printAux();
	int sizePrintMatrix = printMatrix.size();
	cout << seqs[0]->getName() << endl;//Imprime el nombre del cluster
	for(int iSeq = 0; iSeq < numSeq; iSeq++)
	{
		string line = "";	
		string salida = seqs[iSeq]->getSeqId();//El id de la secuencia
		salida.append("\t");
		cout << seqs[iSeq]->getSeqId() << "\t" << seqs[iSeq]->getMapId() << "\t"; //El nombre del mapa
		for(int iCol = 0; iCol < sizePrintMatrix; iCol++)//Columnas
		{			
			int element = printMatrix[iCol][iSeq];
			if(element == 1)//Enzima
			{				
				s = seqs[iSeq]->getEnzymeStr(counts[iSeq]);
				counts[iSeq]++;
			}
			else//Gap
			{
				s = "-.-.-";
			}
			line += s;
			line += ":";
		}
		cout << line.substr(0,line.length() - 1) << endl;
		printSecuencias.push_back(salida.append(line.substr(0,line.length() - 1)));
	}
	cout << "Fitness: " << fitness << endl;
}

/*
	Esta función regresa las secuencias de salida, las que se pintan
	en pantalla, sirve para el AG iterativo, de ahi las toma para agregarlas como siguiente
	entrada.
*/
vector<string> Individuo::getPrintSecuencias(vector<Secuencia*> seqs)
{	
	string s = "";
	printSecuencias.clear();
	vector<int> counts(numSeq,0);
	printAux();
	int sizePrintMatrix = printMatrix.size();
	//cout << seqs[0]->getName() << endl;//Imprime el nombre del cluster
	for(int iSeq = 0; iSeq < numSeq; iSeq++)
	{
		string line = "";	
		string salida = seqs[iSeq]->getSeqId();//El id de la secuencia
		salida.append("\t");
		//cout << seqs[iSeq]->getSeqId() << "\t" << seqs[iSeq]->getMapId() << "\t"; //El nombre del mapa
		for(int iCol = 0; iCol < sizePrintMatrix; iCol++)//Columnas
		{			
			int element = printMatrix[iCol][iSeq];
			if(element == 1)//Enzima
			{				
				s = seqs[iSeq]->getEnzymeStr(counts[iSeq]);
				counts[iSeq]++;
			}
			else//Gap
			{
				s = "-.-.-";
			}
			line += s;
			line += ":";
		}
		//cout << line.substr(0,line.length() - 1) << endl;
		printSecuencias.push_back(salida.append(line.substr(0,line.length() - 1)));
	}
	return printSecuencias;
}

/*
	Se usa para calcular el fitness por bloques
	PENDIENTE!!	
*/
double Individuo::calculateFitnessPerBlock(TFitness* f,vector<Secuencia*> seqs,int p,int bz)
{
	double fitnessBlock = f->calculateBlocks(alignMatrix, p - bz, bz);
	return fitnessBlock;
}
