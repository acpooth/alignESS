#include <fstream>
#include "../h/AG.h"
#include "../h/Selector.h"
#include "../h/FitnessSumaPares.h"
#include "../h/FitnessEntropy.h"
#include "../h/FitnessInfo.h"


/*
Constructor default
*/
AG::AG(void)
{
}

/*
Destructor default
*/
AG::~AG(void)
{
}

/*
Constructor default
*/
AG::AG(int mG, double cR, double mR, Poblacion p)
{
	numSeqSinAlinear = 1;
	crossRate = cR;
	mutRate = mR;
	maxGen = mG;
	population = p;
	individuals = p.individuals;
	sizePop = individuals.size();
	newMaxLength = p.getMaxLength();
}

/*
	Constructor para crear el alineamiento iterativo
	numseqSinalinear es cuantas secuencias ya se alinearon actualmente
*/
AG::AG(int mG, double cR, double mR, Poblacion p,int nSSA)
{
	numSeqSinAlinear = nSSA; //Contando de 1 hasta n
	crossRate = cR;
	mutRate = mR;
	maxGen = mG;
	population = p;
	individuals = p.individuals;
	sizePop = individuals.size();
	newMaxLength = p.getMaxLength();
}

/*
	Devuelve la tasa de cruza
*/
double AG::getCrossRate()
{
	return crossRate;
}

/*
	Devuelve la tasa de mutación
*/
double AG::getMutRate()
{
	return mutRate;
}

/*
	Regresa el maximo de iteraciones sea el maximo de iteraciones sin cambio
*/
int AG::getMaxGen()
{
	return maxGen;
}

/*
	Devuelve el tamaño de la población
*/
int AG::getSizePop()
{
	return sizePop;
}

/*
	Devuelve un vector de individuos
*/
vector<Individuo> AG::getPopulation()
{
	return individuals;
}

////////////////////////////////////////////////////
// Funciones de Selección
//////////////////////////////////////////////////

/*
	Selección con ruleta de varios marcadores
*/
vector<Individuo> AG::sus(vector<Individuo> oldPop, int size,Individuo best)
{
	int   i   = size - 1; //Número de individuos a seleccionar
	size--;
	double sum = 0; //Suma de los fitness relativos
	double mark, inc; //posición del marcador e incremento
	vector<Individuo> newPop(size);
	inc  = 1.0/size;
	mark = inc *(double)rand () / (double) RAND_MAX; //posición del primer marcador
	//////////////Scalar los fitnesses//////////
	vector<double> fitScalado(size);
	double totalFitScalado = 0.0;
	for(int iFit = 0; iFit < size; iFit++)
	{
		fitScalado[iFit] = 10000 - 10*(oldPop[iFit].getFitness());
		totalFitScalado += fitScalado[iFit];
	}
	/////////////////////////////////
	while ((--size >= 0) && (i > 0))
	{
		sum += fitScalado[size] / totalFitScalado;
		while (mark < sum)
		{
			newPop[--i] = oldPop[size];
			mark += inc;
		}
	} //copiar el individuo 
	while (--i >= 0) //Si faltan individuos se llena con el primero..     
		newPop[i] = oldPop[0];
	//Agregamos el mejor individuo
	newPop.push_back(best);
	return newPop;
}

/*
	Selección por torneo
	Cuando el fitness es SumaPares se maximiza, en cualquier otro caso es minimizar..
*/
Individuo AG::tournament(vector<Individuo> oldPop,int tmsize)
{
	Individuo ind;
	Individuo best;  
	int size = oldPop.size();
	best = oldPop[(int)( rand() % (size) )];
	while (--tmsize > 0)
	{ 
		ind = oldPop[(int)(rand() % (size))];
		if (dynamic_cast<FitnessEntropy*>(population.getTipoFitness()) == NULL)
		{
			if (ind.getFitness() > best.getFitness())
				best = ind;
		}
		else
		{
			if (ind.getFitness() < best.getFitness())
				best = ind;

		}
	}

	return best;  
}

/*
	Función principal para la seleción..
*/
vector<Individuo> AG::selection(vector<Individuo> oldPop,int tmsize,int elitist,Individuo best)
{
	int i = sizePop;	
	vector<Individuo> newPop;
	
	if (elitist)
	{
		newPop.push_back(best);
		i--;
	}
	while (--i >= 0)
	{   
		newPop.push_back(tournament(oldPop,tmsize));
	}
	return newPop;
}


/////////////////////////////////////////////////
//	Función de operadores de Mutación...
/////////////////////////////////////////////////

/*
	Este operador se aplica sobre todas las secuencias de la matriz, sobre cada bit del individuo
	Si muta, la operación que se realiza consiste, en tomar 2 bits de la secuencia, e intercambiarlos.
	Los bits deben de ser diferentes, osea, un 1 y un 0.
	El cambio principal esta en el for de las secuencias!
*/
void AG::mutationPermutate(vector<Individuo> oldPop)
{
	int max = population.getMaxLength();
	int nSeqs = population.getNumSeqs();
	//Se lanza el volado por cada bit de cada individuo..
	for(int j = 0 ; j < sizePop; j++)
    {		
		Individuo ind = oldPop[j];		
		for(int iSeq = numSeqSinAlinear ; iSeq <= nSeqs; iSeq++)//Se empieza a contar en la ultima secuencia
		{
			vector<int> sequence = ind.alignMatrix[iSeq];
			max = sequence.size();
			for(int iEnz = 1 ; iEnz <= max; iEnz++)//Hay que acceder al mapa
			{
				double rB = rand() / ((double) RAND_MAX);
				if (rB < mutRate) //Se realiza la mutacion 
				{ //Estos son los bits del individuo..
				  //Hay que elegir primero una secuencia, y luego 2 posiciones dentro de esa secuencia
				  //para permutarlas..
					//Podemos elegir la secuencia tambien o trabajar con la actual..int seq = rand() % nSeqs;
					int rand1 = rand() % max;
					int rand2 = rand() % max;
					while( (rand1 == rand2) || (sequence[rand1] == 0 && sequence[rand2] == 0 )
								|| (sequence[rand1] == 1 && sequence[rand2] == 1) ) //Solamente permutadso 0 y 1
					{ 
						rand2 = rand() % max;
						rand1 = rand() % max;
					}
					int aux = sequence[rand1];
					sequence[rand1] = sequence[rand2];
					sequence[rand2] = aux;
				}
			}
			ind.alignMatrix[iSeq] = sequence; //Actualizamos el individuo..
		}
	}
}

/*
	CHECAR!!! PENDIENTE!!!
	Este operador cambio el tamaño del individuo, osea, aumenta la longitud de la secuencia..
	Se aplica sobre secuencia, es decir, recorremos todas las secuencias de todos los individuos y lanzamos 
	la moneda para ver si se realiza la mutación o no.	
	Este trata de insertar el gap como un bloque, osea, en todas(num aleatorio) las secuencias se inserta el mismo gap
*/
vector<Individuo> AG::mutationGaps2(vector<Individuo> oldPop)
{
	vector<Individuo> newPop;
	vector<Individuo>::iterator it;
	Individuo parent1;
	Individuo parent2;
	int n = getSizePop();
	int numSeqs = population.getNumSeqs();
	int maxLength = population.getMaxLength(); 
	vector<int> lengths = population.getSeqLengths();
	//IMPORTANTE!!! No hay que cambiar el orden de las secuencias, o si se cambian hay que tener
	//en cuenta que el key es consecutivo nos sirve para identificar que secuencias es!!
	int MAYOR = newMaxLength;
	for(int index = 0; index < n; index++)//Individuos
	{
		for(int indSeq = numSeqSinAlinear; indSeq < numSeqs; indSeq++)
		{
			parent1 = oldPop[index]; //El primero		
			map<int, EC> matrix1 = parent1.alignMatrix;
			vector<int> seq1;
			double coin = (double)rand () / (double) RAND_MAX;
			if (coin < mutRate)
			{					
				//Elegimos un numero al azar entre 1 y el número total de columnas..
				int numSeqAMutar = rand() % numSeqs;
				int random = rand() % matrix1[1].size();
				///////////////////////////
				vector<int> auxPermuta;
				for(int iAux = 1; iAux <= numSeqs; iAux++)
				{
					auxPermuta.push_back(iAux);
				}	
				for (int iAux = numSeqs -1 ; iAux  >=0 ; iAux--)
				{
					int aleatorio = rand() % numSeqs;
					int aux = auxPermuta[aleatorio];
					auxPermuta[aleatorio] = auxPermuta[iAux];
					auxPermuta[iAux] = aux;
				}
				//Hacemos mas chico el arreglo...
				auxPermuta.resize(numSeqAMutar);//Aqui estan los indices de las secuencias que van a mutar
				coin = (double)rand () / (double) RAND_MAX;
				for(int iSeqMut = 0; iSeqMut < numSeqAMutar; iSeqMut++)
				{
					int nSeqMuta = auxPermuta[iSeqMut];
					seq1 = matrix1[nSeqMuta];
					vector<int>::iterator it = seq1.begin();
					if (seq1[random] == 0)
					{	//Lanzamos una moneda, si es menor q .5 se aumenta, si es menor se elimina..
						if (coin >= .5) //Aumenta el gap en una unidad
						{							
							//Se inserta desp del random
							seq1.insert(it+random+1,2,0);
						}
						else //Disminuye el gap o se elimina..
						{
							seq1.erase(it+random);	
						}
					}
					else  // Es 1
					{
						//Se inserta un gap desp del simbolo..
						seq1.insert(it+random+1,2,0);
					}
					if ((int)seq1.size() > MAYOR)
					{
						MAYOR = (int)seq1.size();
						maxLength = MAYOR;
					}
					matrix1[nSeqMuta] = seq1;
				}
			}//IF si se muta o no
			parent1.alignMatrix = matrix1;
			newPop.push_back(parent1);		
		}//FOR secuencias
	}//FOR individuos
	//Ahora hay que ver cual es la secuencia mas grande y acompletar las demas..		
	newPop = ajustarMaxLength(newPop, MAYOR);
	newMaxLength = MAYOR;
	return newPop;
}

/*
	Este se va a iterar sobre cada bit del individuo, para usar una probabilidad baja..
	En este se mutan todas las secuencias, con un gap a la vez de tamaño fijo
	El tamaño del vector tambien aumenta
*/
vector<Individuo> AG::mutationGaps(vector<Individuo> oldPop)
{
	vector<Individuo> newPop;
	vector<Individuo>::iterator it;
	Individuo parent1;
	Individuo parent2;
	int n = getSizePop();
	int numSeqs = population.getNumSeqs();
	int maxLength = population.getMaxLength(); 
	vector<int> lengths = population.getSeqLengths();
	//IMPORTANTE!!! No hay que cambiar el orden de las secuencias, o si se cambian hay que tener
	//en cuenta que el key es consecutivo nos sirve para identificar que secuencias es!!
	int MAYOR = newMaxLength;
	for(int index = 0; index < n; index++)
	{
		parent1 = oldPop[index]; //El primero		
		map<int, EC> matrix1 = parent1.alignMatrix;
		vector<int> seq1;
		for(int iSeq = numSeqSinAlinear; iSeq <= numSeqs; iSeq++)
		{
			int lon = matrix1[iSeq].size();
			for(int iLonSeq = 0; iLonSeq < lon; iLonSeq++)
			{
				seq1 = matrix1[iSeq];
				double coin = (double)rand () / (double) RAND_MAX;
				if (coin < mutRate)
				{					
					//Elegimos un numero al azar entre 1 y el número total de columnas..
					int random = rand() % seq1.size();
					vector<int>::iterator it = seq1.begin();
					if (seq1[random] == 0)
					{	//Lanzamos una moneda, si es menor q .5 se aumenta, si es menor se elimina..
						coin = (double)rand () / (double) RAND_MAX;
						if (coin >= .5) //Aumenta el gap en una unidad
						{							
							//Se inserta desp del random
							seq1.insert(it+random+1,3,0);
						}
						else //Disminuye el gap o se elimina..
						{
							seq1.erase(it+random);	
						}
					}
					else  // Es 1
					{
						//Se inserta un gap desp del simbolo..
						seq1.insert(it+random+1,3,0);
					}
					if ((int)seq1.size() > MAYOR)
					{
						MAYOR = (int)seq1.size();
						maxLength = MAYOR;
					}
					
				}				
			
				matrix1[iSeq] = seq1;
			}//FOR enzimas
		} //FOR secuencias
		parent1.alignMatrix = matrix1;
		newPop.push_back(parent1);		
	}//FOR individuos
	//Ahora hay que ver cual es la secuencia mas grande y acompletar las demas..		
	newPop = ajustarMaxLength(newPop, MAYOR);
	newMaxLength = MAYOR;
	return newPop;
}

/*
	Aqui se ajusta la longitud de las secuencias de cada individuo, al tamaño maximo, acompletando con 0's
	Los individuos pueden tener distintas longitudes, pero para cada individuo, las matrices deben de ser
	de las mismas longitudes despues de aplicar esta función.
*/
vector<Individuo> AG::ajustarMaxLength(vector<Individuo> oldPop,int max)
{
	vector<Individuo> newPop;
	int n = getSizePop();
	int numSeqs = population.getNumSeqs();
	//IMPORTANTE!!! No hay que cambiar el orden de las secuencias, o si se cambian hay que tener
	//en cuenta que el key es consecutivo nos sirve para identificar que secuencias es!!
	for(int index = 0; index < n; index++)
	{
		Individuo ind = oldPop[index]; //El primero		
		map<int, EC> matrix = ind.alignMatrix;
		map<int, EC> newMatrix;
		vector<int> seq1;
		for(int iSeq = 1; iSeq <= numSeqs; iSeq++)
		{
			seq1 = matrix[iSeq];
			int sizeActual = seq1.size();
			//Calculamos cuanto le falta..
			int falta = max - sizeActual;
			for(int i = 0; i < falta; i++)
				seq1.push_back(0);
			newMatrix[iSeq] = seq1;
		}
		ind.alignMatrix = newMatrix;
		ind.setMaxLength(max);
		newPop.push_back(ind);	
	}	
	return newPop;
}

/*
	Se aplica por cada bit del individuo.
	Lo que hace este operador, es, después de saber si el bit va a mutar  o no,
	es hacer un corrimiento hacia la derecha de todos los bits de la secuencia actual.
	Aqui no se cambia la longitud de los individuos, pero hay que cuidar que cuando lleguen aqui
	puedan ser de diferentes longitudes
*/
vector<Individuo> AG::mutationCircular(vector<Individuo> oldPop)
{
	vector<Individuo> newPop;
	vector<Individuo>::iterator it;
	Individuo parent1;
	Individuo parent2;
	int n = getSizePop();
	int numSeqs = population.getNumSeqs();
	int maxLength = population.getMaxLength();
	vector<int> lengths = population.getSeqLengths();
	//IMPORTANTE!!! No hay que cambiar el orden de las secuencias, o si se cambian hay que tener
	//en cuenta que el key es consecutivo nos sirve para identificar que secuencias es!!
	for(int index = 0; index < n; index++)
	{
		parent1 = oldPop[index]; //El primero		
		map<int, EC> matrix1 = parent1.alignMatrix;
		for(int iSeq = numSeqSinAlinear; iSeq <= numSeqs; iSeq++)
		{
			//Aqui hay que ver cual es la longitud de la secuencia..
			maxLength = matrix1[iSeq].size();
			for(int iEnz = 0; iEnz < maxLength; iEnz++)
			{
				double coin = (double)rand () / (double) RAND_MAX;
				if (coin < mutRate){					
					vector<int> seq1 = matrix1[iSeq];
					//Creamos el array q sirve para inicializar el valarray..
					int* init = NULL;
					init = new int[maxLength];				
					for(int iAux = 0; iAux < maxLength; iAux++)
					{
						init[iAux] = seq1[iAux];
					}
					valarray<int> circular (init,maxLength);				
					circular = circular.cshift(iEnz); 				
					seq1.clear();
					for(int iAux = 0; iAux < maxLength; iAux++)
					{
						seq1.push_back(circular[iAux]);
					}
					matrix1[iSeq] = seq1;
					delete [] init;
				}//IF mutación
			}//FOR enzimas
		}//FOR secuencias
		parent1.alignMatrix = matrix1;
		newPop.push_back(parent1);		
	}//FOR individuos
	return newPop;
}

//////////////////////////////////////////////////////////////
//	Funciones para el operador de recombinación
//////////////////////////////////////////////////////////////

/*
        Intercambiar secuencias...
        Se toman las segundas mitades (desp del random) y se intercambian
        No se cambia la longitud de las secuencias..
	No aplica para el iterativo, porque solo hay una secuencia!
*/
vector<Individuo> AG::crossoverExchange(vector<Individuo> oldPop)
{
	vector<Individuo> newPop;
	vector<Individuo>::iterator it;
	Individuo parent1;
	Individuo parent2;
	int n = (int)(getSizePop()/2);
	int numSeqs = population.getNumSeqs();
	//IMPORTANTE!!! No hay que cambiar el orden de las secuencias, o si se cambian hay que tener
	//en cuenta que el key es consecutivo nos sirve para identificar que secuencias es!!
	for(int index = 0; index < n; index++)
	{
		parent1 = oldPop[index]; //El primero
		parent2 = oldPop[sizePop-1-index]; //El ultimo..
		map<int, EC> matrix1 = parent1.alignMatrix;
		map<int, EC> matrix2 = parent2.alignMatrix;		
		double coin = (double)rand () / (double) RAND_MAX;
		if (coin < crossRate){
			//Elegimos un numero al azar entre 1 y el número total de secuencias..
			int random = rand() % numSeqs;
			for(int iParent1 = random; iParent1 < numSeqs; iParent1++)
			{
				vector<int> aux = matrix2[iParent1 + 1];
				matrix2[iParent1 + 1] = matrix1[iParent1 + 1];
				matrix1[iParent1 + 1] = aux;
			}
			parent1.alignMatrix = matrix1;
			parent2.alignMatrix = matrix2;
		} //Se hace la cruza...
		newPop.push_back(parent1);
		newPop.push_back(parent2);
	}
	return newPop;
}

/*
        Crossover columna fija modificado
        Se elije una columna, y los hijos conservan la secuencia hasta ese columna de los padres (primera mitad)        
        la segunda mitad, se acompleta con el siguiente padre,contando los 0's y unos que hagan falta..
*/
vector<Individuo> AG::crossoverColumnas2(vector<Individuo> oldPop)
{
	vector<Individuo> newPop;
	vector<Individuo>::iterator it;
	Individuo parent1;
	Individuo parent2;
	int n = (int)(getSizePop()/2);
	int numSeqs = population.getNumSeqs();
	vector<int> lengths = population.getSeqLengths();
	//IMPORTANTE!!! No hay que cambiar el orden de las secuencias, o si se cambian hay que tener
	//en cuenta que el key es consecutivo nos sirve para identificar que secuencias es!!
	for(int index = 0; index < n; index++)
	{
		parent1 = oldPop[index]; //El primero
		parent2 = oldPop[sizePop-1-index]; //El ultimo..
		map<int, EC> matrix1 = parent1.alignMatrix;
		map<int, EC> matrix2 = parent2.alignMatrix;		
		double coin = (double)rand () / (double) RAND_MAX;
		if (coin < crossRate)
		{
			//Aqui elijo la columna al azar...es la misma para todas las secuencias..
			int mL1 = matrix1[1].size();
			int mL2 = matrix2[1].size();
			int random1 = (rand() % mL1) + 1; //genera numeros del 1 al maxLength
			int random2 = (rand() % mL2) + 1; //genera numeros del 1 al maxLength
			for(int iSeq = numSeqSinAlinear; iSeq <= numSeqs; iSeq++)
			{				
				vector<int> seq1 = matrix1[iSeq];
				vector<int> seq2 = matrix2[iSeq];
				//Creamos el array q sirve para inicializar el valarray..
				int* init = NULL;
				int* init2 = NULL;
				init = new int[mL1];
				init2 = new int[mL2];
				for(int iAux = 0; iAux < mL1; iAux++)
				{
					init[iAux] = seq1[iAux];
				}
				for(int iAux = 0; iAux < mL2; iAux++)
				{
					init2[iAux] = seq2[iAux];
				}
					
				//Inicializamos los arreglos circulares..
				valarray<int> circular1 (init,mL1);
				valarray<int> circular2 (init2,mL2);
				valarray<int> aux;
				//Prueba a) Conservando sola la primera mitad de los padres... ESTA FUNCIONA MEJOR
				matrix1[iSeq] = crossoverColumnas2Child(circular1,circular2,lengths[iSeq - 1],mL1,mL2,random1);			
				matrix2[iSeq] = crossoverColumnas2Child(circular2,circular1, lengths[iSeq - 1], mL2,mL1,random2);
				//Prueba b) Tratando de mezclar las 2 mitades de los padres..
				//matrix1[iSeq] = crossoverColumnas3Child(circular1,circular2,lengths[iSeq - 1],maxLength,random);			
				//matrix2[iSeq] = crossoverColumnas3Child(circular2,circular1, lengths[iSeq - 1], maxLength,random);
				delete [] init;
			}//FOR secuencias
			parent1.alignMatrix = matrix1;
			parent2.alignMatrix = matrix2;
		}//IF cruza
		newPop.push_back(parent1);
		newPop.push_back(parent2);
	}//FOR individuos
	return newPop;
}

/*
        Esta es la función auxiliar del crossovercolumnas2, es donde se acompleta la segunda mitad del hijo     
*/
vector<int> AG::crossoverColumnas2Child(valarray<int> circular,valarray<int> circular2,int maxUnos,int maxLength,int maxLength2,int random)
{
	valarray<int> aux;
	vector<int> hijo1;
	//Contadores..
	int maxCeros = maxLength - maxUnos;
	//Contamos cuantos unos hay en la parte que ya cortamos..
	int contaUnosP1 = 0;
	
	if (random == maxLength)
	{
		for(int i = 0; i < (int)circular.size(); i++)
			hijo1.push_back(circular[i]);
	}
	else
	{		
		int iP1 = 0;
		//Ahora empezamos a contar los unos para meterlos al h1
		iP1 = 0;
		while(iP1 < random)
		{
			hijo1.push_back(circular[iP1]);
			if (circular[iP1] == 1)
				contaUnosP1++;
			iP1++;
		}
		//Checamos la longitud actual de H1
		int longActual = hijo1.size();
		int cerosActuales = longActual - contaUnosP1;
		int unosActuales = contaUnosP1;
		//Ahora hay que acompletar la segunda parte del hijo, tomar los elementos del P2

		int iP2 = 0;
		while(cerosActuales < maxCeros || unosActuales < maxUnos)
		{
			if (circular2[iP2] == 0 && cerosActuales < maxCeros)
			{
				hijo1.push_back(0);
				cerosActuales++;
			}
			if (circular2[iP2] == 1 && unosActuales < maxUnos)
			{
				hijo1.push_back(1);
				unosActuales++;
			}
			iP2++;
		}
	}
	return hijo1;
}

/*
        Esta es una prueba para acompletar la segunda mitad del hijo, para el crossovercolumnas2
*/
vector<int> AG::crossoverColumnas3Child(valarray<int> circular,valarray<int> circular2,int maxUnos,int maxLength,int random)
{
	valarray<int> aux;
	vector<int> hijo1;
	//Hacemos un shift en el padre2, para empezar el "acompletamiento" con la segunda mitad del segundo padre..
	//ya que la primera mitad se toma en el segundo hijo..
	circular2 = circular2.cshift(random);
	/////////////////////////////////////////////////////////
	//Contadores..
	int maxCeros = maxLength - maxUnos;
	//Contamos cuantos unos hay en la parte que ya cortamos..
	int contaUnosP1 = 0;
	
	if (random == maxLength)
	{
		for(int i = 0; i < (int)circular.size(); i++)
			hijo1.push_back(circular[i]);
	}
	else
	{		
		int iP1 = 0;
		//Ahora empezamos a contar los unos para meterlos al h1
		iP1 = 0;
		while(iP1 < random)
		{
			hijo1.push_back(circular[iP1]);
			if (circular[iP1] == 1)
				contaUnosP1++;
			iP1++;
		}
		//Checamos la longitud actual de H1
		int longActual = hijo1.size();
		int cerosActuales = longActual - contaUnosP1;
		int unosActuales = contaUnosP1;
		//Ahora hay que acompletar la segunda parte del hijo, tomar los elementos del P2

		int iP2 = 0;
		while(cerosActuales < maxCeros || unosActuales < maxUnos)
		{
			if (circular2[iP2] == 0 && cerosActuales < maxCeros)
			{
				hijo1.push_back(0);
				cerosActuales++;
			}
			if (circular2[iP2] == 1 && unosActuales < maxUnos)
			{
				hijo1.push_back(1);
				unosActuales++;
			}
			iP2++;
		}
	}
	
	return hijo1;
}

/*
	Primer intento de crossover por columnas, cruzando primero, y despues realizando un ajuste en las secuencias
	para conservar el numero de unos, NO FUNCIONO BIEN!
*/
vector<Individuo> AG::crossoverColumns(vector<Individuo> oldPop)
{
	vector<Individuo> newPop;
	vector<Individuo>::iterator it;
	Individuo parent1;
	Individuo parent2;
	int n = (int)(getSizePop()/2);
	int numSeqs = population.getNumSeqs();
	vector<int> seqLengths = population.getSeqLengths();
	int maxLength = population.getMaxLength();
	//IMPORTANTE!!! No hay que cambiar el orden de las secuencias, o si se cambian hay que tener
	//en cuenta que el key es consecutivo nos sirve para identificar que secuencias es!!
	for(int index = 0; index < n; index++)
	{
		parent1 = oldPop[index]; //El primero
		parent2 = oldPop[sizePop-1-index]; //El ultimo..
		map<int, EC> matrix1 = parent1.alignMatrix;
		map<int, EC> matrix2 = parent2.alignMatrix;		
		map<int, EC> childMatrix1;
		map<int, EC> childMatrix2;
		double coin = (double)rand () / (double) RAND_MAX;
		if (coin < crossRate)
		{
			//Elegimos un numero al azar entre 1 y el número total de columnas..
			int random = rand() % maxLength;
			for(int iParent1 = 1; iParent1 <= numSeqs; iParent1++)
			{
				vector<int> child1(maxLength,0);
				vector<int> child2(maxLength,0);
				for(int iCol = 0; iCol < random; iCol++)
				{
					child1[iCol] = matrix1[iParent1][iCol];
					child2[iCol] = matrix2[iParent1][iCol];
				}
				for(int iCol = random; iCol < maxLength; iCol++)
				{
					child1[iCol] = matrix2[iParent1][iCol];
					child2[iCol] = matrix1[iParent1][iCol];
				}
				checkSequences(child1,seqLengths[iParent1-1],maxLength);
				checkSequences(child2,seqLengths[iParent1-1],maxLength);
				childMatrix1[iParent1] = child1;
				childMatrix2[iParent1] = child2;
			}//FOR secuencias
			parent1.alignMatrix = childMatrix1;
			parent2.alignMatrix = childMatrix2;
		}//IF cruza...
		newPop.push_back(parent1);
		newPop.push_back(parent2);
	}//FOR individuos
	return newPop;
}

/*
	Función auxiliar para la cruza por columnas, es la que va haciendo la corrección en las secuencias..
*/
void AG::checkSequences(vector<int>& c1,int numEnz,int maxLength)
{
	int enz = 0;
	for(int iCol = 0; iCol < (int)c1.size(); iCol++)
	{
		enz +=c1[iCol];
	}
	while(enz > numEnz) //Tiene mas enzimas, cambiar por 0's
	{
		int random = rand() % maxLength;
		while(c1[random] == 0)
		{
			random = rand() % maxLength;
		}
		c1[random] = 0;
		enz--;
	}
	while(enz < numEnz) //Tiene mas gaps, cambiar por 1's
	{
		int random = rand() % maxLength;
		while(c1[random] == 1)
		{
			random = rand() % maxLength;
		}
		c1[random] = 1;
		enz++;
	}
}

/*
	Crossover con selector en el crosspoint
	Es igual que el crosssoverColumnas2, pero usando el selector para el punto de corte
*/
vector<Individuo> AG::crossoverColumnsSelectPoint(vector<Individuo> oldPop, Selector select)
{
	vector<Individuo> newPop;
	vector<Individuo>::iterator it;
	Individuo parent1;
	Individuo parent2;
	int n = (int)(getSizePop()/2);
	int numSeqs = population.getNumSeqs();
	int maxLength = population.getMaxLength();
	vector<int> lengths = population.getSeqLengths();
	//IMPORTANTE!!! No hay que cambiar el orden de las secuencias, o si se cambian hay que tener
	//en cuenta que el key es consecutivo nos sirve para identificar que secuencias es!!
	for(int index = 0; index < n; index++)
	{
		parent1 = oldPop[index]; //El primero
		parent2 = oldPop[sizePop-1-index]; //El ultimo..
		map<int, EC> matrix1 = parent1.alignMatrix;
		map<int, EC> matrix2 = parent2.alignMatrix;		
		double coin = (double)rand () / (double) RAND_MAX;
		if (coin < crossRate)
		{			
			//Aqui elijo la columna al azar...es la misma para todas las secuencias..
			int mL1 = matrix1[1].size();
			int mL2 = matrix2[1].size();
			int random = select.selectPunto(getRandom());
			for(int iSeq = numSeqSinAlinear; iSeq <= numSeqs; iSeq++)
			{				
				vector<int> seq1 = matrix1[iSeq];
				vector<int> seq2 = matrix2[iSeq];
				//Creamos el array q sirve para inicializar el valarray..
				int* init = NULL;
				int* init2 = NULL;
				init = new int[maxLength];
				init2 = new int[maxLength];
				for(int iAux = 0; iAux < maxLength; iAux++)
				{
					init[iAux] = seq1[iAux];
					init2[iAux] = seq2[iAux];
				}
				//Inicializamos los arreglos circulares..
				valarray<int> circular1 (init,maxLength);
				valarray<int> circular2 (init2,maxLength);
				valarray<int> aux;
				//Prueba a) Conservando sola la primera mitad de los padres... ESTA FUNCIONA MEJOR
				matrix1[iSeq] = crossoverColumnas2Child(circular1,circular2,lengths[iSeq - 1],mL1,mL2,random);
                                matrix2[iSeq] = crossoverColumnas2Child(circular2,circular1, lengths[iSeq - 1], mL2,mL1,random);
                                //Prueba b) Tratando de mezclar las 2 mitades de los padres..
				//matrix1[iSeq] = crossoverColumnas3Child(circular1,circular2,lengths[iSeq - 1],maxLength,random);			
				//matrix2[iSeq] = crossoverColumnas3Child(circular2,circular1, lengths[iSeq - 1], maxLength,random);
				delete [] init;
			}//FOR secuencias
			parent1.alignMatrix = matrix1;
			parent2.alignMatrix = matrix2;
		} //IF cruza
		newPop.push_back(parent1);
		newPop.push_back(parent2);
	}//FOR individuos
	return newPop;
}

//////////////////////////////////////////////////////
//	Funciones de control del AG
//////////////////////////////////////////////////////

/*
	En esta se minimiza, se usan los operadores normales..
	CHECAR cuales son los operadores que se van a dejar
*/
Individuo AG::start()
{
	//La poblacion inicial se recibe en el constructor
	population.calculateFitness();
	Individuo best = population.getBestIndividual();
	Individuo newBest = best;
	double bestFit = best.getFitness();
	int countGen = 1;
	vector<Individuo> newPopulation;
	int MAX_WITHOUT_CHANGE = 30;
	int countWC = 0;
	/////////////////////////////
	//Selector...
	int l;
	if ((population.getMaxLength() % 2) != 0)
		l = population.getMaxLength() - 1;
	else
		l = population.getMaxLength();
	select = Selector(2,l);
	/////////////////////////////
	//int banAux = 0;
	while ( countWC < MAX_WITHOUT_CHANGE )
	{
		newPopulation = selection(population.individuals, 2, 1, best);
		newPopulation = crossoverColumnas2(newPopulation);
		newPopulation = mutationGaps(newPopulation); population.setMaxLength(newMaxLength);
		//newPopulation = mutationGaps2(newPopulation); population.setMaxLength(newMaxLength);
		//Asignamos a Poblacion la newPopulation
		population.individuals = newPopulation;		
		population.calculateFitness();
		newBest = population.getBestIndividual();
		//Hay que actualizar el maxLength de best..
		best = actualizarMaxLengthBest(best);		
		if (newBest.getFitness() < bestFit)
		{
			bestFit = newBest.getFitness();
			best = newBest;
			countWC = 0;
		}
		//cout << "Gen = " << countGen << " Best Fitness = "  << bestFit << " Average Fitness = "  << population.getAverageFitness() <<endl;    	
		//best.print(population.getSequences());	
		countGen++;
		countWC++;
	}
	return best;
}

/*
	Función auxiliar que se usa cuando se usan operadores que cambian la longitud de las secuencias
*/
Individuo AG::actualizarMaxLengthBest(Individuo b)
{
	int numSeqs = population.getNumSeqs();
	for(int iSeq = 1; iSeq <= numSeqs; iSeq++)
		for(int iBest = b.getMaxLength(); iBest < newMaxLength; iBest++)
				(b.alignMatrix[iSeq]).push_back(0);
	b.setMaxLength(newMaxLength);
	return b;
}

/*
	Start con el crossover de lso blockes, minimiza
*/
Individuo AG::start2()
{
	//La poblacion inicial se recibe en el constructor
	population.calculateFitness();
	Individuo best = population.getBestIndividual();
	Individuo newBest = best;
	double bestFit = best.getFitness();
	int countGen = 1;
	vector<Individuo> newPopulation;
	int MAX_WITHOUT_CHANGE = 20;
	int countWC = 0;
	/////////////////////////////
	//Selector...
	int l;
	if ((population.getMaxLength() % 2) != 0)
		l = population.getMaxLength() - 1;
	else
		l = population.getMaxLength();
	Selector select(2,l);
	/////////////////////////////
	while ( countWC < MAX_WITHOUT_CHANGE )
	{
		newPopulation = selection(population.individuals, 2, 1, best);
		//newPopulation = crossoverExchange(newPopulation);						
		//newPopulation = crossoverCircular(newPopulation);
		//newPopulation = crossoverColumns(newPopulation);
		newPopulation = crossoverColumnsSelectPoint(newPopulation, select);
		mutationPermutate(newPopulation);
		//newPopulation = mutationCircular(newPopulation);
		//newPopulation = mutationGaps(newPopulation);
		//Asignamos a Poblacion la newPopulation
		population.individuals = newPopulation;		
		population.calculateFitness();
		newBest = population.getBestIndividual();
		if (newBest.getFitness() < bestFit)
		{
			bestFit = newBest.getFitness();
			best = newBest;
			countWC = 0;
		}
		//cout << "Gen = " << countGen << " Best Fitness = "  << bestFit << " Average Fitness = "  << population.getAverageFitness() <<endl;    		
		countGen++;
		countWC++;
		select.actualizaFitness(population, 2);
	}
	return best;
}

/*
	Start con el crossover de lso blockes, maximiza!!
*/
Individuo AG::start2Max()
{
	//La poblacion inicial se recibe en el constructor
	population.calculateFitness();
	Individuo best = population.getBestIndividual();
	Individuo newBest = best;
	double bestFit = best.getFitness();
	int countGen = 1;
	vector<Individuo> newPopulation;
	int MAX_WITHOUT_CHANGE = 20;
	int countWC = 0;
	/////////////////////////////
	//Selector...
	int l;
	if ((population.getMaxLength() % 2) != 0)
		l = population.getMaxLength() - 1;
	else
		l = population.getMaxLength();
	Selector select(2,l);
	/////////////////////////////
	while ( countWC < MAX_WITHOUT_CHANGE )
	{
		newPopulation = selection(population.individuals, 2, 1, best);
		//newPopulation = crossoverExchange(newPopulation);						
		//newPopulation = crossoverCircular(newPopulation);
		//newPopulation = crossoverColumns(newPopulation);
		newPopulation = crossoverColumnsSelectPoint(newPopulation, select);
		mutationPermutate(newPopulation);
		//newPopulation = mutationCircular(newPopulation);
		//newPopulation = mutationGaps(newPopulation);
		//Asignamos a Poblacion la newPopulation
		population.individuals = newPopulation;		
		population.calculateFitness();
		newBest = population.getBestIndividual();
		if (newBest.getFitness() > bestFit)
		{
			bestFit = newBest.getFitness();
			best = newBest;
			countWC = 0;
		}
		//cout << "Gen = " << countGen << " Best Fitness = "  << bestFit << " Average Fitness = "  << population.getAverageFitness() <<endl;    		
		countGen++;
		countWC++;
		select.actualizaFitness(population, 2);
	}
	return best;
}

/*
	Control principal del AG , maximizando!
	CHECAR cuales son los operadores que se quedan..
*/
Individuo AG::startMax()
{
	//La poblacion inicial se recibe en el constructor
	population.calculateFitness();
	Individuo best = population.getBestIndividual();
	Individuo newBest = best;
	double bestFit = best.getFitness();
	int countGen = 1;
	vector<Individuo> newPopulation;
	int MAX_WITHOUT_CHANGE = 30;
	int countWC = 0;
	/////////////////////////////
	//Selector...
	int l;
	if ((population.getMaxLength() % 2) != 0)
		l = population.getMaxLength() - 1;
	else
		l = population.getMaxLength();
	select = Selector(2,l);
	/////////////////////////////
	while ( countWC < MAX_WITHOUT_CHANGE )
	{
		newPopulation = selection(population.individuals, 2, 1, best);
		//newPopulation = sus(population.individuals, sizePop,best);
		//newPopulation = crossoverExchange(newPopulation);						
		//mutationPermutate(newPopulation);
		
		//newPopulation = crossoverCircular(newPopulation);
		//newPopulation = crossoverColumns(newPopulation);
		newPopulation = crossoverColumnas2(newPopulation);
		//newPopulation = mutationGaps(newPopulation); population.setMaxLength(newMaxLength);
		//mutationPermutate(newPopulation);
		//if	 (countWC < .23)
		//{
		//	mutRate = .01;
			//mutationPermutate(newPopulation);
			//newPopulation = crossoverExchange(newPopulation);						
			//newPopulation = crossoverCircular(newPopulation);
			//select.actualizaFitness(population, 2);
		//	newPopulation = mutationGaps(newPopulation); population.setMaxLength(newMaxLength);
		//	banAux = 1;						
		//}
		//else
		//{
		//	banAux = 0;
			//newPopulation = crossoverColumnas2(newPopulation);
			newPopulation = mutationCircular(newPopulation);
		//}
		//Asignamos a Poblacion la newPopulation
		population.individuals = newPopulation;		
		population.calculateFitness();
		newBest = population.getBestIndividual();
		//Hay que actualizar el maxLength de best..
		//if (banAux == 1)
		//	best = actualizarMaxLengthBest(best);		
		if (newBest.getFitness() > bestFit)
		{
			bestFit = newBest.getFitness();
			best = newBest;
			countWC = 0;
		}
		//cout << "Gen = " << countGen << " Best Fitness = "  << bestFit << " Average Fitness = "  << population.getAverageFitness() <<endl;    	
		//best.print(population.getSequences());	
		countGen++;
		countWC++;
	}
	return best;
}

