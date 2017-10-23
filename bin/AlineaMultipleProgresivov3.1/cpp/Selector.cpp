#include "../h/Selector.h"


/*
	Esta clase es la que implementa la ruleta para seleccionar el punto de cruza, 
	cuando se esta usando el operador que usa bloques.
	Todavia no esta completamente implementada PENDIENTE!
*/
Selector::Selector(void)
{
}

Selector::~Selector(void)
{
}

//l va a ser el numero par inmediatio inferior si la long no es par
Selector::Selector(int bs,int l){
	blockSize = bs;
	initializePoints(blockSize,l);
	totalPuntos = points.size();
	length = l;
	
	double aux = 1.0 /(double)(length - 1);
	for (int i = 1; i < length; i++)
		fitnessPuntos[i] = aux;
	
	//Calculando la probabilidad acumulada...
	probAcumuladaPuntos.insert(pair<int,double>(1,fitnessPuntos[1]));
	for (int i = 2; i < length; i++)
    {
		probAcumuladaPuntos.insert(pair<int,double>(i,probAcumuladaPuntos[i - 1] + fitnessPuntos[i]));		
    }
	probAcumuladaPuntos[length - 1] = 1.0;
}

void Selector::actualizaFitness(Poblacion& pob, int blockSize){
	//Calcular el número de bloques para cada punto frontera
	//ejem. si la cadena es de 8 bits, entonces tenemos que calcular
	//cuantos bloques hay en 0-2, 2-4, 4-6, 6-8	

	int length = pob.individuals[0].alignMatrix[1].size();
	if ((length % 2) != 0)
		length--;
	
	initializePoints(blockSize,length);
	totalPuntos = points.size();
	
	
	double aux = 1.0 /(double)(length - 1);
	for (int i = 1; i < length; i++)
		fitnessPuntos[i] = aux;
	
	//Calculando la probabilidad acumulada...
	probAcumuladaPuntos.insert(pair<int,double>(1,fitnessPuntos[1]));
	for (int i = 2; i < length; i++)
    {
		probAcumuladaPuntos.insert(pair<int,double>(i,probAcumuladaPuntos[i - 1] + fitnessPuntos[i]));		
    }
	probAcumuladaPuntos[length - 1] = 1.0;

	map<int,double>::iterator it;
	int sizePop = pob.getSizePop();
	//Ahora vamos a contar bloques, en este caso son 4, 2,4,6,8, el último es fuera del for
	double bloques = 0;
	int punto = 0;
	//Para Info double factor = 60;
	///////////////////
	//Factor para block
	int numSeqs = pob.getNumSeqs();
	double factor =  numSeqs*log( 1.0 / (double)numSeqs)*15
			+ numSeqs*log( 1.0 / (double)numSeqs)*10
			+ numSeqs*log( 1.0 / (double)numSeqs)*5 ;
	factor = factor * -1 * 2;
	//////////////////
	for(int i = 0; i < totalPuntos; i++)
	{
		bloques = pob.fitnessPerBlock(points[i],blockSize);
		bloques = bloques / (double) sizePop;
		//if (bloques > .9)
		//	fitnessPuntos[points[i]] = 1.0; 
			fitnessPuntos[points[i]] = (double)bloques / (double)sizePop; 
		//else
		//	fitnessPuntos[points[i]] = 0.0; 
			
		//Ahora hay que agregar los puntos que no son fronteras..
		for(int j = punto + 1; j < points[i]; j++)
			fitnessPuntos[j] = abs(1.0 - fitnessPuntos[points[i]] ); 
		punto = points[i];
	}
	//Hay que agregar el último bloque...
	bloques = pob.fitnessPerBlock(points[points.size()-1]+blockSize,blockSize);
	//if (bloques > .9)
	//	fitnessPuntos[points[points.size()-1]+blockSize] = 1.0; 
		fitnessPuntos[points[points.size()-1]+blockSize] = (double)bloques / (double)sizePop; 
	//else
	//	fitnessPuntos[points[points.size()-1]+blockSize] = 0.0; 
	//Ahora hay que agregar los puntos que no son fronteras..
	for(int j = points[points.size()-1] + 1; j < length; j++)
		fitnessPuntos[j] = abs(1.0 - fitnessPuntos[points[points.size()-1]+blockSize]); 
	//Calculando el fitness relativo de los puntos..
	double fitTotal = 0.0;
	for (int i = 1; i <= length; i++)
		fitTotal += fitnessPuntos[i];

	for (int i = 1; i < length; i++)
		fitnessPuntos[i] = (double) fitnessPuntos[i]/(double) fitTotal;

	//Calculando la probabilidad acumulada...
	probAcumuladaPuntos[1] = fitnessPuntos[1];
	for (int i = 2; i < length; i++)
    {
		probAcumuladaPuntos[i] = probAcumuladaPuntos[i - 1] + fitnessPuntos[i];		
    }
	probAcumuladaPuntos[length - 1] = 1.0;
 }
//RUleta
int Selector::selectPunto(double r){
	int point = 0;
	int j = 0;
	if ( r <= probAcumuladaPuntos[1])
		point = 1;
	else
	    {
		for (j = 2; j < length; j++)
		{
			if (probAcumuladaPuntos[j - 1] < r  &&  r <= probAcumuladaPuntos[j])
				point = j;
			else
				continue;
		}
	    }
	return point;
}


//Para el AGCruzaPuntoFijo
//***Funcion auxiliar para "Selector"
//Points contiene los posibles puntos de cruza
//Blocksize es el tamaño del bloque, si es igual a 1 entonces son consecutivos (no hay bloques)
//Length la longitud de la cadena
void Selector::initializePoints(int blockSize, int length)
{
        int numBloques = length / blockSize;
        int j = 0;
        //cout << "Length: " << length << " BS: " << blockSize <<endl;
        for (int i = blockSize; j < numBloques - 1; i = i + blockSize, j++)
        {
                points.push_back(i);
        }
}

