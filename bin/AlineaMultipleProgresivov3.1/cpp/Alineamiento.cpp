#include "../h/Alineamiento.h"
#include "../h/Secuencia.h"
#include "../h/Individuo.h"
#include "../h/Poblacion.h"
#include "../h/AG.h"

using namespace std;

/*
	Esta es la clase principal, que lleva el control del alineamiento iterativo.
	Aqui se crean las secuencias, y se hace la llamada al AG.
	En la versión 1.0, la función de fitness tiene 2 factores, homogeneidad y penalización por gaps para 
	cualquiera de las 3 opciones de fitness que se elija.
*/

/*
	Auxiliar para dividir las cadenas
*/
void Tokenize(const string& str,
                      vector<string>& tokens,
                      const string& delimiters = " ")
{
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}



Alineamiento::Alineamiento(void)
{
}

Alineamiento::~Alineamiento(void)
{
}

/*
	Constructor
	namefile = Archivo con las secuencias a alinear
	fileCalis = Este es un archivo que contiene la lista de las secuencias del primer archivo y sus calificaciones de alineamientos
			por pares.
*/
Alineamiento::Alineamiento(string nameFile,string fileCalis,double c,double m, int sP, int nG, int fT, int cT, double gap, double inc, double homo, int nExp)
{
	numExp = nExp;
	penGap = gap;
	penCol = inc;
	facHomo = homo;
	readFile(nameFile);//Ahi se crear las secuencias...
	readMatrixCalis(fileCalis);
	ordenarSecuencias();
	setGAParams(c, m, sP, nG);
	createAlineamientoProgresivo(fT, cT);
}

/*
	Setea los parámetros del AG
*/
void Alineamiento::setGAParams(double c,double m, int sP, int nG)
{
	pc = c;
	pm = m;
	sizePob = sP;
	numGen = nG;
}

/*
	El número de secuencias que se van a alinear en este momento
*/
void Alineamiento::setNumSeq(int num)
{
	numSeq = num;
}

/*
	Se calculan las longitudes de todas las secuencias	
*/
void Alineamiento::setSeqLengths()
{
	int MAX  = -1;
	iSeqMasLarga = 0;
	for(int i = 0; i < numSeq; i++)
	{
		if (sequences[i]->getLength() > MAX)
		{
			MAX = sequences[i]->getLength(); 
			iSeqMasLarga = i;
		}
		seqLengths.push_back((sequences[i])->getLength());
	}
	maxLength = MAX + MAX/2;
}

/*
	Número de secuencias totales
*/
int Alineamiento::getNumSeq()
{
	return numSeq;
}

/*
	Regresa probabilidad de cruza
*/
double Alineamiento::getPc()
{
	return pc;
}

/*
	Regresa probabilidad de mutación
*/
double Alineamiento::getPm()
{
	return pm;
}

/*
	Regresa el tamaño de la población
*/
int Alineamiento::getSizePob()
{
	return sizePob;
}


/*
	Regresa el número máximo de iteraciones sin cambio para terminar el AG
*/
int Alineamiento::getNumGen()
{
	return numGen;
}


/*
	Regresa el vector de longitudes de las secuencias
*/
vector<int>  Alineamiento::getSeqLengths()
{
	return seqLengths;
}

/*
	Lee el archivo que contiene las secuencias que se van a alinear, y crea los objetos correspondientes
*/
void Alineamiento::readFile(string nameFile)
{
	int seqTotal = 0;
	int mapTotal = 0;
	
	ifstream read;
	string line;
	
	read.open(nameFile.c_str(),ifstream::in);
	string nameMap;
	getline (read,nameMap);//Leemos la primera línea que es el nombre	
	mapTotal++;
	Secuencia* s;
	while (read.good())//Las demas líneas son las secuencias
	{
		getline (read,line);
		if (line.find(">") == string::npos && !line.empty())
		{ //No es mapa!, empezamos con las secuencias del nameMap
			//Cambio para tener los nombres de los mapas de c/secuencia
			//130510
			vector<string> tokens;
			Tokenize(line, tokens, "\t");	
			seqTotal++;
			s = createSequence(tokens[1], tokens[0], tokens[2],nameMap); //nombre del mapa y secuencia
			//OJOOOOOOOO a lo mejor tenemos q cambiar el orden desp...segun el orden de las calificaciones..
			sequences.push_back(s);
		}
		else
		{ //Es el nombre de un mapa
			mapTotal++;
			nameMap = line;
		}		
	}	
	read.close();	
	setNumSeq(seqTotal);
	setSeqLengths();	
}

/*
	Lee el archivo de las calificaciones de los alineamientos por pares
*/
void Alineamiento::readMatrixCalis(string nameFile)
{
	ifstream read;
	string line;
	read.open(nameFile.c_str(),ifstream::in);		
	while (read.good())//Las demas líneas son las secuencias
	//while(!read.getline(line).eof())
	{
		getline (read,line);
		if (!line.empty()){
			vector<string> tokens;
			Tokenize(line, tokens, "\t");
			string indice = tokens[0].append("-");
			indice.append(tokens[1]);
			matrixCalis[indice] = atof(tokens[2].c_str());
		}
	}	
	read.close();
}

/*
	Crea los objetos Secuencia
*/
Secuencia* Alineamiento::createSequence(string nMap, string seqID, string seqLine,string nomCluster)
{
	Secuencia* seq = new Secuencia(nomCluster, seqLine, nMap, seqID);
	return seq;
}

/*
	Ordena las secuencias de acuerdo a sus calificaciones.
	De menor a mayor, para ir alineando las secuencias de acuerdo a este orden.
	Primero se elige la secuencia mas larga, y de ahi se buscan cuales son las secuencias mejor alineadas a esta
	y se ordenan de mejor a peor, y en ese orden se van alineando.
*/
void Alineamiento::ordenarSecuencias()
{
	//iSeqMasLarga es el indice de la secuencia mas larga..
	string idSeqLarga = (sequences[iSeqMasLarga])->getSeqId();
	//Ahora buscar los valores en la matrizCalis
	map<string, double>::iterator it;
	multimap<double,string> ::iterator itMM;
	vector<string> ids;	
	multimap<double,string> ordenados;
	map<string,double> pares;
	
	for(it =  matrixCalis.begin(); it != matrixCalis.end(); it++)
	{
		string indice = (*it).first;
		size_t found;
		found = indice.find(idSeqLarga);
		if (found!=string::npos)
		{
			vector<string>	tokens;
			Tokenize(indice, tokens, "-");			
			if (tokens[0].compare(idSeqLarga) == 0)
			{
				ids.push_back(tokens[1]);
				ordenados.insert(pair<double,string>((*it).second,tokens[1]));
			}
			else
			{
				ids.push_back(tokens[0]);
				ordenados.insert(pair<double,string>((*it).second,tokens[0]));
			}
		}
	}
	//Pasarlo a un mapa..
	//Multimap calificacion -> ID
	//Ya esta ordenado en "ordenados", ahora hay que cambiar el orden en q aparecen sequences..
	vector<Secuencia*> nuevo;
	vector<int> nuevoLengths;
	//En primer lugar hay que agregar a la secuencia mas larga..la base..
	string buscar = idSeqLarga;
	for(int i = 0; i < (int) sequences.size(); i++)
	{
		if(buscar.compare(sequences[i]->getSeqId()) == 0)
		{
			nuevo.push_back(sequences[i]);
			nuevoLengths.push_back(sequences[i]->getLength());
			i = (int) sequences.size();
		}
	}	
	for(itMM = ordenados.begin(); itMM != ordenados.end(); itMM++)
	{
		buscar = (*itMM).second;
		for(int i = 0; i < (int) sequences.size(); i++)
		{
			if(buscar.compare(sequences[i]->getSeqId()) == 0)
			{
				nuevo.push_back(sequences[i]);
				nuevoLengths.push_back(sequences[i]->getLength());
				i = (int) sequences.size();
			}
		}
	}
	sequences = nuevo;
	seqLengths = nuevoLengths;

}

/*
        Esta es la función que lleva el control del AG, aqui se elige el tipo de ciclo que se va a usar, 
        y recibe el resultado del AG.
	Esta se llama solamente la primera vez para las primeras 2 secuencias.
*/
Individuo Alineamiento::createGA(int fitnessType, int crossType)
{	
	Poblacion pop(sizePob, sequencesParciales, seqLengthsParciales, maxLength, fitnessType, penGap, penCol, facHomo);
	AG ag(numGen, pc, pm, pop);
	if (crossType == 1 && fitnessType == 1)
                best = ag.start();//minimizar sin bloques
        else
                if (crossType == 1 && fitnessType != 1)
                        best = ag.startMax();//maximiza  sin bloques
                else
                        if (crossType != 1 && fitnessType != 1)
                                best = ag.start2Max();//Maximiza con cloques
                        else
                                best = ag.start2(); //Minimiza con bloques
	return best;
}

/*
	Esta función es similar a la anterior, solo que aqui, las poblaciones ya toman en cuenta las secuencias
	que ya estan alineadas, y esas no se pueden mover, los operadores no afectan estas secuencias (numSeqSinAlinear).
	Solo se usan para calcular el fitness de los individuos.
*/
Individuo Alineamiento::createGAAlineados(vector<EC> convertidos, int fitnessType, int crossType, int numSeqSinAlinear)
{	
	Poblacion pop(convertidos, sizePob, sequencesParciales, seqLengthsParciales, maxLength, fitnessType,numSeqSinAlinear, penGap, penCol, facHomo);
	AG ag(numGen, pc, pm, pop, numSeqSinAlinear);
	if (crossType == 1 && fitnessType == 1)
                best = ag.start();//minimizar sin bloques
        else
                if (crossType == 1 && fitnessType != 1)
                        best = ag.startMax();//maximiza  sin bloques
                else
                        if (crossType != 1 && fitnessType != 1)
                                best = ag.start2Max();//Maximiza con cloques
                        else
                                best = ag.start2(); //Minimiza con bloques
	return best;
}

/*
	Esta función es la que lleva el control del alineamiento iterativo, donde se hacen las llamadas a los distintos AG
	El primer AG, es el normal, después el resultado se va agregando a un vector de secuencias alineadas, y estas sirven
	de entrada para la siguiente llamada, donde ya se llama al AGProgresivo que solo afecta a la última secuencia, que
	es la que se esta alineando en este momento.
*/
void Alineamiento::createAlineamientoProgresivo(int fT, int cT)
{
	//El vector de secuencias ya debe de estar en el orden de mayor a menor!!!
	/////////////////////////////////////////
	////////EMPIEZA EL CICLO
	//La primera corrida es el crearAG normal..., hay q ir agregando de uno en uno las secuencias al vector de sequences
	/////////////////////////////////////////
	//Hay que crear un arreglo de sequences y seqLengths parcial, ahi se van a ir agregando las secuencias una a una
	int iSeqActual = 0;
	vector<string> alineados;
	vector<EC> convertidos;
	int total = (int)sequences.size();
	sequencesParciales.clear();
	seqLengthsParciales.clear();
	//Primer alineamiento
	sequencesParciales.push_back(sequences[iSeqActual]);seqLengthsParciales.push_back(seqLengths[iSeqActual]);iSeqActual++;
	sequencesParciales.push_back(sequences[iSeqActual]);seqLengthsParciales.push_back(seqLengths[iSeqActual]);iSeqActual++;
	//AG normal..
	Individuo best;
	double bestFit = 9999999.0;
	vector<string> secs;
	for(int iExp = 0; iExp < numExp; iExp++)
	{
		Individuo resultado = createGA(fT, cT);
		if (resultado.getFitness() < bestFit )
		{
			best = resultado;
			bestFit = resultado.getFitness();
		}	
	}
	//best.print(sequencesParciales);
	secs = best.getPrintSecuencias(sequencesParciales);
	//Agregamos el resultado al vector de secuencias alineadas
	alineados.push_back(secs[0]);alineados.push_back(secs[1]);
	convertidos = convertirAlineadosBin(alineados);
	//////////
	int iSecuencias;
	for(iSecuencias = iSeqActual; iSecuencias < total; iSecuencias++)
	{	//Agregamos la sig secuencia
		sequencesParciales.push_back(sequences[iSecuencias]);seqLengthsParciales.push_back(seqLengths[iSecuencias]);//iSeqActual++;
		bestFit = 9999999.0;
		for(int iExp = 0; iExp < numExp; iExp++)
		{	//creamos el AG..
			Individuo resultado = createGAAlineados(convertidos, fT, cT, iSecuencias + 1);
			if (resultado.getFitness() < bestFit )
			{
				best = resultado;
				bestFit = resultado.getFitness();
			}			
		}
		//best.print(sequencesParciales);
		secs = best.getPrintSecuencias(sequencesParciales);
		int ultima = (int)secs.size();
		ultima--;
		alineados.push_back(secs[ultima]);
		convertidos = convertirAlineadosBin(alineados);
	}
	best.print(sequencesParciales);
}

/*
	Esta función auxiliar toma la salida del AG (la salida a pantalla, las secuencias alineadas), y la convierte a 0 y 1's
	que son los que sirven para los individuos de la siguiente población.
*/
vector<EC> Alineamiento::convertirAlineadosBin(vector<string> alineados)
{
	vector<EC> convertidos;
	vector<int> aux;
	for(int i = 0; i < (int)alineados.size(); i++)
	{
		vector<string> tokens;
		Tokenize(alineados[i], tokens, "\t");
		vector<string> elementos;
		Tokenize(tokens[1], elementos, ":");
		aux.clear();
		for(int iElementos = 0; iElementos < (int)elementos.size(); iElementos++)
		{
			if (elementos[iElementos].compare("-.-.-") != 0)
				aux.push_back(1);
			else
				aux.push_back(0);
		}
		convertidos.push_back(aux);
	}
	return convertidos;
}

/*
	Regresa el mejor alineamiento que encontro el AG
*/
Individuo Alineamiento::getBest()
{
        return best;
}

/*
	Regresa el arreglo de secuencias
*/
vector<Secuencia*> Alineamiento::getSequences()
{
        return sequences;
}

