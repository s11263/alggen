# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <cstring>
# include <chrono>

using namespace std;

# define POPSIZE 25			//  wielkość populacji
# define GENCOUNT 250  		        //  maksymalna ilość pokoleń
# define VARCOUNT 3			//  ilość zmiennych
# define PROPABX 0.5			//  prawdopodobieństwo zajścia crossover
# define PROPABMUT 0.2                  //  prawdopodobieństwo mutacji

struct genotype				//  Każdy GENOTYPE jest osobnikiem populacji
{
  double gene[VARCOUNT]; 	        	//  string z zapisanymi zmiennymi
  double fitness;			//  zdolność prokreacji
  double upper[VARCOUNT];	        	//  górna granica zmiennych
  double lower[VARCOUNT];	         	//  dolna granica zmiennych
  double roundFit;			//  względny fintess
  double cumulFit;			//  kumulacyjny fitness
};

struct genotype population[POPSIZE+1];
struct genotype newpopulation[POPSIZE+1]; 

int main ( );
void crossover ( int &seed );
void elite ( );
void evaluate ( );
int pseudo_I ( int a, int b, int &seed );
void initialize ( int &seed );                       
void top_mem ( );
void mutate ( int &seed );
double pseudo_II ( double a, double b, int &seed );
void report ( int generation );
void selector ( int &seed );
void crossO ( int one, int two, int &seed );

double kingOfTheHill;             //  aktualnie najwyższy wynik
double pretender;                 //  pretender do najwyższego wyniku
int incumbency;                   //  kadencja najwyższego wyniku
double bestRooms[VARCOUNT];
  
struct neighbour 
{
	double rooms[VARCOUNT]; 	     	    //  tablica z wartościami zmiennych
    double rent;			                //  wartość funkcji
    double topBorder[VARCOUNT];	        	//  górna granica zmiennych
    double bottomBorder[VARCOUNT];	        //  dolna granica zmiennych	
};
  
void checkChange ( double pretender );      //  sprawdź czy nastąpiła zmiana najlepszego wyniku
void pickSpot ( int &seed );            	//  wybierz punkt startowy
void visitNeighbours ( );                 	//  odwiedź sąsiadów, sprawdź czy któryś nie jest lepszym wynikiem
void rentCalc ( );							//  obliczenie wartości funkcji


int main ( )
{
  auto t0 = std::chrono::high_resolution_clock::now();
  int generation;
  int i;
  int seed;

  if ( VARCOUNT < 2 )  		// gdy liczba zmiennych jest poniżej 2 nie możliwe jest by zaszedł crossover.
  {
    cout << "\nZ powodu zbyt niskiej ilości zmiennych crossover nigdy nie zajdzie.\n";
  }

  seed = 123456789;

  initialize ( seed );

  evaluate ( );

  top_mem ( );

  for ( generation = 0; generation < GENCOUNT; generation++ )
  {
    selector ( seed );
    crossover ( seed );
    mutate ( seed );
    report ( generation );
    evaluate ( );
    elite ( );
  }

  cout << "\n  Najlepszy osobnik po " << GENCOUNT << " populacji:\n\n";
  
  for ( i = 0; i < VARCOUNT; i++ )
  {
    cout << "  var(" << i << ") = " << population[POPSIZE].gene[i] << "\n";
  }

  cout << "\n";
  cout << "  Najlepszy fitness = " << population[POPSIZE].fitness << "\n";

  cout << "\n  Program wykonał się prawidłowo\n\n";
  auto t1 = std::chrono::high_resolution_clock::now();
  auto dt = 1.e-9*std::chrono::duration_cast<std::chrono::nanoseconds>(t1-t0).count();
  cout << "\n Czas potrzeny na wykonanie algorytmu genetycznego: " << dt;
  
  // ALGORTYM WSPINACZKOWY:
  auto t2 = std::chrono::high_resolution_clock::now();
  
  seed = 123456789;
  incumbency = -1;
  neighbour.topBorder[0] = -1.0;                                //  podajemy przedziały ograniczające zmienne
  neighbour.topBorder[1] = 1.0;
  neighbour.topBorder[2] = 0.0;

  neighbour.bottomBorder[0] =  1.0;
  neighbour.bottomBorder[1] =  4.0;
  neighbour.bottomBorder[2] =  4.0;
  
  while ( incumbency < 1000 )                //  Zapewni nam brak pętli nieskończonej
  {
    pickSpot ( seed );
    checkChange ( neighbour.rent );
	visitNeighbours ( );  
  }
    
  
  cout << "\n Najlepsze wartości zmiennych to: " << bestRooms[0] << ", " << bestRooms[1] << ", " << bestRooms[2] << ". \n" 
  cout << "\n Najwyższy wynik funkcji to: " << kingOfTheHill << ". \n";
   
  auto t3 = std::chrono::high_resolution_clock::now();
  auto dt2 = 1.e-9*std::chrono::duration_cast<std::chrono::nanoseconds>(t3-t2).count();
  cout << "\n Czas potrzeny na wykonanie algorytmu wspinaczkowego: " << dt2;
  
  return 0;
}

void pickSpot ( int &seed ) 
{
  double x;
  int i;
  i = 0;
  
  for ( i = 0; i < VARCOUNT; i++ )
  {
	x = pseudo_II ( neighbour.bottomBorder[i], neighbour.topBorder[i], seed );
	neighbour.rooms[i] = x;
  }
  rentCalc ( );
  return;
}

void rentCalc ( ) {     //  trzeba zmienić w zależności od wzoru funkcji 
	neighbour.rent = 2 * neighbour.rooms[0] * neighbour.rooms[0]   +  neighbour.rooms[1] + 2 * neighbour.rooms[2];
	return; 
}

void checkChange (double pretender) 
{
	int i;
	i = 0;
	
	if ( incumbency == -1 )
	{
		incumbency = 0; 
		kingOfTheHill = pretender;
		
		for ( i; i < VARCOUNT; i++)
		{
			bestRooms[i] = neighbour.rooms[i];
		}
	} else 
	{
		if ( pretender < = kingOfTheHill) 
		{
			incumbency++;
		} else
		{
			incumbency = 0;
			kingOfTheHill = pretender; 
			for ( i; i < VARCOUNT; i++)
			{
				bestRooms[i] = neighbour.rooms[i];
			}
		}
	}
	return;
}

void visitNeighbours ( )						//  kolejno zwiększam każdą zmienną
{
	int i;
	i = 0;
	
	for ( i; i < VARCOUNT; i++ ) {
		if ( neighbour.rooms[i] <= neighbour.topBorder[i] ) 
		{
			neighbour.rooms[i] += 0.01;
			rentCalc ( );
			checkChange( neighbour.rent );
		}
	}
	
	return;
}

void crossover ( int &seed )			//  crossover wybiera dwóch rodziców do pojedynczego punktu crossovera. Seed dla generatora liczb losowych
{
  const double a = 0.0;
  const double b = 1.0;
  int mem;								//  member - osobnik
  int one = 0;
  int first = 0;  						//  ilość wybranych osobników
  double x;

  for ( mem = 0; mem < POPSIZE; ++mem )
  {
    x = pseudo_II ( a, b, seed );
    if ( x < PROPABX )
    {
      ++first;
      if ( first % 2 == 0 )
      {
        crossO ( one, mem, seed );
      }
      else
      {
        one = mem;
      }
    }
  }
  return;
}

void elite ( )					//  ma za zadanie przechowywać najlepszego osobnika poprzedniej populacji
								//  w przypadku, gdy najlepszy z aktualnej populacji jest gorszy od najlepszego z poprzedniej
								//  wtedy ten ostatni stanie się najgorszym osobnikiem obecnej populacji
{
  int i;
  double best; 					//  najlepsza wartość fitness
  int best_mem = 0;
  double worst;					//  najgorsza wartość fitness
  int worst_mem = 0;

  best = population[0].fitness;
  worst = population[0].fitness;

  for ( i = 0; i < POPSIZE - 1; ++i )
  {
    if ( population[i+1].fitness < population[i].fitness )
    {
      if ( best <= population[i].fitness )
      {
        best = population[i].fitness;
        best_mem = i;
      }

      if ( population[i+1].fitness <= worst )
      {
        worst = population[i+1].fitness;
        worst_mem = i + 1;
      }
    }
    else
    {
      if ( population[i].fitness <= worst )
      {
        worst = population[i].fitness;
        worst_mem = i;
      }

      if ( best <= population[i+1].fitness )
      {
        best = population[i+1].fitness;
        best_mem = i + 1;
      }
    }
  }
                   
//  W przypadku gdy najlepszy osobnik z nowej populacji jest lepszy od najlepszego z poprzedniej to kopiujemy wartość z nowej populacji.
//  W przeciwnym razie najgorszy z obecnej populacji zastąpić trzeba najlepszym z poprzedniej. 
  if ( population[POPSIZE].fitness <= best )
  {
    for ( i = 0; i < VARCOUNT; i++ )
    {
      population[POPSIZE].gene[i] = population[best_mem].gene[i];
    }
    population[POPSIZE].fitness = population[best_mem].fitness;
  }
  else
  {
    for ( i = 0; i < VARCOUNT; i++ )
    {
      population[worst_mem].gene[i] = population[POPSIZE].gene[i];
    }
    population[worst_mem].fitness = population[POPSIZE].fitness;
  } 
  return;
}

void evaluate ( )			//  implementuje zdefiniowaną przez użytkownika funkcję skalarną. Przykładowa funkcja na potrzeby tego programu to:
							//  2x[1]^2+(x[2]+x[3])-x[3]
{
  int member;
  int i;
  double x[VARCOUNT+1];

  for ( member = 0; member < POPSIZE; member++ )
  {
    for ( i = 0; i < VARCOUNT; i++ )
    {
      x[i+1] = population[member].gene[i];
    } 
    population[member].fitness = ( 2 * x[1] * x[1] ) + ( x[2] + x[3] ) + x[3];  //  w celu zmiany funkcji, zmienić tutaj. 
  }
  return;
}

int pseudo_I ( int a, int b, int &seed )  //  zwraca przeskalowaną pseudolosową liczbę pomiędzy a i b. Dystrybucja typu uniform.
											   //  a i b to granice przedziału. Seed to waratość seed, która powinna być różna od 0.
{
  int c;
  const int huge = 2147483647;     		   //  minimalny standard, opracowany przez Lewisa, Goodmana i Millera w 1969 r. 
  int k;
  float r;
  int value;

  if ( seed == 0 )
  {
    cerr << "\n Zadana wartość SEED = 0.\n";
    exit ( 1 );
  }

  if ( b < a )  							  //  zapewniamy przedział od a do b
  {
    c = a;
    a = b;
    b = c;
  }

  k = seed / 127773;

  seed = 16807 * ( seed - k * 127773 ) - k * 2836;

  if ( seed < 0 )
  {
    seed = seed + huge;
  }

  r = ( float ) ( seed ) * 4.656612875E-10;
  r = ( 1.0 - r ) * ( ( float ) a - 0.5 )  + r * ( ( float ) b + 0.5 );     //  r leży w przedziale a-0.5 i b+0.5
  value = round ( r ); 														//  zaokgrąglamy do liczby całkowitej z przedziału (a, b)
  if ( value < a )															//  upewniamy, się czy wartość jest w przedziale od a do b.
  {
    value = a;
  }
  if ( b < value )
  {
    value = b;
  }
  return value;
}

void initialize ( int &seed )													//  Losowo generuje liczby w tych zakresach. 															
{															
  int i;
  int j;
  double lbound[VARCOUNT];
  double ubound[VARCOUNT];
  
  lbound[0] = -1.0;                                //  podajemy przedziały ograniczające zmienne
  lbound[1] = 1.0;
  lbound[2] = 0.0;

  ubound[0] =  1.0;
  ubound[1] =  4.0;
  ubound[2] =  4.0;

  for ( i = 0; i < VARCOUNT; i++ )    		   // inicjowanie zmiennych wewnątrz przedziałów
  {
    for ( j = 0; j < POPSIZE; j++ )
    {
      population[j].fitness = 0;
      population[j].roundFit = 0;
      population[j].cumulFit = 0;
      population[j].lower[i] = lbound[i];
      population[j].upper[i]= ubound[i];
      population[j].gene[i] = pseudo_II ( lbound[i], ubound[i], seed );
    }
  }
  return;
}
 
void top_mem ( )      				//  trzyma rękę na pulsie, kto jest najlepszym osobnikiem populacji. Ostatni w tablicy jest najlepszy.
{
  int cur_best;								//  index najlepszego
  int mem;
  int i;

  cur_best = 0;

  for ( mem = 0; mem < POPSIZE; mem++ )
  {
    if ( population[POPSIZE].fitness < population[mem].fitness )
    {
      cur_best = mem;
      population[POPSIZE].fitness = population[mem].fitness;
    }
  }

  for ( i = 0; i < VARCOUNT; i++ )   			//  w przypadku znalezienia najlepszego osobnika, skopiuj jego geny. 
  {
    population[POPSIZE].gene[i] = population[cur_best].gene[i];
  }
  return;
}

void mutate ( int &seed )					//  odpowiada za losową, zunifikowaną mutację. Zmienna zmutowana zostaje zastąpiona nową losową wartością.
{
  const double a = 0.0;
  const double b = 1.0;
  int i;
  int j;
  double lbound;
  double ubound;
  double x;

  for ( i = 0; i < POPSIZE; i++ )
  {
    for ( j = 0; j < VARCOUNT; j++ )
    {
      x = pseudo_II ( a, b, seed );
      if ( x < PROPABMUT )
      {
        lbound = population[i].lower[j];
        ubound = population[i].upper[j];  
        population[i].gene[j] = pseudo_II ( lbound, ubound, seed );
      }
    }
  }

  return;
}

double pseudo_II ( double a, double b, int &seed ) 		//  zwraca liczbę pseudolosową. Dystrybucja typu uniform. Seed różny od 0.
{
  int huge = 2147483647;
  int k;
  double value;

  if ( seed == 0 )
  {
    cerr << "\n Wartość początkowa SEED = 0.\n";
    exit ( 1 );
  }

  k = seed / 127773;
  seed = 16807 * ( seed - k * 127773 ) - k * 2836;

  if ( seed < 0 )
  {
    seed = seed + huge;
  }

  value = ( double ) ( seed ) * 4.656612875E-10;
  value = a + ( b - a ) * value;
  return value;
}

void report ( int generation )					//  raportuje postęp symulacji. 
{
  double avg;							//  średni fitness populacji
  double best_val;						//  najlepszy fitness populacji
  int i;				
  double square_sum;					//  kwadrat sumy dla odchylenia standardowego
  double stddev;						//  odchylenie standardowe fitnessu populacji
  double sum;							//  całkowity fitness populacji
  double sum_square;					//  suma kwadratów dla odchylenia standardowego

  if ( generation == 0 )
  {
    cout << "\n";
    cout << "  Numer            Najlepszy       Średni        Odchylenie \n";
    cout << "  Populacji        wynik           fitness       standardowe \n";
    cout << "\n";
  }

  sum = 0.0;
  sum_square = 0.0;

  for ( i = 0; i < POPSIZE; i++ )
  {
    sum = sum + population[i].fitness;
    sum_square = sum_square + population[i].fitness * population[i].fitness;
  }

  avg = sum / ( double ) POPSIZE;
  square_sum = avg * avg * POPSIZE;
  stddev = sqrt ( ( sum_square - square_sum ) / ( POPSIZE - 1 ) );
  best_val = population[POPSIZE].fitness;

  cout << "  " << setw(8) << generation 
       << "  " << setw(14) << best_val 
       << "  " << setw(14) << avg 
       << "  " << setw(14) << stddev << "\n";

  return;
}

void selector ( int &seed )					//  odpowiada za selekcję oraz zapewnia, że najlepszy osobnik zawsze przetrwa. 
{
  const double a = 0.0;
  const double b = 1.0;
  int i;
  int j;
  int mem;
  double p;
  double sum;

  sum = 0.0;
  for ( mem = 0; mem < POPSIZE; mem++ )		// szukamy całkowitego fitnessu populacji.
  {
    sum = sum + population[mem].fitness;
  }

  for ( mem = 0; mem < POPSIZE; mem++ )		//  obliczamy względny fitness dla każdego osobnika.
  {
    population[mem].roundFit = population[mem].fitness / sum;
  }
  
  population[0].cumulFit = population[0].roundFit;   //  obliczamy kumulacyjny fitness. 
  for ( mem = 1; mem < POPSIZE; mem++ )
  {
    population[mem].cumulFit = population[mem-1].cumulFit +       
      population[mem].roundFit;
  }

  for ( i = 0; i < POPSIZE; i++ )					//  wybieramy ocalałych używając kumulacyjnego fitnessu.
  { 
    p = pseudo_II ( a, b, seed );
    if ( p < population[0].cumulFit )
    {
      newpopulation[i] = population[0];      
    }
    else
    {
      for ( j = 0; j < POPSIZE; j++ )
      { 
        if ( population[j].cumulFit <= p && p < population[j+1].cumulFit )
        {
          newpopulation[i] = population[j+1];
        }
      }
    }
  }

  for ( i = 0; i < POPSIZE; i++ )				//  nadpisujemy starą populację nową.
  {
    population[i] = newpopulation[i]; 
  }

  return;     
}



void crossO ( int one, int two, int &seed )     //  przeprowadza crossover wybranych rodziców. Zmienne oznaczają ich indexy
{
  int i;
  int point;								//  punkty, w których nastąpi crossover
  double t;

  point = pseudo_I ( 0, VARCOUNT - 1, seed );    //  losujemy punkt crossoveru.

  for ( i = 0; i < point; i++ )			//  zamienia geny od punktu 0 do point-1
  {
    t = population[one].gene[i];
    population[one].gene[i] = population[two].gene[i];
    population[two].gene[i] = t;
  }

  return;
}
