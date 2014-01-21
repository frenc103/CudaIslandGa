#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include "mpi.h"

#define N_MESSAGES 200
#define CROSS_RATE 1.0
#define MUT_RATE 0.05
#define N_ROUNDS 1000
#define N_COEFFICIENTS_MAX 10
#define STR_LENGTH 32
#define MAX_VALUE 50
#define MIN_VALUE -50
#define EVAL_MIN -10
#define EVAL_MAX 10
#define EVAL_INC 0.1
#define N_EVALS_MAX 400

double correctTotal[N_EVALS_MAX];
double eval_x[N_EVALS_MAX][N_COEFFICIENTS_MAX];
int N_COEFFICIENTS;
int N_EVALS;

/*
  Got this code from Wikipedia.
  Converts a binary code to binary gray code.
*/
unsigned int grayToBinary(unsigned int num)
{
    unsigned int numBits = 8 * sizeof(num);
    unsigned int shift;
    for (shift = 1; shift < numBits; shift = 2 * shift)
    {
        num = num ^ (num >> shift);
    }
    return num;
}

/*
  Converts an integer coefficient to its double value.
  First converts the integer to the integer equivalent in gray coding.
  It is then normalized between MIN_VALUE and MAX_VALUE
*/
double coefToDouble(unsigned int n)
{
  return (double)grayToBinary(n)/UINT_MAX * (MAX_VALUE-MIN_VALUE)+MIN_VALUE;
}


/*
  Takes as input all the messages and fitnesses.
  Runs a tournament selection algorithm to determine next generation
*/
void reproduction(unsigned int messages[N_MESSAGES][N_COEFFICIENTS_MAX], double fitnesses[], unsigned int *seed)
{
  int i,j,r1,r2;
  unsigned int newMessages[N_MESSAGES][N_COEFFICIENTS];
  double newFitnesses[N_MESSAGES];

  for(i=0;i<N_MESSAGES;i++)
    {
      for(j=0;j<N_COEFFICIENTS;j++)
	{
	  newMessages[i][j]=messages[i][j];
	}
      newFitnesses[i] = fitnesses[i];
    }
  for(i=0;i<N_MESSAGES;i++)
    {
      r1=rand_r(seed)%N_MESSAGES;
      r2=rand_r(seed)%N_MESSAGES;
      if(newFitnesses[r1]>newFitnesses[r2])
	{
	  for(j=0;j<N_COEFFICIENTS;j++)
	    {
	      messages[i][j]=newMessages[r1][j];
	    }
	  fitnesses[i]=newFitnesses[r1];
	}
      else
	{
	  for(j=0;j<N_COEFFICIENTS;j++)
	    {
	      messages[i][j]=newMessages[r2][j];
	    }
	  fitnesses[i]=newFitnesses[r2];
	}
    }
}

/*
  Performs crossover. Takes as input the parents and children.
*/
void crossover(unsigned int parent1[], unsigned int parent2[], 
	       unsigned int child1[], unsigned int child2[], unsigned int *seed)
{
  int r, n;
  unsigned int mask;
  for(n=0;n<N_COEFFICIENTS;n++)
    {
      if(rand_r(seed)<CROSS_RATE*RAND_MAX)
	{
	  r=rand_r(seed)%STR_LENGTH;
	  mask=~0>>r;
	  child1[n]=(parent1[n] & mask) | (parent2[n] & ~mask);
	  child2[n]=(parent2[n] & mask) | (parent1[n] & ~mask);
	}
      else
	{
	  child1[n]=parent1[n];
	  child2[n]=parent2[n];
	}
    }
}

/*
  Performs mutation on a single coefficient.
  Each bit of the coefficient has a MUT_RATE chance of being toggled.
 */
unsigned int mutation(unsigned int *coefficient, unsigned int *seed)
{
  unsigned int numBits = 8 * sizeof(coefficient);
  unsigned int shift;
  for (shift = 0; shift < numBits; shift++)
    {
      if(rand_r(seed)<MUT_RATE*RAND_MAX)
	*coefficient ^= (1 << shift);
    }

  return *coefficient;
}

/*
  Calculates the fitness of the message.
 */
double fitness(unsigned int message[])
{
  int i,n;
  double coefs[N_COEFFICIENTS_MAX];
  double sqrError=0.0;
  double total;

  for(n=0;n<N_COEFFICIENTS;n++)
    {
      coefs[n]=coefToDouble(message[n]);
    }
  for(i=0;i<=N_EVALS;i++)
    {
      total=0.0;
      for(n=0;n<N_COEFFICIENTS;n++)
	{
	  total+=coefs[n]*eval_x[i][n];
	}
      sqrError+=pow(total-correctTotal[i],2.0);
    }
  return -sqrError;
}

/*
  Takes two messages as input.
  Performs crossover to generate two children, then mutates the children.
  Children fitnesses are calculated.
  If a parent's child is superior, the parent is replaced. 
 */
void cme(unsigned int message1[], unsigned int message2[], 
	 double *fitness1, double *fitness2, unsigned int *seed)
{
  unsigned int child1[N_COEFFICIENTS_MAX];
  unsigned int child2[N_COEFFICIENTS_MAX];
  double fitnessC1,fitnessC2;
  int i;
 
  crossover(message1,message2,child1,child2,seed);

  for(i=0;i<N_COEFFICIENTS;i++)
    {
      mutation(&child1[i],seed);
      mutation(&child2[i],seed);
    }

  fitnessC1=fitness(child1);
  fitnessC2=fitness(child2);

  //insert election code here
  if (fitnessC1>*fitness1)
    {
      *fitness1=fitnessC1;
      for(i=0;i<N_COEFFICIENTS;i++)
	{
	  message1[i]=child1[i];
	}
    }
  if (fitnessC2>*fitness2)
    {
      *fitness2=fitnessC2;
      for(i=0;i<N_COEFFICIENTS;i++)
	{
	  message2[i]=child2[i];
	}
    }
}


int main(int argc, char *argv[])
{
  FILE *file; 
  unsigned int messages[N_MESSAGES][N_COEFFICIENTS_MAX];
  double fitnesses[N_MESSAGES];
  int n,i,j;
  unsigned int initTime = time(NULL);
  unsigned int *seed=&initTime;
  double functionCoefs[N_COEFFICIENTS_MAX];
  MPI_Status status;
  MPI_Request messageStatus, fitnessStatus;
  int nprocs;
  int myid;
  int neighborID1,neighborID2;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  *seed+=(myid*100);
  if (myid == 0)
    {
      file = fopen("fitness.txt","a+");
      neighborID1=nprocs-1;
    }
  else
    {
      neighborID1=myid-1;
    }
  if (myid == nprocs-1)
    {
      neighborID2=0;
    }
  else
    {
      neighborID2=myid+1;
    }
  int count=0;
  N_COEFFICIENTS=0;
  if (argc > 1)
    {
      count=1;
      while(count<argc)
	{
	  if (strcmp(argv[count],"-f")==0)
	    {
	      N_COEFFICIENTS=atoi(argv[++count]);
	      for(i=0;i<N_COEFFICIENTS;i++)
		{
		  count++;
		  if(count>=argc)
		    {
		      printf("Error: expected additional inputs\n");
		      MPI_Finalize();
		      exit(1);
		    }
		  functionCoefs[i]=atof(argv[count]);		  
		}
	    }
	  count++;
	}
    }
  else
    {
      printf("Usage: -f N_COEFFICIENTS coef_1 coef_2 ... coef_N\n");
      MPI_Finalize();
      exit(1);
    }

  N_EVALS = (int)((EVAL_MAX-EVAL_MIN)/EVAL_INC + 1);

  for(i=0;i<N_EVALS;i++)
    {
      correctTotal[i]=0;
      for(j=0;j<N_COEFFICIENTS;j++)
      	{
	  eval_x[i][j]=pow(EVAL_MIN*1.0+i*EVAL_INC*1.0,(N_COEFFICIENTS-j-1)*1.0);
	  correctTotal[i]+=functionCoefs[j]*eval_x[i][j];
      	}
    }

  for(i=0;i<N_MESSAGES;i++)
    {
      for(j=0;j<N_COEFFICIENTS;j++)
	{
	  messages[i][j]=rand_r(seed)*2;
	}
      fitnesses[i]=fitness(messages[i]);
    }
    for(n=0;n<N_ROUNDS;n++)
    {
      reproduction(messages,fitnesses,seed);
      if (myid ==0)
        {
	  fprintf(file,"%f\n",fitnesses[0]); 
	}
      if (n % 50 == 0)
	{
	  MPI_Isend(&messages[0],N_COEFFICIENTS,MPI_INT,
		    neighborID1,10,MPI_COMM_WORLD,&messageStatus);
	  MPI_Isend(&fitnesses[0],1,MPI_DOUBLE,
		    neighborID1,11,MPI_COMM_WORLD,&fitnessStatus);
	  
	  MPI_Recv(&messages[N_MESSAGES-1],N_COEFFICIENTS,MPI_INT,
		   neighborID2,10,MPI_COMM_WORLD, &status);
	  MPI_Recv(&fitnesses[N_MESSAGES-1],1,MPI_DOUBLE,
		   neighborID2,11,MPI_COMM_WORLD,&status);
	  
	  MPI_Isend(&messages[0],N_COEFFICIENTS,MPI_INT,
		    neighborID2,10,MPI_COMM_WORLD,&messageStatus);
	  MPI_Isend(&fitnesses[0],1,MPI_DOUBLE,
		    neighborID2,11,MPI_COMM_WORLD,&fitnessStatus);
	  
	  MPI_Recv(&messages[N_MESSAGES-2],N_COEFFICIENTS,MPI_INT,
		   neighborID1,10,MPI_COMM_WORLD,&status);
	  MPI_Recv(&fitnesses[N_MESSAGES-2],1,MPI_DOUBLE,
		   neighborID1,11,MPI_COMM_WORLD,&status);
	}

      #pragma omp parallel for
      for(i=0;i<N_MESSAGES;i+=2)
	{
	  cme(messages[i],messages[i+1],&fitnesses[i],&fitnesses[i+1],seed);
	}
    }
    if (myid == 0)
      {
	fclose(file);
	for(i=0;i<N_MESSAGES;i++)
	  {
	    printf("message %d coefs: ",i);
	    for(j=0;j<N_COEFFICIENTS;j++)
	      {
		printf("%f, ",coefToDouble(messages[i][j]));
	      }
	    printf("fitness: %f\n",fitnesses[i]);
	  }
      }
    MPI_Finalize();
    return 0;
}
