#include "cuda_runtime.h"
#include <stdio.h>
#include <stlib.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include <curand.h>
#include <curand_kernel.h>

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
#define NBLOCKS 8
#define NTHREADS 1

double correctTotal[N_EVALS_MAX];
double eval_x[N_EVALS_MAX][N_COEFFICIENTS_MAX];
int N_COEFFICIENTS;
int N_EVALS;

/*
  Got this code from Wikipedia.
  Converts a binary code to binary gray code.
*/
__device__ unsigned int grayToBinary(unsigned int num)
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
__device__ double coefToDouble(unsigned int n)
{
  return (double)grayToBinary(n)/UINT_MAX * (MAX_VALUE-MIN_VALUE)+MIN_VALUE;
}


/*
  Takes as input all the messages and fitnesses.
  Runs a tournament selection algorithm to determine next generation
*/
__device__ void reproduction(unsigned int messages[N_MESSAGES][N_COEFFICIENTS_MAX], double fitnesses[], unsigned int *seed)
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
__device__ void crossover(unsigned int parent1[], unsigned int parent2[], 
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
__device__  unsigned int mutation(unsigned int *coefficient, unsigned int *seed)
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
__device__ double fitness(unsigned int message[])
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
__device__ void cme(unsigned int message1[], unsigned int message2[], 
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

__global__ void kernel(int *numBlocks, int **messages)
{
  int nprocs;
  int myid = blockIdx.x;
  int neighborID1,neighborID2 = ;
}

__global__ void setup_rng(curandState *state)
{
    int id = threadIdx.x + blockIdx.x * NTHREADS;
    /* Each thread gets same seed, a different sequence 
       number, no offset */
    curand_init(1234, id, 0, &state[id]);
}

int main(int argc, char *argv[])
{
  FILE *file; 
  unsigned int messages[N_MESSAGES][N_COEFFICIENTS_MAX];
  unsigned int **d_messages;
  double fitnesses[N_MESSAGES];
  int n,i,j;
  unsigned int initTime = time(NULL);
  unsigned int *seed=&initTime;
  double functionCoefs[N_COEFFICIENTS_MAX];
  
  // create curand states, allocate memory on device, and initialize rng on device.
  curandState *devStates;
  cudaMalloc((void**)&devStates, NBLOCKS * NTHREADS * sizeof(curandState));
  setup_rng<<<NBLOCKS, NTHREADS>>>(d_nthreads, devStates);

  // create strategies, allocate mem on device, and copy to device.
  for(i=0;i<N_MESSAGES;i++)
  {
    for(j=0;j<N_COEFFICIENTS;j++)
    {
      messages[i][j]=rand_r(seed)*2;
    }
    fitnesses[i]=fitness(messages[i]);
  }

  int size = N_COEFFICIENTS_MAX * sizeof(int);
  cudaMalloc((void**)&d_messages, N_MESSAGES * sizeof(int*));
  for(i=0;i<N_MESSAGES;i++)
    cudaMalloc(&messages[i], size);
  cudaMemcpy2D(&d_messages, size, &messages, size, size, size);



  file = fopen("fitness.txt","a+");

  cudaMalloc( (void **) &d_messages, N_MESSAGES * N_COEFFICIENTS_MAX * sizeof(unsiged int *));


  /* copy messages to the device */
  cudaMemcpy(d_messages, messages, size, cudaMemcpyHostToDevice);

  kernel<<<NBLOCKS,NTHREADS>>>(d_numBlocks, d_messages);
}
