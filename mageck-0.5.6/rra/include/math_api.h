/*
 *  math_api.h
 *  RegulatorPrediction
 *
 *  Created by Han Xu on 3/4/13.
 *  Copyright 2013 Dana Farber Cancer Institute. All rights reserved.
 *
 */

typedef struct
{
	double value;
	int index;
}INDEXED_FLOAT;

//Quicksort an array in real values, in ascending order
void QuicksortF(double *a, int lo, int hi);

//Quicksort an indexed array, in ascending order 
void QuicksortIndexedArray(INDEXED_FLOAT *a, int lo, int hi);

//BTreeSearchingF: Searching value in array, which was organized in ascending order previously
int  bTreeSearchingF(double value, double *a, int lo, int hi);

//Rank the values in a float array and store the rank values in an integer array
void Ranking(int *rank, double *values, int sampleNum);

//Normal score transform
int NormalTransform(double *destA, int *rank, int sampleNum);

//Compute distance correlation. dim: number of samples; inputNum: number of variables in the input; input: the input array with inputNum*dim items; output: the output array
double ComputeDistanceCorrelation(double *input, double *output, int inputNum, int dim);

//Pearson correlation
double PearsonCorrel(double *a, double *b, int dim);

//Partial correlation
double PartialCorrel(double *a, double *b, double *control, int dim);

//Randomly permute an array of float values
void PermuteFloatArrays(double *a, int size);

//Compute CDF of a non-central beta distribution. when lambda is 0.0, it's cpf of beta distribution
double BetaNoncentralCdf(double a, double b, double lambda, double x, double error_max);
