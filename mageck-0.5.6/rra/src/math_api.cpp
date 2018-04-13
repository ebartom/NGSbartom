/*
 *  math_api.c
 *  RegulatorPrediction
 *
 *  Created by Han Xu on 3/4/13.
 *  Copyright 2013 Dana Farber Cancer Institute. All rights reserved.
 *
 */

#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include "math_api.h"
#include "rvgs.h"

// normalInv: from Ziegler's code
double normalInv(double p);

//compute Euclidean distance
double EucliDist(double *a, double *b, int dim);

//Compute logarithm of Gamma function. flag=0, no error; flag=1, x<=0
double LogGamma(double x, int *flag);

//Compute incomplete beta function ratio
double betain (double x, double p, double q, double beta, int *ifault);

//BTreeSearchingF: Searching value in array, which was organized in ascending order previously
int  bTreeSearchingF(double value, double *a, int lo, int hi)
{
	if (value<=a[lo])
	{
		return lo;
	}
	
	if (value>=a[hi])
	{
		return hi;
	}
	
	if (hi-lo<=1)
	{
		if (fabs(a[hi]-value)>fabs(a[lo]-value))
		{
			return lo;
		}
		else
		{
			return hi;
		}
	}
	
	if (value>=a[(lo+hi)/2])
	{
		return bTreeSearchingF(value, a, (lo+hi)/2, hi);
	}
	else
	{
		return bTreeSearchingF(value, a,  lo, (lo+hi)/2);
	}
}

//Quicksort an array in real values, in ascending order
void QuicksortF(double *a, int lo, int hi)
{
	int i=lo, j=hi;
	double x=a[(lo+hi)/2];
	double h;
	
	if (hi<lo)
	{
		return;
	}
	
    //  partition
    while (i<=j)
    {    
		while ((a[i]<x)&&(i<=j))
		{
			i++;
		}
		while ((a[j]>x)&&(i<=j))
		{
			j--;
		}
        if (i<=j)
        {
			h = a[i];
			a[i] = a[j];
			a[j] = h;
            i++; j--;
        }
    } 
	
    //  recursion
    if (lo<j) QuicksortF(a, lo, j);
    if (i<hi) QuicksortF(a, i, hi);
}

//Quicksort an indexed array, in ascending order 
void QuicksortIndexedArray(INDEXED_FLOAT *a, int lo, int hi)
{
	int i=lo, j=hi;
	double x=a[(lo+hi)/2].value;
	INDEXED_FLOAT h;
	
	if (hi<lo)
	{
		return;
	}
	
    //  partition
    while (i<=j)
    {    
		while ((a[i].value<x)&&(i<=j))
		{
			i++;
		}
		while ((a[j].value>x)&&(i<=j))
		{
			j--;
		}
        if (i<=j)
        {
			memcpy(&h, a+i, sizeof(INDEXED_FLOAT));
			memcpy(a+i,a+j,sizeof(INDEXED_FLOAT));
			memcpy(a+j, &h, sizeof(INDEXED_FLOAT));
            i++; j--;
        }
    } 
	
    //  recursion
    if (lo<j) QuicksortIndexedArray(a, lo, j);
    if (i<hi) QuicksortIndexedArray(a, i, hi);
}

//Normal score transform
int NormalTransform(double *destA, int *rank, int sampleNum)
{
	int i;
	
	for (i=0;i<sampleNum;i++)
	{
		destA[i] = normalInv(((double)rank[i]+0.5)/sampleNum);
	}
	
	return 1;
}

//Rank the values in a float array and store the rank values in an integer array
void Ranking(int *rank, double *values, int sampleNum)
{
	double *tmp;
	int i;
	
	tmp = (double *)malloc(sampleNum*sizeof(double));
	
	memcpy(tmp, values, sampleNum*sizeof(double));
	
	QuicksortF(tmp, 0, sampleNum-1);
	
	for (i=0;i<sampleNum;i++)
	{
		rank[i] = (bTreeSearchingF(values[i]-0.0000001, tmp, 0, sampleNum-1)+bTreeSearchingF(values[i]+0.0000001, tmp, 0, sampleNum-1))/2;
	}
	
	free(tmp);
}

// normalInv: from Ziegler's code
double normalInv(double p)
{
	double a1 = -39.69683028665376;
	double a2 = 220.9460984245205;
	double a3 = -275.9285104469687;
	double a4 = 138.3577518672690;
	double a5 =-30.66479806614716;
	double a6 = 2.506628277459239;
	
	double b1 = -54.47609879822406;
	double b2 = 161.5858368580409;
	double b3 = -155.6989798598866;
	double b4 = 66.80131188771972;
	double b5 = -13.28068155288572;
	
	double c1 = -0.007784894002430293;
	double c2 = -0.3223964580411365;
	double c3 = -2.400758277161838;
	double c4 = -2.549732539343734;
	double c5 = 4.374664141464968;
	double c6 = 2.938163982698783;
	
	double d1 = 0.007784695709041462;
	double d2 = 0.3224671290700398;
	double d3 = 2.445134137142996;
	double d4 = 3.754408661907416;
	
	//Define break-points.
	
	double p_low =  0.02425;
	double p_high = 1 - p_low;
	double q, x, r;
	
	//Rational approximation for lower region.
	
	if (0 < p && p < p_low) 
	{
		q = sqrt(-2*log(p));
		x = (((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / ((((d1*q+d2)*q+d3)*q+d4)*q+1);
	}
	
	//Rational approximation for central region.
	
	if (p_low <= p && p <= p_high) 
	{
		q = p - 0.5;
		r = q*q;
		x = (((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q / (((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1);
	}
	
	//Rational approximation for upper region.
	
	if (p_high < p && p < 1) 
	{
		q = sqrt(-2*log(1-p));
		x = -(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / ((((d1*q+d2)*q+d3)*q+d4)*q+1);
	}
	//y = exp(sigma*x+my);
	
	return x;
}

//Compute distance correlation. dim: number of samples; inputNum: number of variables in the input; input: the input array with inputNum*dim items; output: the output array
double ComputeDistanceCorrelation(double *input, double *output, int inputNum, int dim)
{
	int i,j,k;
	double *dist1, *dist2, *vect1,*vect2, *meanRow1, *meanRow2, meanAll1, meanAll2;
	double score, dCov, dVar1, dVar2;
	
	dist1 = (double *)malloc(dim*dim*sizeof(double));
	dist2 = (double *)malloc(dim*dim*sizeof(double));
	vect1 = (double *)malloc(inputNum*sizeof(double));
	vect2 = (double *)malloc(inputNum*sizeof(double));
	meanRow1 = (double *)malloc(dim*sizeof(double));
	meanRow2 = (double *)malloc(dim*sizeof(double));
	meanAll1 = 0;
	meanAll2 = 0;
	
	for (i=0;i<dim;i++)
	{
		meanRow1[i] = 0;
		meanRow2[i] = 0;
		
		for (j=0;j<dim;j++)
		{
			for (k=0;k<inputNum;k++)
			{
				vect1[k] = input[k*dim+i];
				vect2[k] = input[k*dim+j];
			}
			
			dist1[i*dim+j] = EucliDist(vect1,vect2,inputNum);
			dist2[i*dim+j] = sqrt((output[i]-output[j])*(output[i]-output[j]));
			
			meanRow1[i] += dist1[i*dim+j];
			meanRow2[i] += dist2[i*dim+j];
			meanAll1 += dist1[i*dim+j];
			meanAll2 += dist2[i*dim+j];
		}
		
		meanRow1[i] /= dim;
		meanRow2[i] /= dim;
	}
	
	meanAll1 /= (dim*dim);
	meanAll2 /= (dim*dim);
	
	dCov = 0;
	dVar1 = 0;
	dVar2 = 0;
	
	for (i=0;i<dim;i++)
	{
		for (j=0;j<dim;j++)
		{
			dist1[i*dim+j] = dist1[i*dim+j] - meanRow1[i] - meanRow1[j] + meanAll1;
			dist2[i*dim+j] = dist2[i*dim+j] - meanRow2[i] - meanRow2[j] + meanAll2;
			dCov += dist1[i*dim+j]*dist2[i*dim+j];
			dVar1 += dist1[i*dim+j]*dist1[i*dim+j];
			dVar2 += dist2[i*dim+j]*dist2[i*dim+j];
		}
	}
	
	dCov = sqrt(dCov/(dim*dim));
	dVar1 = sqrt(dVar1/(dim*dim));
	dVar2 = sqrt(dVar2/(dim*dim));
	
	score = dCov/sqrt(dVar1*dVar2);
	
	free(dist1);
	free(dist2);
	free(vect1);
	free(vect2);
	free(meanRow1);
	free(meanRow2);
	
	return score;
}

//compute Euclidean distance
double EucliDist(double *a, double *b, int dim)
{
	int i;
	double sum = 0;
	
	for (i=0;i<dim;i++)
	{
		sum += (a[i]-b[i])*(a[i]-b[i]);
	}
	
	return sqrt(sum);
}

//Randomly permute an array of float values
void PermuteFloatArrays(double *a, int size)
{
	int i;
	double tmp;
	double r;
	int index;

	for (i=0;i<size-1;i++)
	{
		r = Uniform(0,1);
		index = i+ (int)(r*(size-i));
		
		if ((index<i)||(index>=size))
		{
			continue;
		}
		
		tmp = a[i];
		a[i] = a[index];
		a[index] = tmp;
	}	
}

//Pearson correlation
double PearsonCorrel(double *a, double *b, int dim)
{
	int i;
	double mean1, mean2, sumAB, sumAA, sumBB;
	
	mean1 = 0;
	mean2 = 0;
	
	for (i=0;i<dim;i++)
	{
		mean1 += a[i];
		mean2 += b[i];
	}
	
	mean1 /= dim;
	mean2 /= dim;
	
	sumAB = 0;
	sumAA = 0;
	sumBB = 0;
	
	for (i=0;i<dim;i++)
	{
		sumAB += (a[i]-mean1)*(b[i]-mean2);
		sumAA += (a[i]-mean1)*(a[i]-mean1);
		sumBB += (b[i]-mean2)*(b[i]-mean2);
	}
	
	return sumAB/sqrt(sumAA*sumBB+0.00000000001);
}

//Partial correlation
double PartialCorrel(double *a, double *b, double *control, int dim)
{
	double corAB, corAC, corBC;
	
	corAB = PearsonCorrel(a, b, dim);
	corAC = PearsonCorrel(a, control, dim);
	corBC = PearsonCorrel(b, control, dim);
	
	return (corAB-corAC*corBC)/(sqrt(1-corAC*corAC)*sqrt(1-corBC*corBC)+0.00000000001);
}

//Compute logarithm of Gamma function. flag=0, no error; flag=1, x<=0
double LogGamma(double x, int *flag)
{
	double f;
	double value;
	double y;
	double z;
	
	if (x<=0.0)
	{
		*flag = 1;
		return 0.0;
	}
	
	*flag = 0;
	y = x;
	
	if (x < 7.0)
	{
		f = 1.0;
		z = y;
		
		while (z < 7.0)
		{
			f = f * z;
			z = z + 1.0;
		}
		y = z;
		f = - log (f);
	}
	else
	{
		f = 0.0;
	}
	
	z = 1.0 / y / y;
	
	value = f + ( y - 0.5 ) * log ( y ) - y
    + 0.918938533204673 +
    (((
	   - 0.000595238095238   * z
	   + 0.000793650793651 ) * z
	  - 0.002777777777778 ) * z
	 + 0.083333333333333 ) / y;
	
	return value;
}

//Compute incomplete beta function ratio
double betain ( double x, double p, double q, double beta, int *ifault )
{
	double acu = 0.1E-14;
	double ai;
	double cx;
	int indx;
	int ns;
	double pp;
	double psq;
	double qq;
	double rx;
	double temp;
	double term;
	double value;
	double xx;
	
	value = x;
	*ifault = 0;
	/*
	 Check the input arguments.
	 */
	if ( p <= 0.0 || q <= 0.0 )
	{
		*ifault = 1;
		return value;
	}
	
	if ( x < 0.0 || 1.0 < x )
	{
		*ifault = 2;
		return value;
	}
	/*
	 Special cases.
	 */
	if ( x == 0.0 || x == 1.0 )
	{
		return value;
	}
	/*
	 Change tail if necessary and determine S.
	 */
	psq = p + q;
	cx = 1.0 - x;
	
	if ( p < psq * x )
	{
		xx = cx;
		cx = x;
		pp = q;
		qq = p;
		indx = 1;
	}
	else
	{
		xx = x;
		pp = p;
		qq = q;
		indx = 0;
	}
	
	term = 1.0;
	ai = 1.0;
	value = 1.0;
	ns = ( int ) ( qq + cx * psq );
	/*
	 Use the Soper reduction formula.
	 */
	rx = xx / cx;
	temp = qq - ai;
	if ( ns == 0 )
	{
		rx = xx;
	}
	
	for ( ; ; )
	{
		term = term * temp * rx / ( pp + ai );
		value = value + term;
		temp = fabs( term );
		
		if ( temp <= acu && temp <= acu * value )
		{
			value = value * exp ( pp * log ( xx )
								 + ( qq - 1.0 ) * log ( cx ) - beta ) / pp;
			
			if ( indx )
			{
				value = 1.0 - value;
			}
			break;
		}
		
		ai = ai + 1.0;
		ns = ns - 1;
		
		if ( 0 <= ns )
		{
			temp = qq - ai;
			if ( ns == 0 )
			{
				rx = xx;
			}
		}
		else
		{
			temp = psq;
			psq = psq + 1.0;
		}
	}
	
	return value;
}

//Compute CDF of a non-central beta distribution. when lambda is 0.0, it's cpf of beta distribution
double BetaNoncentralCdf ( double a, double b, double lambda, double x, double error_max )
{
	double beta_log;
	double bi;
	double bj;
	int i;
	int ifault;
	double p_sum;
	double pb_sum;
	double pi;
	double pj;
	double si;
	double sj;
	double value;
	
	i = 0;
	pi = exp ( - lambda / 2.0 );
	
	beta_log = LogGamma ( a, &ifault )
	+ LogGamma ( b, &ifault )
	- LogGamma ( a + b, &ifault );
	
	bi = betain ( x, a, b, beta_log, &ifault );
	
	si = exp (
			  a * log ( x )
			  + b * log ( 1.0 - x )
			  - beta_log
			  - log ( a ) );
	
	p_sum = pi;
	pb_sum = pi * bi;
	
	while ( p_sum < 1.0 - error_max )
	{
		pj = pi;
		bj = bi;
		sj = si;
		
		i = i + 1;
		pi = 0.5 * lambda * pj / ( double ) ( i );
		bi = bj - sj;
		si = x * ( a + b + i - 1 ) * sj / ( a + i );
		
		p_sum = p_sum + pi;
		pb_sum = pb_sum + pi * bi;
	}
	
	value = pb_sum;
	
	return value;
}

