/***********************************************************************************
// OtherSimilarity.cpp : Other Similairty Measures, e.g. EMD
//            
// 
// Written by Dong-Qing Zhang, 2008
//
// Contact : 
//			Dong-Qing Zhang :	 dongqing@gmail.com, 
//			Prof.Shih-Fu Chang : sfchang@ee.columbia.edu
************************************************************************************/

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <assert.h>

#include "Garg.h"

// download emd.h from http://ai.stanford.edu/~rubner/emd/default.htm

#ifdef EMD_CODE_INCLUDED

	#include "emd.h" // after download, you need to modify feature_t using the following code in #else

#else

	#define FEAT_DIM  17  // Added by Dong-Qing Zhang, Columbia U.
	typedef double feature_t[FEAT_DIM];   // Modified by Dong-Qing Zhang, Columbia U.

	typedef struct
	{
		int n;                /* Number of features in the signature */
		feature_t *Features;  /* Pointer to the features vector */
		float *Weights;       /* Pointer to the weights of the features */
	} signature_t;

#endif

float _fdistance(feature_t * f1, feature_t * f2)
{
	int i;
	double sum=0;

	for(i=0;i<FEAT_DIM;i++)
	{
		sum += ((*f1)[i]-(*f2)[i])*((*f1)[i]-(*f2)[i]);
	}

	return float(sqrt(sum));
}

float EmdSimilarity(GARG *garg1, GARG *garg2)
{
	signature_t g1,g2;
	int i;
	float dist;

// load g1
	g1.n = garg1->node_num;
	g1.Features = (feature_t *)malloc(sizeof(feature_t)*g1.n);
	g1.Weights  = (float *)malloc(sizeof(float)*g1.n);

	for(i=0;i<g1.n;i++)
	{
		garg1->GetNodeFeature(g1.Features[i],i);

		g1.Features[i][0]/=320;
		g1.Features[i][1]/=240;
		g1.Features[i][2]/=256;
		g1.Features[i][3]/=256;
		g1.Features[i][4]/=256;

		g1.Weights[i] = 1.0f/g1.n;
	}

// load g2
	g2.n = garg2->node_num;
	g2.Features = (feature_t *)malloc(sizeof(feature_t)*g2.n);
	g2.Weights  = (float *)malloc(sizeof(float)*g2.n);

	for(i=0;i<g2.n;i++)
	{
		garg2->GetNodeFeature(g2.Features[i],i);

		g2.Features[i][0]/=320;
		g2.Features[i][1]/=240;
		g2.Features[i][2]/=256;
		g2.Features[i][3]/=256;
		g2.Features[i][4]/=256;

		g2.Weights[i] = 1.0f/g2.n;
	}

// Earth Mover's Distance
// download emd.c from http://ai.stanford.edu/~rubner/emd/default.htm
// then change emd.c to emd.cpp

#ifdef EMD_CODE_INCLUDED
	dist = emd(&g1,&g2,_fdistance,NULL,NULL);
#else
	dist = 0;
#endif

	free(g1.Features);
	free(g2.Features);

	return 1-dist;
}
