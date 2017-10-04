/***********************************************************************************
// Sampler.cpp : statistical sampler for Monte Carlo inference
//            
// 
// Written by Dong-Qing Zhang, DVMM lab, Columbia University, 2005
//
// Contact : 
//			Dong-Qing Zhang :	 dongqing@gmail.com, 
//			Prof.Shih-Fu Chang : sfchang@ee.columbia.edu
************************************************************************************/

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "assert.h"
#include "time.h"

#include "BaseClass.h"

#include "Sampler.h"

// statistical sampler for monte carlo analysis

Sampler::Sampler()
{
	m_cumu = NULL;
	m_prob_len = 0;

/*
	srand( (unsigned)time( NULL ) );
	for(i=0;i<10;i++)
		rand(); */
}

Sampler::~Sampler()
{
	if(m_cumu)
		MYFree(m_cumu);
}

void Sampler::InitDiscSampler(double * prob,int prob_len)
{
	int i;

	if(m_cumu)
		MYFree(m_cumu);

//  Generate cumulative distribution

	m_cumu = (double *)MYMalloc(prob_len*sizeof(double));
	for(i=0;i<prob_len;i++)
	{
		if(i==0)
		{
			m_cumu[i] = prob[0];
			continue;
		}
		m_cumu[i] = m_cumu[i-1]+prob[i];
	}
	
//  Initialize random number generator
}

// discrete sampling from a multinomial distribution
int Sampler::DiscSampling()
{
	double rand_num;
	int i;

	srand( (unsigned)time( NULL ) );

	rand_num = (double)rand()/RAND_MAX;
	for(i=0;i<m_prob_len;i++)
		if(rand_num>=m_cumu[i])
			return i;

	return -1;  // error
}

// Discrete sampling in 1 step
int Sampler::DiscSampling1Step(double * prob,int prob_len)
{
	double rand_num,cumu;
	int i,s;

	rand_num = (double)rand()/RAND_MAX;
	if(prob_len==2)  // binary sampling
	{
		if(rand_num>prob[0])
			return 1;
		else
			return 0;
	}

	cumu = 0;
	s = 0;
	for(i=0;i<prob_len;i++)
	{
		if(rand_num>= cumu)
		{
			s++;
			cumu += prob[i];
			continue;
		}

		assert(s>0);
		return s-1;
	}

	assert(s>0);
	return s-1;  // error
}
