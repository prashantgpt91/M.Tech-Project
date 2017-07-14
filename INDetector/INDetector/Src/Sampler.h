/***********************************************************************************
// Sampler.h : statistical sampler for Monte Carlo inference
//            
// 
// Written by Dong-Qing Zhang, DVMM lab, Columbia University, 2005
//
// Contact : 
//			Dong-Qing Zhang :	 dongqing@gmail.com, 
//			Prof.Shih-Fu Chang : sfchang@ee.columbia.edu
************************************************************************************/

#ifndef SAMPLER_H
#define SAMPLER_H

// statistical sampler for monte carlo analysis
class Sampler
{
	double * m_cumu;	  // cumulated probability distribution
	int      m_prob_len;  // the size of a discrete distribution

public:
	int DiscSampling1Step(double * prob,int prob_len);
	int DiscSampling();
	Sampler();
	~Sampler();

// Init a discrete sampler
	void InitDiscSampler(double * prob,int prob_len);  // initialize the sampler
};

#endif
