/***********************************************************************************
// IGraph.cpp : class for performing probabilistic inference, including LBP 
//				and Gibbs Sampling
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
#include "math.h"
#include "assert.h"

#include "Sampler.h"
#include "Igraph.h"

/***************************************************************************/
/*							Inference Graph Construction   				   */
/***************************************************************************/

IGraph::IGraph()
{
	m_mc_table = NULL;
}

IGraph::IGraph(int n_num,int e_num,int state_num,int node_size,int edge_size)
:Graph(n_num,e_num,node_size,edge_size)
{
	  // allocate memory
	InitGraph(n_num,e_num,state_num,node_size,edge_size);
}

IGraph::~IGraph()
{
	if(m_mc_table)
		MYFree(m_mc_table);
}

void IGraph::InitGraph(int n_num, int e_num, int state_num, int node_size, int edge_size)
{
	  // allocate memory
	Graph::InitGraph(n_num,e_num,node_size,edge_size);

// Init the shortcut
	m_state_num = state_num;
	inode_list = (IGraphNode *) node_list;
	iedge_list = (IGraphEdge *) edge_list;
}

void IGraph::ConstructPhi()
{
	int k;
	IGraphNode * inode;

	for(k=0;k<node_num;k++)
	{
		inode = (IGraphNode *)NODE_PTR(k);
		ComputePhi(k,inode->phi);
	}
}

void IGraph::ConstructKsai()
{
	int k;
	IGraphEdge * iedge;

	for(k=0;k<edge_num;k++)
	{
		iedge = (IGraphEdge *)EDGE_PTR(k);
		ComputeKsai(iedge->ni,iedge->nj,k,(double*)(iedge->ksai));
	}
}

void IGraph::CalculateCompatFunc()   
{
	ConstructPhi();
	ConstructKsai();
}

/***************************************************************************/
/*									Inference   						   */
/***************************************************************************/

// do inference in the graph
void IGraph::Inference()
{
	BeliefProp();
//	HOMeanFieldBP();
}

/***************************************************************************/
/*								General BP, LBP   						   */
/***************************************************************************/

// calculate 1 message passing from ni to nj
void IGraph::ComputeMessage(int edge_index, int ni,int nj)
{
	int i,j,k;
	double sum, prod, norm_sum;
	IGraphEdge * edge, *eij;
	int neighbor_n[10000];   // neighborhood nodes
	int neighb_n_num;
	short dir;

	neighb_n_num = GetNeighborNode(ni,neighbor_n);
	eij = IEDGE_PTR(edge_index); //GetEdgeFromNodes(ni,nj);

	norm_sum = 0;
	for(j=0;j<m_state_num;j++)   // mij(xj), loop over xj
	{
		sum = 0.0;
		for(i=0;i<m_state_num;i++)  // loop over xi
		{
			prod = 1.0;
			for(k=0;k<neighb_n_num;k++)   // loop over k E N(i)/j
			{
				if(neighbor_n[k]==nj)
					continue;

				edge = (IGraphEdge *) GetEdgeFromNodes(neighbor_n[k],ni,&dir);
				if(dir==1)
					prod *= edge->mij[i];
				else
					prod *= edge->mji[i];
			}
			sum += INODE_PTR(ni)->phi[i] * eij->ksai[i*m_state_num+j]*prod;
		}

		if(ni==eij->ni && nj==eij->nj)
		{
			eij->mij_p[j] = sum;   // update mij
			norm_sum+=sum;
		}
		else
		{
			assert(ni==eij->nj && nj==eij->ni);
			eij->mji_p[j] = sum;   // update mji*/
			norm_sum+=sum;
		}
	}

// mji
	for(j=0;j<m_state_num;j++)   // mij(xj), loop over xj
	{
		if(ni==eij->ni && nj==eij->nj)
		{
			eij->mij_p[j]/=norm_sum;   // update mij
		}
		else
		{
			assert(ni==eij->nj && nj==eij->ni);
			eij->mji_p[j]/=norm_sum;   // update mji*/
		}
	}
}

// standard belief prorogation
void IGraph::BeliefProp()   // belief prorogation
{
	int iter = 0, max_iter = 7;
	int i,k,ni,nj,j;
	double prod;
	int neighbor_n[10000];   // neighborhood nodes
	int neighb_n_num;
	IGraphEdge *eji;
	short dir;
	double sum,b;

// init all messages
	for(i=0;i<edge_num;i++)
	{
//		iedge_list[i].state_num = state_num;
		for(j=0;j<m_state_num;j++)
		{
//			iedge_list[i].mij[j] = 1;
//			iedge_list[i].mji[j] = 1;
			IEDGE_PTR(i)->mij[j] = 1;
			IEDGE_PTR(i)->mji[j] = 1;
		}
	}

	do
	{
// compute all messages
		for(k=0;k<edge_num;k++)
		{
// mij
			ni = IEDGE_PTR(k)->ni;
			nj = IEDGE_PTR(k)->nj;

			ComputeMessage(k,ni,nj);  // calculate mij
			ComputeMessage(k,nj,ni);  // calculate mji
		}

		iter ++;

// update all messages
		for(k=0;k<edge_num;k++)
		{
			memcpy(IEDGE_PTR(k)->mij,IEDGE_PTR(k)->mij_p,m_state_num*sizeof(double));
			memcpy(IEDGE_PTR(k)->mji,IEDGE_PTR(k)->mji_p,m_state_num*sizeof(double));
		}

	} while(iter < max_iter);

// compute all 1-node belief
	for(k=0;k<node_num;k++)
	{
		neighb_n_num = GetNeighborNode(k,neighbor_n);
		sum = 0;
		for(i=0;i<m_state_num;i++)
		{
			prod = 1.0;
			for(j=0;j<neighb_n_num;j++)
			{
				eji = (IGraphEdge *) GetEdgeFromNodes(neighbor_n[j],k,&dir);
				if(dir==1)
					prod *= eji->mij[i];
				else
					prod *= eji->mji[i];
			}
			INODE_PTR(k)->b1[i] = INODE_PTR(k)->phi[i]*prod;
			sum+= INODE_PTR(k)->b1[i];
		}

		assert(sum>0);
		for(i=0;i<m_state_num;i++)
		{
			INODE_PTR(k)->b1[i] /= sum;
		}
		
		b= INODE_PTR(k)->b1[1];
		INODE_PTR(k)->color = (int)(INODE_PTR(k)->b1[1]*255);
	}
}

// standard belief propagation with 2 node inference
void IGraph::BeliefProp2Node()
{
	int i,k,j,e,p,q;
	double prod;
	int neighbor_n[10000],neighbor_p[10000];   // neighborhood nodes
	int neighb_n_num,neighb_p_num;
	IGraphEdge *eji;
	short dir;
	double sum;

	BeliefProp();
	
// compute all 2-node belief
	for(e=0;e<edge_num;e++)
	{
		k = IEDGE_PTR(e)->nj;
		p = IEDGE_PTR(e)->ni;

		neighb_n_num = GetNeighborNode(k,neighbor_n);
		neighb_p_num = GetNeighborNode(p,neighbor_p);

		sum = 0;
		for(i=0;i<m_state_num;i++)
		{
			for(q=0;q<m_state_num;q++)
			{
				prod = 1.0;
				// k neighbor
				for(j=0;j<neighb_n_num;j++)
				{
					if(neighbor_n[j]==p)
						continue;

					eji = (IGraphEdge *) GetEdgeFromNodes(neighbor_n[j],k,&dir);
					if(dir==1)
						prod *= eji->mij[i];
					else
						prod *= eji->mji[i];
				}

				// p neighbor
				for(j=0;j<neighb_p_num;j++)
				{
					if(neighbor_p[j]==k)
						continue;

					eji = (IGraphEdge *) GetEdgeFromNodes(neighbor_p[j],p,&dir);
					if(dir==1)
						prod *= eji->mij[q];
					else
						prod *= eji->mji[q];
				}

				IEDGE_PTR(e)->b2[i*m_state_num+q] = IEDGE_PTR(e)->ksai[i*m_state_num+q]*INODE_PTR(k)->phi[i]*INODE_PTR(p)->phi[q]*prod;
				sum+= IEDGE_PTR(e)->b2[i*m_state_num+q];
			}
		}

//		assert(sum>0);
		if(sum>0)
		{
			for(i=0;i<m_state_num;i++)
				for(q=0;q<m_state_num;q++)
			{
				IEDGE_PTR(e)->b2[i*m_state_num+q] /= sum;
			}
		}
	}
}

/***************************************************************************/
/*							Gibbs Sampling						       	   */
/***************************************************************************/

void IGraph::GetConditional(int *state_vector, int i, double *prob_table)
{
	int en,* edge,j,k,s,ntmp,e;
	double sum,min;
	int s_debug[500], minus_inf[500];   // tag for minus infinity
	double p_debug[500],ksai_v;

	edge = (int*) MYMalloc(edge_num*sizeof(int));
	en = GetNeighborEdge(i,edge);

	for(j=0;j<node_num;j++)
		s_debug[j] = state_vector[j];

	min = 1.0e300;
	for(j=0;j<m_state_num;j++)
	{
		prob_table[j] = log(((IGraphNode *)NODE_PTR(i))->phi[j]);
		minus_inf[j] = 0;

		for(k=0;k<en;k++)
		{
			e = edge[k];
			ntmp = ((IGraphEdge *)EDGE_PTR(e))->ni;

			assert(i==ntmp || i==((IGraphEdge *)EDGE_PTR(e))->nj);

			if(i==ntmp) // i -- ni
			{
				s = state_vector[((IGraphEdge *)EDGE_PTR(e))->nj];
				ksai_v = ((IGraphEdge *)EDGE_PTR(e))->ksai[j*m_state_num+s];
			}
			else  // i -- nj
			{
				s = state_vector[((IGraphEdge *)EDGE_PTR(e))->ni];
				ksai_v = ((IGraphEdge *)EDGE_PTR(e))->ksai[s*m_state_num+j];
			}

//			prob_table[j] *= ((IGraphEdge *)EDGE_PTR(e))->ksai[j*m_state_num+s];
			
			if(ksai_v<1.0e-300)
			{
				prob_table[j] += -1.0e300;
				minus_inf[j] = 1;
			}
			else
				prob_table[j] += log(ksai_v);
		}

		if(prob_table[j]<min && minus_inf[j]!=1)											// skip zero
			min = prob_table[j];										 // non-zero minimum

		p_debug[j] = prob_table[j];
	}

	assert(min<1.0e300);

	for(j=0;j<m_state_num;j++)
		prob_table[j] = prob_table[j]-min;
	
	sum = 0;
	for(j=0;j<m_state_num;j++)
	{
		if(minus_inf[j]!=1)
		{
			prob_table[j] = exp(prob_table[j]);
			sum += prob_table[j];
		}
	}

	for(j=0;j<m_state_num;j++)
	{
		if(minus_inf[j]!=1)
			prob_table[j] /= sum;
		else
			prob_table[j] = 0;
	}

	MYFree(edge);
}

// Generate Monte Carlo samples : for inference and 
// partition function computation
void IGraph::MonteCarlo(int SampleNum)
{
	int i;
	
	if(m_mc_table)
		MYFree(m_mc_table);

	m_mc_table = (int *)MYMalloc(SampleNum*node_num*sizeof(int));
	m_mc_table_size = SampleNum;

// Use Gibbs Sampling
	for(i=0;i<SampleNum;i++)
		GibbsSampling(i);
}

// sample_index_th round Gibbs sampling
// si = sample index
void IGraph::GibbsSampling(int si)
{
	int i,j,state,s_sum;
	int s_debug[250];
	Sampler sampler;
	double * prob,p_debug[250];

	if(si==0)
	{
		GibbsSamplingInit(m_mc_table);  // Init GibbsSampling: finding the mode of the distribution
		return;
	}

//	prob = new double[m_state_num];
	prob = (double*)MYMalloc(m_state_num*sizeof(double));

// copy the previous sample
	for(i=0;i<node_num;i++)
		m_mc_table[si*node_num+i] = m_mc_table[(si-1)*node_num+i];

// iterate over nodes
	s_sum = 0;
	for(i=0;i<node_num;i++)
	{
	//  Get Conditional Distribution(prob,2000);
		GetConditional(&(m_mc_table[si*node_num+0]),i,prob); // get p(x_i|X/x_i);

		for(j=0;j<m_state_num;j++)  
			p_debug[j] = prob[j];

		if(m_state_num!=2)
		{
			state = sampler.DiscSampling1Step(prob,m_state_num);
			assert(state>=0);
		}
		else
		{
			if((double)rand()/RAND_MAX>prob[0])
				state = 1;
			else
				state = 0;
		}
		s_debug[i] = m_mc_table[si*node_num+i] = state;
		s_sum += s_debug[i];
	}

	MYFree(prob);
}

// Monte Carlo Inference
void IGraph::MCInference(int samplex)
{
	int sample_num,i,s,j,s2,s1,ni,nj;
	int begin = 0,k,a;
	double p,sai,b;

	sample_num = node_num*m_state_num*samplex;

	MonteCarlo(sample_num);

	for(j=0;j<node_num;j++)
		for(i=0;i<m_state_num;i++)
			((IGraphNode *)NODE_PTR(j))->b1[i] = 0;

	for(j=0;j<edge_num;j++)
		for(i=0;i<m_state_num;i++)
			for(k=0;k<m_state_num;k++)
				((IGraphEdge *)EDGE_PTR(j))->b2[i*m_state_num+k] = 0;

// 1-node belief
	begin =0;

	for(i=begin;i<sample_num;i++)
	{
		for(j=0;j<node_num;j++)
		{
			s = m_mc_table[i*node_num+j];
			((IGraphNode *)NODE_PTR(j))->b1[s]++;
		}
	}

	for(j=0;j<node_num;j++)
	{
		for(i=0;i<m_state_num;i++)
		{
			((IGraphNode *)NODE_PTR(j))->b1[i] /= (sample_num-begin);
		}
	}

#ifdef USE_ANA_2NODE_BLF
	Analytic2NodeBelief();
#else

// 2-node belief
	for(i=begin;i<sample_num;i++)
	{
		for(j=0;j<edge_num;j++)
		{
			ni = ((IGraphEdge *)EDGE_PTR(j))->ni;
			nj = ((IGraphEdge *)EDGE_PTR(j))->nj;

			s1 = m_mc_table[i*node_num+ni];
			s2 = m_mc_table[i*node_num+nj];

			((IGraphEdge *)EDGE_PTR(j))->b2[s1*m_state_num+s2]++;
		}
	}

	for(j=0;j<edge_num;j++)
	{
		p = 0;
		for(i=0;i<m_state_num;i++)
			for(k=0;k<m_state_num;k++)
			{
				((IGraphEdge *)EDGE_PTR(j))->b2[i*m_state_num+k]/=(sample_num-begin);
				b = ((IGraphEdge *)EDGE_PTR(j))->b2[i*m_state_num+k];
				if((sai = ((IGraphEdge *)EDGE_PTR(j))->ksai[i*m_state_num+k])<1.0e-300)
					a = 1;

				p += ((IGraphEdge *)EDGE_PTR(j))->b2[i*m_state_num+k];
			}

		assert(fabs(p-1.0)<1.0e-2);
	}
#endif

}

// compute two-node belief using one-node belief, according to Max Welling's Belief Optimization paper
// only for binary graph with Ising form
void IGraph::Analytic2NodeBelief()
{
	int j,ni,nj;
	double qi,qj,alpha_ij,Qij,ita_ij;
	double b1,b2,b3,b4;

	assert(m_state_num==2);

	for(j=0;j<edge_num;j++)
	{
		ni = EDGE_PTR(j)->ni;
		nj = EDGE_PTR(j)->nj;

		qi = ((IGraphNode *)NODE_PTR(ni))->b1[1];
		qj = ((IGraphNode *)NODE_PTR(nj))->b1[1];
		alpha_ij = ((IGraphEdge *)EDGE_PTR(j))->ksai[1*2+1];
		alpha_ij -= 1;

//		assert(((IGraphEdge *)EDGE_PTR(j))->ksai[0*2+0]==1.0);

		if(fabs(alpha_ij)>1.0e-300)
		{
			Qij = 1+alpha_ij*(qi+qj);
			ita_ij = (Qij-sqrt(Qij*Qij-4*alpha_ij*(1+alpha_ij)*qi*qj))/(2*alpha_ij);
		}
		else
		{
			ita_ij = qi*qj;
		}

		((IGraphEdge *)EDGE_PTR(j))->b2[3] = ita_ij;
//		assert(ita_ij>=0 && ita_ij<=1.0);

		((IGraphEdge *)EDGE_PTR(j))->b2[1*2+0] = qi-ita_ij>=0?qi-ita_ij:0;
		((IGraphEdge *)EDGE_PTR(j))->b2[0*2+1] = qj-ita_ij>=0?qj-ita_ij:0;
		((IGraphEdge *)EDGE_PTR(j))->b2[0] = 1-((IGraphEdge *)EDGE_PTR(j))->b2[0*2+1]
										-((IGraphEdge *)EDGE_PTR(j))->b2[1*2+0]-((IGraphEdge *)EDGE_PTR(j))->b2[3];

		if(((IGraphEdge *)EDGE_PTR(j))->b2[0]<0.0)
			((IGraphEdge *)EDGE_PTR(j))->b2[0] = 0.0;

		b1 = ((IGraphEdge *)EDGE_PTR(j))->b2[0];
		b2 = ((IGraphEdge *)EDGE_PTR(j))->b2[1];
		b3 = ((IGraphEdge *)EDGE_PTR(j))->b2[2];
		b4 = ((IGraphEdge *)EDGE_PTR(j))->b2[3];

//		assert(((IGraphEdge *)EDGE_PTR(j))->b2[0]>=0 && ((IGraphEdge *)EDGE_PTR(j))->b2[0]<=1.0);
	}
}

/***************************************************************************/
/*								MC max likelihood state approximation\	   	*/
/***************************************************************************/

// input -- state vector
// output -- log likelihood - log partition (constant)
double IGraph::GetLogLikelihood(int *state_vector)
{
	double loglike=0;
	int i,si,sj,ni,nj;

// phi
	for(i=0;i<node_num;i++)
	{
		si = state_vector[i];
		loglike += log(((IGraphNode *)NODE_PTR(i))->phi[si]);
	}

// ksai
	for(i=0;i<edge_num;i++)
	{
		ni = ((IGraphEdge *)EDGE_PTR(i))->ni;
		nj = ((IGraphEdge *)EDGE_PTR(i))->nj;

		si = state_vector[ni];
		sj = state_vector[nj];

		loglike += log(((IGraphEdge *)EDGE_PTR(i))->ksai[si*m_state_num+sj]);
	}

	return loglike;
}

// Monte Carlo max likelihood state approximation
// i.e. get maximum likelihood configuration/assignment
// input -- monte carlo table
// output -- state vector
void IGraph::MCMaxState(int samplex)
{
	int sample_num,i,s,j,s1,s2,k,ni,nj;
	double loglike,loglike_max=-1.0e38;
	int    loglike_ind=-1;
	int  * max_state;

	sample_num = node_num*m_state_num*samplex;
	sample_num = 40;

	MonteCarlo(sample_num);

	for(i=0;i<sample_num;i++)
	{
		loglike = GetLogLikelihood(&(m_mc_table[i*node_num]));

		if(loglike>loglike_max)
		{
			loglike_max = loglike;
			loglike_ind = i;
		}
	}

	max_state = &(m_mc_table[loglike_ind*node_num]);

// 1-node viterbi
	for(j=0;j<node_num;j++)
		for(i=0;i<m_state_num;i++)
			((IGraphNode *)NODE_PTR(j))->b1[i] = 0;

	for(j=0;j<node_num;j++)
	{
		s = max_state[j];
		((IGraphNode *)NODE_PTR(j))->b1[s] = 1.0;
	}

// 2-node viterbi
	for(j=0;j<edge_num;j++)
		for(i=0;i<m_state_num;i++)
			for(k=0;k<m_state_num;k++)
				((IGraphEdge *)EDGE_PTR(j))->b2[i*m_state_num+k] = 0;

	for(j=0;j<edge_num;j++)
	{
		ni = ((IGraphEdge *)EDGE_PTR(j))->ni;
		nj = ((IGraphEdge *)EDGE_PTR(j))->nj;

		s1 = max_state[ni];
		s2 = max_state[nj];

		((IGraphEdge *)EDGE_PTR(j))->b2[s1*m_state_num+s2] = 1.0;
	}
}
