/*****************************************************************************
// GraphMatcher.cpp : ARG or Bag-of-Parts similarity by MRF inference
// 
// Written by Dong-Qing Zhang, DVMM lab, Columbia University, 2005
//
// Contact : 
//			Dong-Qing Zhang :	 dongqing@gmail.com, 
//			Prof.Shih-Fu Chang : sfchang@ee.columbia.edu
*****************************************************************************/

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <assert.h>

#define min(a,b)  ((a)<(b)? (a):(b))
#define max(a,b)  ((a)>=(b)? (a):(b))

#include "GraphMatcher.h"
#include "Matrix.h"

#include "Configuration.h"

/***************************************************************************/
/*								1. Support Functions		   			   */
/***************************************************************************/
/***************************************************************************/
/*			2. Construction of the product MRF					    	   */
/***************************************************************************/
/***************************************************************************/
/*			3. Approximate Inference or Optimization				       */
/***************************************************************************/
/***************************************************************************/
/*							4. Similarity Measure calculation			   */
/***************************************************************************/
/***************************************************************************/
/*							   5. Learning								   */
/***************************************************************************/



/***************************************************************************/
/*							1.	Support Functions		   				   */
/***************************************************************************/

// Evaluate normal distribution with diagonal covariance
double MatchMRFGraph::SimpleNormal(double *x, double *m, double *v, int dim)
{
	int i;
	double ret = 1,a;

	for(i=0;i<dim;i++)
	{
		a = exp(-(x[i]-m[i])*(x[i]-m[i])/(2*v[i]))/(sqrt(2*3.1415926*v[i]));
		ret *= a;
	}

	return ret;
}

// switch two graphs, if necessary
void MatchMRFGraph::SwitchTwoGraphs(GARG **garg1, GARG **garg2)
{
	GARG * tmp;

	tmp = *garg1;
	*garg1 = *garg2;
	*garg2 = tmp;
}

void MatchMRFGraph::FreeMemory()
{
	if(edge_list)
	{
		MYFree(edge_list);  edge_list = NULL;
	}
	if(node_list)
	{
		MYFree(node_list);	node_list = NULL;
	}
}

MatchMRFGraph::MatchMRFGraph()
{
	m_LearningMode = 0;
	m_LearningSign = 1;
	m_InitFlag = 1;

	mnode_list = NULL;
	medge_list = NULL;
	edge_list = NULL;
	node_list = NULL;

	m_MatchPara.ForegroundVar1n = NULL;
	m_MatchPara.BackgroundMean1n = NULL;
	m_MatchPara.BackgroundVar1n  = NULL;

	m_MatchPara.ForegroundVar2n  = NULL;
	m_MatchPara.BackgroundMean2n  = NULL;
	m_MatchPara.BackgroundVar2n  = NULL;
}

MatchMRFGraph::~MatchMRFGraph()
{
	FreeMemory();

	MYFree(m_MatchPara.ForegroundVar1n);
	MYFree(m_MatchPara.BackgroundMean1n);
	MYFree(m_MatchPara.BackgroundVar1n);

	MYFree(m_MatchPara.ForegroundVar2n);
	MYFree(m_MatchPara.BackgroundMean2n);
	MYFree(m_MatchPara.BackgroundVar2n);

	m_MatchPara.ForegroundVar1n = 0;
	m_MatchPara.BackgroundMean1n = 0;
	m_MatchPara.BackgroundVar1n = 0;

	m_MatchPara.ForegroundVar2n = 0;
	m_MatchPara.BackgroundMean2n = 0;
	m_MatchPara.BackgroundVar2n = 0;
}

/***************************************************************************/
/*		2. Construction of the association graph (product graph)   		   */
/***************************************************************************/

// form product graph for graph matching
// Garg1, Garg2: source and target graph
double MatchMRFGraph::FormAdjGraph(GARG * garg1,GARG * garg2,int NodeBound)
{
	int i,j,n_num,e_num,n,n_num1,n_num2,be_num;
	int part_size[2],i1,i2,j1,j2,e1,e2;
	short dir, * bip_aff=NULL;
//	GARG * tmp;

/*
	if(garg2->node_num < garg1->node_num)
	{
		tmp = garg2;
		garg2 = garg1;
		garg1 = tmp;
	}
*/
	m_Garg1 = garg1;  m_Garg2 = garg2;
	n_num1 = m_Garg1->node_num; n_num2 = m_Garg2->node_num;

	be_num = NodeBound;  // the limit of the node number
	bip_aff = Prefiltering(garg1,garg2,&be_num); // Prune MRF node to speed up

	if(be_num>NodeBound)
	{
		if(bip_aff)
			MYFree(bip_aff);
		return 1.0;
	}

	part_size[0] = n_num1;	part_size[1] = n_num2;

	n_num = be_num; 
	e_num = n_num*(n_num-1)/2;

	InitGraph(n_num,e_num,2,sizeof(MRFNode),sizeof(MRFEdge));

	mnode_list = (MRFNode *) node_list;
	medge_list = (MRFEdge *) edge_list;

	n = 0;
	for(i=0;i<n_num1;i++)
	{
		for(j=0;j<n_num2;j++)
		{
			if(bip_aff[i*n_num2+j]==0)
				continue;

			mnode_list[n].feat_vec_int[0] = i;	// ARG node i
			mnode_list[n].feat_vec_int[1] = j;	// ARG node u
			n++;
		}
	}

	n = 0;
	for(i=0;i<n_num;i++)
	{
		for(j=i+1;j<n_num;j++)
		{
			medge_list[n].ni = i;
			medge_list[n].nj = j;
			medge_list[n].weight = 1;

			i1 = mnode_list[i].feat_vec_int[0];  // node i in the arg1
			j1 = mnode_list[i].feat_vec_int[1];  // node u in the arg2
			i2 = mnode_list[j].feat_vec_int[0];  // node j in the arg1
			j2 = mnode_list[j].feat_vec_int[1];  // node v in the arg2

// to assign the two graph edges
			e1 = garg1->GetEdgeFromNodesInt(i1,i2,&dir);

			medge_list[n].feat_vec_int[0] = e1;	  // ARG edge 1
			medge_list[n].feat_vec_int[1] = dir;
			assert(dir>=0);   // dir1 should be always larger than 0

			e2 = garg2->GetEdgeFromNodesInt(j1,j2,&dir);

			medge_list[n].feat_vec_int[2] = e2;	  // ARG edge 2
			medge_list[n].feat_vec_int[3] = dir;

			medge_list[n].feat_vec_int[4] = i1;	  // node i in arg1
			medge_list[n].feat_vec_int[5] = i2;	  // node j in arg1
			medge_list[n].feat_vec_int[6] = j1;	  // node u in arg2
			medge_list[n].feat_vec_int[7] = j2;   // node v in arg2

			n++;
		}
	}

	if(bip_aff)
		MYFree(bip_aff);

	return 1;
}

// prefitlering using spatial locations
// need to change for other applications
int PrefilterCompare(double *feat1, double *feat2, double relax)
{
	int i;
	double dist_pos=0;
	
	for(i=0;i<2;i++)   // spatial distance
		dist_pos += (feat1[i]-feat2[i])*(feat1[i]-feat2[i]);
	dist_pos = sqrt(dist_pos);

	if(dist_pos > relax*g_SpatialPrefilterThreshold)   // use distance threshold 80, needs to change if images are of different sizes
		return 0;
	else
		return 1;
}

// prefilter the Adj graph node by sorting to speed up MRF inference
short * MatchMRFGraph::Prefiltering(GARG *garg1, GARG *garg2,int * edge_num)
{
	int nn1,nn2;
	short * aff;
	double *like,like_max; // likelihood
	double LikeThreshold = -15;
	double * garg_like;
	int i,j,max_ind,n=0;
	int node_bound; 

	int node_feat_dim;
	double c1[50],c2[50],phi[2];

	node_bound = *edge_num;
	nn1 = garg1->node_num;
	nn2 = garg2->node_num;

	like = (double *)MYMalloc(nn1*nn2*sizeof(double));
	aff = (short *)MYMalloc(nn1*nn2*sizeof(short));
	memset(aff,0,nn1*nn2*sizeof(short));

// compute 1-node likelihood ratio
	n = 0;
	for(i=0;i<nn1;i++)
	{
		node_feat_dim = m_Garg1->GetNodeFeature(c1,i);

		for(j=0;j<nn2;j++)
		{
			node_feat_dim = m_Garg2->GetNodeFeature(c2,j);

// phi 0
			phi[0] = SimpleNormal(c2,m_MatchPara.BackgroundMean1n,
					m_MatchPara.BackgroundVar1n,node_feat_dim);

// phi 1
			phi[1] = SimpleNormal(c2,c1,
					m_MatchPara.ForegroundVar1n,node_feat_dim);

			if(PrefilterCompare(c1,c2,1)==1)			 // spatial filtering
				like[i*nn2+j] = log(phi[1])-log(phi[0]);
			else
				like[i*nn2+j] = -1.0e30;

			if(like[i*nn2+j]<LikeThreshold)	// thresholding the likelihood bias
				like[i*nn2+j] = -1.0e30;

			if(like[i*nn2+j]>-1.0e29)		
				n++;						
		}
	}

	node_bound = (n>node_bound)? node_bound:n;

	garg_like = (double *)MYMalloc(nn2*sizeof(double));
	for(i=0;i<nn2;i++)
		garg_like[i] = -1.0e100;

// naive sorting
	for(i=0;i<node_bound;i++)
	{
		like_max = -1.0e30;
		max_ind = -1;

		for(j=0;j<nn1*nn2;j++)
		{
			if(like[j]>like_max)
			{
				like_max = like[j];
				max_ind = j;
			}
		}

		if(max_ind<0)
			break;
//		assert(max_ind>=0);

		aff[max_ind] = 1;
		like[max_ind] = -1.0e30;

		if(garg_like[max_ind%nn2]<-1.0e99)
			garg_like[max_ind%nn2] = like_max;
	}

	*edge_num = node_bound;

	MYFree(like);
	MYFree(garg_like);

	return aff;
}

/*
// Very rough approximation of z
double MatchMRFGraph::ComputePrior_z()
{
	double r,C=1;  // Using low bound, can be set to other values, would not affect results much

	r = m_MatchPara.r;
	return r/((1-r)*C);
}
*/

// approximation of z according to r
double MatchMRFGraph::ComputePrior_z()
{
	double C = max(m_Garg1->node_num,m_Garg2->node_num);  // using the upper bound
	double r = m_MatchPara.r;

	return r/((1-r)*C);
}

void MatchMRFGraph::ComputePhi(int nodei,double * phi)
{
//	PRGaussian g0,g1;
	int node_feat_dim,ri,rj;
	double c1[50],c2[50],z;

	ri = mnode_list[nodei].feat_vec_int[0];  // GARG node i
	rj = mnode_list[nodei].feat_vec_int[1];  // GARG node j

	node_feat_dim = m_Garg1->GetNodeFeature(c1,ri);
	node_feat_dim = m_Garg2->GetNodeFeature(c2,rj);

// phi 0
	phi[0] = SimpleNormal(c2,m_MatchPara.BackgroundMean1n,
					m_MatchPara.BackgroundVar1n,node_feat_dim);

// phi 1
	z = ComputePrior_z();
	phi[1] = SimpleNormal(c2,c1,
					m_MatchPara.ForegroundVar1n,node_feat_dim)*z;

	if(phi[0]<1.0e-300)
		assert(0);

	phi[1] /= phi[0];
	phi[0] = 1;
}

// compute 2-node potential
void MatchMRFGraph::ComputeKsai(int nodei,int nodej,int edgei, double * ksai)
{
	int * feat_int1,*feat_int2;
	double /*feat,*/e_ksai[4];  // emission ksai
	int ri,rj,dir1,dir2,edge_feat_dim,ni,nj,nu,nv;
	double c1[50],c2[50];

	ri = medge_list[edgei].feat_vec_int[0];
	dir1 = medge_list[edgei].feat_vec_int[1];
	rj = medge_list[edgei].feat_vec_int[2];
	dir2 = medge_list[edgei].feat_vec_int[3];

	ni = medge_list[edgei].feat_vec_int[4];
	nj = medge_list[edgei].feat_vec_int[5];
	nu = medge_list[edgei].feat_vec_int[6];
	nv = medge_list[edgei].feat_vec_int[7];

	feat_int1 = mnode_list[nodei].feat_vec_int;
	feat_int2 = mnode_list[nodej].feat_vec_int;

	if(ri==-1 || rj==-1)   // the edge in G1 is mapped to one node in G2, there's no definition of 2-node feature.
	{
		e_ksai[3] = 1.0;
		e_ksai[2] = e_ksai[1] = e_ksai[0] = 1.0;
	}
	else
	{
		if(g_Representation!=0)
		{
			edge_feat_dim = m_Garg1->GetEdgeFeature(c1,ri,dir1);
			edge_feat_dim = m_Garg2->GetEdgeFeature(c2,rj,dir2);

			e_ksai[3] = SimpleNormal(c2,c1,m_MatchPara.ForegroundVar2n,edge_feat_dim);
			e_ksai[0] = SimpleNormal(c2,m_MatchPara.BackgroundMean2n,
													m_MatchPara.BackgroundVar2n,edge_feat_dim);

			e_ksai[3] /= e_ksai[0];
			e_ksai[0] = e_ksai[2] = e_ksai[1] = 1.0;
		}
	}

	if(feat_int1[0]==feat_int2[0] || feat_int1[1]==feat_int2[1])
	{
		ksai[0] = 1;		ksai[1] = 1;
		ksai[2] = 1;		

		if(g_InferenceAlgorithm!=0)
			ksai[3] = 0.00000001;   // for BP stability
		else
			ksai[3] = 0;
	}
	else
	{
		ksai[0] =1;		ksai[1] =1;
		ksai[2] =1;		ksai[3] =1;
	}

	if(g_Representation!=0)   // use bag-of-parts or ARG 
	{
		ksai[0]*= e_ksai[0];	ksai[1]*= e_ksai[1];
		ksai[2]*= e_ksai[2];	ksai[3]*= e_ksai[3];
	}
}

/***************************************************************************/
/*				3. Approximate Inference or Optimization				   */
/***************************************************************************/

void MatchMRFGraph::Inference(int samplex)
{
	if(g_InferenceAlgorithm==0)
		MCInference(samplex);
	else if(g_InferenceAlgorithm==1)
		BeliefProp();
	else
		MeanFiledBP();
}

// Depth-First-Search to get the maximum likelihood assignment of matching
// used by Gibbs sampling
void MatchMRFGraph::DFSSearch(int * state_seq)
{
	SortMachine sm;
	int *node_fill, *state_fill,i,j,n,nid,sid;

	memset(state_seq,0,sizeof(int)*node_num);
	sm.Init(node_num*m_state_num);
	for(i=0;i<node_num;i++)
		for(j=0;j<m_state_num;j++)
	{
		sm.m_Arr[i*m_state_num+j].ind = i;
		sm.m_Arr[i*m_state_num+j].ind2 = j;
		sm.m_Arr[i*m_state_num+j].d = ((IGraphNode *)NODE_PTR(i))->phi[j];
	}
	sm.Sort();

// occupation array
	node_fill = new int[node_num];
	memset(node_fill,0,node_num*sizeof(int));
	state_fill = new int[m_state_num];
	memset(state_fill,0,m_state_num*sizeof(int));

	state_seq[sm.m_Arr[0].ind] = sm.m_Arr[0].ind2;
	node_fill[sm.m_Arr[0].ind] = 1; state_fill[sm.m_Arr[0].ind2] = 1;

	n = 1;
	while(n<node_num*m_state_num)
	{
		nid = sm.m_Arr[n].ind;
		sid = sm.m_Arr[n++].ind2;
		if(node_fill[nid]>0 || state_fill[sid]>0)
			continue;

		state_seq[nid] = sid;
		node_fill[nid] = 1;
		state_fill[sid] = 1;
	}

	delete state_fill;
	delete node_fill;
}

// this is a virtual function
void MatchMRFGraph::GibbsSamplingInit(int * state_seq)
{
	DFSSearch(state_seq);
}

// calculate 1 message passing from ni to nj
void MatchMRFGraph::ComputeMessage(int edge_index, int ni,int nj)
{
	int i,j;
	double sum, prod, norm_sum;
	IGraphEdge *eij;
	int neighbor_n[10000];   // neighborhood nodes
	int neighb_n_num;

	neighb_n_num = GetNeighborNode(ni,neighbor_n);
	eij = IEDGE_PTR(edge_index); //GetEdgeFromNodes(ni,nj);

	norm_sum = 0;
	for(j=0;j<m_state_num;j++)   // mij(xj), loop over xj
	{
		sum = 0.0;
		for(i=0;i<m_state_num;i++)  // loop over xi
		{
			if(ni==eij->ni && nj==eij->nj)
				prod = eij->mji[i];
			else
				prod = eij->mij[i];

			assert(fabs(prod)>1.0e-100);
			sum += INODE_PTR(ni)->phi[i] *
					eij->ksai[i*m_state_num+j] * m_M[ni*m_state_num+i]/prod;
		}

		assert(sum<1.0e300);

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

	assert(norm_sum<1.0e300);

// mji
	for(j=0;j<m_state_num;j++)   // mij(xj), loop over xj
	{
		if(ni==eij->ni && nj==eij->nj)
		{
			assert(fabs(norm_sum)>1.0e-400);
			eij->mij_p[j]/=norm_sum;   // update mij
		}
		else
		{
			assert(fabs(norm_sum)>1.0e-400);
			assert(ni==eij->nj && nj==eij->ni);
			eij->mji_p[j]/=norm_sum;   // update mji*/
		}
	}
}

// Modified Belief Propagation for complete-graph MRF
void MatchMRFGraph::BeliefProp()
{
	int iter = 0, max_iter = 7;
	// max_iter is typically from 3 to 7
	int i,k,ni,nj,j;
	int neighbor_n[10000];   // neighborhood nodes
	int neighb_n_num;
	double sum,b,a,c;

// init all messages
	for(i=0;i<edge_num;i++)
	{
		for(j=0;j<m_state_num;j++)
		{
			IEDGE_PTR(i)->mij[j] = 1;
			IEDGE_PTR(i)->mji[j] = 1;
		}
	}

	m_M = (double *) MYMalloc(node_num*m_state_num*sizeof(double));
	for(i=0;i<node_num*m_state_num;i++)  m_M[i] = 1.0;

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
			for(i=0;i<m_state_num;i++)
			{
				a = IEDGE_PTR(k)->mij[i];
				c = IEDGE_PTR(k)->mij_p[i];
				IEDGE_PTR(k)->mij[i] = IEDGE_PTR(k)->mij_p[i];  
				IEDGE_PTR(k)->mji[i] = IEDGE_PTR(k)->mji_p[i];
			}
		}

// update M
		for(i=0;i<node_num*m_state_num;i++)  
			m_M[i] = 0;

		for(k=0;k<edge_num;k++)
		{
			ni = IEDGE_PTR(k)->ni;
			nj = IEDGE_PTR(k)->nj;

			for(i=0;i<m_state_num;i++)
			{
				// to nj
				m_M[nj*m_state_num+i] += log(IEDGE_PTR(k)->mij[i]);
				// to ni
				m_M[ni*m_state_num+i] += log(IEDGE_PTR(k)->mji[i]);
			}
		}

		for(i=0;i<node_num*m_state_num;i++)  
			m_M[i] = exp(m_M[i]);
	} while(iter < max_iter);

// compute all 1-node belief
	for(k=0;k<node_num;k++)
	{
		neighb_n_num = GetNeighborNode(k,neighbor_n);

		sum = 0;
		for(i=0;i<m_state_num;i++)
		{
			INODE_PTR(k)->b1[i] = INODE_PTR(k)->phi[i]*m_M[k*m_state_num+i];
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

// Disable 2Node belief
	Compute2NodeBelief();
	MYFree(m_M);
}

void MatchMRFGraph::Compute2NodeBelief()
{
	int i,k,e,p,q;
	double prod,div; //,prod2;
	int neighbor_n[10000],neighbor_p[10000];   // neighborhood nodes
	int neighb_n_num,neighb_p_num;
	IGraphEdge *eji;
	short dir;
	double sum,b2;

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

				eji = (IGraphEdge *)GetEdgeFromNodes(p,k,&dir);
				if(dir==1)
					div = eji->mij[i];
				else
					div = eji->mji[i];

				prod *= m_M[k*m_state_num+i]/div;

				eji = (IGraphEdge *)GetEdgeFromNodes(k,p,&dir);
				if(dir==1)
					div = eji->mij[q];
				else
					div = eji->mji[q];
				prod *= m_M[p*m_state_num+q]/div;

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
				b2 = IEDGE_PTR(e)->b2[i*m_state_num+q];
			}
		}
	}	
}

// Mean-Field Inference
void MatchMRFGraph::MeanFiledBP()
{
	int iter=0, max_iter = 4,n,nk,ei,e,i,j;
	int k, result_e[10000];
	double bsum,esum,a;	

// initialize all belief	
	for(k=0;k<node_num;k++)
		for(i=0;i<m_state_num;i++)
			INODE_PTR(k)->b1[i] = 1.0/m_state_num;
	do
	{
// compute all belief
		for(k=0;k<node_num;k++)
		{
			bsum = 0.0;
			for(i=0;i<m_state_num;i++)
			{
// E_phi
				esum = 0.0;
				esum += log(INODE_PTR(k)->phi[i]);

// E_ksai, edge
				n = GetNeighborEdge(k,result_e);
				for(ei=0;ei<n;ei++)
				{
					e = result_e[ei];
					if(IEDGE_PTR(e)->ni == k)
						nk = IEDGE_PTR(e)->nj;
					else
						nk = IEDGE_PTR(e)->ni;
					for(j=0;j<m_state_num;j++)
					{
						esum += INODE_PTR(nk)->b1[j]*log(IEDGE_PTR(e)->ksai[i*m_state_num+j]);
					}
				}

				a = INODE_PTR(k)->b1_new[i] = exp(esum);
				bsum += INODE_PTR(k)->b1_new[i];
			} 

			assert(bsum>0);

// normalize all belief
			for(i=0;i<m_state_num;i++)
				a = INODE_PTR(k)->b1_new[i] /= bsum;
		}
		iter++;

// update all belief
		for(k=0;k<node_num;k++)
		{
			for(i=0;i<m_state_num;i++)
			{
				a = INODE_PTR(k)->b1[i] = INODE_PTR(k)->b1_new[i];
			}
		}
	} while(iter < max_iter);

	for(k=0;k<node_num;k++)
		a = INODE_PTR(k)->color = (int)(INODE_PTR(k)->b1[1]*255);
}

/***************************************************************************/
/*						4. Similarity Measure Calculation				   */
/***************************************************************************/

// compute entropy for likelihood computation with approximation
double MatchMRFGraph::ComputeEntropy()
{
	int neighb_n_num, neighbor_n[10000];   // neighborhood nodes
	int k,e,i,j;
	double b, entropy,en;

	entropy = 0;
	for(k=0;k<node_num;k++)
	{
		neighb_n_num = GetNeighborNode(k,neighbor_n);

		en = 0;
		for(j=0;j<m_state_num;j++)
		{
			b= INODE_PTR(k)->b1[j];
			if(b>0)
				en+= b*log(b);
		}
		entropy += en*(neighb_n_num-1);
	}

	for(e=0;e<edge_num;e++)
	{
		en = 0;
		for(i=0;i<m_state_num;i++)
			for(j=0;j<m_state_num;j++)
		{
			b = IEDGE_PTR(e)->b2[i*m_state_num+j];
			if(b>0)
				en+= b*log(b);
		}

		entropy -= en;
	}

//	assert(entropy>=0);

	return entropy;
}

// compute the average energy
double MatchMRFGraph::ComputeAverageEnergy()
{
	int k,e,i,j;
	double b, aver_en1,aver_en2,ae,a,c,tmp,like;

	aver_en1 = 0;
	for(k=0;k<node_num;k++)
	{
		ae = 0;

		if((b=INODE_PTR(k)->b1[1])>1.0e-10)  // for debug
		{
			a = log(INODE_PTR(k)->phi[1]);
			c = log(INODE_PTR(k)->phi[0]);
		}

		for(j=0;j<m_state_num;j++)
		{
			b = INODE_PTR(k)->b1[j];
			if(b<1.0e-300)
				continue;

			like = log(INODE_PTR(k)->phi[j]);
			ae+= b*log(INODE_PTR(k)->phi[j]);
		}

		aver_en1 += ae;
	}

	aver_en2 = 0;

	for(e=0;e<edge_num;e++)
	{
		ae = 0;
		for(i=0;i<m_state_num;i++)
			for(j=0;j<m_state_num;j++)
		{
			b = IEDGE_PTR(e)->b2[i*m_state_num+j];
			tmp = IEDGE_PTR(e)->ksai[i*m_state_num+j];  // for debug

			if(b<1.0e-300)
				continue;

			if(IEDGE_PTR(e)->ksai[i*m_state_num+j]>0)
				ae+= b*log(IEDGE_PTR(e)->ksai[i*m_state_num+j]);
			else
			{
				tmp = IEDGE_PTR(e)->ksai[i*m_state_num+j];
			}
		}
		aver_en2 += ae;
	}

	return aver_en1+aver_en2;
}

#define PI 3.1415926
float MatchMRFGraph::StructureDiff()
{
	double Spos = m_MatchPara.Spos, Mneg = m_MatchPara.Mneg, Sneg = m_MatchPara.Sneg;
	int N = m_Garg1->node_num, M = m_Garg2->node_num;

	float struct_diff = log(Sneg*sqrt(2*PI))-log(Spos*sqrt(2*PI))-(N-M)*(N-M)/(2*Spos*Spos)+(M-Mneg)*(M-Mneg)/(2*Sneg*Sneg);

	return struct_diff;
}

// rough approximation of LogZ
double MatchMRFGraph::ComputeLogZ()
{
	int N = min(m_Garg1->node_num,m_Garg1->node_num);
	double r = m_MatchPara.r;

	return -log(1-r)*N;
}

// LogZp
double MatchMRFGraph::ComputeLogZp()
{
	double ae,en;

	ae = ComputeAverageEnergy();
	en = ComputeEntropy();

	return ae+en;
}

// Similarity from log likelihood ratio 
double MatchMRFGraph::LogLikelihoodRatio()
{
	double LogZp,LogZ;

	LogZp = ComputeLogZp();
	LogZ = ComputeLogZ();

	float sdif = StructureDiff();

// normalization to define size-invariant similarity
//   return (LogZp-LogZ+sdif)/sqrt(m_Garg1->node_num*m_Garg2->node_num);
	if(g_Normalization==0)
		return (LogZp-LogZ+sdif)/min(m_Garg1->node_num,m_Garg2->node_num);
	else if(g_Normalization==1)
		return (LogZp-LogZ+sdif)/sqrt(m_Garg1->node_num*m_Garg2->node_num);
	else if(g_Normalization==2)
		return (LogZp-LogZ+sdif)/(m_Garg1->node_num*m_Garg2->node_num);
	else 
		return (LogZp-LogZ+sdif)/max(m_Garg1->node_num,m_Garg2->node_num);
}

/***************************************************************************/
/*							   5. Learning								   */
/***************************************************************************/

void MatchMRFGraph::LoadParameters(const char * ParameterDir)
{
	char filename[2000];
	FILE * fp;
	int n_feat_d,e_feat_d,j,id;
	int ret;
	char str[1000];
	float f,f2,f3;

// read topological parameters
	
	sprintf(filename,"%s%s",ParameterDir,"para.txt");
	fp = fopen(filename,"r");
	fscanf(fp,"%d",&n_feat_d);
	fscanf(fp,"%d",&e_feat_d);
	fclose(fp);

	m_MatchPara.nFeatDim = n_feat_d;
	m_MatchPara.eFeatDim = e_feat_d;

// transfer prob
	sprintf(filename,"%s%s",ParameterDir,"r.txt");
	fp = fopen(filename,"r");
	fscanf(fp,"%f",&(f));
	m_MatchPara.r = f;
	fclose(fp);

// structure change prob
	sprintf(filename,"%s%s",ParameterDir,"struct.txt");
	fp = fopen(filename,"r");
	fscanf(fp,"%f %f %f",&f,&f2,&f3);
	m_MatchPara.Spos = f;
	m_MatchPara.Mneg = f2;
	m_MatchPara.Sneg = f3;
	fclose(fp);

// foreground var 1n
	if(m_MatchPara.ForegroundVar1n==NULL)
		m_MatchPara.ForegroundVar1n = (double *)MYMalloc(n_feat_d*sizeof(double));

	sprintf(filename,"%sV0000.var.txt",ParameterDir);
	fp = fopen(filename,"r");
	fscanf(fp,"%d",&id);
	for(j=0;j<n_feat_d;j++)
	{
		ret = fscanf(fp,"%s",str);
		m_MatchPara.ForegroundVar1n[j] = atof(str);
	}
	fclose(fp);

// background mean 1n
	if(m_MatchPara.BackgroundMean1n==NULL)
		m_MatchPara.BackgroundMean1n = (double *)MYMalloc(n_feat_d*sizeof(double));

	sprintf(filename,"%sV9999.mean.txt",ParameterDir);
	fp = fopen(filename,"r");
	fscanf(fp,"%d",&id);
	for(j=0;j<n_feat_d;j++)
	{
		ret = fscanf(fp,"%s",str);
		m_MatchPara.BackgroundMean1n[j] = atof(str);
	}
	fclose(fp);

// background var 1n
	if(m_MatchPara.BackgroundVar1n==NULL)
		m_MatchPara.BackgroundVar1n = (double *)MYMalloc(n_feat_d*sizeof(double));

	sprintf(filename,"%sV9999.var.txt",ParameterDir);
	fp = fopen(filename,"r");
	fscanf(fp,"%d",&id);
	for(j=0;j<n_feat_d;j++)
	{
		ret = fscanf(fp,"%s",str);
		m_MatchPara.BackgroundVar1n[j] = atof(str);
	}
	fclose(fp);

//  read the attribute parameters
//  of the edges
// foreground var 2n
	if(g_Representation!=0)
	{
		if(m_MatchPara.ForegroundVar2n==NULL)
			m_MatchPara.ForegroundVar2n = (double *)MYMalloc(n_feat_d*sizeof(double));

		sprintf(filename,"%sE0000.var.txt",ParameterDir);
		fp = fopen(filename,"r");
		fscanf(fp,"%d",&id);
		fscanf(fp,"%d",&id);
		for(j=0;j<e_feat_d;j++)
		{
			ret = fscanf(fp,"%s",str);
			m_MatchPara.ForegroundVar2n[j] = atof(str);
		}
		fclose(fp);

	// background mean 2n
		if(m_MatchPara.BackgroundMean2n==NULL)
			m_MatchPara.BackgroundMean2n = (double *)MYMalloc(n_feat_d*sizeof(double));

		sprintf(filename,"%sE9999.mean.txt",ParameterDir);
		fp = fopen(filename,"r");
		fscanf(fp,"%d",&id);
		fscanf(fp,"%d",&id);

		for(j=0;j<e_feat_d;j++)
		{
			ret = fscanf(fp,"%s",str);
			m_MatchPara.BackgroundMean2n[j] = atof(str);
		}
		fclose(fp);

	// background var 2n
		if(m_MatchPara.BackgroundVar2n==NULL)
			m_MatchPara.BackgroundVar2n = (double *)MYMalloc(n_feat_d*sizeof(double));

		sprintf(filename,"%sE9999.var.txt",ParameterDir);
		fp = fopen(filename,"r");
		fscanf(fp,"%d",&id);
		fscanf(fp,"%d",&id);

		for(j=0;j<e_feat_d;j++)
		{
			ret = fscanf(fp,"%s",str);
			m_MatchPara.BackgroundVar2n[j] = atof(str);
		}
		fclose(fp);
	}
}

void MatchMRFGraph::SaveBeliefInit(const char *NodeBlfFile1, const char *NodeBlfFile2,
											GARG *garg1, GARG *garg2)
{
	int i,j;
	FILE * fp_xiu;
//	GARG * tmp;

/*	if(garg2->node_num < garg1->node_num)
	{
		tmp = garg2;
		garg2 = garg1;
		garg1 = tmp;
	}
*/

	fp_xiu = fopen(NodeBlfFile1,"w");

// 1-node belief
	for(i=0;i<garg1->node_num;i++)
		for(j=0;j<garg2->node_num;j++)
	{
		if((float)rand()/RAND_MAX<0.02)
			fprintf(fp_xiu,"%d %d %f\n",i,j,1.0);
	}
	
	fclose(fp_xiu);
}

// save the correspondence results
void MatchMRFGraph::SaveBelief(const char *NodeBlfFile1,const char *NodeBlfFile2)
{
	int k,e,i1,i2,j1,j2,ri,rj,dir1,dir2;
	FILE * fp_xiu,*fp_xiujv;

	fp_xiu = fopen(NodeBlfFile1,"w");
	fp_xiujv = fopen(NodeBlfFile2,"w");

// 1-node belief
	for(k=0;k<node_num;k++)
	{
		fprintf(fp_xiu,"%d %d %f\n",mnode_list[k].feat_vec_int[0],
							mnode_list[k].feat_vec_int[1],mnode_list[k].b1[1]);
	}
	
// 2-node belief
	for(e=0;e<edge_num;e++)
	{
		ri = medge_list[e].feat_vec_int[0];
		dir1 = medge_list[e].feat_vec_int[1];
		rj = medge_list[e].feat_vec_int[2];
		dir2 = medge_list[e].feat_vec_int[3];
		i1 = medge_list[e].feat_vec_int[4];  i2 = medge_list[e].feat_vec_int[5];
		j1 = medge_list[e].feat_vec_int[6];  j2 = medge_list[e].feat_vec_int[7];

		fprintf(fp_xiujv,"%d %d %d %d %d %d %d %d %f\n",ri,dir1,rj,dir2,i1,i2,j1,j2,medge_list[e].b2[1*m_state_num+1]);
	}

	fclose(fp_xiujv);
	fclose(fp_xiu);
}
