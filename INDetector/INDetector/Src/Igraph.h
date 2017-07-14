/***********************************************************************************
// IGraph.h :   class for performing probabilistic inference, including LBP 
//				and Gibbs Sampling
//            
// 
// Written by Dong-Qing Zhang, DVMM lab, Columbia University, 2005
//
// Contact : 
//			Dong-Qing Zhang :	 dongqing@gmail.com, 
//			Prof.Shih-Fu Chang : sfchang@ee.columbia.edu
************************************************************************************/

#ifndef IGRAPH
#define IGRAPH

#include "Graph.h"

/********************************************************************************/
// inference graph state
#define MAX_INST	2  // maximum state number for binary MRF

// IGraphNode
class IGraphNode : public Node
{
public:
	double   phi[MAX_INST];   
	double   b1[MAX_INST];   
	double   b1_new[MAX_INST];

	short  * n_e_n;        // neighboring edges and nodes index
};

// IGraphEdge , used to perform message passing
class IGraphEdge : public Edge
{
public:
	double    ksai[MAX_INST*MAX_INST];  

	double    mij[MAX_INST];   // mij(xi)  , message from i to j
	double    mji[MAX_INST];   // mji(xi)  , message from j to i

	double    mij_p[MAX_INST];   // mij(xi)  plus, the updated message from i to j
	double    mji_p[MAX_INST];   // mji(xi)  plus, the updated message from j to i
	double    b2[MAX_INST*MAX_INST]; // two node belief
};

class IGraph : public Graph
{
	void ComputeMessage(int edge_index, int ni,int nj); // compute message passing from i to j

public:
	int m_state_num;
	int * m_mc_table;		 // Monte Carlo table
	int   m_mc_table_size;   // the row number of the state vector
	int   m_mc_table_begin;  // the beginning of the table

	IGraphNode * inode_list;      // short-cut to the node list
	IGraphEdge * iedge_list;      // short-cut to the edge list

public:
// inference
	void Inference();

//  General BP
	virtual void BeliefProp();   // belief prorogation
	void BeliefProp2Node();

//  Monte Carlo methods
	double GetLogLikelihood(int * state_vector);
	void MCMaxState(int samplex=40); // max likelihood state approximation

	void MCInference(int sample=40);

	void MonteCarlo(int SampleNum);
	void GibbsSampling(int sample_index);
	virtual void GibbsSamplingInit(int * state_seq) {};
	void GetConditional(int * state_vector,int i,double * prob_table);
	void Analytic2NodeBelief();  // analytic 2-node belief, not used
	void InitSampler();

//  Inference Graph construction
	virtual void ConstructPhi();
	virtual void ConstructKsai();
	virtual void ComputePhi(int nodei,double * phi) {};
	virtual void ComputeKsai(int nodei,int nodej,int edgei, double * Ksai) {};
	void CalculateCompatFunc();   // Calculate compatibility functions

	void InitGraph(int n_num,int e_num,int state_num,int node_size,int edge_size);

public:
	IGraph();
	IGraph(int n_num,int e_num,int state_num,int node_size, int edge_size);
	virtual ~IGraph();
};

#define INODE_PTR(i)	((IGraphNode *)NODE_PTR(i))
#define IEDGE_PTR(i)	((IGraphEdge *)EDGE_PTR(i))

#endif
