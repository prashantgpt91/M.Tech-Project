/*****************************************************************************
// GraphMatcher.h : ARG or Bag-of-Parts similarity by MRF inference
// 
// Written by Dong-Qing Zhang, DVMM lab, Columbia University, 2005
//
// Contact : 
//			Dong-Qing Zhang :	 dongqing@gmail.com, 
//			Prof.Shih-Fu Chang : sfchang@ee.columbia.edu
*****************************************************************************/

#ifndef GRAPH_MATCHER
#define GRAPH_MATCHER

#include "Graph.h"
#include "Garg.h"
#include "Igraph.h"

// all graph matching parameters
typedef struct tag_MatchPara
{
	double r;
	double Spos,Mneg,Sneg;

	int		nFeatDim;
	int		eFeatDim;

	double * ForegroundVar1n;
	double * BackgroundMean1n;
	double * BackgroundVar1n;

	double * ForegroundVar2n;
	double * BackgroundMean2n;
	double * BackgroundVar2n;

} MatchPara;

// MRF graph node, used to hold the one node features
class MRFNode : public IGraphNode
{
public:
	double * feat_vec;
	int		 feat_dim;
	float	 label;   // the label of the node
	int      feat_vec_int[4];   // integer feature vector 
};

// MRF graph edge, used to hold the two-node features
class MRFEdge : public IGraphEdge
{
public:
	double * feat_vec;   // feature vector in edge
	int		 feat_dim;
	int      feat_vec_int[10];   // integer feature vector 
};

// MRF model for relational graph matching
class MatchMRFGraph : public IGraph
{
	double * m_M;   // the message cache for BP

	int  m_LearningMode;	   // Learning mode, will save the correspondence data
	int  m_LearningSign;	   // 1 -- positive, 0 -- negative
	int  m_InitFlag;		   // 1 -- Initialization, 0 -- follow-up training
	char m_LearningDir[500];   // the directory to save the training data for learning

	MatchPara m_MatchPara;     // graph matching parameters

public:
	GARG  *m_Garg1,*m_Garg2;

	MRFNode * mnode_list;   // short-cut to the nodes list
	MRFEdge * medge_list;   // short-cut to the nodes list

public:
// form adj graph
	virtual short * Prefiltering(GARG *garg1, GARG *garg2,int * edge_num);
	double FormAdjGraph(GARG * img_rgns1,GARG * img_rgns2,int node_bound=40);  // form adjacency graph
	double ComputePrior_z();
	virtual void ComputePhi(int nodei,double * phi);
	virtual void ComputeKsai(int nodei,int nodej,int edgei, double * Ksai);

// learning-related
	void SaveBeliefInit(const char *NodeBlfFile1,const char *NodeBlfFile2,GARG *garg1, GARG *garg2);
	void SaveBelief(const char * NodeBlfFile1,const char * NodeBlfFile2);
	void LoadParameters(const char * TrainingDataDir);

// inference
	void Inference(int samplex=40);

// simplified BP for complete graph
	void ComputeMessage(int edge_index, int ni,int nj);
	virtual void BeliefProp();   // belief prorogation
	void   Compute2NodeBelief();

// Gibbs sampling, implemented in IGraph
	void DFSSearch(int * state_seq);
	virtual void GibbsSamplingInit(int * state_seq);

// inference others
	void MeanFiledBP();

// similarity-related
	double LogLikelihoodRatio();

	double ComputeAverageEnergy();
	double ComputeEntropy();
	double ComputeLogZp();
	double ComputeLogZ();
	float  StructureDiff();

// other functions
	void FreeMemory();
	double SimpleNormal(double * feature,double * mean,double * var,int dim);
	void SwitchTwoGraphs(GARG **garg1, GARG **garg2);

public:
	MatchMRFGraph();
	~MatchMRFGraph();
};

#endif
