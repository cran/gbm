//  GBM by Greg Ridgeway  Copyright (C) 2003

#include "node.h"

CNode::CNode()
{
    dPrediction = 0.0;
    dTrainW = 0.0;
    isTerminal = false;
}


CNode::~CNode()
{
    // the nodes get deleted by deleting the node factory
}


HRESULT CNode::Adjust
(
    unsigned long cMinObsInNode
)
{
    HRESULT hr = E_NOTIMPL;
    ErrorTrace(hr);
    return hr;
}


HRESULT CNode::Predict
(
    CDataset *pData, 
    unsigned long iRow, 
    double &dFadj
)
{
    HRESULT hr = E_NOTIMPL;
    ErrorTrace(hr);
    return hr;
}


double CNode::TotalError()
{
    HRESULT hr = E_NOTIMPL;
    ErrorTrace(hr);
    return hr;
}


HRESULT CNode::PrintSubtree
(
    unsigned long cIndent
)
{
    HRESULT hr = E_NOTIMPL;
    ErrorTrace(hr);
    return hr;
}


HRESULT CNode::GetVarRelativeInfluence
(
    double *adRelInf
)
{
    HRESULT hr = E_NOTIMPL;
    ErrorTrace(hr);
    return hr;
}


HRESULT CNode::TransferTreeToRList
(
    int &iNodeID,
    CDataset *pData,
    int *aiSplitVar,
    double *adSplitPoint,
    int *aiLeftNode,
    int *aiRightNode,
    int *aiMissingNode,
    double *adErrorReduction,
    double *adWeight,
    VEC_VEC_CATEGORIES &vecSplitCodes,
    int cCatSplitsOld,
    double dShrinkage
)
{
    return E_NOTIMPL;
}
