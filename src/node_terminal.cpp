//------------------------------------------------------------------------------
//  GBM by Greg Ridgeway  Copyright (C) 2003
//
//  File:       node_terminal.cpp
//
//------------------------------------------------------------------------------

#include "node_terminal.h"
#include "node_factory.h"

CNodeTerminal::CNodeTerminal()
{
    isTerminal = true;
}


CNodeTerminal::~CNodeTerminal()
{
    #ifdef NOISY_DEBUG
    Rprintf("terminal destructor\n");
    #endif
}

HRESULT CNodeTerminal::Adjust
(
    unsigned long cMinObsInNode
)
{
    return S_OK;
}

HRESULT CNodeTerminal::ApplyShrinkage
(
    double dLambda
)
{
    HRESULT hr = S_OK;
    dPrediction *= dLambda;

    return hr;
}



HRESULT CNodeTerminal::Predict
(
    CDataset *pData, 
    unsigned long iRow, 
    double &dFadj
)
{
    dFadj = dPrediction;

    return S_OK;
}




HRESULT CNodeTerminal::Predict
(
    double *adX,
    unsigned long cRow,
    unsigned long cCol,
    unsigned long iRow,
    double &dFadj
)
{
    dFadj = dPrediction;

    return S_OK;
}



HRESULT CNodeTerminal::PrintSubtree
(
    unsigned long cIndent
)
{
    unsigned long i = 0;
    
    for(i=0; i< cIndent; i++) Rprintf("  ");
    Rprintf("N=%f, Prediction=%f *\n",
           dTrainW,
           dPrediction);

    return S_OK;
}


HRESULT CNodeTerminal::GetVarRelativeInfluence
(
    double *adRelInf
)
{
    return S_OK;
}


HRESULT CNodeTerminal::RecycleSelf
(
    CNodeFactory *pNodeFactory
)
{
    pNodeFactory->RecycleNode(this);
    return S_OK;
};



HRESULT CNodeTerminal::TransferTreeToRList
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
    HRESULT hr = S_OK;

    aiSplitVar[iNodeID] = -1;
    adSplitPoint[iNodeID] = dShrinkage*dPrediction;
    aiLeftNode[iNodeID] = -1;
    aiRightNode[iNodeID] = -1;
    aiMissingNode[iNodeID] = -1;
    adErrorReduction[iNodeID] = 0.0;
    adWeight[iNodeID] = dTrainW;

    iNodeID++;

    return hr;
}
