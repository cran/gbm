//------------------------------------------------------------------------------
//  GBM by Greg Ridgeway  Copyright (C) 2003
//
//  File:       node_search.cpp
//
//------------------------------------------------------------------------------

#include "node_search.h"

CNodeSearch::CNodeSearch()
    :k_cMaxClasses(256)
{
    iBestSplitVar = 0;
    dBestSplitValue = 0.0;
    fIsSplit = false;

    dBestMissingTotalW = 0.0;
    dCurrentMissingTotalW = 0.0;
    dBestMissingSumZ = 0.0;
    dCurrentMissingSumZ = 0.0;

    adGroupSumZ = NULL;
    adGroupW = NULL;
    acGroupN = NULL;
    adGroupMean = NULL;
    aiCurrentCategory = NULL;
    aiBestCategory = NULL;

    iRank = UINT_MAX;
}


CNodeSearch::~CNodeSearch()
{
    if(adGroupSumZ != NULL)
    {
        delete [] adGroupSumZ; delete_item(adGroupSumZ);
        adGroupSumZ = NULL;
    }
    if(adGroupW != NULL)
    {
        delete [] adGroupW; delete_item(adGroupW);
        adGroupW = NULL;
    }
    if(acGroupN != NULL)
    {
        delete [] acGroupN; delete_item(acGroupN);
        acGroupN = NULL;
    }
    if(adGroupMean != NULL)
    {
        delete [] adGroupMean; delete_item(adGroupMean);
        adGroupMean = NULL;
    }
    if(aiCurrentCategory != NULL)
    {
        delete [] aiCurrentCategory; delete_item(aiCurrentCategory);
        aiCurrentCategory = NULL;
    }
    if(aiBestCategory != NULL)
    {
        delete [] aiBestCategory; delete_item(aiBestCategory);
        aiBestCategory = NULL;
    }
}


HRESULT CNodeSearch::Initialize
(
    unsigned long cMinObsInNode
)
{
    HRESULT hr = S_OK;

    adGroupSumZ = new double[k_cMaxClasses]; new_item(adGroupSumZ);
    if(adGroupSumZ == NULL)
    {
        hr = E_OUTOFMEMORY;
        ErrorTrace(hr);
        goto Error;
    }
    adGroupW = new double[k_cMaxClasses]; new_item(adGroupW);
    if(adGroupW == NULL)
    {
        hr = E_OUTOFMEMORY;
        ErrorTrace(hr);
        goto Error;
    }
    acGroupN = new ULONG[k_cMaxClasses]; new_item(acGroupN);
    if(acGroupN == NULL)
    {
        hr = E_OUTOFMEMORY;
        ErrorTrace(hr);
        goto Error;
    }
    adGroupMean = new double[k_cMaxClasses]; new_item(adGroupMean);
    if(adGroupMean == NULL)
    {
        hr = E_OUTOFMEMORY;
        ErrorTrace(hr);
        goto Error;
    }
    aiCurrentCategory = new int[k_cMaxClasses]; new_item(aiCurrentCategory);
    if(aiCurrentCategory == NULL)
    {
        hr = E_OUTOFMEMORY;
        ErrorTrace(hr);
        goto Error;
    }
    aiBestCategory = new ULONG[k_cMaxClasses]; new_item(aiBestCategory);
    if(aiBestCategory == NULL)
    {
        hr = E_OUTOFMEMORY;
        ErrorTrace(hr);
        goto Error;
    }

    this->cMinObsInNode = cMinObsInNode;

Cleanup:
    return hr;
Error:
    goto Cleanup;
}


HRESULT CNodeSearch::Set
(
    double dSumZ,
    double dTotalW,
    unsigned long cTotalN,
    CNodeTerminal *pThisNode,
    CNode **ppParentPointerToThisNode,
    CNodeFactory *pNodeFactory
)
{
    HRESULT hr = S_OK;

    dInitSumZ = dSumZ;
    dInitTotalW = dTotalW;
    cInitN = cTotalN;

    dBestLeftSumZ       = 0.0;
    dBestLeftTotalW     = 0.0;
    cBestLeftN          = 0;
    dCurrentLeftSumZ    = 0.0;
    dCurrentLeftTotalW  = 0.0;
    cCurrentLeftN       = 0;

    dBestRightSumZ      = dSumZ;
    dBestRightTotalW    = dTotalW;
    cBestRightN         = cTotalN;
    dCurrentRightSumZ   = 0.0;
    dCurrentRightTotalW = dTotalW;
    cCurrentRightN      = cTotalN;

    dBestMissingSumZ      = 0.0;
    dBestMissingTotalW    = 0.0;
    cBestMissingN         = 0;
    dCurrentMissingSumZ   = 0.0;
    dCurrentMissingTotalW = 0.0;
    cCurrentMissingN      = 0;

    dBestImprovement    = 0.0;
    iBestSplitVar       = UINT_MAX;

    dCurrentImprovement = 0.0;
    iCurrentSplitVar    = UINT_MAX;
    dCurrentSplitValue  = -HUGE_VAL;

    fIsSplit = false;

    this->pThisNode = pThisNode;
    this->ppParentPointerToThisNode = ppParentPointerToThisNode;
    this->pNodeFactory = pNodeFactory;

    return hr;
}


HRESULT CNodeSearch::ResetForNewVar
(
    unsigned long iWhichVar,
    long cCurrentVarClasses
)
{
    HRESULT hr = S_OK;
    long i=0;

    if(fIsSplit) goto Cleanup;

    for(i=0; i<cCurrentVarClasses; i++)
    {
        adGroupSumZ[i] = 0.0;
        adGroupW[i] = 0.0;
        acGroupN[i] = 0;
    }

    iCurrentSplitVar = iWhichVar;
    this->cCurrentVarClasses = cCurrentVarClasses;

    dCurrentLeftSumZ      = 0.0;
    dCurrentLeftTotalW    = 0.0;
    cCurrentLeftN         = 0;
    dCurrentRightSumZ     = dInitSumZ;
    dCurrentRightTotalW   = dInitTotalW;
    cCurrentRightN        = cInitN;
    dCurrentMissingSumZ   = 0.0;
    dCurrentMissingTotalW = 0.0;
    cCurrentMissingN      = 0;

    dCurrentImprovement = 0.0;

    dLastXValue = -HUGE_VAL;

Cleanup:
    return hr;
}



HRESULT CNodeSearch::WrapUpCurrentVariable()
{
    HRESULT hr = S_OK;
    if(iCurrentSplitVar == iBestSplitVar)
    {
        if(cCurrentMissingN > 0)
        {
            dBestMissingSumZ   = dCurrentMissingSumZ;
            dBestMissingTotalW = dCurrentMissingTotalW;
            cBestMissingN      = cCurrentMissingN;
        }
        else // DEBUG: consider a weighted average with parent node?
        {
            dBestMissingSumZ   = dInitSumZ;
            dBestMissingTotalW = dInitTotalW;
            cBestMissingN      = 0;
        }
    }

    return hr;
}



HRESULT CNodeSearch::EvaluateCategoricalSplit()
{
    HRESULT hr = S_OK;
    long i=0;
    long j=0;
    unsigned long cFiniteMeans = 0;

    if(fIsSplit) goto Cleanup;

    if(cCurrentVarClasses == 0)
    {
        hr = E_INVALIDARG;
        ErrorTrace(hr);
        goto Error;
    }

    cFiniteMeans = 0;
    for(i=0; i<cCurrentVarClasses; i++)
    {
        aiCurrentCategory[i] = i;
        if(adGroupW[i] != 0.0)
        {
            adGroupMean[i] = adGroupSumZ[i]/adGroupW[i];
            cFiniteMeans++;
        }
        else
        {
            adGroupMean[i] = HUGE_VAL;
        }
    }

    rsort_with_index(adGroupMean,aiCurrentCategory,cCurrentVarClasses);

    // if only one group has a finite mean it will not consider
    // might be all are missing so no categories enter here
    for(i=0; (cFiniteMeans>1) && ((ULONG)i<cFiniteMeans-1); i++)
    {
        dCurrentSplitValue = (double)i;

        dCurrentLeftSumZ    += adGroupSumZ[aiCurrentCategory[i]];
        dCurrentLeftTotalW  += adGroupW[aiCurrentCategory[i]];
        cCurrentLeftN       += acGroupN[aiCurrentCategory[i]];
        dCurrentRightSumZ   -= adGroupSumZ[aiCurrentCategory[i]];
        dCurrentRightTotalW -= adGroupW[aiCurrentCategory[i]];
        cCurrentRightN      -= acGroupN[aiCurrentCategory[i]];

        dCurrentImprovement = 
            CNode::Improvement(dCurrentLeftTotalW,dCurrentRightTotalW,
                               dCurrentMissingTotalW,
                               dCurrentLeftSumZ,dCurrentRightSumZ,
                               dCurrentMissingSumZ);
        if((cCurrentLeftN >= cMinObsInNode) && 
           (cCurrentRightN >= cMinObsInNode) &&
           (dCurrentImprovement > dBestImprovement))
        {
            dBestSplitValue = dCurrentSplitValue;
            if(iBestSplitVar != iCurrentSplitVar)
            {
                iBestSplitVar = iCurrentSplitVar;
                cBestVarClasses = cCurrentVarClasses;
                for(j=0; j<cCurrentVarClasses; j++)
                {
                    aiBestCategory[j] = aiCurrentCategory[j];
                }
            }

            dBestLeftSumZ      = dCurrentLeftSumZ;
            dBestLeftTotalW    = dCurrentLeftTotalW;
            cBestLeftN         = cCurrentLeftN;
            dBestRightSumZ     = dCurrentRightSumZ;
            dBestRightTotalW   = dCurrentRightTotalW;
            cBestRightN        = cCurrentRightN;
            dBestImprovement   = dCurrentImprovement;
        }
    }

Cleanup:
    return hr;
Error:
    goto Cleanup;
}




HRESULT CNodeSearch::SetupNewNodes
(
    PCNodeNonterminal &pNewSplitNode,
    PCNodeTerminal &pNewLeftNode,
    PCNodeTerminal &pNewRightNode,
    PCNodeTerminal &pNewMissingNode
)
{
    HRESULT hr = S_OK;
    CNodeContinuous *pNewNodeContinuous = NULL;
    CNodeCategorical *pNewNodeCategorical = NULL;
    unsigned long i=0;

    pNewLeftNode    = pNodeFactory->GetNewNodeTerminal();
    pNewRightNode   = pNodeFactory->GetNewNodeTerminal();
    pNewMissingNode = pNodeFactory->GetNewNodeTerminal();

    // set up a continuous split
    if(cBestVarClasses==0)
    {
        pNewNodeContinuous = pNodeFactory->GetNewNodeContinuous();

        pNewNodeContinuous->dSplitValue = dBestSplitValue;
        pNewNodeContinuous->iSplitVar = iBestSplitVar;

        pNewSplitNode = pNewNodeContinuous;
    }
    else
    {
        // get a new categorical node and its branches
        pNewNodeCategorical = pNodeFactory->GetNewNodeCategorical();

        // set up the categorical split
        pNewNodeCategorical->iSplitVar = iBestSplitVar;
        pNewNodeCategorical->cLeftCategory = (ULONG)dBestSplitValue + 1;
        pNewNodeCategorical->aiLeftCategory = 
            new ULONG[pNewNodeCategorical->cLeftCategory]; new_item(pNewNodeCategorical->aiLeftCategory);
        for(i=0; i<pNewNodeCategorical->cLeftCategory; i++)
        {
            pNewNodeCategorical->aiLeftCategory[i] = aiBestCategory[i];
        }

        pNewSplitNode = pNewNodeCategorical;
    }

    *ppParentPointerToThisNode = pNewSplitNode;

    pNewSplitNode->dPrediction  = pThisNode->dPrediction;
    pNewSplitNode->dImprovement = dBestImprovement;
    pNewSplitNode->dTrainW      = pThisNode->dTrainW;
    pNewSplitNode->pLeftNode    = pNewLeftNode;
    pNewSplitNode->pRightNode   = pNewRightNode;
    pNewSplitNode->pMissingNode = pNewMissingNode;

    pNewLeftNode->dPrediction    = dBestLeftSumZ/dBestLeftTotalW;
    pNewLeftNode->dTrainW        = dBestLeftTotalW;
    pNewLeftNode->cN             = cBestLeftN;
    pNewRightNode->dPrediction   = dBestRightSumZ/dBestRightTotalW;
    pNewRightNode->dTrainW       = dBestRightTotalW;
    pNewRightNode->cN            = cBestRightN;
    pNewMissingNode->dPrediction = dBestMissingSumZ/dBestMissingTotalW;
    pNewMissingNode->dTrainW     = dBestMissingTotalW;
    pNewMissingNode->cN          = cBestMissingN;

    pThisNode->RecycleSelf(pNodeFactory);

    return hr;
}
