//------------------------------------------------------------------------------
//
//  GBM by Greg Ridgeway  Copyright (C) 2003
//  File:       gbm.cpp
//
//------------------------------------------------------------------------------

#include<vector>
#include "gbm.h"

unsigned long gbm_setup
(
    double *adY,
    double *adOffset,
    double *adX,
    int *aiXOrder,
    double *adWeight,
    double *adMisc,
    int cRows,
    int cCols,
    int *acVarClasses,
    int *alMonotoneVar,
    char *pszFamily,
    int cTrees,
    int cDepth,
    int cMinObsInNode,
    double dShrinkage,
    double dBagFraction,
    int cTrain,
    CDataset *pData,
    PCDistribution &pDist
)
{
    unsigned long hr = 0;

    hr = pData->SetData(adX,aiXOrder,adY,adOffset,adWeight,adMisc,
                        cRows,cCols,acVarClasses,alMonotoneVar);
    if(FAILED(hr))
    {
        goto Error;
    }

    // set the distribution
    if(strncmp(pszFamily,"bernoulli",2) == 0)
    {
        pDist = new CBernoulli(); new_item(pDist);
        if(pDist==NULL)
        {
            hr = E_OUTOFMEMORY;
            goto Error;
        }
    }
    else if(strncmp(pszFamily,"gaussian",2) == 0)
    {
        pDist = new CGaussian(); new_item(pDist);
        if(pDist==NULL)
        {
            hr = E_OUTOFMEMORY;
            goto Error;
        }
    }
    else if(strncmp(pszFamily,"poisson",2) == 0)
    {
        pDist = new CPoisson(); new_item(pDist);
        if(pDist==NULL)
        {
            hr = E_OUTOFMEMORY;
            goto Error;
        }
    }
    else if(strncmp(pszFamily,"adaboost",2) == 0)
    {
        pDist = new CAdaBoost(); new_item(pDist);
        if(pDist==NULL)
        {
            hr = E_OUTOFMEMORY;
            goto Error;
        }
    }
    else if(strncmp(pszFamily,"coxph",2) == 0)
    {
        pDist = new CCoxPH(); new_item(pDist);
        if(pDist==NULL)
        {
            hr = E_OUTOFMEMORY;
            goto Error;
        }
    }
    else if(strncmp(pszFamily,"laplace",2) == 0)
    {
        pDist = new CLaplace(); new_item(pDist);
        if(pDist==NULL)
        {
            hr = E_OUTOFMEMORY;
            goto Error;
        }
    }
    if(pDist==NULL)
    {
        hr = E_INVALIDARG;
        goto Error;
    }

Cleanup:
    return hr;
Error:
    goto Cleanup;
}



HRESULT gbm_transfer_to_R
(
    CGBM *pGBM,
    VEC_VEC_CATEGORIES &vecSplitCodes,
    int *aiSplitVar,
    double *adSplitPoint,
    int *aiLeftNode,
    int *aiRightNode,
    int *aiMissingNode,
    double *adErrorReduction,
    double *adWeight,
    int cCatSplitsOld
)
{
    HRESULT hr = S_OK;

    hr = pGBM->TransferTreeToRList(aiSplitVar,
                                   adSplitPoint,
                                   aiLeftNode,
                                   aiRightNode,
                                   aiMissingNode,
                                   adErrorReduction,
                                   adWeight,
                                   vecSplitCodes,
                                   cCatSplitsOld);
    if(FAILED(hr)) goto Error;

Cleanup:
    return hr;
Error:
    goto Cleanup;
}


HRESULT gbm_transfer_catsplits_to_R
(
    int iCatSplit,
    VEC_VEC_CATEGORIES &vecSplitCodes,
    int *aiSplitCodes
)
{
    unsigned long i=0;

    for(i=0; i<vecSplitCodes[iCatSplit].size(); i++)
    {
        aiSplitCodes[i] = vecSplitCodes[iCatSplit][i];
    }

    return S_OK;
}


int size_of_vector
(
    VEC_VEC_CATEGORIES &vec,
    int i
)
{
    return vec[i].size();
}



