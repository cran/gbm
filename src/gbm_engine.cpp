//  GBM by Greg Ridgeway  Copyright (C) 2003

#include "gbm_engine.h"

CGBM::CGBM()
{
    adFadj = NULL;
    adZ = NULL;
    afInBag = NULL;
    aiNodeAssign = NULL;
    aNodeSearch = NULL;
    
    cDepth = 0;
    cMinObsInNode = 0;
    dBagFraction = 0.0;
    dLambda = 0.0;
    fInitialized = false;
    cTotalInBag = 0;
    cTrain = 0;
    cValid = 0;

    pData = NULL;
    pDist = NULL;
    pNodeFactory = NULL;
    ptreeTemp = NULL;
}


CGBM::~CGBM()
{
    if(adFadj != NULL)
    {
        delete [] adFadj; delete_item(adFadj);
        adFadj = NULL;
    }
    if(adZ != NULL)
    {
        delete [] adZ; delete_item(adZ);
        adZ = NULL;
    }
    if(afInBag != NULL)
    {
        delete [] afInBag; delete_item(afInBag);
        afInBag = NULL;
    }
    if(aiNodeAssign != NULL)
    {
        delete [] aiNodeAssign; delete_item(aiNodeAssign);
        aiNodeAssign = NULL;
    }
    if(aNodeSearch != NULL)
    {
        delete [] aNodeSearch; delete_item(aNodeSearch);
        aNodeSearch = NULL;
    }
    if(ptreeTemp != NULL)
    {
        delete ptreeTemp; delete_item(ptreeTemp);
        ptreeTemp = NULL;
    }
    // must delete the node factory last!!! at least after deleting trees
    if(pNodeFactory != NULL)
    {
        delete pNodeFactory; delete_item(pNodeFactory);
        pNodeFactory = NULL;
    }
}


HRESULT CGBM::Initialize
(
    CDataset *pData,
    CDistribution *pDist,
    double dLambda,
    unsigned long cTrain,
    double dBagFraction,
    unsigned long cDepth,
    unsigned long cMinObsInNode
)
{
    HRESULT hr = S_OK;
    unsigned long i=0;

    if(pData == NULL)
    {
        hr = E_INVALIDARG;
        ErrorTrace(hr);
        goto Error;
    }
    if(pDist == NULL)
    {
        hr = E_INVALIDARG;
        ErrorTrace(hr);
        goto Error;
    }

    this->pData = pData;
    this->pDist = pDist;
    this->dLambda = dLambda;
    this->cTrain = cTrain;
    this->dBagFraction = dBagFraction;
    this->cDepth = cDepth;
    this->cMinObsInNode = cMinObsInNode;

    // allocate the tree structure
    ptreeTemp = new CCARTTree; new_item(ptreeTemp);
    if(ptreeTemp == NULL)
    {
        hr = E_OUTOFMEMORY;
        ErrorTrace(hr);
        goto Error;
    }

    cValid = pData->cRows - cTrain;
    cTotalInBag = (unsigned long)(dBagFraction*cTrain);

    adZ = new double[cTrain]; new_item(adZ);
    if(adZ == NULL)
    {
        hr = E_OUTOFMEMORY;
        ErrorTrace(hr);
        goto Error;
    }
    adFadj = new double[pData->cRows]; new_item(adFadj);
    if(adFadj == NULL)
    {
        hr = E_OUTOFMEMORY;
        ErrorTrace(hr);
        goto Error;
    }

    pNodeFactory = new CNodeFactory(); new_item(pNodeFactory);
    if(pNodeFactory == NULL)
    {
        hr = E_OUTOFMEMORY;
        ErrorTrace(hr);
        goto Error;
    }
    hr = pNodeFactory->Initialize(cDepth);
    if(FAILED(hr))
    {
        goto Error;
    }
    ptreeTemp->Initialize(pNodeFactory);

    // array for flagging those observations in the bag
    afInBag = new bool[cTrain]; new_item(afInBag);
    if(afInBag==NULL)
    {
        hr = E_OUTOFMEMORY;
        ErrorTrace(hr);
        goto Error;
    }
    // aiNodeAssign tracks to which node each training obs belongs
    aiNodeAssign = new ULONG[cTrain]; new_item(aiNodeAssign);
    if(aiNodeAssign==NULL)
    {
        hr = E_OUTOFMEMORY;
        ErrorTrace(hr);
        goto Error;
    }
    // NodeSearch objects help decide which nodes to split
    aNodeSearch = new CNodeSearch[2*cDepth+1]; new_item(aNodeSearch);
    if(aNodeSearch==NULL)
    {
        hr = E_OUTOFMEMORY;
        ErrorTrace(hr);
        goto Error;
    }
    for(i=0; i<2*cDepth+1; i++)
    {
        aNodeSearch[i].Initialize(cMinObsInNode);
    }
    vecpTermNodes.resize(2*cDepth+1,NULL);

    fInitialized = true;

Cleanup:
    return hr;
Error:
    goto Cleanup;
}




HRESULT CGBM::Predict
(
    unsigned long iVar,
    unsigned long cTrees,
    double *adF,
    double *adX,
    unsigned long cLength
)
{
    HRESULT hr = S_OK;


    return hr;
}


HRESULT CGBM::Predict
(
    double *adX,
    unsigned long cRow,
    unsigned long cCol,
    unsigned long cTrees,
    double *adF
)
{
    HRESULT hr = S_OK;


    return hr;
}



HRESULT CGBM::GetVarRelativeInfluence
(
    double *adRelInf,
    unsigned long cTrees
)
{
    HRESULT hr = S_OK;
    int iVar=0;

    for(iVar=0; iVar<pData->cCols; iVar++)
    {
        adRelInf[iVar] = 0.0;
    }

    return hr;
}


HRESULT CGBM::PrintTree()
{
    HRESULT hr = S_OK;

    hr = ptreeTemp->Print();
    if(FAILED(hr)) goto Error;

Cleanup:
    return hr;
Error:
    goto Cleanup;
}




HRESULT CGBM::iterate
(
    double *adF,
    double &dTrainError,
    double &dValidError,
    double &dOOBagImprove,
    int &cNodes
)
{
    HRESULT hr = S_OK;
    unsigned long i = 0;
    unsigned long cBagged = 0;

    if(!fInitialized)
    {
        hr = E_FAIL;
        ErrorTrace(hr);
        goto Error;
    }

    dTrainError = 0.0;
    dValidError = 0.0;
    dOOBagImprove = 0.0;

    vecpTermNodes.assign(2*cDepth+1,NULL);

    // randomly assign observations to the Bag
    cBagged = 0;
    for(i=0; i<cTrain; i++)
    {
        if(unif_rand()*(cTrain-i) < cTotalInBag-cBagged)
        {
            afInBag[i] = true;
            cBagged++;
        }
        else
        {
            afInBag[i] = false;
        }
    }

    #ifdef NOISY_DEBUG
    Rprintf("Compute working response\n");
    #endif
    hr = pDist->ComputeWorkingResponse(pData->adY, 
                                       pData->adMisc,
                                       pData->adOffset,
                                       adF, 
                                       adZ,
                                       pData->adWeight,
                                       afInBag,
                                       cTrain);
    if(FAILED(hr))
    {
        goto Error;
    }

    #ifdef NOISY_DEBUG
    Rprintf("Reset tree\n");
    #endif
    hr = ptreeTemp->Reset();
    #ifdef NOISY_DEBUG
    Rprintf("grow tree\n");
    #endif
    hr = ptreeTemp->grow(adZ,pData,pData->adWeight,adFadj,
                         cTrain,cTotalInBag,dLambda,cDepth,
                         cMinObsInNode,
                         afInBag,
                         aiNodeAssign,aNodeSearch,vecpTermNodes);
    if(FAILED(hr))
    {
        goto Error;
    }
    #ifdef NOISY_DEBUG
    Rprintf("get node count\n");
    #endif
    hr = ptreeTemp->GetNodeCount(cNodes);
    if(FAILED(hr))
    {
        goto Error;
    }

    // Now I have adF, adZ, and vecpTermNodes (new node assignments)
    // Fit the best constant within each terminal node
    #ifdef NOISY_DEBUG
    Rprintf("fit best constant\n");
    #endif
    hr = pDist->FitBestConstant(pData->adY,
                                pData->adMisc,
                                pData->adOffset,
                                pData->adWeight,
                                adF,
                                adZ,
                                aiNodeAssign,
                                cTrain,
                                vecpTermNodes,
                                (2*cNodes+1)/3, // number of terminal nodes
                                cMinObsInNode,
                                afInBag,
                                adFadj);
    if(FAILED(hr))
    {
        goto Error;
    }

    // update training predictions
    // fill in missing nodes where N < cMinObsInNode
    hr = ptreeTemp->Adjust(aiNodeAssign,adFadj,cTrain,
                           vecpTermNodes,cMinObsInNode);
    if(FAILED(hr))
    {
        goto Error;
    }

    ptreeTemp->SetShrinkage(dLambda);

    dOOBagImprove = pDist->BagImprovement(pData->adY,
                                          pData->adMisc,
                                          pData->adOffset,
                                          pData->adWeight,
                                          adF,
                                          adFadj,
                                          afInBag,
                                          dLambda,
                                          cTrain);

    // update the training predictions
    for(i=0; i < cTrain; i++)
    {
        adF[i] += dLambda*adFadj[i];
    }
    dTrainError = pDist->LogLikelihood(pData->adY,
                                       pData->adMisc,
                                       pData->adOffset,
                                       pData->adWeight,
                                       adF,
                                       cTrain);
    // update the validation predictions
    hr = ptreeTemp->PredictValid(pData,cValid,adFadj);
    for(i=cTrain; i < cTrain+cValid; i++)
    {
        adF[i] += adFadj[i];
    }
    if(pData->fHasOffset)
    {
        dValidError = 
            pDist->LogLikelihood(&(pData->adY[cTrain]),
                                 &(pData->adMisc[cTrain]),
                                 &(pData->adOffset[cTrain]),
                                 &(pData->adWeight[cTrain]),
                                 &(adF[cTrain]),
                                 cValid);
    }
    else
    {
        dValidError = pDist->LogLikelihood(&(pData->adY[cTrain]),
                                           &(pData->adMisc[cTrain]),
                                           NULL,
                                           &(pData->adWeight[cTrain]),
                                           &(adF[cTrain]),
                                           cValid);
    }

Cleanup:
    return hr;

Error:
    goto Cleanup;
}


HRESULT CGBM::TransferTreeToRList
(
    int *aiSplitVar,
    double *adSplitPoint,
    int *aiLeftNode,
    int *aiRightNode,
    int *aiMissingNode,
    double *adErrorReduction,
    double *adWeight,
    VEC_VEC_CATEGORIES &vecSplitCodes,
    int cCatSplitsOld
)
{
    HRESULT hr = S_OK;

    hr = ptreeTemp->TransferTreeToRList(pData,
                                        aiSplitVar,
                                        adSplitPoint,
                                        aiLeftNode,
                                        aiRightNode,
                                        aiMissingNode,
                                        adErrorReduction,
                                        adWeight,
                                        vecSplitCodes,
                                        cCatSplitsOld,
                                        dLambda);

    return hr;
}


