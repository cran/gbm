//  GBM by Greg Ridgeway  Copyright (C) 2003

#include "poisson.h"

CPoisson::CPoisson()
{
}

CPoisson::~CPoisson()
{
}


HRESULT CPoisson::ComputeWorkingResponse
(
    double *adY,
    double *adMisc,
    double *adOffset,
    double *adF, 
    double *adZ, 
    double *adWeight,
    bool *afInBag,
    unsigned long nTrain
)
{
    unsigned long i = 0;
    double dF = 0.0;

    assert(adY != NULL);
    assert(adF != NULL);
    assert(adZ != NULL);
    assert(adWeight != NULL);

    // compute working response
    for(i=0; i < nTrain; i++)
    {
        dF = adF[i] + ((adOffset==NULL) ? 0.0 : adOffset[i]);
        adZ[i] = adY[i] - exp(dF);
    }

    return S_OK;
}



HRESULT CPoisson::InitF
(
    double *adY,
    double *adMisc,
    double *adOffset,
    double *adWeight,
    double &dInitF, 
    unsigned long cLength
)
{
    HRESULT hr = S_OK;

    double dSum = 0.0;
    double dDenom = 0.0;
    unsigned long i = 0;

    assert(adY != NULL);
    assert(adWeight != NULL);

    if(adOffset == NULL)
    {
        for(i=0; i<cLength; i++)
        {
            dSum += adWeight[i]*adY[i];
            dDenom += adWeight[i];
        }
    }
    else
    {
        for(i=0; i<cLength; i++)
        {
            dSum += adWeight[i]*adY[i];
            dDenom += adWeight[i]*exp(adOffset[i]);
        }
    }

    dInitF = log(dSum/dDenom);

    return hr;
}


double CPoisson::LogLikelihood
(
    double *adY,
    double *adMisc,
    double *adOffset,
    double *adWeight,
    double *adF,
    unsigned long cLength
)
{
    unsigned long i=0;
    double dL = 0.0;

    if(adOffset == NULL)
    {
        for(i=0; i<cLength; i++)
        {
            dL += adWeight[i]*(adY[i]*adF[i] - exp(adF[i]));
        }
    }
    else
    {
        for(i=0; i<cLength; i++)
        {
            dL += adWeight[i]*(adY[i]*(adOffset[i]+adF[i]) - 
                               exp(adOffset[i]+adF[i]));
        }
    }

    return dL;
}


HRESULT CPoisson::FitBestConstant
(
    double *adY,
    double *adMisc,
    double *adOffset,
    double *adW,
    double *adF,
    double *adZ,
    unsigned long *aiNodeAssign,
    unsigned long nTrain,
    VEC_P_NODETERMINAL vecpTermNodes,
    unsigned long cTermNodes,
    unsigned long cMinObsInNode,
    bool *afInBag,
    double *adFadj
)
{
    HRESULT hr = S_OK;

    unsigned long iObs = 0;
    unsigned long iNode = 0;
    vecdNum.resize(cTermNodes);
    vecdNum.assign(vecdNum.size(),0.0);
    vecdDen.resize(cTermNodes);
    vecdDen.assign(vecdDen.size(),0.0);

    if(adOffset == NULL)
    {
        for(iObs=0; iObs<nTrain; iObs++)
        {
            if(afInBag[iObs])
            {
                vecdNum[aiNodeAssign[iObs]] += adW[iObs]*adY[iObs];
                vecdDen[aiNodeAssign[iObs]] += adW[iObs]*exp(adF[iObs]);
            }
        }
    }
    else
    {
        for(iObs=0; iObs<nTrain; iObs++)
        {
            if(afInBag[iObs])
            {
                vecdNum[aiNodeAssign[iObs]] += adW[iObs]*adY[iObs];
                vecdDen[aiNodeAssign[iObs]] += 
                    adW[iObs]*exp(adOffset[iObs]+adF[iObs]);
            }
        }        
    }
    for(iNode=0; iNode<cTermNodes; iNode++)
    {
        if(vecpTermNodes[iNode]!=NULL)
        {
            if(vecdDen[iNode] == 0.0)
            {
                vecpTermNodes[iNode]->dPrediction = 0.0;
            }
            else
            {
                vecpTermNodes[iNode]->dPrediction = 
                    log(vecdNum[iNode]/vecdDen[iNode]);
            }
        }
    }

    return hr;
}


double CPoisson::BagImprovement
(
    double *adY,
    double *adMisc,
    double *adOffset,
    double *adWeight,
    double *adF,
    double *adFadj,
    bool *afInBag,
    double dStepSize,
    unsigned long nTrain
)
{
    double dReturnValue = 0.0;
    double dF = 0.0;
    unsigned long i = 0;

    for(i=0; i<nTrain; i++)
    {
        if(!afInBag[i])
        {
            dF = adF[i] + ((adOffset==NULL) ? 0.0 : adOffset[i]);

            dReturnValue = adWeight[i]*
                           (adY[i]*dStepSize*adFadj[i] - 
                            exp(dF+dStepSize*adFadj[i]) + 
                            exp(dF));
        }
    }

    return dReturnValue;
}


