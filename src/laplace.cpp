//  GBM by Greg Ridgeway  Copyright (C) 2003

#include "laplace.h"

CLaplace::CLaplace()
{

}

CLaplace::~CLaplace()
{

}


HRESULT CLaplace::ComputeWorkingResponse
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

    assert(adY != NULL);
    assert(adF != NULL);
    assert(adZ != NULL);
    assert(adWeight != NULL);

    if(adOffset == NULL)
    {
        for(i=0; i<nTrain; i++)
        {
            adZ[i] = (adY[i] - adF[i]) > 0.0 ? 1.0 : -1.0;
        }
    }
    else
    {
        for(i=0; i<nTrain; i++)
        {
            adZ[i] = (adY[i] - adOffset[i] - adF[i]) > 0.0 ? 1.0 : -1.0;
        }
    }

    return S_OK;
}



// DEBUG: needs weighted median
HRESULT CLaplace::InitF
(
    double *adY,
    double *adMisc,
    double *adOffset, 
    double *adWeight,
    double &dInitF, 
    unsigned long cLength
)
{
    double dOffset=0.0;
    unsigned long i=0;

    assert(adY != NULL);
    assert(adWeight != NULL);

    vecd.resize(cLength);
    for(i=0; i<cLength; i++)
    {
        dOffset = (adOffset==NULL) ? 0.0 : adOffset[i];
        vecd[i] = adY[i] - dOffset;
    }

    itMedian = vecd.begin() + (vecd.end() - vecd.begin())/2;
    nth_element(vecd.begin(), itMedian, vecd.end()) ;

    dInitF = *itMedian;

    return S_OK;
}


double CLaplace::LogLikelihood
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

    assert(adY != NULL);
    assert(adF != NULL);

    if(adOffset == NULL)
    {
        for(i=0; i<cLength; i++)
        {
            dL += adWeight[i]*fabs(adY[i]-adF[i]);
        }
    }
    else
    {
        for(i=0; i<cLength; i++)
        {
            dL += adWeight[i]*fabs(adY[i]-adOffset[i]-adF[i]);
        }
    }

    return dL;
}


// DEBUG: needs weighted median
HRESULT CLaplace::FitBestConstant
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

    unsigned long iNode = 0;
    unsigned long iObs = 0;
    unsigned long iVecd = 0;
    double dOffset;

    vecd.resize(nTrain); // should already be this size from InitF
    for(iNode=0; iNode<cTermNodes; iNode++)
    {
        if(vecpTermNodes[iNode]->cN >= cMinObsInNode)
        {
            iVecd = 0;
            for(iObs=0; iObs<nTrain; iObs++)
            {
                if(afInBag[iObs] && (aiNodeAssign[iObs] == iNode))
                {
                    dOffset = (adOffset==NULL) ? 0.0 : adOffset[iObs];
                    vecd[iVecd] = adY[iObs] - dOffset - adF[iObs];
                    iVecd++;
                }
            }
            itMedian = vecd.begin() + ((vecd.begin()+iVecd) - vecd.begin())/2;
            nth_element(vecd.begin(), itMedian, (vecd.begin()+iVecd)) ;
            vecpTermNodes[iNode]->dPrediction = *itMedian;
        }
    }

    return hr;
}



double CLaplace::BagImprovement
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
            
            dReturnValue += 
                adWeight[i]*(fabs(adY[i]-dF) - fabs(adY[i]-dF-dStepSize*adFadj[i]));
        }
    }

    return dReturnValue;
}


