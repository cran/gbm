//  GBM by Greg Ridgeway  Copyright (C) 2003

#include "gaussian.h"

CGaussian::CGaussian()
{
}

CGaussian::~CGaussian()
{
}


HRESULT CGaussian::ComputeWorkingResponse
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
    HRESULT hr = S_OK;
    unsigned long i = 0;

    if((adY == NULL) || (adF == NULL) || (adZ == NULL) || (adWeight == NULL))
    {
        hr = E_INVALIDARG;
        goto Error;
    }

    if(adOffset == NULL)
    {
        for(i=0; i<nTrain; i++)
        {
            adZ[i] = adY[i] - adF[i];
        }
    }
    else
    {
        for(i=0; i<nTrain; i++)
        {
            adZ[i] = adY[i] - adOffset[i] - adF[i];
        }
    }

Cleanup:
    return hr;
Error:
    goto Cleanup;
}



HRESULT CGaussian::InitF
(
    double *adY,
    double *adMisc,
    double *adOffset, 
    double *adWeight,
    double &dInitF, 
    unsigned long cLength
)
{
    double dSum=0.0;
    double dTotalWeight = 0.0;
    unsigned long i=0;

    assert(adY != NULL);
    assert(adWeight != NULL);

    // compute the mean
    if(adOffset==NULL)
    {
        for(i=0; i<cLength; i++)
        {
            dSum += adWeight[i]*adY[i];
            dTotalWeight += adWeight[i];
        }
    }
    else
    {
        for(i=0; i<cLength; i++)
        {
            dSum += adWeight[i]*(adY[i] - adOffset[i]);
            dTotalWeight += adWeight[i];
        }
    }
    dInitF = dSum/dTotalWeight;

    return S_OK;
}


double CGaussian::LogLikelihood
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
    assert(adWeight != NULL);
    assert(adF != NULL);

    if(adOffset == NULL)
    {
        for(i=0; i<cLength; i++)
        {
            dL += adWeight[i]*(adY[i]-adF[i])*(adY[i]-adF[i]);
        }
    }
    else
    {
        for(i=0; i<cLength; i++)
        {
            dL += adWeight[i]*(adY[i]-adOffset[i]-adF[i])*
                              (adY[i]-adOffset[i]-adF[i]);
        }
    }

    return dL;
}


HRESULT CGaussian::FitBestConstant
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
    // the tree aready stores the mean prediction
    // no refitting necessary

    return S_OK;
}

double CGaussian::BagImprovement
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
            
            dReturnValue += adWeight[i]*dStepSize*adFadj[i]*
                            (2.0*(adY[i]-dF) - dStepSize*adFadj[i]);
        }
    }

    return dReturnValue;
}

