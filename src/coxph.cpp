//  GBM by Greg Ridgeway  Copyright (C) 2003

#include "coxph.h"

CCoxPH::CCoxPH()
{
}

CCoxPH::~CCoxPH()
{
}


HRESULT CCoxPH::ComputeWorkingResponse
(
    double *adT,
    double *adDelta,
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
    double dTot = 0.0;
    double dRiskTot = 0.0;

    assert(adT != NULL);
    assert(adDelta != NULL);
    assert(adF != NULL);
    assert(adZ != NULL);
    assert(adWeight != NULL);

    vecdRiskTot.resize(nTrain);
    dRiskTot = 0.0;
    for(i=0; i<nTrain; i++)
    {
        if(afInBag[i])
        {
            dF = adF[i] + ((adOffset==NULL) ? 0.0 : adOffset[i]);
            dRiskTot += adWeight[i]*exp(dF);
            vecdRiskTot[i] = dRiskTot;
        }
    }

    dTot = 0.0;
    for(i=nTrain-1; (i>=0) && (i<ULONG_MAX); i--)
    {
        if(afInBag[i])
        {
            if(adDelta[i]==1.0)
            {
                dTot += adWeight[i]/vecdRiskTot[i];
            }
            dF = adF[i] + ((adOffset==NULL) ? 0.0 : adOffset[i]);
            adZ[i] = adDelta[i] - exp(dF)*dTot;
        }
    }
    return S_OK;
}



HRESULT CCoxPH::InitF
(
    double *adY,
    double *adMisc,
    double *adOffset, 
    double *adWeight,
    double &dInitF, 
    unsigned long cLength
)
{
    dInitF = 0.0;

    return S_OK;
}


double CCoxPH::LogLikelihood
(
    double *adT,
    double *adDelta,
    double *adOffset, 
    double *adWeight,
    double *adF,
    unsigned long cLength
)
{
    unsigned long i=0;
    double dL = 0.0;
    double dF = 0.0;
    double dTotalAtRisk = 0.0;

    assert(adT != NULL);
    assert(adDelta != NULL);
    assert(adWeight != NULL);
    assert(adF != NULL);

    dTotalAtRisk = 0.0;
    for(i=0; i<cLength; i++)
    {
        dF = adF[i] + ((adOffset==NULL) ? 0.0 : adOffset[i]);
        dTotalAtRisk += adWeight[i]*exp(dF);
        if(adDelta[i]==1.0)
        {
            dL += adWeight[i]*(dF - log(dTotalAtRisk));
        }
    }

    return dL;
}


HRESULT CCoxPH::FitBestConstant
(
    double *adT,
    double *adDelta,
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

    double dF = 0.0;
    double dRiskTot = 0.0;
    unsigned long i = 0;
    unsigned long k = 0;
    unsigned long m = 0;

    double dTemp = 0.0;
    bool fTemp = false;

    unsigned long K = 0;
    veciK2Node.resize(cTermNodes);
    veciNode2K.resize(cTermNodes);

    for(i=0; i<cTermNodes; i++)
    {
        veciNode2K[i] = 0;
        if(vecpTermNodes[i]->cN > cMinObsInNode)
        {
            veciK2Node[K] = i;
            veciNode2K[i] = K;
            K++;
        }            
    }

    vecdP.resize(K);

    matH.setactualsize(K-1);
    vecdG.resize(K-1);
    vecdG.assign(K-1,0.0);

    // zero the Hessian
    for(k=0; k<K-1; k++)
    {
        for(m=0; m<K-1; m++)
        {
            matH.setvalue(k,m,0.0);
        }
    }

    // get the gradient & Hessian, Ridgeway (1999) pp. 100-101
    // correction from Ridgeway (1999): fix terminal node K-1 to be 0
    //      for identifiability
    dRiskTot = 0.0;
    vecdP.assign(K,0.0);
    for(i=0; i<nTrain; i++)
    {
        if(afInBag[i] && (vecpTermNodes[aiNodeAssign[i]]->cN > cMinObsInNode))
        {
            dF = adF[i] + ((adOffset==NULL) ? 0.0 : adOffset[i]);
            vecdP[veciNode2K[aiNodeAssign[i]]] += adW[i]*exp(dF);
            dRiskTot += adW[i]*exp(dF);

            if(adDelta[i]==1.0)
            {
                // compute g and H
                for(k=0; k<K-1; k++)
                {
                    vecdG[k] += 
                        adW[i]*((aiNodeAssign[i]==veciK2Node[k]) - vecdP[k]/dRiskTot);

                    matH.getvalue(k,k,dTemp,fTemp);
                    matH.setvalue(k,k,dTemp - 
                        adW[i]*vecdP[k]/dRiskTot*(1-vecdP[k]/dRiskTot));
                    for(m=0; m<k; m++)
                    {
                        matH.getvalue(k,m,dTemp,fTemp);
                        dTemp += adW[i]*vecdP[k]/dRiskTot*vecdP[m]/dRiskTot;
                        matH.setvalue(k,m,dTemp);
                        matH.setvalue(m,k,dTemp);
                    }
                }
            }
        }
    }

    /*
    for(k=0; k<K-1; k++)
    {
        for(m=0; m<K-1; m++)
        {
            matH.getvalue(k,m,dTemp,fTemp);
            Rprintf("%f ",dTemp);
        }
        Rprintf("\n");
    }
    */

    // one step to get leaf predictions
    matH.invert();

    for(k=0; k<cTermNodes; k++)
    {
        vecpTermNodes[k]->dPrediction = 0.0;
    }
    for(k=0; k<K-1; k++)
    {
        for(m=0; m<K-1; m++)
        {
            matH.getvalue(k,m,dTemp,fTemp);
            vecpTermNodes[veciK2Node[k]]->dPrediction -= dTemp*vecdG[m];
            dTemp = vecdG[m];
        }
        dTemp = vecpTermNodes[veciK2Node[k]]->dPrediction;
    }
    // vecpTermNodes[veciK2Node[K-1]]->dPrediction = 0.0;

    return hr;
}


double CCoxPH::BagImprovement
(
    double *adT,
    double *adDelta,
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
    double dNum = 0.0;
    double dDen = 0.0;
    double dF = 0.0;
    unsigned long i = 0;

    dNum = 0.0;
    dDen = 0.0;
    for(i=0; i<nTrain; i++)
    {
        if(!afInBag[i])
        {
            dNum += adWeight[i]*exp(dF + dStepSize*adFadj[i]);
            dDen += adWeight[i]*exp(dF);
            if(adDelta[i]==1.0)
            {
                dReturnValue += 
                    adWeight[i]*(dStepSize*adFadj[i] - log(dNum) + log(dDen));
            }
        }
    }

    return dReturnValue;
}

