//------------------------------------------------------------------------------
//  GBM by Greg Ridgeway  Copyright (C) 2003
//
//  File:       node_search.h
//
//  License:    GNU GPL (version 2 or later)
//
//  Contents:   does the searching for where to split a node
//        	  
//  Owner:      gregr@rand.org
//
//  History:    3/26/2001   gregr created
//              2/14/2003   gregr: adapted for R implementation
//
//------------------------------------------------------------------------------

#ifndef NODESEARCH_H
#define NODESEARCH_H

#include "dataset.h"
#include "node_factory.h"

using namespace std;

class CNodeSearch
{
public:

    CNodeSearch();
    ~CNodeSearch();
    HRESULT Initialize(unsigned long cMinObsInNode);

    HRESULT IncorporateObs
    (
        double dX,
        double dZ,
        double dW,
        long lMonotone
    )
    {
        HRESULT hr = S_OK;
        static double dWZ = 0.0;

        if(fIsSplit) goto Cleanup;

        dWZ = dW*dZ;

        if(ISNAN(dX))
        {
            dCurrentMissingSumZ += dWZ;
            dCurrentMissingTotalW += dW;
            cCurrentMissingN++;
            dCurrentRightSumZ -= dWZ;
            dCurrentRightTotalW -= dW;
            cCurrentRightN--;
        }
        else if(cCurrentVarClasses == 0)   // variable is continuous
        {
            assert(dLastXValue <= dX);

            // Evaluate the current split
            // the newest observation is still in the right child
            dCurrentSplitValue = 0.5*(dLastXValue + dX);
            if((dLastXValue != dX) && 
               (cCurrentLeftN >= cMinObsInNode) && 
               (cCurrentRightN >= cMinObsInNode) &&
               ((lMonotone==0) ||
                (lMonotone*(dCurrentRightSumZ*dCurrentLeftTotalW - 
                            dCurrentLeftSumZ*dCurrentRightTotalW) > 0)))
            {
                dCurrentImprovement = 
                    CNode::Improvement(dCurrentLeftTotalW,dCurrentRightTotalW,
                                       dCurrentMissingTotalW,
                                       dCurrentLeftSumZ,dCurrentRightSumZ,
                                       dCurrentMissingSumZ);
                if(dCurrentImprovement > dBestImprovement)
                {
                    iBestSplitVar = iCurrentSplitVar;
                    dBestSplitValue = dCurrentSplitValue;
                    cBestVarClasses = 0;

                    dBestLeftSumZ    = dCurrentLeftSumZ;
                    dBestLeftTotalW  = dCurrentLeftTotalW;
                    cBestLeftN       = cCurrentLeftN;
                    dBestRightSumZ   = dCurrentRightSumZ;
                    dBestRightTotalW = dCurrentRightTotalW;
                    cBestRightN      = cCurrentRightN;
                    dBestImprovement = dCurrentImprovement;
                }
            }

            // now move the new observation to the left
            // if another observation arrives we will evaluate this
            dCurrentLeftSumZ += dWZ;
            dCurrentLeftTotalW += dW;
            cCurrentLeftN++;
            dCurrentRightSumZ -= dWZ;
            dCurrentRightTotalW -= dW;
            cCurrentRightN--;

            dLastXValue = dX;
        }
        else // variable is categorical, evaluates later
        {
            adGroupSumZ[(unsigned long)dX] += dWZ;
            adGroupW[(unsigned long)dX] += dW;
            acGroupN[(unsigned long)dX] ++;
        }

    Cleanup:
        return hr;
    }

    HRESULT Set(double dSumZ,
                double dTotalW,
                unsigned long cTotalN,
                CNodeTerminal *pThisNode,
                CNode **ppParentPointerToThisNode,
                CNodeFactory *pNodeFactory);
    HRESULT ResetForNewVar(unsigned long iWhichVar,
                           long cVarClasses);

    double BestImprovement() { return dBestImprovement; }
    HRESULT SetToSplit() 
    {    
        fIsSplit = true;
        return S_OK;
    };
    HRESULT SetupNewNodes(PCNodeNonterminal &pNewSplitNode,
                          PCNodeTerminal &pNewLeftNode,
                          PCNodeTerminal &pNewRightNode,
                          PCNodeTerminal &pNewMissingNode);

    HRESULT EvaluateCategoricalSplit();
    HRESULT WrapUpCurrentVariable();
    double ThisNodePrediction() {return pThisNode->dPrediction;}
    bool operator<(const CNodeSearch &ns) {return dBestImprovement<ns.dBestImprovement;}

    unsigned long iBestSplitVar;
    double dBestSplitValue;

    double dBestLeftSumZ;
    double dBestLeftTotalW;
    unsigned long cBestLeftN;

    double dBestRightSumZ;
    double dBestRightTotalW;
    unsigned long cBestRightN;

    double dBestMissingSumZ;
    double dBestMissingTotalW;
    unsigned long cBestMissingN;

    double dCurrentMissingSumZ;
    double dCurrentMissingTotalW;
    unsigned long cCurrentMissingN;

    long cCurrentVarClasses;

    unsigned long iRank;
    double dInitTotalW;
    double dInitSumZ;
    unsigned long cInitN;
    double dBestImprovement;

private:
    bool fIsSplit;

    unsigned long cMinObsInNode;

    long cBestVarClasses;

    double dCurrentLeftSumZ;
    double dCurrentLeftTotalW;
    unsigned long cCurrentLeftN;
    double dCurrentRightSumZ;
    double dCurrentRightTotalW;
    unsigned long cCurrentRightN;
    double dCurrentImprovement;
    unsigned long iCurrentSplitVar;
    double dCurrentSplitValue;

    double dLastXValue;

    double *adGroupSumZ;
    double *adGroupW;
    unsigned long *acGroupN;
    double *adGroupMean;
    int *aiCurrentCategory;
    unsigned long *aiBestCategory;
    const unsigned long k_cMaxClasses;

    CNodeTerminal *pThisNode;
    CNode **ppParentPointerToThisNode;
    CNodeFactory *pNodeFactory;
};

typedef CNodeSearch *PCNodeSearch;

#endif // NODESEARCH_H



