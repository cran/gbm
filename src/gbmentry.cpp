// GBM by Greg Ridgeway  Copyright (C) 2003

#include "gbm.h"

extern "C" {

#include <R.h>
#include <Rinternals.h>


SEXP gbm
(
    SEXP radY,       // outcome or response
    SEXP radOffset,  // offset for f(x), NA for no offset
    SEXP radX,        
    SEXP raiXOrder,        
    SEXP radWeight,
    SEXP radMisc,   // other row specific data (eg failure time), NA=no Misc
    SEXP rcRows,
    SEXP rcCols,
    SEXP racVarClasses,
    SEXP ralMonotoneVar,
    SEXP rszFamily, 
    SEXP rcTrees,
    SEXP rcDepth,       // interaction depth
    SEXP rcMinObsInNode,
    SEXP rdShrinkage,
    SEXP rdBagFraction,
    SEXP rcTrain,
    SEXP radFOld,
    SEXP rcCatSplitsOld,
    SEXP rcTreesOld
)
{
    unsigned long hr = 0;

    SEXP rAns = NULL;
    SEXP rNewTree = NULL;
    SEXP riSplitVar = NULL;
    SEXP rdSplitPoint = NULL;
    SEXP riLeftNode = NULL;
    SEXP riRightNode = NULL;
    SEXP riMissingNode = NULL;
    SEXP rdErrorReduction = NULL;
    SEXP rdWeight = NULL;

    SEXP rdInitF = NULL;
    SEXP radF = NULL;
    SEXP radTrainError = NULL;
    SEXP radValidError = NULL;
    SEXP radOOBagImprove = NULL;

    SEXP rSetOfTrees = NULL;
    SEXP rSetSplitCodes = NULL;
    SEXP rSplitCode = NULL;

    VEC_VEC_CATEGORIES vecSplitCodes;

    int i = 0;
    int iT = 0;
    int cTrees = INTEGER(rcTrees)[0];
    const int cResultComponents = 7;
    // rdInitF, radF, radTrainError, radValidError, radOOBagImprove
    // rSetOfTrees, rSetSplitCodes
    const int cTreeComponents = 7;
    // riSplitVar, rdSplitPoint, riLeftNode,
    // riRightNode, riMissingNode, rdErrorReduction, rdWeight
    int cNodes = 0;
    int cTrain = INTEGER(rcTrain)[0];

    double dTrainError = 0.0;
    double dValidError = 0.0;
    double dOOBagImprove = 0.0;

    CGBM *pGBM = NULL;
    CDataset *pData = NULL;
    CDistribution *pDist = NULL;

    // set up the dataset
    pData = new CDataset(); new_item(pData);
    if(pData==NULL)
    {
        hr = E_OUTOFMEMORY;
        goto Error;
    }

    // initialize R's random number generator
    GetRNGstate();

    // initialize some things
    hr = gbm_setup(REAL(radY),
                   REAL(radOffset),
                   REAL(radX),
                   INTEGER(raiXOrder),
                   REAL(radWeight),
                   REAL(radMisc),
                   INTEGER(rcRows)[0],
                   INTEGER(rcCols)[0],
                   INTEGER(racVarClasses),
                   INTEGER(ralMonotoneVar),
                   CHAR(STRING_ELT(rszFamily,0)),
                   INTEGER(rcTrees)[0],
                   INTEGER(rcDepth)[0],
                   INTEGER(rcMinObsInNode)[0],
                   REAL(rdShrinkage)[0],
                   REAL(rdBagFraction)[0],
                   INTEGER(rcTrain)[0],
                   pData,
                   pDist);
    if(FAILED(hr))
    {
        goto Error;
    }
        
    // allocate the GBM
    pGBM = new CGBM(); new_item(pGBM);
    if(pGBM==NULL)
    {
        hr = E_OUTOFMEMORY;
        goto Error;
    }

    // initialize the GBM
    hr = pGBM->Initialize(pData,
                          pDist,
                          REAL(rdShrinkage)[0], 
                          cTrain, 
                          REAL(rdBagFraction)[0],
                          INTEGER(rcDepth)[0],
                          INTEGER(rcMinObsInNode)[0]);
    if(FAILED(hr))
    {
        goto Error;
    }

    // allocate the main return object
    PROTECT(rAns = allocVector(VECSXP, cResultComponents));

    // allocate the initial value
    PROTECT(rdInitF = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(rAns,0,rdInitF);
    UNPROTECT(1); // rdInitF

    // allocate the predictions
    PROTECT(radF = allocVector(REALSXP, pData->cRows));
    SET_VECTOR_ELT(rAns,1,radF);
    UNPROTECT(1); // radF

    if(ISNAN(REAL(radFOld)[0])) // check for old predictions
    {
        // set the initial value of F as a constant
        hr = pDist->InitF(pData->adY,
                          pData->adMisc,
                          pData->adOffset,
                          pData->adWeight,
                          REAL(rdInitF)[0], 
                          cTrain);
        for(i=0; i < pData->cRows; i++)
        {
            REAL(radF)[i] = REAL(rdInitF)[0];
        }
    }
    else
    {
        for(i=0; i < pData->cRows; i++)
        {
            REAL(radF)[i] = REAL(radFOld)[i];
        }
    }

    // allocate space for the performance measures
    PROTECT(radTrainError = allocVector(REALSXP, cTrees));
    PROTECT(radValidError = allocVector(REALSXP, cTrees));
    PROTECT(radOOBagImprove = allocVector(REALSXP, cTrees));
    SET_VECTOR_ELT(rAns,2,radTrainError);
    SET_VECTOR_ELT(rAns,3,radValidError);
    SET_VECTOR_ELT(rAns,4,radOOBagImprove);
    UNPROTECT(3); // radTrainError , radValidError, radOOBagImprove

    // allocate the component for the tree structures
    PROTECT(rSetOfTrees = allocVector(VECSXP, cTrees));
    SET_VECTOR_ELT(rAns,5,rSetOfTrees);
    UNPROTECT(1); // rSetOfTrees

    Rprintf("Iter        TrainLL      ValidLL     StepSize      Improve\n");
    for(iT=0; iT<cTrees; iT++)
    {
        hr = pGBM->iterate(REAL(radF),
                           dTrainError,dValidError,dOOBagImprove,
                           cNodes); 
        // store the performance measures
        REAL(radTrainError)[iT] = dTrainError;
        REAL(radValidError)[iT] = dValidError;
        REAL(radOOBagImprove)[iT] = dOOBagImprove;

        // allocate the new tree component for the R list structure
        PROTECT(rNewTree = allocVector(VECSXP, cTreeComponents));

        PROTECT(riSplitVar = allocVector(INTSXP, cNodes));
        PROTECT(rdSplitPoint = allocVector(REALSXP, cNodes));
        PROTECT(riLeftNode = allocVector(INTSXP, cNodes));
        PROTECT(riRightNode = allocVector(INTSXP, cNodes));
        PROTECT(riMissingNode = allocVector(INTSXP, cNodes));
        PROTECT(rdErrorReduction = allocVector(REALSXP, cNodes));
        PROTECT(rdWeight = allocVector(REALSXP, cNodes));
        SET_VECTOR_ELT(rNewTree,0,riSplitVar);
        SET_VECTOR_ELT(rNewTree,1,rdSplitPoint);
        SET_VECTOR_ELT(rNewTree,2,riLeftNode);
        SET_VECTOR_ELT(rNewTree,3,riRightNode);
        SET_VECTOR_ELT(rNewTree,4,riMissingNode);
        SET_VECTOR_ELT(rNewTree,5,rdErrorReduction);
        SET_VECTOR_ELT(rNewTree,6,rdWeight);
        UNPROTECT(cTreeComponents); // riNodeID,riSplitVar,rdSplitPoint,riLeftNode,
                                    // riRightNode,riMissingNode,rdErrorReduction,rdWeight
        SET_VECTOR_ELT(rSetOfTrees,iT,rNewTree);
        UNPROTECT(1); // rNewTree

        hr = gbm_transfer_to_R(pGBM,
                               vecSplitCodes,
                               INTEGER(riSplitVar),
                               REAL(rdSplitPoint),
                               INTEGER(riLeftNode),
                               INTEGER(riRightNode),
                               INTEGER(riMissingNode),
                               REAL(rdErrorReduction),
                               REAL(rdWeight),
                               INTEGER(rcCatSplitsOld)[0]);

        if((iT <= 9) || 
           ((iT+1+INTEGER(rcTreesOld)[0])/100 == 
            (iT+1+INTEGER(rcTreesOld)[0])/100.0) || 
            (iT==cTrees-1))
        {
            Rprintf("%6d %12.4f %12.4f %12.4f %12.4f\n",
                    iT+1+INTEGER(rcTreesOld)[0],
                    REAL(radTrainError)[iT],
                    REAL(radValidError)[iT],
                    REAL(rdShrinkage)[0],
                    REAL(radOOBagImprove)[iT]);
        }
    }
    Rprintf("\n");

    // transfer categorical splits to R
    PROTECT(rSetSplitCodes = allocVector(VECSXP, vecSplitCodes.size()));
    SET_VECTOR_ELT(rAns,6,rSetSplitCodes);
    UNPROTECT(1); // rSetSplitCodes

    for(i=0; i<(int)vecSplitCodes.size(); i++)
    {
        PROTECT(rSplitCode = 
                    allocVector(INTSXP, size_of_vector(vecSplitCodes,i)));
        SET_VECTOR_ELT(rSetSplitCodes,i,rSplitCode);
        UNPROTECT(1); // rSplitCode

        hr = gbm_transfer_catsplits_to_R(i,
                                         vecSplitCodes,
                                         INTEGER(rSplitCode));
    }
    // dump random number generator seed
    #ifdef NOISY_DEBUG
    Rprintf("PutRNGstate\n");
    #endif
    PutRNGstate();

Cleanup:
    UNPROTECT(1); // rAns
    #ifdef NOISY_DEBUG
    Rprintf("destructing\n");
    #endif

    if(pGBM != NULL)
    {
        delete pGBM;
        pGBM = NULL;
    }
    if(pDist != NULL)
    {
        delete pDist;
        pDist = NULL;
    }
    if(pData != NULL)
    {
        delete pData;
        pData = NULL;
    }

    return rAns;
Error:
    goto Cleanup;
}


SEXP gbm_pred
(
    SEXP radX,
    SEXP rcRows,
    SEXP rcCols,
    SEXP rcTrees,
    SEXP rdInitF,
    SEXP rTrees,
    SEXP rCSplits,
    SEXP raiVarType
)
{
    unsigned long hr = 0;
    int iTree = 0;
    int iObs = 0;
    int cRows = INTEGER(rcRows)[0];
    int cTrees = INTEGER(rcTrees)[0];

    SEXP rThisTree = NULL;
    int *aiSplitVar = NULL;
    double *adSplitCode = NULL;
    int *aiLeftNode = NULL;
    int *aiRightNode = NULL;
    int *aiMissingNode = NULL;
    int iCurrentNode = 0;
    double dX = 0.0;
    int iCatSplitIndicator = 0;

    SEXP radPredF = NULL;

    // allocate the predictions to return
    PROTECT(radPredF = allocVector(REALSXP, cRows));
    if(radPredF == NULL)
    {
        hr = E_OUTOFMEMORY;
        goto Error;
    }

    for(iObs=0; iObs<cRows; iObs++)
    {
        REAL(radPredF)[iObs] = REAL(rdInitF)[0];
    }

    for(iTree=0; iTree<cTrees; iTree++)
    {
        rThisTree     = VECTOR_ELT(rTrees,iTree);
        aiSplitVar    = INTEGER(VECTOR_ELT(rThisTree,0));
        adSplitCode   = REAL   (VECTOR_ELT(rThisTree,1));
        aiLeftNode    = INTEGER(VECTOR_ELT(rThisTree,2));
        aiRightNode   = INTEGER(VECTOR_ELT(rThisTree,3));
        aiMissingNode = INTEGER(VECTOR_ELT(rThisTree,4));
        for(iObs=0; iObs<cRows; iObs++)
        {
            iCurrentNode = 0;
            while(aiSplitVar[iCurrentNode] != -1)
            {
                dX = REAL(radX)[aiSplitVar[iCurrentNode]*cRows + iObs];
                // missing?
                if(ISNAN(dX))
                {
                    iCurrentNode = aiMissingNode[iCurrentNode];
                }
                // continuous?
                else if(INTEGER(raiVarType)[aiSplitVar[iCurrentNode]] == 0)
                {
                    if(dX < adSplitCode[iCurrentNode])
                    {
                        iCurrentNode = aiLeftNode[iCurrentNode];
                    }
                    else
                    {
                        iCurrentNode = aiRightNode[iCurrentNode];
                    }
                }
                else // categorical
                {
                    iCatSplitIndicator = INTEGER(
                        VECTOR_ELT(rCSplits,
                                   (int)adSplitCode[iCurrentNode]))[(int)dX];
                    if(iCatSplitIndicator==-1)
                    {
                        iCurrentNode = aiLeftNode[iCurrentNode];
                    }
                    else if(iCatSplitIndicator==1)
                    {
                        iCurrentNode = aiRightNode[iCurrentNode];
                    }
                    else
                    {
                        // handle unused level
                    }
                }
            }
            REAL(radPredF)[iObs] += adSplitCode[iCurrentNode]; // add the prediction
        } // iObs
    } // iTree

Cleanup:
    UNPROTECT(1); // radPredF
    return radPredF;
Error:
    goto Cleanup;
}



SEXP gbm_plot
(
    SEXP radX,        // vector or matrix of points to make predictions
    SEXP rcRows,      // number of rows in X
    SEXP rcCols,      // number of columns in X
    SEXP raiWhichVar, // length=cCols, index of which var cols of X are
    SEXP rcTrees,     // number of trees to use
    SEXP rdInitF,     // initial value
    SEXP rTrees,      // tree list object
    SEXP rCSplits,    // categorical split list object
    SEXP raiVarType   // vector of variable types
)
{
    unsigned long hr = 0;
    int i = 0;
    int iTree = 0;
    int iObs = 0;
    int cRows = INTEGER(rcRows)[0];
    int cCols = INTEGER(rcCols)[0];
    int cTrees = INTEGER(rcTrees)[0];

    SEXP rThisTree = NULL;
    int *aiSplitVar = NULL;
    double *adSplitCode = NULL;
    int *aiLeftNode = NULL;
    int *aiRightNode = NULL;
    double *adW = NULL;
    int iCurrentNode = 0;
    double dCurrentW = 0.0;
    double dX = 0.0;
    int iCatSplitIndicator = 0;

    SEXP radPredF = NULL;
    int aiNodeStack[40];
    double adWeightStack[40];
    int cStackNodes = 0;
    int iPredVar = 0;

    // allocate the predictions to return
    PROTECT(radPredF = allocVector(REALSXP, cRows));
    if(radPredF == NULL)
    {
        hr = E_OUTOFMEMORY;
        goto Error;
    }
    for(iObs=0; iObs<cRows; iObs++)
    {
        REAL(radPredF)[iObs] = REAL(rdInitF)[0];
    }
    for(iTree=0; iTree<cTrees; iTree++)
    {
        rThisTree   = VECTOR_ELT(rTrees,iTree);
        aiSplitVar  = INTEGER(VECTOR_ELT(rThisTree,0));
        adSplitCode = REAL   (VECTOR_ELT(rThisTree,1));
        aiLeftNode  = INTEGER(VECTOR_ELT(rThisTree,2));
        aiRightNode = INTEGER(VECTOR_ELT(rThisTree,3));
        adW         = REAL   (VECTOR_ELT(rThisTree,6));
        for(iObs=0; iObs<cRows; iObs++)
        {
            aiNodeStack[0] = 0;
            adWeightStack[0] = 1.0;
            cStackNodes = 1;
            while(cStackNodes > 0)
            {
                cStackNodes--;
                iCurrentNode = aiNodeStack[cStackNodes];

                if(aiSplitVar[iCurrentNode] == -1) // terminal node
                {
                    REAL(radPredF)[iObs] += 
                        adWeightStack[cStackNodes]*adSplitCode[iCurrentNode];
                }
                else // non-terminal node
                {
                    // which split variable and am I interested in it?
                    iPredVar = -1;
                    for(i=0; (iPredVar == -1) && (i < cCols); i++)
                    {
                        if(INTEGER(raiWhichVar)[i] == aiSplitVar[iCurrentNode])
                        {
                            iPredVar = i; // split is on one that interests me
                        }
                    }

                    if(iPredVar != -1)
                    {
                        dX = REAL(radX)[iPredVar*cRows + iObs];
                        // continuous?
                        if(INTEGER(raiVarType)[aiSplitVar[iCurrentNode]] == 0)
                        {
                            if(dX < adSplitCode[iCurrentNode])
                            {
                                aiNodeStack[cStackNodes] = aiLeftNode[iCurrentNode];
                                cStackNodes++;
                            }
                            else
                            {
                                aiNodeStack[cStackNodes] = aiRightNode[iCurrentNode];
                                cStackNodes++;
                            }
                        }
                        else // categorical
                        {
                            iCatSplitIndicator = INTEGER(
                                VECTOR_ELT(rCSplits,
                                           (int)adSplitCode[iCurrentNode]))[(int)dX];
                            if(iCatSplitIndicator==-1)
                            {
                                aiNodeStack[cStackNodes] = aiLeftNode[iCurrentNode];
                                cStackNodes++;
                            }
                            else if(iCatSplitIndicator==1)
                            {
                                aiNodeStack[cStackNodes] = aiRightNode[iCurrentNode];
                                cStackNodes++;
                            }
                            else
                            {
                                // handle unused level... maybe nothing
                            }
                        }
                    } // iPredVar != -1
                    else // not interested in this split, average left and right 
                    {
                        aiNodeStack[cStackNodes] = aiRightNode[iCurrentNode];
                        dCurrentW = adWeightStack[cStackNodes];
                        adWeightStack[cStackNodes] = dCurrentW *
                            adW[aiRightNode[iCurrentNode]]/
                            (adW[aiLeftNode[iCurrentNode]]+
                             adW[aiRightNode[iCurrentNode]]);
                        cStackNodes++;
                        aiNodeStack[cStackNodes] = aiLeftNode[iCurrentNode];
                        adWeightStack[cStackNodes] = 
                            dCurrentW-adWeightStack[cStackNodes-1];
                        cStackNodes++;
                    }
                } // non-terminal node
            } // while(cStackNodes > 0)
        } // iObs
    } // iTree

Cleanup:
    UNPROTECT(1); // radPredF
    return radPredF;
Error:
    goto Cleanup;
} // gbm_plot


} // end extern "C"