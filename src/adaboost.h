//------------------------------------------------------------------------------
//  GBM by Greg Ridgeway  Copyright (C) 2003
//
//  File:       adaboost.h
//
//  License:    GNU GPL (version 2 or later)
//
//  Contents:   Object for fitting for the AdaBoost loss function
//            
//  Owner:      gregr@rand.org
//
//  History:    3/26/2001   gregr created
//              2/14/2003   gregr: adapted for R implementation
//
//------------------------------------------------------------------------------

#ifndef ADABOOST_H
#define ADABOOST_H

#include "distribution.h"

class CAdaBoost : public CDistribution
{

public:

    CAdaBoost();

    virtual ~CAdaBoost();

    HRESULT ComputeWorkingResponse(double *adY,
                                   double *adMisc,
                                   double *adOffset,
                                   double *adWeight,
                                   double *adF, 
                                   double *adZ,
                                   bool *afInBag,
                                   unsigned long nTrain);

    HRESULT InitF(double *adY, 
                  double *adMisc,
                  double *adOffset,
                  double *adWeight,
                  double &dInitF, 
                  unsigned long cLength);

    HRESULT FitBestConstant(double *adY,
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
                            double *adFadj);

    double LogLikelihood(double *adY,
                         double *adMisc,
                         double *adOffset,
                         double *adWeight,
                         double *adF,
                         unsigned long cLength);

    double BagImprovement(double *adY,
                          double *adMisc,
                          double *adOffset,
                          double *adWeight,
                          double *adF,
                          double *adFadj,
                          bool *afInBag,
                          double dStepSize,
                          unsigned long nTrain);

private:
    vector<double> vecdNum;
    vector<double> vecdDen;
};

#endif // ADABOOST_H


