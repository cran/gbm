//------------------------------------------------------------------------------
//  GBM by Greg Ridgeway  Copyright (C) 2003
//  File:       poisson.h
//
//  License:    GNU GPL (version 2 or later)
//
//  Contents:   poisson object
//        	  
//  Owner:      gregr@rand.org
//
//  History:    3/26/2001   gregr created
//              2/14/2003   gregr: adapted for R implementation
//
//------------------------------------------------------------------------------

#ifndef POISSON_H
#define POISSON_H

#include "distribution.h"

class CPoisson : public CDistribution
{

public:

    CPoisson();

    virtual ~CPoisson();

    HRESULT ComputeWorkingResponse(double *adY,
                                   double *adMisc,
                                   double *adOffset,
                                   double *adWeight,
                                   double *adF,
                                   double *adZ,
                                   bool *afInBag,
                                   unsigned long nTrain);

    double LogLikelihood(double *adY,
                         double *adMisc,
                         double *adOffset,
                         double *adWeight,
                         double *adF,
                         unsigned long cLength);

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

#endif // POISSON_H