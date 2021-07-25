#ifndef LEVENBERGMARQUARDT_H_INCLUDED
#define LEVENBERGMARQUARDT_H_INCLUDED

/**
    Header file for Levenberg Marquardt.c
**/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

struct PARAMSCORRECTION{
    double *NewParams;

    double NewLoss;
    double LossChange;
    double Gain;
    double GradientNorm;
    double RelativeParamsChange;
};


void LevenbergMarquardt(
    double (*ModelFunc)(double*, double* , double*),
    double *Data,
    size_t DataSize,
    int VariableSize,
    double *Params,
    int ParamsSize);

struct PARAMSCORRECTION ParamsCorrection(
    double (*ModelFunc)(double*, double*, double*),
    double *Data,
    size_t DataSize,
    int VariableSize,
    double *Params,
    int ParamsSize,
    double Hessian[][ParamsSize],
    double Gradient[],
    double Lambda[],
    double MinimalLoss,
    int DerivativeFlag);

double *QRDecomposition(
    int ParamsSize,
    double M[][ParamsSize],
    double Gradient[]);


#endif // LEVENBERGMARQUARDT_H_INCLUDED
