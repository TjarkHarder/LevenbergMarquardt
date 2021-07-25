#ifndef PARAMSCORRECTION_H_INCLUDED
#define PARAMSCORRECTION_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stddef.h>

double *ParamsCorrection(
    double (*ModelFunc)(double*, double*, double*),
    double *Data,
    size_t DataSize,
    int VariableSize,
    double *Params,
    int ParamsSize,
    double Hessian[][ParamsSize],
    double Gradient[],
    struct LEVENBERGMARQUARDT Args);

double *QRDecomposition(
    int ParamsSize,
    double M[][ParamsSize],
    double Gradient[]);

#endif // PARAMSCORRECTION_H_INCLUDED
