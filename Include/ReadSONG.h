#ifndef READSONG_H_INCLUDED
#define READSONG_H_INCLUDED

/**
    Header file for "ReadSONG.c"
**/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

struct DATA {
    double *Data;
    size_t DataSize;
    int VariableSize;
};

struct DATA ReadBinaryCDMKernelSONG(
    size_t *ArraySizes,
    size_t *PositionByte,
    char *FileName,
    int *RandomSampleSize);

#endif // READSONG_H_INCLUDED
