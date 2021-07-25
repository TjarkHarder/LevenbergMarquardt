#include "../Include/LevenbergMarquardt.h"
#include "../Include/ReadSONG.h"

double KernelGR(double *Data, double *Params, double *Jacobian){

    double F = Data[0], k1 = pow(10, 2)*Data[1], k2 = pow(10, 2)*Data[2], k3 = pow(10, 2)*Data[3];
    double f[4], k1k2;

    k1k2 = (pow(k3, 2) - (pow(k1, 2) + pow(k2, 2)))/(2*k1*k2);

    f[0] = 1.0;
    f[1] = k1k2 * (k1/k2 + k2/k1);
    f[2] = pow(k1k2, 2);
    f[3] = pow((k1/k2 - k2/k1), 2);


    double Residue = F - ((Params[1]/pow(k3, 0)+Params[2]/pow(k3, 2)+Params[3]/pow(k3, 4))*f[0]
                        + (Params[4]/pow(k3, 0)+Params[5]/pow(k3, 2)+Params[6]/pow(k3, 4))*f[1]
                        + (Params[7]/pow(k3, 0)+Params[8]/pow(k3, 2)+Params[9]/pow(k3, 4))*f[2]
                        + (Params[10]/pow(k3, 0)+Params[11]/pow(k3, 2)+Params[12]/pow(k3, 4))*f[3])
                        / ((1+Params[0]/pow(k1, 2)) * (1+Params[0]/pow(k2, 2)));

    if (Jacobian != NULL){

        Jacobian[0] = (F - Residue) * (1/(pow(k1, 2) + Params[0]) + 1/(pow(k2, 2) + Params[0]));

        Jacobian[1] = -f[0] / pow(k3, 0) / ((1+Params[0]/pow(k1, 2)) * (1+Params[0]/pow(k2, 2)));
        Jacobian[2] = -f[0] / pow(k3, 2) / ((1+Params[0]/pow(k1, 2)) * (1+Params[0]/pow(k2, 2)));
        Jacobian[3] = -f[0] / pow(k3, 4) / ((1+Params[0]/pow(k1, 2)) * (1+Params[0]/pow(k2, 2)));

        Jacobian[4] = -f[1] / pow(k3, 0) / ((1+Params[0]/pow(k1, 2)) * (1+Params[0]/pow(k2, 2)));
        Jacobian[5] = -f[1] / pow(k3, 2) / ((1+Params[0]/pow(k1, 2)) * (1+Params[0]/pow(k2, 2)));
        Jacobian[6] = -f[1] / pow(k3, 4) / ((1+Params[0]/pow(k1, 2)) * (1+Params[0]/pow(k2, 2)));

        Jacobian[7] = -f[2] / pow(k3, 0) / ((1+Params[0]/pow(k1, 2)) * (1+Params[0]/pow(k2, 2)));
        Jacobian[8] = -f[2] / pow(k3, 2) / ((1+Params[0]/pow(k1, 2)) * (1+Params[0]/pow(k2, 2)));
        Jacobian[9] = -f[2] / pow(k3, 4) / ((1+Params[0]/pow(k1, 2)) * (1+Params[0]/pow(k2, 2)));

        Jacobian[10] = -f[3] / pow(k3, 0) / ((1+Params[0]/pow(k1, 2)) * (1+Params[0]/pow(k2, 2)));
        Jacobian[11] = -f[3] / pow(k3, 2) / ((1+Params[0]/pow(k1, 2)) * (1+Params[0]/pow(k2, 2)));
        Jacobian[12] = -f[3] / pow(k3, 4) / ((1+Params[0]/pow(k1, 2)) * (1+Params[0]/pow(k2, 2)));
    }

    return Residue;
}




int main()
{
    size_t ArraySizes[6] = {125, 7875, 157500, 1375, 11, 6}, PositionByte[5] = {17436, 17936, 49436, 680144, 685712};
    char *FileName = "Data/sources_song_z000.dat";

    struct DATA Data;

    Data = ReadBinaryCDMKernelSONG(ArraySizes, PositionByte, FileName, NULL);

    double Params[] = {0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0};
    int ParamsSize = 13;

    LevenbergMarquardt(KernelGR, Data.Data, Data.DataSize, Data.VariableSize, Params, ParamsSize);

    return 0;
}
