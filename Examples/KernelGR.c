#include "../Include/LevenbergMarquardt.h"
#include "../Include/ReadSONG.h"

/**
    Example file using the Levenberg-Marquardt algorithm to fit the GR kernel
    to the numerical kernel from SONG.
**/

double KernelGRResidue(double *Data, double *Params, double *Jacobian){

    double F = Data[0]; // Data point
    double k1 = pow(10, 2)*Data[1], k2 = pow(10, 2)*Data[2], k3 = pow(10, 2)*Data[3];  // Rescale the corresponding variables to obtain convergence; remember to adjust the optimal parameters as they will also be scaled
    double f[4], k1k2;

    k1k2 = (pow(k3, 2) - (pow(k1, 2) + pow(k2, 2)))/(2*k1*k2);

    f[0] = 1.0;
    f[1] = k1k2 * (k1/k2 + k2/k1);
    f[2] = pow(k1k2, 2);
    f[3] = pow((k1/k2 - k2/k1), 2);

    // Calculate the residue
    double Residue = F - ((Params[1]/pow(k3, 0)+Params[2]/pow(k3, 2)+Params[3]/pow(k3, 4))*f[0]
                        + (Params[4]/pow(k3, 0)+Params[5]/pow(k3, 2)+Params[6]/pow(k3, 4))*f[1]
                        + (Params[7]/pow(k3, 0)+Params[8]/pow(k3, 2)+Params[9]/pow(k3, 4))*f[2]
                        + (Params[10]/pow(k3, 0)+Params[11]/pow(k3, 2)+Params[12]/pow(k3, 4))*f[3])
                         / ((1+Params[0]/pow(k1, 2)) * (1+Params[0]/pow(k2, 2)));

    // Calculate the derivatives if needed
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
    // The following arrays contain informations on the location and size of the data arrays for the SONG .dat file

    size_t ArraySizes[6] = {125, 7875, 157500, 1375, 11, 6};            /* Sizes of the arrays from the .dat file from SONG, where the index in "ArraySizes" corresponds to:
                                                                            0 -> "k" array
                                                                            1 -> "k3_size" array
                                                                            2 -> "k3" and "sources" arrays
                                                                            3 -> "pvec_sources" (first order perturbation) array
                                                                            4 -> number of first order perturbations
                                                                            5 -> every 6'th element in "pvec_sources" corresponds to the "delta_cdm" first order perturbations
                                                                            */
    size_t PositionByte[5] = {17436, 17936, 49436, 680144, 685712};     // Position (in bytes) of the first element of the arrays in the .dat file
    char *FileName = "Data/sources_song_z001.dat";                      // Path of the .dat file

    // Declare DATA structure; contains SONG's data, as well as the number of data points and number of variables
    struct DATA Data;

    // Read the data from the binary .dat file
    Data = ReadBinaryCDMKernelSONG(ArraySizes, PositionByte, FileName, NULL);

    /*      // Uncomment and comment the previous "Data = ..." line to use a small subset of given size "Size" instead of the large set
    int Size = 1000;
    Data = ReadBinaryCDMKernelSONG(ArraySizes, PositionByte, FileName, &Size);
    */

    // Initial parameters and number of parameters
    double Params[] = {0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0};
    int ParamsSize = 13;

    // Run the algorithm
    LevenbergMarquardt(KernelGRResidue, Data.Data, Data.DataSize, Data.VariableSize, Params, ParamsSize);
    
    return 0;
}
