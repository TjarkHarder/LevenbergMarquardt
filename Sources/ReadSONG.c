#include "../Include/ReadSONG.h"

int *__k3_size__IndexTo__k__Index(int __k3_size__Index){

    /**
        This function transforms the index of the "k3_size" array __k3_size__Index from SONG into the
        corresponding indices for the "k" array __k__Index[0] and __k__Index[1] (for k1 and k2), such
        that __k__Index[0] >= __k__Index[1].

        To understand this function, consider the first few indices of the "k3_size" array
        and their corresponding "k" indices:

            __k3_size__Index:       0   1   2   3   4   5   6   7   8   9  10 ...
            __k__Index[1]:          0   1   1   2   2   2   3   3   3   3   4 ...
            __k__Index[2]:          0   0   1   0   1   2   0   1   2   3   0 ...

        We see that __k3_size__Index is given by:

            __k3_size__Index = 1/2 __k__Index[0] * (__k__Index[0] + 1) + __k__Index[1] ................ (1)


        We can obtain a lower bound for __k__Index[0] by setting __k__Index[0]_Bound = __k__Index[1]
        and inverting equation (1):

            __k__Index[0]_Bound = 1/2 * (-3 + sqrt(8 * __k3_size__Index + 9)) ......................... (2)

        Note that:

            8 * __k3_size__Index + 9 = 4 * __k__Index[0] * (__k__Index[0] + 1) + 8 * __k__Index[1] + 9

        is only a perfect square (and equation (2) an integer) if "__k__Index[0] = __k__Index[1]" since
        "0 <= __k__Index[1] <= __k__Index[0]" and therefore:

            4 * ((__k__Index[0] + 1/2)^2 + 2) <=  8 * __k3_size__Index + 9 <= 4 * (__k__Index[0] + 3/2)^2
                <=> 4 * (__k__Index[0] + 1/2)^2 < 8 * __k3_size__Index + 9 <= 4 * (__k__Index[0] + 3/2)^2

        We then have with equation (2) and the inequalities from above:

            __k__Index[0] - 1 < __k__Index[0]_Bound <= __k__Index[0]

        We therefore can obtain __k__Index[0] by taking the ceiling of equation (2). We then obtain
        __k__Index[1] using equation (1).
    **/

    // Declare the storage variable for __k__Index
    int *__k__Index = malloc(sizeof(int)*2);

    // Solve for the indices as described
    __k__Index[0] = (int) ceil((-3 + pow(8.0*__k3_size__Index+9.0, 0.5))/2);
    __k__Index[1] = __k3_size__Index - __k__Index[0]*(__k__Index[0]+1)/2;

    return __k__Index;
}


int __k3__IndexTo__k3_size__Index(int k3_size[], size_t __k3__Index){

    /**
        This function takes the current index of the "k3" or "sources" array and returns the
        corresponding index __k3_size__Index in the "k3_size" array, such that

            __k3__Index = k3_size[0] + ... + k3_size[__k3_size__Index - 1] + __k3_size__I ............. (1)

        where __k3_size__I is the current value in the range from 0 to k3_size[__k3_size__Index] - 1.


        To obtain __k3_size__Index[0], we first note that:

            __k3_size__I < k3_size[__k3_size_Index]

        we therefore obtain an upper bound for __k3__Index by setting the inequality equal:

            __k3__Index < k3_size[0] + ... + k3_size[__k3_size__Index] ................................. (2)

        By adding the "k3_size" values iteratively will at some point fulfill the inequality (2). The
        index of the last added "k3_size" value then is  __k3_size__Index.
    **/

    // Declare the storage variable for  __k3_size__Index
    int __k3_size__Index;

    // Iterate through the indices of "k3_size"
    int IterationIndex = 0;
    size_t __k3_size__Sum = k3_size[0];

    while (1 == 1){
        if (__k3__Index < __k3_size__Sum){
            // __k3_size__Index[0] is found
            __k3_size__Index = IterationIndex;
            break;
        }

        // Update IterationIndex and __k3_size__Sum
        IterationIndex++;
        __k3_size__Sum += k3_size[IterationIndex];
    }

    // Calculate __k3_size__Index[1]
    // Note that the sum contains one additional term "k3_size[__k3_size__Index[0]]", which is not included in equation (1). We therefore have to correct the sum for this term.
    //__k3_size__Index[1] = __k3__Index - __k3_size__Sum + k3_size[__k3_size__Index[0]];

    return __k3_size__Index;
}


struct DATA ReadBinaryCDMKernelSONG(size_t *ArraySizes, size_t *PositionByte, char *FileName, int *RandomSampleSize){

    /**
        This function reads the binary data file from SONG and outputs the data in an array, such that
        every data point is followed by its corresponding variables "k1", "k2" and "k3". Since SONG
        does not output the CDM kernel directly, but multiplied with the first order perturbations:

            sources(k1, k2, k3) = kernel_cdm(k1, k2, k3) * delta_cdm(k1) * delta_cdm(k2)

        this function also also corrects for this.
        If one wishes to only use a small subset of the large data set, the function also chooses a
        random subset of given size, instead of the full set.

        Input parameters:

            ArraySizes                          -  Sizes of the data arrays from SONG. The indices of
                                                   "ArraySizes" correspond to the following arrays from SONG:

                                                      0 -> "k" array
                                                      1 -> "k3_size" array
                                                      2 -> "k3" and "sources" array (both are of the same size)
                                                      3 -> "pvec_sources" array (first order perturbations)
                                                      4 -> number of first order perturbations
                                                      5 -> every 6'th element in "pvec_sources" corresponds to
                                                           the "delta_cdm" first order perturbations

            PositionByte                        - Position of the arrays in the data file in bytes. The indices
                                                  of "PositionByte" correspond to the following arrays from SONG:

                                                      0 -> "k" array
                                                      1 -> "k3_size" array
                                                      2 -> "k3" array
                                                      3 -> "pvec_sources"
                                                      4 -> "sources"

            FileName                            -  Path to to the file

            RandomSampleSize                    -  Size of the random subset. If the full set is to be used, set
                                                   "RandomSampleSize == NULL".


        Output parameters (as structure DATA):

            Data                                -  Data array of the kernel. Every data point in followed by its
                                                   corresponding variables "k1", "k2" and "k3".

            DataSize                            -  Number of data points

            VariableSize                        -  Number of variables (in this case always 3)
    **/


    // Declare storage variables
    struct DATA DataSONG;
    size_t *RandomIndices;                                 // Sorted array of random indices for the "sources" array of length "RandomSampleSize"

    // Three variables k1, k2 and k3
    DataSONG.VariableSize = 3;

    // Check if a subset should be used and set the size of the data array
    if (RandomSampleSize == NULL){
        DataSONG.DataSize = ArraySizes[2];
        RandomIndices = NULL;
    }
    else {
        DataSONG.DataSize = *RandomSampleSize;
        RandomIndices = malloc(sizeof(size_t) * DataSONG.DataSize);
    }

    // Allocate memory to the "Data" array for the output
    DataSONG.Data = malloc(sizeof(double) * DataSONG.DataSize * (DataSONG.VariableSize + 1));

    // Declare and allocate memory for temporary storage variables for the data from SONG
    float *k = malloc(sizeof(float) * ArraySizes[0]);
    float *k3 = malloc(sizeof(float) * ArraySizes[2]);
    float *DataTemp = malloc(sizeof(float) * ArraySizes[2]);
    float *PerturbTemp = malloc(sizeof(float) * ArraySizes[3]);

    int *k3_size = malloc(sizeof(int) * ArraySizes[1]);

    double *Perturb = malloc(sizeof(double) * ArraySizes[0]);

    // Read in the data from the binary file
    FILE *fptr;

    fptr = fopen(FileName, "rb");

    fseek(fptr, PositionByte[0], SEEK_SET);
    fread(k, sizeof(float), ArraySizes[0], fptr);

    fseek(fptr, PositionByte[1], SEEK_SET);
    fread(k3_size, sizeof(int), ArraySizes[1], fptr);

    fseek(fptr, PositionByte[2], SEEK_SET);
    fread(k3, sizeof(float), ArraySizes[2], fptr);

    fseek(fptr, PositionByte[3], SEEK_SET);
    fread(PerturbTemp, sizeof(float), ArraySizes[3], fptr);

    fseek(fptr, PositionByte[4], SEEK_SET);
    fread(DataTemp, sizeof(float), ArraySizes[2], fptr);

    fclose(fptr);

    // Create temporary indices for iterations
    size_t TempIndex1, TempIndex2;

    // Store the CDM first order perturbations from the "pvec_sources" array
    TempIndex1 = 0;
    for (int i = 0; i < ArraySizes[3]; i++){
            if (i % ArraySizes[4] == ArraySizes[5]){
                Perturb[TempIndex1] = (double) PerturbTemp[i];
                TempIndex1++;
            }
    }

    // If a subset is to be used, choose unique random indices in the range from 0 to ArraysSizes[2] - 1 (where ArraySizes[2] is the number of data points for the large set).
    if (RandomSampleSize != NULL){
        // Choose a seed for rand()
        //srand(1627786677);            // This seed has been used to generate the subsets of the large set in my bachelor thesis
        time_t t = time(NULL);
        srand(t);

        // Choose "DataSize" many unique random indices
        for (int i = 0; i < DataSONG.DataSize; i++){
            // Use a flag to break out of the while loop
            int Flag = 0;

            // Use a while loop to choose unique values
            while (1 == 1){
                TempIndex1 = 0;

                // Since "RAND_MAX = 32767", we need to enlarge the range in which random numbers can be generated
                for (int m = 0; m < 3*(ArraySizes[2] - ArraySizes[2]%32767)/32767; m++)
                    TempIndex1 += rand();

                TempIndex1 = TempIndex1 % ArraySizes[2];    // Random value might be larger than ArraySizes[2]

                // No other values are in "RandomIndices", "TempIndex1" is therefore unique and can be included
                if (i == 0){
                    RandomIndices[0] = TempIndex1;
                    break;
                }

                // Check if "TempIndex1" is already in "RandomIndices" and sort "RandomIndices" from lowest to largest
                for (int m = 0; m < i; m++){    // There are a maximum of "i" values in "RandomIndices"

                    // If the value is in the array, return to the while loop
                    if (TempIndex1 == RandomIndices[m])
                        break;

                    /* If the value is smaller than the m'th value of "RandomIndices", shift every large element of
                       "RandomIndices" one space up and insert "TempIndex1". We automatically sort the array this way. */
                    if (TempIndex1 < RandomIndices[m]){

                        for (int n = i; n > m; n--){
                            int Temp = RandomIndices[n-1];
                            RandomIndices[n] = Temp;
                        }

                        RandomIndices[m] = TempIndex1;
                        Flag = 1;
                        break;
                    }

                    // If the value is larger than all elements, place it at the end of the array "RandomIndices"
                    if (m == i - 1){
                        RandomIndices[i] = TempIndex1;
                        Flag = 1;
                    }
                }

                // Break out of the while loop if a value has been accepted
                if (Flag)
                    break;
            }
        }
    }

    // Store the data in the "Data" array
    TempIndex2 = 0;

    for (int i = 0; i < DataSONG.DataSize; i++){

        // Check if a subset is to be used. If true, use the indices stored in "RandomSampleSize"
        if (RandomSampleSize == NULL)
            TempIndex1 = i;
        else
            TempIndex1 = RandomIndices[i];

        // Get the indices __k3_size__Index, __k__Index[0] and __k__Index[1]
        int __k3_size__Index = __k3__IndexTo__k3_size__Index(k3_size, TempIndex1);
        int *__k__Index = __k3_size__IndexTo__k__Index(__k3_size__Index);

        // Store the "sources" data and correct for the first order perturbations
        DataSONG.Data[TempIndex2] = (double) DataTemp[TempIndex1];
        DataSONG.Data[TempIndex2] *= 1/(Perturb[__k__Index[0]]*Perturb[__k__Index[1]]);

        // Store the variables
        DataSONG.Data[TempIndex2+1] = (double) k[__k__Index[0]];
        DataSONG.Data[TempIndex2+2] = (double) k[__k__Index[1]];
        DataSONG.Data[TempIndex2+3] = (double) k3[TempIndex1];

        TempIndex2 += DataSONG.VariableSize + 1;
    }

    // Free allocated memory
    free(k3);
    free(k3_size);
    free(DataTemp);
    free(PerturbTemp);
    free(Perturb);
    free(k);
    free(RandomIndices);

    return DataSONG;
};
