#include "../Include/ReadSONG.h"

int *InverseIndexTransform1(int index){

    int *index_k = malloc(sizeof(int)*2);

    double index_;

    index_ = (double) index;

    while (1 == 1){
        double a = (-3 + pow(8*index_+9, 0.5))/2;
        if (ceil(a) == a){
            index_k[0] = (int) a;
            break;
        }
        index_ += 1;
    }

    index_k[1] = index - index_k[0]*(index_k[0]+1)/2;

    return index_k;
}


int *InverseIndexTransform2(int k3_size[], int index){

    int *index_ = malloc(sizeof(int)*2);

    int i = 0, n = k3_size[0];

    while (1 == 1){
        if (index < n){
            index_[0] = i;
            break;
        }
        i++;
        n += k3_size[i];
    }

    index_[1] = index - n + k3_size[i];

    return index_;
}


struct DATA ReadBinaryCDMKernelSONG(size_t *ArraySizes, size_t *PositionByte, char *FileName, int *RandomSampleSize){

    struct DATA DataSONG;
    int *RandomIndices;                  // Sorted array of random indices of length "RandomSampleSize"

    DataSONG.VariableSize = 3;                              // three variables k1, k2 and k3

    if (RandomSampleSize == NULL){
        DataSONG.DataSize = ArraySizes[2];
        RandomIndices = NULL;
    }
    else {
        DataSONG.DataSize = *RandomSampleSize;
        RandomIndices = malloc(sizeof(int) * DataSONG.DataSize);
    }

    DataSONG.Data = malloc(sizeof(double) * DataSONG.DataSize * (DataSONG.VariableSize + 1));

    float *k = malloc(sizeof(float) * ArraySizes[0]);
    float *k3 = malloc(sizeof(float) * ArraySizes[2]);
    float *DataTemp = malloc(sizeof(float) * ArraySizes[2]);
    float *TransferTemp = malloc(sizeof(float) * ArraySizes[3]);

    int *k3_size = malloc(sizeof(int) * ArraySizes[1]);

    double *Transfer = malloc(sizeof(double) * ArraySizes[0]);


    FILE *fptr;

    fptr = fopen(FileName, "rb");

    fseek(fptr, PositionByte[0], SEEK_SET);
    fread(k, sizeof(float), ArraySizes[0], fptr);

    fseek(fptr, PositionByte[1], SEEK_SET);
    fread(k3_size, sizeof(int), ArraySizes[1], fptr);

    fseek(fptr, PositionByte[2], SEEK_SET);
    fread(k3, sizeof(float), ArraySizes[2], fptr);

    fseek(fptr, PositionByte[3], SEEK_SET);
    fread(TransferTemp, sizeof(float), ArraySizes[3], fptr);

    fseek(fptr, PositionByte[4], SEEK_SET);
    fread(DataTemp, sizeof(float), ArraySizes[2], fptr);

    fclose(fptr);

    int I, j = 0;

    for (int i = 0; i < ArraySizes[3]; i++){
            if (i % ArraySizes[4] == ArraySizes[5]){
                Transfer[j] = (double) TransferTemp[i];
                j++;
            }
    }

    if (RandomSampleSize != NULL){
        time_t t = time(NULL);
        srand(t);

        for (int i = 0; i < DataSONG.DataSize; i++){
            int Flag = 0;

            while (1 == 1){
                I = 0;

                for (int m = 0; m < 3*(ArraySizes[2] - ArraySizes[2]%32767)/32767; m++)
                    I += rand();

                I = I % ArraySizes[2];

                if (i == 0)
                    break;

                for (int m = 0; m < i; m++){
                    if (I == RandomIndices[m])
                        break;
                    if (I < RandomIndices[m]){

                        for (int n = i; n > m; n--){
                            int Temp = RandomIndices[n-1];
                            RandomIndices[n] = Temp;
                        }

                        RandomIndices[m] = I;
                        Flag = 1;
                        break;
                    }

                    if (m == i - 1){
                        RandomIndices[i] = I;
                        Flag = 1;
                    }
                }

                if (Flag)
                    break;
            }
        }
    }

    j = 0;

    for (int i = 0; i < DataSONG.DataSize; i++){

        if (RandomSampleSize == NULL)
            I = i;
        else {
            I = RandomIndices[i];

            if (RandomIndices[i] > ArraySizes[2]-1)
                break;
        }

        int *index2 = InverseIndexTransform2(k3_size, I);
        int *index1 = InverseIndexTransform1(index2[0]);

        DataSONG.Data[j] = (double) DataTemp[I];
        DataSONG.Data[j] *= 1/(Transfer[index1[0]]*Transfer[index1[1]]);

        DataSONG.Data[j+1] = (double) k[index1[0]];
        DataSONG.Data[j+2] = (double) k[index1[1]];
        DataSONG.Data[j+3] = (double) k3[I];

        j += DataSONG.VariableSize + 1;
    }

    free(k3);
    free(k3_size);
    free(DataTemp);
    free(TransferTemp);
    free(Transfer);
    free(k);
    free(RandomIndices);

    return DataSONG;
};
