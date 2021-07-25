#include "Include/LevenbergMarquardt.h"
#include "Include/ReadSONG.h"
#include "InputFunction/FuncDataTest.h"

int main()
{
    size_t lengths[6] = {125, 7875, 157500, 1375, 11, 6}, offset[5] = {17436, 17936, 49436, 680144, 685712};
    char *file_name = "InputData/sources_song_z001.dat";

    struct DATA Data;

    Data = ReadBinaryCDMKernelSONG(lengths, offset, file_name, NULL);

    double params[] = {0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0};
    int len_par = 13;

    LevenbergMarquardt(KernelGR, Data.Data, Data.DataSize, Data.VariableSize, params, len_par);

    return 0;
}
