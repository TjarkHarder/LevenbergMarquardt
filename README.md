# LevenbergMarquardt

I have written this Levenberg-Marquardt algorithm in C during my bachelor thesis, to approximate the numerical kernel of the CDM intrinsic bispectrum in second order perturbation theory from SONG (see https://github.com/coccoinomane/song), by using separable source functions. 

I have also included the file "ReadSONG.c", which reads the output file from SONG and outputs the needed arrays for the Levenberg-Marquardt algorithm "LevenbergMarquardt.c". 

The example files "KernelGR.c" and "KernelExp.c" are provided, which portray the use of the "LevenbergMarquardt.c" file, as well as obtaining the data from the "sources_song_z000.dat" file using the "ReadSONG.c" file. 
