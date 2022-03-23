#include "../Include/LevenbergMarquardt.h"

void LevenbergMarquardt(double (*ModelFuncResidue)(double*, double* , double*), double* Data, size_t DataSize, int VariableSize, double *Params, int ParamsSize){

    /**
        This function takes a model function "f(x1, x2, ...; p1, p2, ...)" and a set of data points
        "(F_1, ..., F_n)", and attempts to find the optimal parameters, such that the loss
        "S = 1/2 (r_1^2 + ... + r_n^2)", where "r_i = F_i - f(x1_i, x2_i, ...; p1, p2, ...)" is the
        i'th residual, is minimal. The algorithm is based on the Levenberg-Marquardt method.

        Notes on notation:
            - "residual vector" -> Residuals in vector notation "r = (r_1, ..., r_n)"
            - "Jacobian" -> Jacobian of the residual vector "J = (dr/dp_1, ..., dr/dp_m)", where "m"
                            is the number of parameters
            - "step" or "correction" -> Correction to the parameters "dp = (dp_1, ..., dp_m)"
            - "squared loss" -> "L = 1/2 ||r - J dp||^2"
            - "Hessian" -> Hessian of the squared loss "H = J^T J"
            - "gradient" -> (negative) gradient "g" of the loss

        Input parameters:

            ModelFuncResidue(Data_, Params, Jacobian)   -  Model function returns the i'th residue, where:
                                                            - "Data_" contains the i'th data point from "Data"
                                                              and its corresponding variables x1, x2, ...
                                                            - "Params" are the parameters
                                                            - "Jacobian" is used to store the i'th row of the
                                                              Jacobian matrix if the input Jacobian != NULL

            Data                                        -  Contains the data, where the i'th data point
                                                           in the array is followed by its corresponding
                                                           variables x1, x2, ...

            DataSize                                    -  Number of data points

            VariableSize                                -  Number of variables

            Params                                      -  Initial guess for the parameters, will be used to
                                                           store the corrected parameters

            ParamsSize                                  -  Number of parameters
    **/

        // Declare some storage variables

    double Lambda[1];                                   // Damping parameter
    double Nu;                                          // Parameter to increase Lambda if the step "dp" was unsuccessful
    double LambdaDown;                                  // Parameter to decrease Lambda if the step "dp" was successful

    double MinimalLoss;                                 // Minimal loss "S_min" yet encountered

    double StoppingGradient;                            // Stopping condition for the norm for the gradient "||g||"
    double StoppingParams;                              // Stopping condition for the relative change of the parameters "||dp|| / ||p||"
    double StoppingLoss;                                // Stopping condition for the change of the loss "S_new - S_min"

    int MaxIterations;

    struct PARAMSCORRECTION Args;                       // Output parameters of the function "ParamsCorrection"; for more see this function

    double Hessian[ParamsSize][ParamsSize];             // Hessian H of the quadratic loss
    double Gradient[ParamsSize];                        // (negative) Gradient "g" of the loss

    int DerivativeFlag;                                 // Decides whether "ModelFunc" should calculate the Jacobian or not


        // Initialise some values

    Lambda[0] = 0;                                      // Will be used as a flag for "ParamsCorrection" to calculate the actual initial value
    Nu = 2.0;
    LambdaDown = 1.0/3.0;

    // If the algorithm stops too early (first few iterations), the following values might be too large. Decrease them if needed.
    StoppingGradient = 1e-8 * DataSize * ParamsSize;
    StoppingParams = 1e-8;
    StoppingLoss = 1e-8;

    // For more iterations increase the value
    MaxIterations = 500;


        // Calculate the initial loss

    MinimalLoss = 0;
    size_t I = 0;

    for (size_t i = 0; i < DataSize; i++){

        // Store the i'th data point with its variables temporarely
        double TempData[VariableSize + 1];

        for (int j = 0; j < VariableSize + 1; j++)
            TempData[j] = Data[I+j];

        I += VariableSize + 1;

        // Caclulate the loss
        MinimalLoss += 0.5*pow(ModelFuncResidue(TempData, Params, NULL), 2);       // Set Jacobian = NULL to not calculate any derivatives
    }


        // Print initial informations to screen

    printf("  Iteration     Minimal Loss     Loss Change       Lambda       ||p|| / ||dp||        ||g||     \n");
    printf(" -----------   --------------   -------------   ------------   ----------------   --------------\n");
    printf("%8d %18.3e        -------         ------           ------           -------    \n\n", 0, MinimalLoss);


        // Start the algorithm

    // Begin by calculating the initial Jacobian, as it has not yet been calculated
    DerivativeFlag = 1;

    for (int n = 1; n <= MaxIterations; n++){

        // Calculate the new parameters
        Args = ParamsCorrection(ModelFuncResidue, Data, DataSize, VariableSize, Params, ParamsSize, Hessian, Gradient, Lambda, MinimalLoss, DerivativeFlag);

        // Determine whether the calculated step was successful
        if (Args.Gain > 0){

           // Update Lambda depending on the gain of the step
           double LambdaDown_ = fabs(1-pow(2*Args.Gain-1, 3));

            if (LambdaDown > LambdaDown_ && LambdaDown_ != 0)
                Lambda[0] *= LambdaDown_;
            else
                Lambda[0] *= LambdaDown;

            // Reset Nu to the initial value
            Nu = 2.0;

            // Store new parameters in Params
            for (int i = 0; i < ParamsSize; i++)
                Params[i] = Args.NewParams[i];

            // Update the minimal loss
            MinimalLoss = Args.NewLoss;

            // Derivatives must be calculated in the next iteration
            DerivativeFlag = 1;
        }

        else {

            // Update Lambda
            Lambda[0] *= Nu;

            // Increase Nu
            Nu *= 2;

            // Changes are not accepted; set parameters to be printed to screen to zero
            Args.LossChange = 0;
            Args.RelativeParamsChange = 0;

            // The derivatives do not have to be calculated again
            DerivativeFlag = 0;
        }

        // Print informations of the current iteration to screen
        printf("%8d %18.3e %15.2e %14.1e %17.3e %16.2e\n", n, MinimalLoss, Args.LossChange, Lambda[0], Args.RelativeParamsChange, Args.GradientNorm);


        // Check if the stopping conditions are met. If they are, break the loop and print a short text.
        if (Args.GradientNorm < StoppingGradient){
            printf("\n\nConvergence with ||g|| = %e < %e\n\nFinal loss: %.6e\n\n", Args.GradientNorm, StoppingGradient, MinimalLoss);
            break;
        }

        if (Args.RelativeParamsChange < StoppingParams && Args.RelativeParamsChange != 0.0){            // Exclude unsuccessful iterations, i.e. where the change of the loss was zero
            printf("\n\nConvergence with ||p|| / ||dp|| = %e < %e\n\nFinal loss: %.6e\n\n", Args.RelativeParamsChange, StoppingParams, MinimalLoss);
            break;
        }

        if (Args.LossChange < StoppingLoss && Args.LossChange != 0.0){                                  // Exclude unsuccessful iterations, i.e. where the relative change of the parameters was zero
            printf("\n\nConvergence with Loss Change = %e < %e\n\nFinal loss: %.6e\n\n", Args.LossChange, StoppingLoss, MinimalLoss);
            break;
        }

        if (Lambda[0] > pow(10, 50)){
            printf("\n\nConvergence with Lambda = %e > %.e\n\nFinal loss: %.6e\n\n", Lambda[0], pow(10, 50), MinimalLoss);
            break;
        }

        if (n == MaxIterations)
            printf("\n\nMaximal number of iterations reached: n_max = %d\n\nFinal loss: %.6e\n\n", MaxIterations, MinimalLoss);
    }

    // Print the optimal parameters
    for (int i = 0; i < ParamsSize; i++)
        printf("P%d = %.10e\n", i+1, Params[i]);
}


struct PARAMSCORRECTION ParamsCorrection(double (*ModelFuncResidue)(double*, double*, double*), double *Data, size_t DataSize, int VariableSize, double *Params, int ParamsSize, double Hessian[][ParamsSize], double Gradient[], double Lambda[], double MinimalLoss, int DerivativeFlag){

    /**
        This function computes the updated parameters for given damping parameter "Lambda", based
        on the Levenberg-Marquardt method. In the process of doing so, it also computes the Hessian of
        the quadratic loss "L(dp) = 1/2 ||r + J dp||^2", where "r = (r_1, ..., r_n)^T" is the residual
        vector, "J = (dr/dp_1, ..., dr/dp_m)" is the Jacobian of the residual vector and
        "dp = (dp_1, ..., dp_m)" is the correction to the parameters.
        and the gradient of the loss, if needed.

        Inpute parameters:

            ModelFuncResidue(Data_, Params, Jacobian)   -  Model function returns the i'th residue, where:
                                                            - "Data_" contains the i'th data point from "Data"
                                                              and its corresponding variables x1, x2, ...
                                                            - "Params" are the parameters
                                                            - "Jacobian" is used to store the i'th row of the
                                                              Jacobian matrix if the initial Jacobian != NULL

            Data                                        -  Contains the data, where the i'th data point
                                                           in the array is followed by its corresponding
                                                           variables x1, x2, ...

            DataSize                                    -  Number of data points

            VariableSize                                -  Number of variables

            Params                                      -  Initial guess for the parameters, will be used to
                                                           store the corrected parameters

            ParamsSize                                  -  Number of parameters

            Hessian                                     -  Hessian matrix of the quadratic loss

            Gradient                                    -  Gradient of the loss

            Lambda                                      -  Damping parameter

            MinimalLoss                                 -  Smallest loss yet encountered

            DerivativeFlag                              -  Decides whether the Jacobian must be calculated or not


        Output parameters:

            NewParams                                   -  Updated parameters "p_new = p_old + dp"

            NewLoss                                     -  Loss for the updated parameters

            LossChange                                  -  Difference between "MinimalLoss" and "NewLoss"

            Gain                                        -  Gain of the computed corrections "dp": "Gain = (S(p_old) - S(p_new)) / (L(0) - L(dp))"

            GradientNorm                                -  Norm of the gradient

            RelativeParamsChange                        -  Relative change of the updated parameters "||dp|| / ||p_new||"
    **/


    // Calculate the Hessian and the gradient if needed
    if (DerivativeFlag){

        // Initialise the Hessian and the gradient as 0
        for (int i = 0; i < ParamsSize; i++){
            Gradient[i] = 0;
            for (int j = i; j < ParamsSize; j++)
                Hessian[i][j] = 0;                  // It is sufficient to initialise only the upper triangular elements of the Hessian, since the Hessian is symmetric
        }

        size_t I = 0;

        // Loop through the "Data" array
        for (size_t i = 0; i < DataSize; i++){

            double Jacobian[ParamsSize];                // Storage for the i'th row of the Jacobian
            double TempData[VariableSize + 1];          // Storage for the i'th data point, and its corresponding variables

            for (int j = 0; j < VariableSize + 1; j++)
                TempData[j] = Data[I+j];

            double Residue = ModelFuncResidue(TempData, Params, Jacobian);         // Calculate the i'th residue, as well as the i'th row of the Jacobian


            // Calculate the gradient and the Hessian
            for (int j = 0; j < ParamsSize; j++){
                Gradient[j] -= Residue*Jacobian[j];

                for (int k = j; k < ParamsSize; k++){
                    Hessian[j][k] += Jacobian[j]*Jacobian[k];
                    Hessian[k][j] = Hessian[j][k];                          // Lower triangular matrix is equal to the upper triangular matrix
                }
            }

            I += VariableSize + 1;                  // Increase I to the index of the next data point in the Data array
        }
    }


    // If Lambda has not yet been initialised(i.e. Lambda = 0), choose Lambda proportional to the largest square root of the diagonal elements of the Hessian
    if (Lambda[0] == 0){
        Lambda[0] = pow(Hessian[0][0], 0.5);

        for (int i = 1; i < ParamsSize; i++){

            if (pow(Hessian[i][i], 0.5) > Lambda[0])
                Lambda[0] = pow(Hessian[i][i], 0.5);
        }

        Lambda[0] *= 0.001;             // Choose proportionality constant as 0.001
    }


            /** Calculate the correction to the parameters based on the Levenberg-Marquardt method **/

        // Declare new storage variables

    struct PARAMSCORRECTION Args;                                       // Output parameters of this function

    double M[ParamsSize][ParamsSize];                                   // where "M = Hessian + Lambda*Id"
    double TempData[VariableSize + 1];                                  // TempData will be used to calculate the loss of the step: LastLoss
    double *ParamsCorrection;                                           // Pointer to the corrections of the parameters
    double TempGain;                                                    // Temporary variables to calculate the denominator of the gain


    // Calculate M = H + Lambda*Id

    for (int i = 0; i < ParamsSize; i++){

        for (int j = 0; j < ParamsSize; j++){
            if (i == j)
                M[i][j] = Hessian[i][j] + Lambda[0];
            else
                M[i][j] = Hessian[i][j];
        }
    }


    /* ---------------------------------------------------------------------------------------------- */
    // Print the matrix M and the gradient to screen if necessary

    /*
    printf("\n\n");

    for (int i = 0; i < ParamsSize; i++){
        for (int j = 0; j < ParamsSize; j++)
            printf("%.1e ", M[i][j]);

        printf(" \t\t ||  %.1e\n", Gradient[i]);
    }

    printf("\n");
    */
    /* ---------------------------------------------------------------------------------------------- */


    // Calculate the corrections to the parameters
    ParamsCorrection = QRDecomposition(ParamsSize, M, Gradient);


    // Calculate the new parameters, the gain and some of the printed variables
    Args.Gain = 0;
    Args.GradientNorm = 0;
    Args.RelativeParamsChange = 0;
    Args.NewParams = malloc(sizeof(double) * ParamsSize);
    TempGain = 0;

    for (int i = 0; i < ParamsSize; i++){

        Args.NewParams[i] = Params[i] + ParamsCorrection[i];

        Args.Gain += Lambda[0] * pow(ParamsCorrection[i], 2);
        Args.GradientNorm += pow(Gradient[i], 2);
        Args.RelativeParamsChange += pow(Params[i], 2);

        double TempGain_ = 0;

        for (int j = 0; j < ParamsSize; j++)
            TempGain_ += Hessian[i][j]*ParamsCorrection[j];             // Calculate the i'th element of H*dp

        TempGain += 0.5 * TempGain_ * ParamsCorrection[i];                   // Scalar product of dp^T*H*dp
    }

    Args.RelativeParamsChange = pow(Args.Gain/(Lambda[0] * Args.RelativeParamsChange), 0.5);
    Args.GradientNorm = pow(Args.GradientNorm, 0.5);
    Args.Gain += fabs(TempGain);

    size_t I = 0;
    Args.NewLoss = 0;

    for (int i = 0; i < DataSize; i++){

        for (int j = 0; j < VariableSize + 1; j++)
            TempData[j] = Data[I+j];

        I += VariableSize + 1;

        Args.NewLoss += 0.5*pow(ModelFuncResidue(TempData, Args.NewParams, NULL), 2);
    }

    Args.LossChange = MinimalLoss - Args.NewLoss;

    Args.Gain = Args.LossChange/Args.Gain;

    /* ---------------------------------------------------------------------------------------------- */
    // Print the upper triangular matrix R = Q^t * M and the gradient to screen if necessary
    /*
    printf("\n\n");
    for (int i = 0; i < ParamsSize; i++){
        for (int j = 0; j < ParamsSize; j++)
            printf("%.1e ", M[i][j]);

        printf(" \t\t ||  %.1e\n", Gradient[i]);
    }
    printf("\n");
    */
    /* ---------------------------------------------------------------------------------------------- */

    return Args;
}


double *QRDecomposition(int ParamsSize, double M[][ParamsSize], double Gradient[]){

    /**
        This function performs a QR decomposition for a positive definite matrix M, using Householder
        transformations. It then solves the set of linear equations R * dp = Q^T * g for a given
        vector g and outputs the solution dp.


        Input parameters:

            ParamsSize                          -  Number of parameters / Dimension of M

            M                                   -  M = H + Lambda*Id / positive definite matrix

            Gradient                            -  Gradient of the loss / some vector


        Output parameters:

            ParamsCorrection                    -  Solution to the equation M * dp = g
    **/

    double v[ParamsSize];
    double *ParamsCorrection = malloc(sizeof(double) * ParamsSize);
    double Temp1, Temp2;

    //double Temp3[ParamsSize], Q[ParamsSize][ParamsSize];      // Declare variables to calculate the orthogonal matrix Q

    for (int i = 0; i < ParamsSize; i++){
        ParamsCorrection[i] = Gradient[i];

        /*
        Q[i][i] = 1;
        for (int j = i+1; j < ParamsSize; j++){
            Q[i][j] = 0;
            Q[j][i] = 0;
        }
        */
    }

    for (int i = 0; i < ParamsSize-1; i++){
        Temp1 = 0;

        for (int j = i; j < ParamsSize; j++){
            Temp1 += pow(M[j][i], 2);
            v[j] = M[j][i];
        }

        if (M[i][i] > 0)
            Temp1 = -pow(Temp1, 0.5);
        else
            Temp1 = pow(Temp1, 0.5);

        v[i] -= Temp1;

        Temp1 = pow(Temp1 * (Temp1 - M[i][i]), 0.5);

        for (int j = i; j < ParamsSize; j++)
            v[j] /= Temp1;

        Temp2 = 0;

        /*
        for (int j = 0; j < ParamsSize; j++)
            Temp3[j] = 0;
        */

        for (int j = i; j < ParamsSize; j++){
            Temp1 = 0;

            for (int k = i; k < ParamsSize; k++)
                Temp1 += M[k][j]*v[k];

            for (int k = i; k < ParamsSize; k++)
                M[k][j] -= v[k]*Temp1;

            Temp2 += ParamsCorrection[j]*v[j];

            /* //

            for (int k = 0; k < ParamsSize; k++)
                Temp3[k] += Q[j][k]*v[j];
            */
        }

        for (int j = i; j < ParamsSize; j++){
            ParamsCorrection[j] -= v[j]*Temp2;

            /* //

            for (int k = 0; k < ParamsSize; k++)
                Q[j][k] -= v[j]*Temp3[k];
            */
        }
    }

    /* ---------------------------------------------------------------------------------------------- */
    // Print the matrix product Q^T*Q = Id if necessary. Uncomment the prior lines to do so.
    /*
    printf("\n\n");
    for (int i = 0; i < ParamsSize; i++){
        for (int j = 0; j < ParamsSize; j++){
            Temp1 = 0;

            for (int k = 0; k < ParamsSize; k++)
                Temp1 += Q[i][k]*Q[j][k];

            printf("%.1e ", Temp1);

            if (j == ParamsSize-1)
                printf("\n");
        }
    }
    printf("\n");
    */
    /* ---------------------------------------------------------------------------------------------- */


    // Solve R*dp = Q^T*g
    for (int i = ParamsSize-1; i > -1; i--){

        for (int j = i+1; j < ParamsSize; j++)
            ParamsCorrection[i] -= ParamsCorrection[j]*M[i][j];

        ParamsCorrection[i] /= M[i][i];
    }

    return ParamsCorrection;
}
