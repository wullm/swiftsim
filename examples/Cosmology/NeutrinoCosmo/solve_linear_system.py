# Solve the linear system for the control variate optimal coefficients
# By Willem Elbers
# 23 June, 2020

import numpy as np;

X = np.genfromtxt("control_vars_1.txt", names = True);

N = 3;
M = int(X["Step"][-1]);

for i in np.arange(M):
    A = np.zeros(N*N).reshape(N, N);
    b = np.zeros(N);

    # Retrieve the values form the dataset
    A[0,0] = X["Covar_11"][i];
    A[0,1] = X["Covar_12"][i];
    A[0,2] = X["Covar_13"][i];
    A[1,1] = X["Covar_22"][i];
    A[1,2] = X["Covar_23"][i];
    A[2,2] = X["Covar_33"][i];

    b[0] = X["Covar_01"][i];
    b[1] = X["Covar_02"][i];
    b[2] = X["Covar_03"][i];

    Time = X["Time"][i];
    Scale = X["Scalefactor"][i];


    # Make matrix symmetric
    A = A + np.tril(A.T, k = -1);

    if (np.linalg.det(A) != 0):

        # Solve the system
        c = np.linalg.inv(A).dot(b);

        print(Time, Scale, c[0], c[1], c[2]);
