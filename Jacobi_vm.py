import numpy as np

A = [[4, -1, -1],               # Coefficients Matrix
     [-2, 6, 1],
     [-1, 1, 7]]
b = [3, 9, -6]                  # RHS of system of equations
x_old = [0, 0, 0]               # initial guess
x_new = np.zeros_like(x_old)    # initialize array to store new x values
D = np.zeros_like(A)            # diagonal matrix
L = np.zeros_like(A)            # lower triangular matrix
U = np.zeros_like(A)            # upper triangular matrix
for i in range(len(A)):         # in these loops we create D, L, and U
    for j in range(len(A)):
        if i == j:
            D[i][j] = A[i][j]
        elif i < j:
            U[i][j] = A[i][j]
        else:
            L[i][j] = A[i][j]
for k in range(50):             # now we solve fpr a specified number of iterations
    parenthesis = np.matmul((-L-U), x_old) + b
    invD = np.linalg.inv(D)
    x_new = np.matmul(invD, parenthesis)
    x_old = x_new
print(x_new)


