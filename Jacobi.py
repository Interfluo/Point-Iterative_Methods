
def jacobi(A, b, x_old, max_error, max_its):
    """
    Jacobi Method
    Assumes that the coefficient matrix, A, has no zeros along it's main diagonal
    Convergence conditions
        (1) spectral radius of iteration matrix < 1
        (2) matrix A is strictly or irreducibly diagonally dominant (abs(aii) > sum(abs(aij)))
    """
    x_new = x_old.copy()
    its = 0
    error = 1e8
    while its < max_its and error > max_error:
        for i in range(len(A)):
            sum_val = 0
            for j in range(len(A)):
                if j != i:
                    sum_val += (A[i][j] * x_old[j])
            x_new[i] = (1 / A[i][i]) * (b[i] - sum_val)
        error = 0
        for k in range(len(A)):
            error += abs(x_old[k]-x_new[k])
            x_old[k] = x_new[k]
        its += 1
    print("max_its to convergence:", its)
    print("error:", error)
    return x_new


max_its = 50
max_error = 1e-8

A = [[2, 1],
     [5, 7]]
b = [11, 13]
x_old = [0, 0]
print("================================================")
print("A:", A)
x_new = jacobi(A, b, x_old, max_error, max_its)
print("x =", x_new)

A = [[2, 1, 1],
     [-1, 3, -1],
     [1, -1, 2]]
b = [7, 2, 5]
x_old = [0, 0, 0]
print("================================================")
print("A:", A)
x_new = jacobi(A, b, x_old, max_error, max_its)
print("x =", x_new)

A = [[10, -1, 2, 0],
     [-1, 11, -1, 3],
     [2, -1, 10, -1],
     [0, 3, -1, 8]]
b = [6, 25, -11, 15]
x_old = [0, 0, 0, 0]
print("================================================")
print("A:", A)
x_new = jacobi(A, b, x_old, max_error, max_its)
print("x =", x_new)