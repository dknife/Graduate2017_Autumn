import numpy as np

def convertLongVector(x) :
    longX = np.zeros(shape=(len(x)*3,))
    for i in range(0, len(x)) :
        for j in range(0,3) :
            longX[i*3+j] = x[i,j]
    return longX

def set33SubMatrix(M, i, j, m) :
    for row in range(0,3) :
        for col in range(0,3) :
            M[i*3+row, j*3+col] = m[row,col]

def set33SubMatrixSymetric(M, i, j, m):
    for row in range(0, 3):
        for col in range(0, 3):
            M[i * 3 + row, j * 3 + col] = m[row, col]
            M[j * 3 + row, i * 3 + col] = m[row, col]