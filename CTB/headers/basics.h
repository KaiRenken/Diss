#ifndef BASICS_H
#define BASICS_H

// Returns an integer array whose entries are the digits of the binary representation of a given integer
int* intToBinary(int n);

// Returns the binomial coefficient n choose k
int nChooseK(int n, int k);

// Returns a matrix with n lines and nChooseK(n,k) columns whose columns represent all possibilities how one can arrange k ones on n places (the other entries are zeros)
Matrix* binaryCombinations(int n, int k);

// Returns a matrix with k lines and nChooseK(v->getLines(), k) columns whose columns represent all possible combinations of k entries from v, where v is a column matrix
Matrix* vChooseK(Matrix* v, int k);

// Returns the product of matrix mat1 with matrix mat2
Matrix* multMats(Matrix* mat1, Matrix* mat2);

// Returns the matrix representing the boundary operator from the (k+1)-th chain group to the k-th chain group of the n-simplex.
Matrix* boundaryMatrix(int n, int k);

// Returns the matrix representing the coboundary operator from the k-th cochain group to the (k+1)-th cochain group of the n-simplex.
Matrix* coboundaryMatrix(int n, int k);

// Returns the transposed of the given matrix M
Matrix* transposeMatrix(Matrix* M);

// Returns the sum of matrix mat1 with matrix mat2
Matrix* addMats(Matrix* mat1, Matrix* mat2);

// Returns the k-th Cheeger constant of the n-simplex.
float cheegerConstant(int n, int k);

#endif
