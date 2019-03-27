#ifndef MATRIX_H
#define MATRIX_H

typedef std::vector<int> line;

class Matrix
{
    public:

        Matrix(int pLines, int pColumns);

        // Returns the number of lines of the matrix
        int getLines();

        // Returns the number of columns of the matrix
        int getColumns();

        // Returns a certain entry of the matrix
        int getEntry(int pLine, int pColumn);

        // Sets a certain entry of the matrix
        void setEntry(int pLine, int pColumn, int pEntry);

        // Returns a certain column of the matrix as a column matrix
        Matrix* getColumn(int j);

        // Checks if a given column matrix is contained in the column matrix as a subset
        bool containsColumn(Matrix* column);

        // Returns the norm of the column matrix
        int getNorm();

        // Returns the cosystolic norm of the column matrix representing a k-dimensional cochain of the n-simplex (requires the coboundaryMatrix(n,k-1) as coboundaryMat)
        int getCosystolicNorm(int n, int k, Matrix* coboundaryMat);

        // Returns the coboundary expansion of the column matrix representing a k-dimensional cochain of the n-simplex (requires the coboundaryMatrix(n,k) as coboundaryMat and the coboundaryMatrix(n,k-1) as coboundaryMat1)
        float getCoboundaryExpansion(int n, int k, Matrix* coboundaryMat, Matrix* coboundaryMat1);

        // Checks if the column matrix is a k-dimensional cosystole of the n-simplex (requires the coboundaryMatrix(n,k-1) as coboundaryMat)
        bool isCosystole(int n, int k, Matrix* coboundaryMat);

    private:
        int lines;
        int columns;
        std::vector<line> matrix;
};

#endif
