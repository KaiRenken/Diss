#pragma once
#include <vector>

typedef std::vector<int> line;

class Matrix
{

    public:
        Matrix(int pLines, int pColumns);

        // Get number of lines
        int getLines();

        // Get size of a line
        int getColumns();

        // Get a certain entry of the matrix
        int getEntry(int pLine, int pColumn);

        // Set a certain entry of the matrix
        void setEntry(int pLine, int pColumn, int pEntry);

    private:
        int lines;
        int columns;
        std::vector<line> matrix;
};
