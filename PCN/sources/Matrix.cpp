#include "../headers/Matrix.h"
#include "../headers/basics.h"

Matrix::Matrix(int pLines, int pColumns)
{
    lines = pLines;
    columns = pColumns;

    matrix.resize(lines);

    for (int i = 0; i < lines; i++)
    {
        matrix[i].resize(columns);
    }
}

int Matrix::getLines()
{
    return lines;
}

int Matrix::getColumns()
{
    return columns;
}

int Matrix::getEntry(int pLine, int pColumn)
{
    return matrix[pLine][pColumn];
}

void Matrix::setEntry(int pLine, int pColumn, int pEntry)
{
    matrix[pLine][pColumn] = pEntry;
}
