#include <vector>
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

Matrix* Matrix::getColumn(int j)
{
	Matrix* result = new Matrix(lines, 1);

	for (int i = 0; i < lines; i++)
	{
		result->setEntry(i, 0, matrix[i][j]);
	}

	return result;
}

bool Matrix::containsColumn(Matrix* column)
{
	int checker;

	for (int i = 0; i < column->getLines(); i++)
	{
		checker = 0;

		for (int l = 0; l < lines; l++)
		{
			if (matrix[l][0] == column->getEntry(i,0))
			{
				checker = 1;
			}
		}

		if (checker == 0)
		{
			return false;
		}
	}

	return true;
}

int Matrix::getNorm()
{
	int result = 0;

	for (int i = 0; i < lines; i++)
	{
		result = result + matrix[i][0];
	}

	return result;
}

int Matrix::getCosystolicNorm(int n, int k, Matrix* coboundaryMat)
{
	int length = nChooseK(n, k);
	int result = this->getNorm();
	int tempNorm;
	Matrix* simplices;
	Matrix* column;
	Matrix* cobound;
	Matrix* tempVec;

	for (int i = 1; i <= length; i++)
    {
        simplices = binaryCombinations(length, i);

        for (int j = 0; j < simplices->getColumns(); j++)
        {
            column = simplices->getColumn(j);
            cobound = multMats(coboundaryMat, column);
            tempVec = addMats(cobound, this);

            tempNorm = tempVec->getNorm();

            if (tempNorm < result)
            {
                result = tempNorm;
            }

            delete column;
            delete cobound;
            delete tempVec;
        }

        delete simplices;
    }

    return result;
}

float Matrix::getCoboundaryExpansion(int n, int k, Matrix* coboundaryMat, Matrix* coboundaryMat1)
{

	Matrix* coboundary = multMats(coboundaryMat, this);

	int cosNorm = this->getCosystolicNorm(n, k, coboundaryMat1);

	if (cosNorm == 0)
	{
		return 0;
	}

	float result = (float)coboundary->getNorm() / (float)cosNorm;

	return result;
}

bool Matrix::isCosystole(int n, int k, Matrix* coboundaryMat)
{
    int length = nChooseK(n, k);
	int norm = this->getNorm();
	int tempNorm;
	Matrix* simplices;
	Matrix* column;
	Matrix* cobound;
	Matrix* tempVec;

	for (int i = 1; i <= length; i++)
    {
        simplices = binaryCombinations(length, i);

        for (int j = 0; j < simplices->getColumns(); j++)
        {
            column = simplices->getColumn(j);
            cobound = multMats(coboundaryMat, column);
            tempVec = addMats(cobound, this);

            tempNorm = tempVec->getNorm();

            if (tempNorm < norm)
            {
                delete column;
                delete cobound;
                delete tempVec;
                delete simplices;
                return false;
            }

            delete column;
            delete cobound;
            delete tempVec;
        }

        delete simplices;
    }

    return true;
}
