#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <vector>
#include "../headers/Matrix.h"
#include "../headers/basics.h"

int* intToBinary(int n)
{
	int digits = floor(log2(n)) + 1;
	int* result = new int[digits];
	int temp = n;

	for (int i = digits - 1; i >= 0; i--)
    {
		result[i] = temp % 2;
		temp = floor(temp / 2);
	}

	return result;
}

int nChooseK(int n, int k)
{

	if (k == 0)
	{
		return 1;
	}

	if (2*k > n)
	{
		k = n-k;
	}

	int result = 1;

	for (int i = 1; i <= k; i++)
	{
		result = result * (n - k + i) / i;
	}

	return result;
}

Matrix* binaryCombinations(int n, int k)
{
    if (n == k)
    {
        Matrix* result = new Matrix(n,1);

        for (int i = 0; i < n; i++)
        {
            result->setEntry(i, 0, 1);
        }
        return result;
    }

    if (k == 1)
    {
        Matrix* result = new Matrix(n,n);

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                if (i == j)
                {
                    result->setEntry(i, j, 1);
                } else {
                    result->setEntry(i, j, 0);
                }
            }
        }
        return result;
    }

    Matrix* result = new Matrix(n, nChooseK(n, k));
    Matrix* temp1 = binaryCombinations(n - 1, k - 1);
    Matrix* temp2 = binaryCombinations(n - 1, k);

    for (int c = 0; c < nChooseK(n - 1, k - 1); c++)
    {
        result->setEntry(0, c, 1);

        for(int i = 1; i <= n-1; i++)
        {
            for (int j = 0; j < nChooseK(n - 1, k - 1); j++)
            {
                result->setEntry(i, j, temp1->getEntry(i-1, j));
            }
        }
    }

    for (int c = nChooseK(n - 1, k - 1); c < nChooseK(n, k); c++)
    {
        result->setEntry(0, c, 0);

        for(int i = 1; i <= n-1; i++)
        {
            for (int j = nChooseK(n - 1, k - 1); j < nChooseK(n, k); j++)
            {
                result->setEntry(i, j, temp2->getEntry(i-1, j - nChooseK(n - 1, k - 1)));
            }
        }
    }

    delete temp1;
    delete temp2;

    return result;
}

Matrix* vChooseK(Matrix* v, int k)
{
    Matrix* result = new Matrix(k, nChooseK(v->getLines(), k));
    Matrix* binCombs = binaryCombinations(v->getLines(), k);
    int counter;

    for (int i = 0; i < binCombs->getColumns(); i++)
    {
        counter = 0;

        for (int j = 0; j < binCombs->getLines(); j++)
        {
            if (binCombs->getEntry(j, i) == 1)
            {
                result->setEntry(counter, i, v->getEntry(j,0));
                counter++;
            }
        }
    }

    delete binCombs;
    return result;
}

Matrix* multMats(Matrix* mat1, Matrix* mat2)
{
	int temp;
    int prod;

	Matrix* result = new Matrix(mat1->getLines(), mat2->getColumns());

	for (int i = 0; i < mat1->getLines(); i++)
	{
		for (int j = 0; j < mat2->getColumns(); j++)
		{
			temp = 0;

			for (int n = 0; n < mat1->getColumns(); n++)
			{
			    if (mat1->getEntry(i,n) == 1 && mat2->getEntry(n,j) == 1)
                {
                    if (temp == 1)
                    {
                        temp = 0;
                    } else {
                        temp = 1;
                    }
                }
			}
			result->setEntry(i,j,temp);
		}
	}

	return result;
}

Matrix* transposeMatrix(Matrix* M)
{
    Matrix* result = new Matrix(M->getColumns(), M->getLines());

    for (int i = 0; i < M->getLines(); i++)
    {
        for (int j = 0; j < M->getColumns(); j++)
        {
            result->setEntry(j, i, M->getEntry(i,j));
        }
    }

    return result;
}

Matrix* boundaryMatrix(int n, int k)
{
	Matrix* result = new Matrix(nChooseK(n, k+1), nChooseK(n, k+2));
    Matrix* vertices = new Matrix(n, 1);
	Matrix* columnJ;
	Matrix* columnI;

	for (int i = 0; i < n; i++)
	{
		vertices->setEntry(i, 0, i+1);
	}

	Matrix* s1 = vChooseK(vertices, k+2);
	Matrix* s0 = vChooseK(vertices, k+1);

	for (int j = 0; j < result->getColumns(); j++)
	{
		columnJ = s1->getColumn(j);

		for (int i = 0; i < result->getLines(); i++)
		{
			columnI = s0->getColumn(i);

			if (columnJ->containsColumn(columnI) == true)
			{
				result->setEntry(i, j, 1);
			} else {
				result->setEntry(i, j, 0);
			}

			delete columnI;
		}

		delete columnJ;
	}

	delete vertices;
	delete s1;
	delete s0;

	return result;
}

Matrix* coboundaryMatrix(int n, int k)
{
	Matrix* result = transposeMatrix(boundaryMatrix(n, k));

	return result;
}

Matrix* addMats(Matrix* mat1, Matrix* mat2)
{
    Matrix* result = new Matrix(mat1->getLines(), mat1->getColumns());

    for (int i = 0; i < result->getLines(); i++)
    {
        for (int j = 0; j < result->getColumns(); j++)
        {
            if (mat1->getEntry(i, j) == mat2->getEntry(i, j))
            {
                result->setEntry(i, j, 0);
            } else {
                result->setEntry(i, j, 1);
            }
        }
    }

    return result;
}

float cheegerConstant(int n, int k)
{
	float result = ( n / (k + 2) ) + 1;
	float tempRes;
	Matrix* column;
	Matrix* cobound;
	Matrix* simplices;
	Matrix* coboundaryMat = coboundaryMatrix(n, k);
	Matrix* coboundaryMat1 = coboundaryMatrix(n, k-1);

	int length = nChooseK(n, k+1);

	for (int i = 1; i <= length; i++)
    {
        simplices = binaryCombinations(length, i);

        for (int j = 0; j < simplices->getColumns(); j++)
        {
            column = simplices->getColumn(j);

            if (column->isCosystole(n, k, coboundaryMat1) == true)
            {
                cobound = multMats(coboundaryMat, column);
                tempRes = (float)cobound->getNorm() / (float)column->getNorm();

                if (tempRes < result)
                {
                    result = tempRes;
                }
                delete cobound;
            }

            delete column;
        }

        delete simplices;
    }

    delete coboundaryMat;
    delete coboundaryMat1;

	return result;
}
