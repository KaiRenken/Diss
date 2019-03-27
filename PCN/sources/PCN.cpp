#include "../headers/Matrix.h"
#include "../headers/findPartitions.h"
#include "../headers/basics.h"
#include <iostream>
#include <math.h>
#include <cstdlib>

int main(int argn, char *argv[])
{
	// Check number of input arguments
    if (argn != 4)
    {
		displayHowTo();
		return 0;
	}

	// Convert input arguments to integers
    int pN = strtol(argv[1], NULL, 0);
    int pA = strtol(argv[2], NULL, 0);
    int pB = strtol(argv[3], NULL, 0);

	// Check if conversion succeeded
    if (pN == 0 || pA == 0 || pB == 0)
    {
		displayHowTo();
		return 0;
	}

	// Check if input arguments satisfy the assumptions of the statement
    double n = (double) pN;
    double a = (double) pA;
    double b = (double) pB;

    if ((n * (n + 1)) / 2 != ((a + b) / 2) * (b - a + 1))
    {
		displayHowTo();
		return 0;
	}

	// Create the output matrix
    Matrix* mat = new Matrix(b - a + 1, n);

	// Compute the partitioning
    findPartitions(n, a, b, mat, 0, 0);

	// Display the output
    for (int i = 0; i < b - a + 1; i++)
    {
        for (int j = 0; j < n; j++)
        {
           std::cout << mat->getEntry(i,j) << " ";
        }
        std::cout << std::endl;
    }

    delete mat;

    return 0;
}
