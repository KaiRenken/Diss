#include "../headers/basics.h"
#include "../headers/Matrix.h"
#include "../headers/findPartitions.h"
#include <map>

void findPartitions(int pN, int pA, int pB, Matrix* pMat, int pIter, int pFirstLine)
{
	// Set entries for all numbers from 1 to pN, which are between pA and pB and call next iteration step
    if (pA <= pN)
    {
        for (int i = pA; i <= pN; i++)
        {
            pMat->setEntry(i - pA + pFirstLine, pIter, i);

            for (int j = pIter + 1; j < pMat->getColumns(); j++)
            {
                pMat->setEntry(i - pA + pFirstLine, j, 0);
            }
        }

        for (int i = pN + 1; i <= pB; i++)
        {
            pMat->setEntry(i - pA + pFirstLine, pIter, 0);
        }

        if (pA != 1)
        {
            int nextLine = pN - pA + 1 + pFirstLine;
            findPartitions(pA - 1, pN + 1, pB, pMat, pIter + 1, nextLine);
        }

        return;
    }

	// Define some variables needed for the computation
    int s = pB - pA + 1;
    int c = (2 * pN) - (2 * s) + 1;
    int m;
    int steps = c - pA;
    int posX;
    int posY;
    int temp;
    int nextLine;

    if (steps < 0)
    {
		steps = 0;
	}

    std::map<int,int> result;

    for (int i = pN; i > pN - s; i--)
    {
        result[i] = c - i;
    }

    int i = pA;
    m = c - i;

	// Find pairs for negative summing and swap second components
    while (m > 0)
    {
        posX = findPairsByStep(m, steps, pN - s + 1, 0);
        posY = findPairsByStep(m, steps, pN - s + 1, 1);
        temp = result[posX];
        result[posX] = result[posY];
        result[posY] = temp;
        pMat->setEntry(i - pA + pFirstLine, pIter, posX);
        pMat->setEntry(i - pA + pFirstLine, pIter + 1, result[posX]);
        pMat->setEntry(i - pA + (2 * m) + pFirstLine, pIter, posY);
        pMat->setEntry(i - pA + (2 * m) + pFirstLine, pIter + 1, result[posY]);
        i++;
        m = c - i;
    }

	// Find pair for zero summing
    if (steps > 0)
    {
        posX = findFirstFreeSpace(steps, pN - s + 1);
        pMat->setEntry(i - pA + pFirstLine, pIter, posX);
        pMat->setEntry(i - pA + pFirstLine, pIter + 1, result[posX]);
    }

    posX = findSecondFreeSpace(steps, pN - s + 1);

	// Find pairs for non-zero summing
    if (steps == 0)
    {
        for (i = 0; i < s; i++)
        {
            pMat->setEntry(i + pFirstLine, pIter, posX);
            pMat->setEntry(i + pFirstLine, pIter + 1, result[posX]);
            posX++;
        }

    } else {
        for (i = (2 * steps) + 1; i < s; i++)
        {
            pMat->setEntry(i + pFirstLine, pIter, posX);
            pMat->setEntry(i + pFirstLine, pIter + 1, result[posX]);
            posX++;
        }
    }

	// Call next iteration step
    if (pN - (2 * s) > 0)
    {
        for (int i = pFirstLine; i < s + pFirstLine; i++)
        {
            posX = pMat->getEntry(i, pIter);
            posY = pMat->getEntry(i, pIter + 1);

            if (posX + posY != pA + i - pFirstLine)
            {
                nextLine = i;
                break;
            }
        }

        for (int i = pFirstLine; i < nextLine; i++)
        {
            for (int j = pIter + 2; j < pMat->getColumns(); j++)
            {
                pMat->setEntry(i, j, 0);
            }
        }

        if (c - pA >= 0)
        {
            findPartitions(pN - (2 * s), pA + (2 * steps) + 1 - c, pB - c, pMat, pIter + 2, nextLine);
        } else {
            findPartitions(pN - (2 * s), pA - c, pB - c, pMat, pIter + 2, nextLine);
        }

    // Set remaining entries to zero
    } else {
        for (int i = pIter + 2; i < pMat->getColumns(); i++)
        {
            pMat->setEntry(pMat->getLines() - 1, i, 0);
        }
    }
}
