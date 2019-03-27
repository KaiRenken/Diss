#include "../headers/basics.h"
#include <iostream>
#include <math.h>
#include <vector>

int findPairsByStep(int pStep, int pSteps, int pShift, bool pComp)
{
    if (pStep > pSteps)
    {
		return 0;
	}

    double step = (double) pStep;
    double steps = (double) pSteps;
    double result;

    if (pStep % 2 == 1)
    {
        result = ceil(steps / 2) - ceil(step / 2);
    } else {
        result = (ceil(steps / 2) * 2) + floor(steps / 2) - (step / 2);
    }

    if (pComp == 0)
    {
        return (int) result + pShift;
    } else {
        return (int) result + pStep + pShift;
    }
}

int findFirstFreeSpace(int pSteps, int pShift)
{
    double steps = (double) pSteps;
    double result;
    result = (2 * ceil(steps / 2)) + floor(steps / 2);
    return (int) result + pShift;
}

int findSecondFreeSpace(int pSteps, int pShift)
{
    if (pSteps == 1)
    {
        return 2 + pShift;
    } else if (pSteps == 0) {
        return pShift;
    } else {
        double steps = (double) pSteps;
        double result;
        result = (2 * ceil(steps / 2)) + (2 * floor(steps / 2)) + 1;
        return (int) result + pShift;
    }
}

void displayHowTo()
{
    std::cout << "PCN - Partitioning Consecutive Numbers (written by Kai Renken, 2018)" << std::endl;
    std::cout << std::endl;
    std::cout << "Type 'PCN [Parameter n] [Parameter a] [Parameter b]'," << std::endl;
    std::cout << "such that the sum of the numbers 1 to n equals the sum of the numbers a to b" << std::endl;
}
