#pragma once

// Find pairs of natural numbers with distance 'pStep' within the most narrow range
int findPairsByStep(int pStep, int pSteps, int pShift, bool pComp);

// Find first non-used number within the range of findPairsByStep
int findFirstFreeSpace(int pSteps, int pShift);

// Find second non-used number within the range of findPairsByStep
int findSecondFreeSpace(int pSteps, int pShift);

// Display how to use PCN
void displayHowTo();
