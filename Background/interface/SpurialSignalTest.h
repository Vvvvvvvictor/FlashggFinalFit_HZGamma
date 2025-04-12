#ifndef SPURIALSIGNALTEST_H
#define SPURIALSIGNALTEST_H

#include "TString.h"

/**
 * Spurious signal test function that uses template signal models 
 * @param cat           - Category name (e.g. "VBF0")
 * @param funcName      - Name of function to test (e.g. "bern3", "exp2", etc.)
 * @param sig           - Signal strength multiplier
 * @param runPeriod     - Run period ("run2" or "run3")
 * @return true if test passes, false otherwise
 */
bool SpurialSignalTest(TString cat, TString funcType, int order, int sig, TString runPeriod, TString outDir);

#endif // SPURIALSIGNALTEST_H
