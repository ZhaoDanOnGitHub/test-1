
/*
 * sample.h for MSIsensor1.1
 * Copyright (c) 2013 Beifang Niu && Kai Ye CNIC && Xi’an Jiaotong University All Rights Reserved.
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#ifndef _SAMPLE_H_
#define _SAMPLE_H_

#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>

#include "somatic.h"
#include "structs.h"

// sample
class Sample {
public:
    Sample();
    ~Sample();

    std::string outputPrefix;

    std::ofstream output;
    //std::ofstream outputPSomatic;
    std::ofstream outputSomatic;
    std::ofstream outputGermline;
    std::ofstream outputDistribution;
    std::ofstream outputFeature;

    unsigned numberOfSites;

    unsigned precisionNumS;
    unsigned precisionNumL; 

    unsigned numberOfDataPoints; 
    unsigned numberOfMsiDataPoints;
    unsigned numberOftotalSites;

    // container for FDR
    std::vector< SomaticSite > totalSomaticSites;

    void iniOutput( const std::string &gavePrefix );
    void iniTumorDisOutput( const std::string &gavePrefix );
    void pourOutMsiScore();
    void closeOutStream();
    void calculateFDR();
    void pourOutSomaticFDR();
    void VerboseInfo();
    protected:
        // xxx
};

#endif //_SAMPLE_H_

