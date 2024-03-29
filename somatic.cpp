/*
 * somatic.cpp for MSIsensor1.1 
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

#include <iostream>
#include <sstream>
#include <bitset>
#include <omp.h>
#include <cmath>
#include "somatic.h"

SomaticSite::SomaticSite()
    : chr("")
    , length(0)
    , location(0)
    , bases("")
    , fbases("")
    , ebases("")
    , diff( 0.0 )
    , pValue( 1.0 )
    , somatic( false )
    , FDR( 1.0 )
    , rank( 1 )
{
    //xxxxxxxx
};


SomaticSite::~SomaticSite() {
    // xxxxx
};

// PourOut values
void SomaticSite::PourOut() {
    std::cerr << chr << "\t"
              << location << "\t"
              << bases << "\t"
              << length << "\t"
              << fbases << "\t"
              << bases << "\t"
              << ebases << "\t"
              << diff << "\t"
              << pValue << "\n";

};

