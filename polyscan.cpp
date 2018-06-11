/*
 * polyscan.cpp for MSIsensor1.1 
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
#include <map>
#include <omp.h>

#include "utilities.h"
#include "polyscan.h"
#include "bamreader.h"
#include "param.h"
#include "sample.h"

extern Param paramd;
extern bit8_t alphabet[];
extern bit8_t rev_alphabet[];
extern char uhomo_code[];
extern char homo_code[];
extern Sample sample;

PolyScan::PolyScan() { 
    homosBuffer.reserve(paramd.bufSize);
    totalSites.reserve(paramd.bufSize);
}

PolyScan::~PolyScan() { 
    totalSites.clear();
}

// eliminates a character from the input string
void PolyScan::eliminate(const char ch, std::string & str){
    size_t eliminateCharPos = str.find(ch);
    while (eliminateCharPos != std::string::npos) {
        str.erase(eliminateCharPos, 1);
        eliminateCharPos = str.find(ch);
    }
}

// Parse one region 
bool PolyScan::ParseOneRegion(const std::string & regionString) {
    size_t separatorPos = regionString.find(":");
    bool correctParse   = false;
    bool m_endDefined   = false;
    bool m_startDefined = false;
    int m_start = -1;
    int m_end   = -1;
    std::string m_targetChromosomeName;
    // found a separator
    if (separatorPos != std::string::npos) {
        m_targetChromosomeName = regionString.substr(0, separatorPos);
        std::string coordinates = regionString.substr(separatorPos + 1);
        // removes the ',' in 1,000 or 1,000,000 that users may add 
        // for readability but wreak havoc with atoi
        eliminate(',', coordinates); 
        size_t startEndSeparatorPos = coordinates.find("-");
        // there are two coordinates
        if (startEndSeparatorPos != std::string::npos) {
            std::string secondPositionStr = coordinates.substr(startEndSeparatorPos + 1);
            m_end = atoi(secondPositionStr.c_str());
            m_endDefined = true;
        }

        m_start = atoi(coordinates.c_str());
        m_startDefined = true;
        if (m_start < 0 || (m_endDefined && (m_end < m_start))) {
            correctParse = false;
        } else { correctParse = true; }
    }
    // no separator found
    else {
        m_targetChromosomeName = regionString;
        correctParse = true;
    }
    // assign values
    region_one.chr   = m_targetChromosomeName;
    region_one.start = m_start;
    region_one.end   = m_end;

    return correctParse;
}

// Loading bed regions
void PolyScan::LoadBeds(std::ifstream &fin) {
    std::string chr;
    std::string line;
    std::string tempChr = "";
    int start;
    int stop;
    int i = -1;
    while (getline(fin, line)){
        std::stringstream linestream(line);
        linestream >> chr;
        linestream >> start;
        linestream >> stop;
        /*
        std::cout<< chr <<"\t"
                 << start <<"\t"
                 << stop << "\n";
        */
        if (chr == tempChr) {
            BedRegion tempBedRegion;
            tempBedRegion.start = start;
            tempBedRegion.end = stop;
            beds[i].regions_list.push_back(tempBedRegion);
        } else {
            ++i;
            BedChr tempBedChr;
            tempBedChr.chr = chr;
            beds.push_back(tempBedChr);
            // load mapping
            chrMaptoIndex.insert(std::pair<std::string, bit16_t>(chr, i));
            BedRegion tempBedRegion;
            tempBedRegion.start = start;
            tempBedRegion.end = stop;
            beds[i].regions_list.push_back(tempBedRegion);
            tempChr = chr;
        } 
        linestream.clear();
        linestream.str("");
    } 
}

// loading bam list
// only load one bam file
void PolyScan::LoadBams(const std::string &bam1, const std::string &bam2) {
    BamPairs t_bampair;
    t_bampair.sName = "sample_name";

    if (bam1.find(".bam") != std::string::npos) { 
        t_bampair.normal_bam = bam1;
    } else { std::cerr << "please provide valid format normal bam file ! \n"; exit(0); }
    if (bam2.find(".bam") != std::string::npos) {
        t_bampair.tumor_bam = bam2; 
    } else { std::cerr << "please provide valid format tumor bam file ! \n"; exit(0); }

    // loading
    totalBamPairs.push_back(t_bampair);
    totalBamPairsNum++;
}


void PolyScan::LoadBam(const std::string &bam) {
    BamTumors t_bamtumor;
    t_bamtumor.sName = "sample_name";

    if (bam.find(".bam") != std::string::npos) { 
        t_bamtumor.tumor_bam = bam;
    } else { std::cerr << "please provide valid format normal bam file ! \n"; exit(0); }

    // loading
    totalBamTumors.push_back(t_bamtumor);
    totalBamTumorsNum++;
}
// read and load sites
void PolyScan::LoadHomosAndMicrosates(std::ifstream &fin) {
    std::string chr;
    std::string bases;
    std::string fbases;
    std::string ebases;
    std::string line;
    std::string tChr = "";
    // count total loading sites
    //
    totalHomosites = 0;
    int loc;
    bit8_t  siteLength;
    bit16_t tsiteLength;
    bit16_t siteBinary;
    bit16_t siteRepeats;
    bit16_t frontF;
    bit16_t tailF;

    int j = 0;
    BedChr tbedChr;
    BedRegion tbedRegion;
    bit16_t tIndex;
    
    // skip title
    getline(fin, line);
    while (getline(fin, line)) {
        std::stringstream linestream(line);
        linestream >> chr;
        linestream >> loc;
        linestream >> tsiteLength;
        linestream >> siteBinary;
        linestream >> siteRepeats;
        linestream >> frontF;
        linestream >> tailF;
        //xxx
        linestream >> bases;
        linestream >> fbases;
        linestream >> ebases;

        // filtering
        if (tsiteLength > 1 && paramd.HomoOnly == 1) continue;
        if (tsiteLength == 1 && paramd.MicrosateOnly == 1) continue;
        if (tsiteLength == 1 && ((siteRepeats < paramd.MininalHomoForDis) || (siteRepeats > paramd.MaxHomoSize)) ) continue;
        if (tsiteLength > 1 && ((siteRepeats < paramd.MinMicrosateForDis) || (siteRepeats > paramd.MaxMicrosateForDis)) ) continue;

        siteLength = tsiteLength & 255;
        // defined one region
        if (ifUserDefinedRegion) {
            if (chr != region_one.chr) {
                continue;
            } else {
                if ( 
                     (loc < region_one.start) 
                     || 
                     ((loc + siteLength * siteRepeats) > region_one.end)
                   ) { continue; }
            }
        }
        // bed filtering
        if (ifUserDefinedBed) {
            // new chr 
            if (tChr != chr) {
                j = 0;
                if (chrMaptoIndex.count(chr) > 0) {
                    tbedChr = beds[chrMaptoIndex[chr]];
                    tChr = tbedChr.chr;
                    tbedRegion = tbedChr.regions_list[j];
                } else { continue; }
            }
            // filtering 
            if (loc < tbedRegion.start) continue;
            if (loc > tbedRegion.end) {
                for (j; j<tbedChr.regions_list.size(); j++) {
                    tbedRegion = tbedChr.regions_list[j];
                    if (loc < tbedRegion.end) { break; }
                }
                if (j >= tbedChr.regions_list.size()) continue;
            }
            if ((loc + siteLength * siteRepeats) > tbedRegion.end) continue;
        }

        // load sites 
        //HomoSite *toneSite = new HomoSite;
        HomoSite toneSite;
        toneSite.chr = chr;
        toneSite.location = loc;
        toneSite.typeLen = siteLength;
        toneSite.homoType = siteBinary;
        toneSite.length = siteRepeats;
        toneSite.frontKmer = frontF;
        toneSite.endKmer = tailF;
        toneSite.bases = bases;
        toneSite.fbases = fbases;
        toneSite.ebases = ebases;

        toneSite.lowcut = ((loc - MAX_READ_LENGTH) > 0) ? (loc - MAX_READ_LENGTH) : 0;
        toneSite.highcut = loc + MAX_READ_LENGTH;

        totalSites.push_back(toneSite);
        totalHomosites++;

        linestream.clear();
        linestream.str("");

    } // end while

}

// bed regions ?
void PolyScan::BedFilterorNot() {
    if (beds.size() > 0) ifUserDefinedBed = true;
}

// losd repeat file
void PolyScan::LoadRepeats(std::ifstream &fin) {
    std::string line;
    std::string chr;
    std::string tempChr = "";
    int start;
    int end;
    int repeat_len;
    int repeat_bin;
    int repeat_times;
    int i = -1;
    while (getline(fin, line)){
	std::stringstream linestream(line);
	linestream >> chr;
	linestream >> start;
	linestream >> repeat_len;
	linestream >> repeat_bin;
	linestream >> repeat_times;
	end = start + repeat_len * repeat_times -1;
	if (chr == tempChr) {
	    RepeatRegion tempRepeatRegion;
	    tempRepeatRegion.start = start;
	    tempRepeatRegion.end = end;
	    repeats[i].repeatregions_list.push_back(tempRepeatRegion);
	} else {
	    ++i;
	    RepeatChr tempRepeatChr;
	    //change chr1 to 1
	    if(chr.find("chr")!= std::string::npos) {
		tempRepeatChr.chr = chr.substr(3);
	    } else {
		tempRepeatChr.chr = chr;
	    }
	    repeats.push_back(tempRepeatChr);
	    repeatChrMaptoIndex.insert(std::pair<std::string, bit16_t>(chr.substr(3), i));
	    RepeatRegion tempRepeatRegion;
	    tempRepeatRegion.start = start;
	    tempRepeatRegion.end = end;
	    repeats[i].repeatregions_list.push_back(tempRepeatRegion);
	    tempChr = chr;
	}
	linestream.clear();
	linestream.str("");
    }
}

//load maf file
void PolyScan::LoadMaffile(std::ifstream &fin, const std::string &prefix) {
    std::ifstream finM;
    std::ofstream output;
    output.open( prefix.c_str() );
    std::string line;
    std::string maffile;
    std::string status; 
    int MSI = 0;
    int len = 0;
    while (getline(fin, line)){
	std::vector<std::string> p = split(line, " ");
	maffile = p[0];
    std::cout << maffile << std::endl;
	if (p.size() == 2)  {
        if(isAllDigit(p[1])){
            len = atoi(p[1].c_str());
        }else{
            status = p[1];
	        if (status == "MSS") MSI = 0;
	        if (status == "MSI-L") MSI = 0;
	        if (status == "MSI-H") MSI = 1;
        }
    }
	if(p.size() == 3) {
        status = p[1];
	    if (status == "MSS") MSI = 0;
	    if (status == "MSI-L") MSI = 0;
	    if (status == "MSI-H") MSI = 1;
        len = atoi(p[2].c_str());

    }
	finM.open(maffile.c_str());
	if (!finM) {
	    std::cerr << "fatal error: failed to open maf file\n";
	    exit(1);
	}
        std::cout << "begin to comput input variables ..." << std::endl;
    double arr[15] = {0};
	pourOUtFeature(finM, arr, len);
	for (int m = 0; m <= 15; m++){
	    output << arr[m] << "\t";
	}
	if (p.size() >= 2) output << MSI << std::endl;
	finM.close();
   }
}

// test sites loading
void PolyScan::TestHomos() {
    for (unsigned long i=0; i<totalHomosites; i++) {
        HomoSite *toneSite = &totalSites[i];
        std::cout << toneSite->chr<<"\t"
                  << toneSite->location<<"\t"
                  << int(toneSite->typeLen)<<"\t"
                  << toneSite->homoType<<"\t"
                  << toneSite->length<<"\t"
                  << toneSite->frontKmer<<"\t"
                  << toneSite->endKmer<<"\t"
                  << sizeof(*toneSite) <<"\n";
    }
}

// split windows
void PolyScan::SplitWindows() {
    Window oneW;
    HomoSite *second;
    HomoSite *first = &totalSites[0];

    oneW._start = first->location;
    oneW._end = oneW._start;
    oneW._chr = first->chr;
    oneW._startSite = oneW._endSite = &totalSites[0];
    for (int i=1; i< totalHomosites; i++) {
        first = &totalSites[i];
        if ( (first->chr == oneW._chr) 
             && 
             (first->location - oneW._start) < paramd.windowSize) {
            continue;
        }
        oneW._end = totalSites[i-1].location + MAX_SPAN_SIZE;
        oneW._endSite = &totalSites[i-1];
        oneW._siteCount = oneW._endSite - oneW._startSite + 1;
        // record one window
        oneW.ChangeStart();
        totalWindows.push_back(oneW);
        totalWindowsNum++;
        oneW._start = first->location;
        oneW._end = oneW._start;
        oneW._chr = first->chr;
        oneW._startSite = oneW._endSite = &totalSites[i];
    }
    oneW._end = first->location + MAX_SPAN_SIZE;
    oneW._endSite = first;
    oneW._siteCount = oneW._endSite - oneW._startSite + 1;
    // record this window
    oneW.ChangeStart();
    totalWindows.push_back(oneW);
    totalWindowsNum++;
}

//test whether a string consist of numbers or not.
bool PolyScan::isAllDigit(const std::string & str) {
    for (int i=0; i < str.size(); i++){
        int tmp = (int)str[i];
        if (tmp >= 48 && tmp <= 57){
            continue;
        }else{
            return false;
        }
    }
    return true;
}
    

// test windows
void PolyScan::TestWindows() {
    Window *oneW;
    for (int i=0; i< totalWindowsNum; i++) {
        oneW = &totalWindows[i];
        std::cout << oneW->_chr <<"\t"
                  << oneW->_start <<"\t"
                  << oneW->_siteCount <<"\t"
                  << oneW->_startSite->chr<<"\t"
                  << oneW->_startSite->location<<"\t"
                  << oneW->_endSite->chr<<"\t"
                  << oneW->_endSite->location<<"\n";
    }
}

// initial distribution
void PolyScan::InithializeDistributions() {
    for (int i=0; i< totalWindowsNum; i++) {
        totalWindows[i].InitialDisW();
    }
}

// release distribution
void PolyScan::releaseDistributions() {
    for (int i=0; i< totalWindowsNum; i++) {
        totalWindows[i].ClearDis();
    }
}

// output distribution
void PolyScan::outputDistributions() {
    for (int i=0; i< totalWindowsNum; i++) {
        totalWindows[i].OutputDisW();
    }
}

// get distribution 
//void PolyScan::GetHomoDistribution( std::ofstream &fout ) {
void PolyScan::GetHomoDistribution( Sample &oneSample, const std::string &prefix ) {
    oneSample.iniOutput(prefix);
    std::vector< SPLIT_READ > readsInWindow;
    for (int i=0; i< totalWindowsNum; i++) {
        totalWindows[i].InitialDisW();
        totalWindows[i].GetDistribution(readsInWindow);
        totalWindows[i].PouroutDisW(oneSample);
        totalWindows[i].DisGenotypingW(oneSample);
        totalWindows[i].ClearDis();
        readsInWindow.clear();
        std::cout << "window: " << i << " done...:" <<  totalWindows[i]._chr << ":" << totalWindows[i]._start << "-" << totalWindows[i]._end << std::endl;
    }
    // FDR
    oneSample.calculateFDR();
    oneSample.pourOutSomaticFDR();
    // MSI score
    oneSample.pourOutMsiScore();
    oneSample.closeOutStream();
    oneSample.VerboseInfo();

}

void PolyScan::GetHomoTumorDistribution( Sample &oneSample, const std::string &prefix ) {
    oneSample.iniTumorDisOutput(prefix);
    std::vector< SPLIT_READ > readsInWindow;
    for (int i=0; i< totalWindowsNum; i++) {
        totalWindows[i].InitialTumorDisW();
        totalWindows[i].GetTumorDistribution(readsInWindow);
        totalWindows[i].PourTumoroutDisW(oneSample);
        totalWindows[i].PouroutTumorSomatic(oneSample);
        totalWindows[i].ClearTumorDis();
        readsInWindow.clear();
        std::cout << "window: " << i << " done...:" <<  totalWindows[i]._chr << ":" << totalWindows[i]._start << "-" << totalWindows[i]._end << std::endl;
    }
    // FDR
    // MSI score
    oneSample.pourOutMsiScore();
    oneSample.closeOutStream();
    oneSample.VerboseInfo();

}

std::vector<std::string> PolyScan::split(const std::string &s, const std::string &seperator){
    std::vector<std::string> result;
    typedef std::string::size_type string_size;
    string_size i = 0;
    while(i != s.size()){
	int flag = 0;
	while(i != s.size() && flag == 0){
	    flag = 1;
	    for(string_size x = 0; x < seperator.size(); ++x)
	        if(s[i] == seperator[x]){
		    ++i;
		    flag = 0;
		    break;
	        } 
	}
	flag = 0;
	string_size j = i;
	while(j != s.size() && flag == 0){
	    for(string_size x = 0; x < seperator.size(); ++x)
		if(s[j] == seperator[x]){
		    flag = 1;
		    break;
		}
	    if(flag == 0)
		++j;
	}
	if(i != j){
	    result.push_back(s.substr(i, j-i));
	    i = j;
	}
    }
    return result;
}

//Get and output features
void PolyScan::pourOUtFeature(std::ifstream &maffile, double arr[],int actuallen, int defaultlen) {
    RepeatChr trepeatChr;
    RepeatRegion trepeatRegion;
    std::string tchr = "";
    std::string Chrom;
    int Start_Position;
    int End_Position;
    std::string Variant_Type;
    std::string Tumor_Seq_Allele2;
    std::string Tumor_Sample_Barcode;
    std::string line;

    //features
    double T_sns = 0.0;
    double S_sns = 0.0;
    double T_ins = 0.0;
    double S_ins = 0.0;
    double T_del = 0.0;
    double S_del = 0.0;
    double T_ind = 0.0;
    double S_ind = 0.0;
    double T = 0.0;
    double S = 0.0;
    double ratio_sns = 0.0;
    double ratio_ind = 0.0;
    double ratio = 0.0;
    double seq = 0.0;
    double PI = 0.0;
    double PD = 0.0;

    int j;
    int i = 0;
    while (getline(maffile, line)){
        i++;
        if(i >= 5) {
            std::vector<std::string> p = split(line, "\t");
            Chrom = p[4];
            Start_Position =  atoi(p[5].c_str());
            End_Position = atoi(p[6].c_str());
            Variant_Type = p[9];
            Tumor_Seq_Allele2 = p[12];
            Tumor_Sample_Barcode = p[15];
            if ( Variant_Type == "SNP" ) T_sns += 1;
            else if ( Variant_Type == "DEL" ) {
                T_del += (End_Position - Start_Position + 1) ;
		 T_ind += (End_Position - Start_Position + 1) ;
            }
            else if ( Variant_Type == "INS" ){
                T_ins += Tumor_Seq_Allele2.size();
                T_ind += Tumor_Seq_Allele2.size();
            }

            if(tchr != Chrom) {
                j = 0;
                if (repeatChrMaptoIndex.count(Chrom) > 0) {
                    trepeatChr = repeats[repeatChrMaptoIndex[Chrom]];
                    tchr = trepeatChr.chr;
                    trepeatRegion = trepeatChr.repeatregions_list[j++];
                } else { continue; }
            }
            if (End_Position < trepeatRegion.start) continue;
            if (Start_Position > trepeatRegion.end){
                for (j; j < trepeatChr.repeatregions_list.size() && Start_Position > trepeatRegion.end; j++){
                    trepeatRegion = trepeatChr.repeatregions_list[j];
                }
                if ( j >= trepeatChr.repeatregions_list.size()) continue;
            }
            if ( Variant_Type == "SNP" ){
                if ( (trepeatRegion.start <= Start_Position) && (Start_Position <= trepeatRegion.end)) {
                    S_sns += 1;
                }
            }
            else if ( Variant_Type == "DEL" ){
                if ((Start_Position <= trepeatRegion.start) && (End_Position >= trepeatRegion.start) && (End_Position <= trepeatRegion.end)){
                    S_del += (End_Position - trepeatRegion.start + 1);
                    S_ind += (End_Position - trepeatRegion.start + 1);
                }
                else if ((Start_Position <= trepeatRegion.start)&&(End_Position >= trepeatRegion.end)) {
                    S_del += (trepeatRegion.end - trepeatRegion.start + 1);
                    S_ind += (trepeatRegion.end - trepeatRegion.start + 1);
                }
                else if ((Start_Position >= trepeatRegion.start) && (End_Position <= trepeatRegion.end)){
                    S_del += (End_Position - Start_Position +1);
                    S_ind += (End_Position - Start_Position +1);
                }
                else if ((Start_Position >= trepeatRegion.start)&&(Start_Position <= trepeatRegion.end)&&(End_Position >= trepeatRegion.end)){
                    S_del += (trepeatRegion.end - Start_Position + 1);
                    S_ind += (trepeatRegion.end - Start_Position + 1);
                }
            }
            else if ( Variant_Type == "INS" ){
                if((trepeatRegion.start <= Start_Position) && (Start_Position < trepeatRegion.end)){
                    S_ins += Tumor_Seq_Allele2.size();
		    S_ind += Tumor_Seq_Allele2.size();
                }
            }
        }
    }
    if( actuallen > 0){
        T_sns = T_sns/actuallen;
        T_ind = T_ind/actuallen;
        S_sns = S_sns/actuallen;
        S_ind = S_ind/actuallen;
        S_ins = S_ins/actuallen;
        S_del = S_del/actuallen;
        T_ins = T_ins/actuallen;
        T_del = T_del/actuallen;
    }else{
        T_sns = T_sns/defaultlen;
        T_ind = T_ind/defaultlen;
        S_sns = S_sns/defaultlen;
        S_ind = S_ind/defaultlen;
        S_ins = S_ins/defaultlen;
        S_del = S_del/defaultlen;
        T_ins = T_ins/defaultlen;
        T_del = T_del/defaultlen;
    }
    T = T_sns + T_ind;
    S = S_sns + S_ind;
    if(T_sns != 0.0) ratio_sns = S_sns/T_sns;
    if(T_ind != 0.0) ratio_ind = S_ind/T_ind;
    if(T != 0) ratio = S/T;
    if(T_ins != 0) PI = S_ins/T_ins;
    if(T_del != 0) PD = S_del/T_del;
    if(PD != 0) seq = PI/PD;

    arr[0] = T_sns;
    arr[1] = T_ind;
    arr[2] = S_sns;
    arr[3] = S_ind;
    arr[4] = T_ins;
    arr[5] = T_del;
    arr[6] = S_ins;
    arr[7] = S_del;
    arr[8] = T;
    arr[9] = S;
    arr[10] = ratio_sns;
    arr[11] = ratio_ind;
    arr[12] = ratio;
    arr[13] = PI;
    arr[14]= PD;
    arr[15] = seq;
}

