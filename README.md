MSIsensor1.1
======================
MSIsensor is a program written in C++ and python to detect replication slippage variants at microsatellite regions, and differentiate them as somatic or germline. When given paired tumor and normal sequence data, it builds a distribution for expected (normal) and observed (tumor) lengths of repeated sequence per microsatellite, and compares them using Pearson's Chi-Squared Test. MSIsensor is a C++ program for automatically detecting somatic and germline variants at microsatellite regions. When given tumor only, it computes comentropy for the distribution of each site. Sites whose information entropy exceeds the threshold are marked as somatic sites. Finally, the ratio of the number of somatic sites to the total number of microsatellite points is calculated as the MSI score. when given maf file,it will build data frame for train or test. Comprehensive testing indicates MSIsensor is an efficient and effective tool for deriving MSI status from standard tumor-normal paired sequence data, tumor only and mutation data. 

Usage
-----

        Version 2.0
        Usage:  msisensor <command> [options]

Key commands:

        scan            scan homopolymers and miscrosatelites
        msi             msi scoring

msisensor scan [options]:
       
       -d   <string>   reference genome sequences file, *.fasta format
       -o   <string>   output homopolymer and microsatelites file

       -l   <int>      minimal homopolymer size, default=5
       -c   <int>      context length, default=5
       -m   <int>      maximal homopolymer size, default=50
       -s   <int>      maximal length of microsatellite, default=5
       -r   <int>      minimal repeat times of microsatellite, default=3
       -p   <int>      output homopolymer only, 0: no; 1: yes, default=0
       
       -h   help
 
msisensor msi [options]:

       -d   <string>   homopolymer and microsatellites file
       -n   <string>   normal bam file ( bam index file is needed )
       -t   <string>   tumor  bam file ( bam index file is needed. If you have tumor-normal paired sequence data, you need to specify them with -n and -t. If you have tumor only, you just need to specify tumor file with -t. )
       -o   <string>   output distribution file

       -e   <string>   bed file, to select a few resions
       -f   <double>   When using paired tumor-normal sequence data, FDR threshold for somatic sites detection, default=0.05 
       -v   <double>   When using tumor sequence data, comentropy threshold for somatic sites detection, default=0.5
       -c   <int>      coverage threshold for msi analysis, WXS: 20; WGS: 15, default=20
       -r   <string>   choose one region, format: 1:10000000-20000000
       -l   <int>      mininal homopolymer size, default=5
       -p   <int>      mininal homopolymer size for distribution analysis, default=10
       -m   <int>      maximal homopolymer size for distribution analysis, default=50
       -q   <int>      mininal microsatellites size, default=3
       -s   <int>      mininal number of repeats in microsatellites for distribution analysis, default=5
       -w   <int>      maximal microsatellites size for distribution analysis, default=40
       -u   <int>      span size around window for extracting reads, default=500
       -b   <int>      threads number for parallel computing, default=1
       -x   <int>      output homopolymer only, 0: no; 1: yes, default=0
       -y   <int>      output microsatellite only, 0: no; 1: yes, default=0
       -P   <string>   the path of maf files(If you want to use mutation data, you just need to specify homopolymer and microsatellites file and maf files path with -d and -P respectively)
       
       -h   help

Install
-------
The Makefile assumes that you have the samtools source code in an environment variable `$SAMTOOLS_ROOT`. 

you don't know what that means, then simply follow these steps from any directory that you have permissions to write into:
Install some prerequisite packages if you are using Debian or Ubuntu:

    sudo apt-get install git libbam-dev zlib1g-dev

If you are using Fedora, CentOS or RHEL, you'll need these packages instead:

    sudo yum install git samtools-devel zlib-devel

Download the samtools-0.1.19 from SOURCEFORGE (http://sourceforge.net/projects/samtools/files/samtools/0.1.19):

    tar jxf samtools-0.1.19.tar.bz2
    cd samtools-0.1.19
    make
    export SAMTOOLS_ROOT=$PWD

Clone the msisensor repos, and build the `msisensor` binary:

    git clone https://github.com/ding-lab/msisensor.git
    cd msisensor
    make

Now you can put the resulting binary where your `$PATH` can find it. If you have su permissions, then
I recommend dumping it in the system directory for locally compiled packages:

    sudo mv msisensor /usr/local/bin/

We have also provided pre-build binary distributions for Linux x86_64 and Mac OS X in directory: ./binary

    msisensor_Linux_x86_64: for Linux x86_64
    msisensor_Mac_OS_X    : for Mac OS X

If you want to use mutation data to train or test, you must install python and sklearn model

Example
-------
1. Scan microsatellites from reference genome:
  
        msisensor scan -d referen.fa -o microsatellites.list

2. Msi scorring: 

   If you have paired tumor-normal sequence data, you can use -t and -n to specify the tumor and normal files, respectively:
        msisensor msi -d microsatellites.list -n normal.bam -t tumor.bam -e bed.file -o output.prefix -l 1 -q 1 -b 2
   
   If you only have the tumor file, just use the -t option:
        msisensor msi -d microsatellites.list -t tumor.bam -e bed.file -o output.prefix -l 1 -q 1 -b 2 

   If you have maf file, you can use -P option:
        msisensor msi -d microsatellites.list -P /Path/to/maf -o data_frame
   ( /Path/to/maf: a file of the list of maf specimens, with each line being the absolute path to one sample. If you want to train model, being the absolute path should have microsatellite status. 
    Microsatellite status must be MSS/MSI-L/MSS-H. You can specify the gene length(/Mb) after microsatellite status, default 38.
    If you want to test, each line just contain path and gene length. Path, microsatellite status and gene length are separated by sapce.
    A File of the list of maf specimens as follow:
    /home/sczhaod/TCGA/mutation_data/gdac.broadinstitute.org_COAD.Mutation_Packager_Oncotated_Calls.Level_3.2016012800.0.0/TCGA-A6-2672-01.hg19.oncotator.hugo_entrez_remapped.maf.txt MSI-H
    /home/sczhaod/TCGA/mutation_data/gdac.broadinstitute.org_COAD.Mutation_Packager_Oncotated_Calls.Level_3.2016012800.0.0/TCGA-A6-2674-01.hg19.oncotator.hugo_entrez_remapped.maf.txt MSS
    /home/sczhaod/TCGA/mutation_data/gdac.broadinstitute.org_COAD.Mutation_Packager_Oncotated_Calls.Level_3.2016012800.0.0/TCGA-A6-2676-01.hg19.oncotator.hugo_entrez_remapped.maf.txt MSI-H
    /home/sczhaod/TCGA/mutation_data/gdac.broadinstitute.org_COAD.Mutation_Packager_Oncotated_Calls.Level_3.2016012800.0.0/TCGA-A6-2677-01.hg19.oncotator.hugo_entrez_remapped.maf.txt MSS
    /home/sczhaod/TCGA/mutation_data/gdac.broadinstitute.org_COAD.Mutation_Packager_Oncotated_Calls.Level_3.2016012800.0.0/TCGA-A6-2678-01.hg19.oncotator.hugo_entrez_remapped.maf.txt MSS
    
    You can refer to the format of the maf file in test.

    After getting data_frame, you can use "python train.py" to train or "python test.py" to test )


   Note: normal and tumor bam index files are needed in the same directory as bam files 

Output
-------
There will be one microsatellite list output in "scan" step :
 
   microsatellites.list: microsatellite list output ( columns with *_binary means: binary conversion of DNA bases based on A=00, C=01, G=10, and T=11 )

        chromosome      location        repeat_unit_length     repeat_unit_binary    repeat_times    left_flank_binary     right_flank_binary      repeat_unit_bases      left_flank_bases       right_flank_bases
        1       10485   4       149     3       150     685     GCCC    AGCCG   GGGTC
        1       10629   2       9       3       258     409     GC      CAAAG   CGCGC
        1       10652   2       2       3       665     614     AG      GGCGC   GCGCG
        1       10658   2       9       3       546     409     GC      GAGAG   CGCGC
        1       10681   2       2       3       665     614     AG      GGCGC   GCGCG

Msi scorring step will give 4 or 3 output files based on given output prefix acording to whether you use paired data or tumor data:
When using paired tumor-normal sequence data, output files are :
        output.prefix
        output.prefix_dis
        output.prefix_germline
        output.prefix_somatic

1. output.prefix: msi score output

        Total_Number_of_Sites   Number_of_Somatic_Sites %
        640     75      11.72

2. output.prefix_dis: read count distribution (N: normal; T: tumor)

        1 10529896 CTTTC 15[T] GAGAC
        N: 0 0 0 0 0 0 0 1 0 0 8 9 1 7 17 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
        T: 0 0 0 0 0 0 0 0 0 1 19 14 17 9 32 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 

3. output.prefix_somatic: somatic sites detected ( FDR: false discovery rate ) 

        chromosome   location        left_flank     repeat_times    repeat_unit_bases    right_flank      difference      P_value    FDR     rank
        1       16200729        TAAGA   10      T       CTTGT   0.55652 2.8973e-15      1.8542e-12      1
        1       75614380        TTTAC   14      T       AAGGT   0.82764 5.1515e-15      1.6485e-12      2
        1       70654981        CCAGG   21      A       GATGA   0.80556 1e-14   2.1333e-12      3
        1       65138787        GTTTG   13      A       CAGCT   0.8653  1e-14   1.6e-12 4
        1       35885046        TTCTC   11      T       CCCCT   0.84682 1e-14   1.28e-12        5
        1       75172756        GTGGT   14      A       GAAAA   0.57471 1e-14   1.0667e-12      6
        1       76257074        TGGAA   14      T       GAGTC   0.66023 1e-14   9.1429e-13      7
        1       33087567        TAGAG   16      A       GGAAA   0.53141 1e-14   8e-13   8
        1       41456808        CTAAC   14      T       CTTTT   0.76286 1e-14   7.1111e-13      9

4. output.prefix_germline: germline sites detected
    
        chromosome   location        left_flank     repeat_times    repeat_unit_bases    right_flank      genotype
        1       1192105 AATAC   11      A       TTAGC   5|5
        1       1330899 CTGCC   5       AG      CACAG   5|5
        1       1598690 AATAC   12      A       TTAGC   5|5
        1       1605407 AAAAG   14      A       GAAAA   1|1
        1       2118724 TTTTC   11      T       CTTTT   1|1

When using tumor sequence data, output files are :
        output.prefix
        output.prefix_dis
        output.prefix_somatic

1. output.prefix: msi score output

        Total_Number_of_Sites   Number_of_Somatic_Sites %
        640     75      11.72

2. output.prefix_dis: read count distribution (T: tumor)

        1 10529896 CTTTC 15[T] GAGAC
        T: 0 0 0 0 0 0 0 0 0 1 19 14 17 9 32 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0

3. output.prefix_somatic: somatic sites detected  

        chromosome   location        left_flank     repeat_times    repeat_unit_bases   comentropy 
        1	16248728	ACCTC	11	T	AAAGG	0.68491
	1	16890814	GAGGT	12	A	TTATT	1.09340

when using mutation data, output files are:
        data_frame
for examplei(train data):

20.6053	2.44737	0.263158	0.368421	1.31579	1.13158	0.315789	0.0526316	23.0526	0.631579	0.0127714	0.150538	0.0273973	0.24	0.0465116	5.16	1
2.55263	0.131579	0.0263158	0.0263158	0.105263	0.0263158	0.0263158	0	2.68421	0.0526316	0.0103093	0.2	0.0196078	0.25	0	0	0
3.15789	0.105263	0.105263	0.0263158	0.0526316	0.0526316	0.0263158	0	3.26316	0.131579	0.0333333	0.25	0.0403226	0.5	0	0	0
23.8684	1.07895	0.368421	0.105263	0.368421	0.710526	0.105263	0	24.9474	0.473684	0.0154355	0.097561	0.0189873	0.285714	0	0	1

If you don't have maf file and have mutation data, you can provide mutation data as follow to train or test:
T_sns T_ind S_sns S_ind T_ins T_del S_ins S_del T S ration_sns tation_ind ratio PI PD seq microsatellite_status
(If for testing, you don't need the last column "microsatellite_status")

Test sample
-------
We provided one small sized sample data (tumor and matched normal bam files) for user to try msi scoring step.
It is very simple to run this test using sample data:

        cd ./test
        bash run.sh

Contact
-------
Please contact Beifang Niu by bniu@genome.wustl.edu and Kai Ye by kye@genome.wustl.edu if you have any questions.

