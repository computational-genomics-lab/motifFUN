[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![minimal R version](https://img.shields.io/badge/R%3E%3D-3.5.0-6666ff.svg)](https://cran.r-project.org/)

# motifFUN : A tool for searching user-defined motif along with functional annotation
An R package for searching user defined motif with functional annotation features
The **motifFUN** package is an R package designed to search for the motifs of interest using user-defined regular expression searches and hidden Markov models (HMM) along with a prediction of some functional criteria such as signal peptides, subcellular location of Eukaryotes and transmembrane helices.

Several command line tools and web-interfaces exist to perform predictions of individual motifs and domains (SignalP, TargetP, TMHMM) however the interface that combines the outputs in a single flexible workflow is lacking, So that developed a motifFUN package to fulfill that gap.

## Installation motifFUN

To install **motifFUN** package:

``` r
library("devtools")
install_github("computational-genomics-lab/motifFUN")
library("motifFUN")
```

## 1. OUTLINE

The **motifFUN** package provides a platform to search motifs and build automated multi-step secretome prediction pipelines that can be applied to largescale protein datasets. The features of this package is described as below.

         
**ORF extraction** 

In this package get_orf function extracts six-frame translation of all ORF(open reading frame) of the input nucleotide sequence fasta file. In order to run this, we must have EMBOSS program installed in the client machine that has the getorf function.

**Pattern Search**

This package searching user-defined motif of interest using a regular expression for both genome and protein fasta sequences.

+ This function gives output as fasta file having motif sequences.
+ This file will be used for HMM search and Functional domains.

**HMMER**

+ The motifFUN package aligns the PATT_REG search results using **MAFFT** and builds an HMM profile based on the multiple sequence alignment result using the hmmbuild program from HMMER. 

+ The HMM profile is used to search across ORF of the genome of interest using the hmmsearch binary from HMMER. 

+ The search step will retain sequences with significant hits to the profile of interest. motifFUN also combines the redundant sequences found in both PATT_REG and HMM searches into a single data set and generate fasta file. 

+ motifFUN reads and returns the HMM profile to the user and allows for the creation of a **MOTIF plot** using ggplot2.

**Functional domains**

Secretome prediction often involves multiple steps. 

+ The motifFUN package starts with the prediction of short signal peptides(signalP) at the N-terminal end of a protein.
+ Predicts the sub cellular location of eukaryotes(TargetP), by using get_targetp function.
+ Predicts trans membrane helices(TMHMM) by using get_tmhmm function.

![motif](https://user-images.githubusercontent.com/43873570/56793098-cbe0d200-6828-11e9-9df9-77b6ea027616.png)
 
 
 **Summary of Functionality**:
 
 |S/No |Function Name       |    Function                                            |Input parameters           |Dependencies      |
|-----|--------------------|----------------------------------------------------|---------------------------|----------------------|
|  1 |get_orf()        |Extraction of ORF from genome sequence              |Genome fasta file          | EMBOSS getorf                        
| 2  |orf_discard()     |Discards residues by lower limit to upper limit(length)|protein fasta file, upper limit, lower limit | No dependencies |                                                                                                                
| 3  |parse_file()      |To covert single fasta file into multiple fasta files |Input file, Output file name, Number of proteins to be parse|No dependencies |
| 4  |pattern.search()  |For searching user_defined motif| Fasta sequence file, Regular expression pattern|No dependencies |
| 5  |hmm.search()      |For searching sequence motifs from sequence database|Original input fasta file, Pattern search candidates, MAFFT path, HMMER path| HMMER|                
| 6  |get_signalp()    |To predict secretory proteins|version, Complete Path Signalp, Input file, organism type|signalp-3.0/signalp-5.0|  
|  7 |get_targetp()   |To predict subcellular location of eukaryotes      |Path of targetp, organism group, input file|TargetP 1.1|
|  8 |get_tmhmm()       |To predict transmembrane helices                   |complete Path of tmhmm, Input file| Tmhmm |
| 9  |hmm.plot()        |The plot shows the bits (amino acid scores) of each amino acid and its position in the HMM profile|motif candidate data frame|No dependencies |
|10  |summary_motifs()  |To extract all the non-redundant sequences & a summary table with the information about the motifs|motif candidates, motif pattern without range, signalp vesrion, input filename used in signalp|No dependencies|
|    |                |     |  
 
 
Due to limitations imposed by the external dependencies, some of the motifFUN wrapper(get_signalp, get_targetp, get_tmhmm) functions won't work in Windows or Mac, however, are fully functional on Linux.

## 2. REQUIREMENTS

 R packages:

|ID|Name|Function|
|------|---|-----|
|1|seqinr<!-- .element: style="text-align:center;" -->|Reading fasta file<!-- .element: style="text-align:right;" -->|
|2|ggplot2<!-- .element: style="text-align:center;" -->|used for producing Plots<!-- .element: style="text-align:right;" -->|

External Tools:

|ID|Name|Function|
|------|---|-----|
|3|[signalp 3.0](http://www.cbs.dtu.dk/cgi-bin/sw_request?signalp+3.0),[signalp 5.0](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?signalp)<!-- .element: style="text-align:center;" -->|To predict secretary proteins<!-- .element: style="text-align:right;" -->|
|4|[targetp 1.1](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?targetp)<!-- .element: style="text-align:center;" -->|To prediction of sub cellular location of eukaryotes<!-- .element: style="text-align:right;" -->|
|5|[tmhmm](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?tmhmm)<!-- .element: style="text-align:center;" -->|To predicts trans membrane helices<!-- .element: style="text-align:right;" -->|
|6|[mafft](https://mafft.cbrc.jp/alignment/software/)<!-- .element: style="text-align:center;" -->|For alignment of sequences<!-- .element: style="text-align:right;" -->|
|7|[HMMER](http://hmmer.org/)<!-- .element: style="text-align:center;" -->|For searching motifs<!-- .element: style="text-align:right;" -->|
|8|[EMBOSS](http://emboss.sourceforge.net/download/#Stable/)<!-- .element: style="text-align:center;" -->|For extraction of six-frame tranlastion of ORF<!-- .element: style="text-align:right;" -->|


## 3. EXTERNAL SOFTWARES

The **motifFUN** package uses signalp, targetp, tmhmm for prediction of extracellular proteins that are secreted via classical pathways, getorf for extraction of ORF.

MAFFT and HMMER3 used to perform the hidden Markov model search across the results from the REGEX step. 

These packages should be installed before running any of the motifFUN functions.

**3.1.Downloading EMBOSS**

* EMBOSS getorf program finds and outputs the sequences of open reading frames (ORFs) in one or more nucleotide sequences. 

* For installing in Linux wget ftp://emboss.open-bio.org/pub/EMBOSS/EMBOSS-6.6.0.tar.gz

* For installing through Cygwin http://emboss.sourceforge.net/download/cygwin.html

* For windows ftp://emboss.open-bio.org/pub/EMBOSS/windows/

Read instructions and install

**3.2.Downloading signalP**

**3.2.1.signalp-3.0**

* This version will run on the most common UNIX platforms.

* Download stand-alone [signalp3.0](http://www.cbs.dtu.dk/cgi-bin/sw_request?signalp+3.0)

* Unpack the archive

<font size="2">tar -zxvf signalp-3.0.Linux.tar.Z</font>

 cd signalp-3.0

Edit"General settings" at the top of the signalp file. Set the value of 'SIGNALP' variable to be a path to your signalp-3.0 directory. Other variables usually do not require changes. For more details please check signalp-3.0.readme.

**3.2.2.signalp 5.0**

* This version will run on the most common UNIX platforms.

* Download stand-alone [signalp5.0](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?signalp)

* Unpack the archive

<font size="2">tar -zxvf signalp-5.0.Linux.tar.gz</font>

 cd signalp-5.0

**3.3.Downloading targetp-1.1**

**3.3.targetp-1.1**

* This version will run on the most common UNIX platforms.

* Download stand-alone [targetp-1.1](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?targetp)

* Unpack the archive

<font size="2">tar -zxvf targetp-1.1b.Linux.tar.Z</font>

cd targetp-1.1

Edit the paragraph labeled "GENERAL SETTINGS, customize" at the top of the targetp file. Set values for 'TARGETP' and 'TMP' variables. Ensure, that the path to targetp does not exceed 60 characters, otherwise targetp-1.1 might fail.

**3.4.Downloading tmhmm,v.2.0**

**3.4.tmhmm,v.2.0**

* This version will run on the most common UNIX platforms.

* Download stand-alone [tmhmm v,2.0](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?tmhmm)

* Unpack the archive

<font size="2">tar -zxvf tmhmm-2.0c.Linux.tar.gz</font>

 cd tmhmm-2.0c

**3.5.Downloading and installing MAFFT**

**3.5.MAFFT**

MAFFT is a multiple sequence alignment program that uses Fourier-transform algorithms to align multiple sequences[@Katoh2002]. We recommend downloading and installing MAFFT by following the instructions and steps in the [MAFFT installation](https://mafft.cbrc.jp/alignment/software/) web site.

     wget https://mafft.cbrc.jp/alignment/software/mafft_7.427-1_amd64.deb

**Linux/OS X Users**

Make sure that you remember the directory in which MAFFT is installed, or if the installation is successful, make sure to obtain the path via bash/tsh/console:

On the Ubuntu window, run the following command to download MAFFT package.

    $ wget https://mafft.cbrc.jp/alignment/software/mafft_7.427-1_amd64.deb

For extraction

    $ sudo dpkg -i mafft_7.427-1_amd64.deb
    [sudo] password for username: (Type the password that was set in step 2)

Check location and version number of MAFFT.

    which mafft
    /usr/local/bin/mafft

For more information about MAFFT go to the MAFFT website: http://mafft.cbrc.jp/

**Windows Users**

MAFFT comes in two main distributions for windows:

+ [A version that requires a UNIX-like interface called Cygwin](https://mafft.cbrc.jp/alignment/software/windows_cygwin.html)

+ [An “all-in-one” version“](https://mafft.cbrc.jp/alignment/software/windows_without_cygwin.html)

Please, download and install the all-in-one version.

**3.6.Downloading and installing HMMER**

**3.6.HMMER**

HMMER is used for searching sequence databases for sequence homologs. It uses hidden Markov models[@Finn2011] (profile HMMs) to search for sequences with hits to similar patterns than the profile. We use three main HMMER tools:

+ **hmmbuild** to create the HMM database from the sequences contained in the PATT_REG step of motifFUN
+ **hmmpress** converts the HMM database into a format usable by other HMMER programs
+ **hmmsearch** to execute the HMM search in our sequence queries based on the HMM profile

The motifFUN package requires all of these tools. A correct HMMER installation will install all three programs.

**Linux/OS X users**

We recommend downloading and installing HMMER by following the instructions and steps in the HMMER installation web site. Make sure that you remember the directory in which HMMER is installed, or if the installation is successful, make sure to obtain the path via bash/tsh/console:

    which hmmbuild
    which hmmpress
    which hmmsearch
    
    /usr/local/bin/hmmbuild
    /usr/local/bin/hmmpress
    /usr/local/bin/hmmsearch

For more information about HMMER go to the HMMER website: http://hmmer.org/

**Windows users**

To use the motifFUN package in Windows, the user must download the Windows binaries of HMMER. motifFUN will not work with any other version of HMMER.

## 4.WORK FLOW

**4.1.Input Data**

The motifFUN package design to predicts sequence motifs for both nucleotide and amino acid sequence by using a user-defined regular expression. This package supports both Gene FASTA and protein FASTA file as input. 

The **getorf** function can be used to translate gene fasta to protien fasta. 

INPUT

+ Protein FASTA
+ Genome FASTA
+ Regular expression

```{r}
library(motifFUN)
library(seqinr)
library(ggplot2)
#orf_fasta = system.file("tests","testfile.fasta", package = "motifFUN")
orf_fasta = "/home/cglab/motifFUN/tests/testfile.fasta"
```

**4.2.ORF extraction**

Emboss getorf is a software tool to finds and outputs the sequences of open reading frames (ORFs) in one or more nucleotide sequences. An ORF is a part of the reading frame has the ability to be translated. ORF having the continuous stretch of codons begin with a start codon(AUG) and stop codon(UAA, UAG or UGA). 

+ get_orf function use the EMBOSS getorf and gives sequences of open reading frames (ORFs) in one or more nucleotide sequences.
+ The get_orf function requires to provides a complete path of the getorf, input genome file, output.file.
+ If user did not provide output.file then by default function will set outputfile name .
+ If user did not provide getorf path then by default function will take getorf path.

```{r}
ORF_filename <- get_orf(getorf.path= NULL, input.file= orf_fasta, output.file = NULL)
```

**4.3.ORF discard**

Bhattacharjee et al.[@Bhattacharjee2006] noted that P. infestans–candidate effectors that contain at least 100 residues after the predicted SS cleavage site, highlighting the conservation of the RxLR motif.

+ Based on that study I have included discard function to discard residues. 
+ user can use this function discards residues by his requirement from **lower limit** to **upper limit**. 

```{r}
ORF_disacrd_filename <- orf_discard(orf_file= ORF_filename, upper.limit= 1800, lower.limit= 100)
```

**4.4.Motif Pattern Search**

motifFUN package has the function pattern.search to perform the search of the motif of interest. 

+ The function pattern. search requires the list of SeqFastadna objects and the gene family of interest.
+ The read.fasta function from the seqinr package reads the FASTA amino acid file into R, creating a list of SeqFastadna objects that represent each of the translated ORF’s from the original FASTA file

Example with sample data:

Here we show an example to search for sequences with RxLR-EER motifs from the 63 ORF subset of testfile.fasta This ORF example data set contains 45 sequences with RxLR-EER motifs.

```{r}
pattern <- "^\\w{10,40}\\w{1,96}R\\wLR\\w{1,40}eer"
PATT_REG <- pattern.search(fasta.file = ORF_disacrd_filename, reg.pat = pattern)
head( PATT_REG, n = 2)
```

This function generates one fasta file having motif sequences, We observe that the PATT_REG object has 24 sequences with the RxLR motif. These sequences will be aligned using  MAFFT and used to build an HMM profile to search for similar sequences. 

This fasta file will become input to signalp, targetp, tmhmm functional programs.

**4.5.MAFFT & HMMER Search**

To perform the HMM search and obtain all possible motif from a proteome, motifFUN uses the PATT_REG results as a template to create an HMM profile and perform a search across the proteome of interest.  MotifFUN package have the hmm.search function in order to perform this search. The hmm.search function requires a local installation of MAFFT and HMMER in order to perform the searches.

The absolute paths of the binaries must be specified in the mafft.path and hmmer.path options of the hmm.search function.

    Note for Windows users: Please use the ABSOLUTE PATH for HMMER and MAFFT or motifFUN will not work (e.g. mafft.path ="C:/User/Banana/Desktop/mafft/")

In addition, the hmm.function requires the path of the original FASTA file containing the translated ORF in the original.seq parameter of the function. hmm.search will use this file as a query in the hmm search software from HMMER, and search for all sequences with hits against the HMM profile created with the PATT_REG results.

We are performing a hmmsearch in our example data set. This function requires original FASTA file location (stored in the filepath object), the location of the MAFFT binary and the location of the HMMER binaries:

If user did not provide MAFFT & HMMER path then by default function will take MAFFT & HMMER paths.

```{r}
motif_candidates <- hmm.search(original.seq = ORF_disacrd_filename, regex.seq = PATT_REG, mafft.path = NULL, hmm.path = NULL)
```

The hmm.search function has resulted in 41 motif candidates. As a reminder, we used the PATT_REG results of an RxLR motif search, so we can consider this hmm.search results as RxLR candidate effectors.

The hmm.search object returns a list of 3 elements:

+ The PATT_REG sequences used to build the HMM profile in a SeqFastadna class
+ The sequences from the original translated ORF files with hits to the HMM profile in a SeqFastadna class
+ The HMM profile table created by HMMER’s hmmbuild as a data frame

This function combines and returns  PATT_REG and HMM search resulted in fasta file, User can use this file as input to functional domains or can use PATT_REG resulted fasta file as input to functional domains depends on the user.

**4.6.signalp**

The signal peptide is a short peptide usually 16-30 amino acids long [@Hemminger1998]present at the N-terminus of the majority of newly synthesized proteins that are bound towards the secretory pathway. These proteins incorporate those that reside either inside certain organelles, secreted from the cell, or inserted into most cellular membranes.

signalP is a software tool to predicts the presence and location of signal peptide cleavage sites[@Nielsen1997] in amino acid sequences from different organisms: Gram-positive prokaryotes, Gram-negative prokaryotes, and eukaryotes. 

**get_signalp** function from motifFUN package provides an interface for two versions of signalp-3.0 and signalp-5.0.

get_signalp() requires to provide version, signalp path, input file, organism type argument: 

For signalp-3.0

+ ‘euk’ - for eukaryotes;
+ ‘gram+’ - for gram-positive bacteria;
+ ‘gram-’ - for gram-negative bacteria.

For signalp-5.0

+ -org euk - for eukaryotes(default "euk")
+ -org arch - for Archaea
+ -org gram+ - for Gram-positive:
+ -org gram- - for Gram-negative

If the input fasta files containing more than 300 sequences, get_signalp will automatically switch to parallel mode. It will split the input into smaller chunks and run prediction as a massive parallel process using a specified number of CPUs. If user did not provide signalp path then by default function takes the latest version of signalp-5 path. For signalp-3 version user needs to provide signalp-3 path.

```{r}
input_file <- "testfile_orf_descard_REGEX.fasta"
signalp5 <- get_signalp(signalp.path = NULL, signalp.version = 5, input_file = input_file, org.type = "-org euk")
head(signalp5)
```

This function generates signalp resulted text file and returns summary dataframe.

**4.7.targetp**

TargetP 1.1 is a software tool to predicts the subcellular location of eukaryotic proteins.TargetP provides a potential cleavage site for sequences predicted to contain a cTP, mTP or SP. The get_targetp function requires to provide a targetp path, organism group, input file.

+ -P use plant networks (mandatory)
+ -N use non-plant (mandatory)

If user did not provide targetp path then by default function will take targetp path.

```{r}
targetp <- get_targetp(targetp.path = NULL, organismgroup = "-P", file = input_file)
head(targetp)
```

+ This function generates text file having Targetp result
+ This function returns data frame contains NAME and its LOCALIZATION columns.

LOCALIZATION column:

+ C Chloroplast, i.e. the sequence contains cTP, a chloroplast transit peptide;
+ M Mitochondrion, i.e. the sequence contains mTP, a mitochondrial targeting peptide;
+ S Secretory pathway, i.e. the sequence contains SP, a signal peptide;
+ _ Any other location;
+ &#42; "don't know"; indicates that cutoff restrictions were set (see instructions) and the winning network output score was below the requested cutoff for that category.

**4.8.tmhmm**

TMHMM predicts transmembrane α-helices and identifies integral membrane proteins based on HMMs [@Krogh2001a].
The get_tmhmm function requires to provides a tmhmm path, input fasta file, If user did not provide tmhmm path then by default function will take tmhmm path.

```{r}
tmhmm <- get_tmhmm(tmhmm.path = NULL, file = "testfile_orf_descard_REGEX.fasta")
head(tmhmm)
```

+ This function generates text file having tmhmm result
+ This function returns data frame contains the ID and NUMBER_OF_PREDICTED_TMHS columns

**4.9.motif summarys and motif sequences**

The user can extract all of the non-redundant sequences and a summary table with the information about the motifs using the summary_motifs function. This function uses the results from either hmm.search or pattern. search functions to generate a table that includes the name of the candidate motif sequence, the number of motifs of interest per sequence and its location within the sequence. 

```{r}
motif_summary <- summary_motifs(hmm.result = motif_candidates, reg.pat= pattern, signalp_version = 5, input_file = input_file)
head(motif_summary$consensus.sequences, n = 2)
head(motif_summary$motif.table, n=5)
```

**4.10.Visuvalizing HMM profile**

To determine if the HMM profile includes the motifs of interest, MotifFUN have The function hmm.plot reads the HMM profile (obtained from the hmm.search step) and uses ggplot2 to create a point plot. The plot will illustrate the bits (amino acid scores) of each amino acid used to construct the HMM profile according to its position in the HMM profile.

```{r, fig.align="center", fig.width=6, fig.height=4, fig.cap="Figure: Sequence logo plot of the motif candidates "}
hmm.plot(hmm_data = motif_candidates$HMM_Table)
```
### Reporting bugs

Please raise an issue (preferred option) or email <bairupravi999@gmail.com> about bugs and strange things.

