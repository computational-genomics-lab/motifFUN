# motifFUN
An R package for searching user defined motif with functional annotation features
The **motifFUN** package is an R package designed to search for the motifs of interest using user-defined regular expression searches and hidden Markov models (HMM) along with a prediction of some functional criteria such as signal peptides, subcellular location of Eukaryotes and transmembrane helices.

Several command line tools and web-interfaces exist to perform predictions of individual motifs and domains (SignalP, TargetP, TMHMM) however the interface that combines the outputs in a single flexible workflow is lacking, So that developed a motifFUN package to fulfill that gap.

<font size="4">**1.OUTLINE**</font>

The **motifFUN** package provides a platform to search motifs and build automated multi-step secretome prediction pipelines that can be applied to large protein. The features of this package is described below.

         
**ORF extraction** 

In this package get_orf function extract six-frame translation of all ORF(open reading frame) in the genome sequence, To perform this we should recommend installing EMBOSS getorf.

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

+ The motifFUN package starts with the prediction of short **signal peptides(signalP)** at the N-terminal end of a protein.
+ Predicts the sub cellular location of eukaryotes (**TargetP**), by using get_targetp function.
+ Predicts **trans membrane helices(TMHMM)** by using get_tmhmm function.
