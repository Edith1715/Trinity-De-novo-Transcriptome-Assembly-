# Genome-based RNA-Seq analysis using the Tuxedo package

The following details the steps involved in:

The following details the steps involved in:
*	Aligning RNA-Seq reads to a genome using Tophat
*	Assembling transcript structures from read alignments using Cufflinks
*	Visualizing reads and transcript structures using IGV
*	Performing differential expression analysis using Cuffdiff
*	Expression analysis using CummeRbund

All required software and data are provided pre-installed on a VirtualBox image. Be sure to [run the workshop VM and open a terminal](Import-and-run-workshop-VM).

After installing the VM, be sure to quickly update the contents of the rnaseq_workshop_data directory by:

        %   cd RNASeq_Trinity_Tuxedo_Workshop

        %   git pull

This way, you’ll have the latest content, including any recent bugfixes.

### Data Content:

This demo uses RNA-Seq data corresponding to Schizosaccharomyces pombe (fission yeast), involving paired-end 76 base strand-specific RNA-Seq reads corresponding  to four samples:  Sp_log (logarithmic growth), Sp_plat (plateau phase), Sp_hs (heat shock), and Sp_ds (diauxic shift). 

There are 'left.fq' and 'right.fq' FASTQ formatted Illlumina read files for each of the four samples.  All RNA-Seq data sets can be found in the **RNASEQ_data/** subdirectory.

Also included is a 'genome.fa' file corresponding to a genome sequence, and annotations for reference genes ('genes.bed' or 'genes.gff3'). These resources can be found in the **GENOME_data/** subdirectory.  

>Note, although the genes, annotations, and reads represent genuine sequence data, they were artificially selected and organized for use in this tutorial, so as to provide varied levels of expression in a very small data set, which could be processed and analyzed within an approximately one hour time session and with minimal computing resources.


### Automated and Interactive Execution of Activities 

To avoid having to cut/paste the numerous commands shown below into a unix terminal, the VM includes a script ‘runTrinityDemo.pl’ that enables you to run each of the steps interactively.  To begin, simply run:

      %  ./runTuxedoDemo.pl 

Note, by default and for convenience, the demo will show you the commands that are to be executed. This way, you don’t need to type them in yourself.

The protocol followed below is that described in 

[Trapnell C, Roberts A, Goff L, Pertea G, Kim D, Kelley DR, Pimentel H, Salzberg SL, Rinn JL, Pachter L.  Differential gene and transcript expression analysis of RNA-seq experiments with TopHat and Cufflinks. Nat Protoc. 2012 Mar 1;7(3):562-78. doi: 10.1038/nprot.2012.016.](http://www.nature.com/nprot/journal/v7/n3/full/nprot.2012.016.html)

and follows this general framework as illustrated in the above publication:
<img src="https://raw.githubusercontent.com/wiki/trinityrnaseq/RNASeq_Trinity_Tuxedo_Workshop/images/TuxedoWorkshop/nprot.2012.016-F2.jpg" width=500 />
