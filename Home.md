# RNA-Seq Analysis Workshop

Welcome to the Trinity RNA-Seq analysis workshop! Here we will cover RNA-Seq analysis using genome-guided and/or genome-free methods:

* For the genome-free methods, we'll be using Trinity (of course!), and cover de novo transcript assembly followed by transcript quantitation and differential expression analysis.

* For the genome-guided RNA-Seq workshop, we'll use the highly popular Tuxedo toolkit.

An overview of these methodologies is presented in the [accompanying workshop slides](https://github.com/trinityrnaseq/RNASeq_Trinity_Tuxedo_Workshop/blob/master/docs/rnaseq_workshop_slides.pdf). 

The workshop involves hands-on learning in applying the above computational methods to sample RNA-Seq data.  The RNA-Seq data involve paired-end 76 base strand-specific Illumina RNA-Seq reads corresponding  to Schizosaccharoymyces pombe (fission yeast) being grown under 4 different conditions: logarithmic growth (Sp_log), plateau phase (Sp_plat), heat shock (Sp_hs), and diauxic shift (Sp_ds). These data were generated as part of our previous publications ( [Rhind et al.]( http://www.ncbi.nlm.nih.gov/pubmed/21511999) and [Grabherr et al.]( http://www.ncbi.nlm.nih.gov/pubmed/21572440) ).

The raw data and all the software required to complete the workshop are built into a VirtualBox image as [Trinity2015.ova](ftp://ftp.broadinstitute.org/pub/Trinity/RNASEQ_WORKSHOP/Trinity2015.ova).  To install it, you'll first need [VirtualBox](https://www.virtualbox.org/wiki/Downloads) installed - which is free and easy to do.  Then simply follow the [instructions to import and run the workshop VM](Import-and-run-workshop-VM).

The Trinity and Tuxedo workshops leverage the same input data, but are otherwise independent analysis trajectories.  Simply choose which one you'd like to run through and begin:

* [Trinity de novo transcriptome assembly workshop](Trinity-De-novo-Transcriptome-Assembly-Workshop)
* [Tuxedo genome-guided transcriptome assembly workshop](Tuxedo-Genome-Guided-Transcriptome-Assembly-Workshop)

