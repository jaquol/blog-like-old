# 03. Choice of Trimmomatic parameters

[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) is a tool for trimming low-quality bases and technical sequences (e.g. Illumina adapters) from sequencing reads. When I started using Trimmomatic routinely I was unsure about some choices:

### What trimming mode should be used to remove low-quality positions?

According to [Trimmomatic's paper](http://bioinformatics.oxfordjournals.org/content/30/15/2114.full), the Maximum Information mode seems to be the recommended option.


### What technical sequences need to be removed?

A common practice is removing Illumina sequencing adapters as these are amongst the most common contaminants in the reads. The trimmomatic package comes with a directory that contains the sequence of the adapters used by the different Illumina machines (in FASTA format so that they are ready to use by Trimmomatic). I noticed that for most machines there are two similar FASTA files. For instance:

- TruSeq3-PE.fa
- TruSeq3-PE-2.fa (which contains the same FASTA sequences as in the first file but with another header and reverse format)

I got the following rule of thumb from Trimmomatic's developer [Anthony Bolger](http://www.usadellab.org/cms/index.php?page=BolgerAnthony):

> "For most libraries, TruSeq3-PE is sufficient. This will detect adapter read-throughs only (which are by far the most common kind of contamination), and it should have a low false positive rate. TruSeq3-PE-2 will find less common sources of contamination, but it is necessary only for libraries created from degraded kits or via
non-official methods. As a general rule, i try TruSeq3-PE first, and if there is still non-trivial amounts of contamination (e.g. visible in FastQC), i switch to TruSeq3-PE-2." 


### What values should be specified for the different parametes?

I normally use the parameters used in the Trimmomatic paper for paired-end 100-bp reads, as provided in the [supplementary data](http://bioinformatics.oxfordjournals.org/content/suppl/2014/03/29/btu170.DC1/commands.txt):

- seedMismatches=2
- palindromeClipThreshold=30
- simpleClipThreshold=12
- minAdapterLength=1
- keepBothReads=true
- minQual=3
- targetLength=40
- strictness=0.5 (the authors used a range of values for this parameter; since it can range from 0 to 1, I use an intermediate value of 0.5)
- minLength=36
