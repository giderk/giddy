```bashscript

READ MAPPING for the LyzMcre transgene (cre recombinase):

OVERVIEW
======
To validate that the Cpt2KO genotype is indeed set-up correctly we can look at read mapping and expression of the LyzM gene (synonomous to Lys2). This is straight forward to achieve from the standard "ExonCoverage" analysis and can also be assessed from the featureCounts output, as the Lys2 gene is present in the mouse reference assembly. However, the transgene 'cre recombinase' must be checked separately, as it is a foreign gene from enterobacteria phage 1, that is inserted after LyzM to enable conditional gene knockouts in myeloid tissue (the mouse lineage with this transgene is referred to as LyzMcre, and it can be crossed with floxed mice lineages to achieve tissue specific gene deletion). To check for expression of cre recombinase, we need to create a new dB (reference) file to map to, or append the cre gene to an existing dB file. Rather than mapping the entire RNAseq dataset to the new/altered dB, we can extract the portion of reads that failed to map to the mouse assembly and attempt to remap them with the new dB. This will drastically reduce processing time, as only a fraction of the reads did not map.

ANTICIPATED WORKFLOW:

Extract unmapped reads from .bam and/or .sam files (samtools view) >>
Convert .bam/.sam to fasta or fastq, or whatever file type that BLASTn accepts >>
Download and install BLAST+ and magicBLAST software >>
Download cre recombinase fasta file >>
Create BLAST dB from cre fasta file with NCBI tool makeblastDb >>
map reads to cre gene per sample with magicBLAST tool (which accepts fastq or fasta)

RESOURCES:
============
samtools - /gscmnt/gc2732/mitrevalab/TOOLS/SOFTWARE/Samtools1.3/SAMTOOLS_1.3/samtools-1.3/samtools

		See these pages for help configuring the command:
		http://www.htslib.org/doc/samtools.html # commands and parameters
		http://broadinstitute.github.io/picard/explain-flags.html # use FLAGs to subset mapped and unmapped data
		

BLAST software - https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download

Cre recombinase sequence - https://www.ncbi.nlm.nih.gov/nuccore/NC_005856.1?report=fasta&from=436&to=1467



RESULTS
=======
Used magicblast to output .tsv files (can read easily with excel) that indicate all of the unmapped reads that are mapping to somewhere on the cre Recombinase gene. The mapping looks really good, 100% identity scores almost every time. See files here....

</Users/gideon/Box/RNAseq/creRecombinase/Cpt2_fastq>


CONCLUSION
===========

I find some strangeness. There is expression of cre in most of the samples, ranging from zero up to 100 reads. Not sure if this represents a problem that could result from bad crossing, bad genotyping of mice before the experiment, contamination during library prep, etc .This magicblast method to explore cre expression is not a standard part of any RNAseq analysis workflow, so it may be good to redo this with another dataset (in which we have complete confidence in the KO genotype) to see if the results look cleaner, and, if not, maybe there is something flawed about what I am doing here.


DETAILS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


20192507
========
 - creating cre recombinase dB on local environment (since magicBLAST software is not on server) 
 
 cd /Users/gideon/Box/RNAseq/creRecombinase
 
makeblastdb -in cre.fa -input_type fasta -parse_seqids -title creRecombinase -logfile cre_makeblsdB_log -dbtype nucl



20192507
========
  - extracting unmapped RNAseq reads (on server)
  
find . -type f -name *.sorted.bam | xargs -I % bsub -q research-hpc -M 61000000 -R 'rusage[mem=63000]' -n 8 -a 'docker(jmartin7/ubuntu:mitrevalab_env2)' -oo %.unmapped.bsub.log '/gscmnt/gc2732/mitrevalab/TOOLS/SOFTWARE/Samtools1.3/SAMTOOLS_1.3/samtools-1.3/samtools view -b -f 4 -o %.unmapped.bam %'

bjobs | cut -b -7 | xargs echo -n
JOBID 2045615 2045606 2045616 2045603 2045617 2045619 2045601 2045609 2045602 2045599 2045611 2045608 2045600 2045604 2045605 2045613 2045598 2045610 2045621 2045614 2045612 2045620 2045607 2045618


  - converting .bam back to fastq format
  
find . -type f -name '*.unmapped.bam' | xargs -I % bsub -q research-hpc -M 61000000 -R 'rusage[mem=63000]' -n 8 -a 'docker(jmartin7/ubuntu:mitrevalab_env2)' -o %.fastq '/gscmnt/gc2732/mitrevalab/TOOLS/SOFTWARE/Samtools1.3/SAMTOOLS_1.3/samtools-1.3/samtools fastq %'

NOTE: for some reason bsub does not handle this command like others. For example -oo option for a bsub.log output does not produce the stder file. If you make the option -o then it will produce the actual .fastq.  The command 'samtools fastq' itself has no standard output option such as -o, instead they say to us the '>' to write to a file, and there are lots of options for splitting up or merging .bam files. So the above command works fine, but no log files are produced. LATER learned, the stder information is included at the beginning, end, and a little bit within the body of the fastq file. This information needs to be removed for downstream processing of the fastq file.

20192507
========
  - sending just the .fastq files to local environment
  
rsync -zt --progress gideon.e@virtual-workstation1.gsc.wustl.edu:'/gscmnt/gc2665/philips/Cpt2_seq_data/ExonCoverage/for_IGV/*.fastq' /Users/gideon/Box/RNAseq/creRecombinase/fastq/
  

20192507
========
  - testing magicBLAST with a fastq
  
magicblast -query reads.fastq -db my_reference -infmt fastq -outfmt tabular -no_unaligned -out 

NOTE: The resulting fastq is almost ready for magicblast but needs to be cleaned in three ways
		1) remove all header lines associated with the bsub run
		2) somewhere near the end of the file, in a different place each time, an insertion of [M::bam2fq_mainloop] processed XXX reads] occurs. This insertion breaks up the fastq formating and causes a parsing error. Remove just this insertion and the read will run perfectly fine
		
		EXAMPLE OF INSERTION:
		
		DDDCDIIIIHHHIHIIIIHIDHGFHIIIIIIIIIIIHHIIIIIIHIHIIH
		@SN1063:704:H5FJHBCX2:1:2216:9907:22278
		TCCTAGGTTTTTTGTT[M::bam2fq_mainloop] processed 549952 reads
		ATTCCAGATGAATTTGCAAATTGCTCCTTCTAAT
		+
		
		EXAMPLE OF APPROPRIATE CORRECTION:
		
		DDDCDIIIIHHHIHIIIIHIDHGFHIIIIIIIIIIIHHIIIIIIHIHIIH
		@SN1063:704:H5FJHBCX2:1:2216:9907:22278
		TCCTAGGTTTTTTGTTATTCCAGATGAATTTGCAAATTGCTCCTTCTAAT
		+
		
		3) remove the last lines of the file pertaining to bsub run reporting
		


END END END
```