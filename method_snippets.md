# Method snippets


## Expression Microarrays

### Illumina arrays

Arrays were processed using the 'oligo' [Carvalho B. S., and Irizarry, R. A. (2010). A Framework for Oligonucleotide Microarray Preprocessing. Bioinformatics, 26(19):2363-7.] BioConductor package, quality-controlled with arrayQualityMetrics [Kauffmann, A., Gentleman, R.,, Huber, W. (2009) arrayQualityMetrics--a bioconductor package for quality assessment of microarray data. Bioinformatics, 25(3):415-6.]
 and RMA [Rafael. A. Irizarry, Benjamin M. Bolstad, Francois Collin, Leslie M. Cope, Bridget Hobbs and Terence P. Speed (2003), Summaries of Affymetrix GeneChip probe level data Nucleic Acids Research 31(4):e15]
 normalized. Differentially expressed genes (as defined by a minimun 2 fold change in expression and  a  Benjamini-Hochberg FDR adjusted pvalue of less than 0.1) were identified using limma  [Smyth, GK (2005). Limma: linear models for microarray data. In:  'Bioinformatics and Computational Biology Solutions using R and  Bioconductor'. R. Gentleman, V. Carey, S. Dudoit, R. Irizarry, W. Huber  (eds), Springer, New York, pages 397-420.].”

### Affymetrix arrays

Arrays were processed using the 'oligo' [Carvalho B. S., and Irizarry, R. A. (2010). A Framework for Oligonucleotide Microarray Preprocessing. Bioinformatics, 26(19):2363-7.] BioConductor package, quality-controlled with arrayQualityMetrics [Kauffmann, A., Gentleman, R.,, Huber, W. (2009) arrayQualityMetrics--a bioconductor package for quality assessment of microarray data. Bioinformatics, 25(3):415-6.]
 and RMA [Rafael. A. Irizarry, Benjamin M. Bolstad, Francois Collin, Leslie M. Cope, Bridget Hobbs and Terence P. Speed (2003), Summaries of Affymetrix GeneChip probe level data Nucleic Acids Research 31(4):e15]
 normalized. Differentially expressed genes (as defined by a minimun 2 fold change in expression and  a  Benjamini-Hochberg FDR adjusted pvalue of less than 0.1) were identified using limma  [Smyth, GK (2005). Limma: linear models for microarray data. In:  'Bioinformatics and Computational Biology Solutions using R and  Bioconductor'. R. Gentleman, V. Carey, S. Dudoit, R. Irizarry, W. Huber  (eds), Springer, New York, pages 397-420.].”



## Next-gen sequencing

## Variant calling methods

Sequencing data will be obtained, tested for quality and trimmed as before. Filtered reads will be aligned with Novoalign (http://www.novocraft.com/main/index.php) before running a standard GATK best practice workflow (http://gatkforums.broadinstitute.org/discussion/1186/best-practice-variant-detection-with-the-gatk-v4-for-release-2-0) including deduplication, base recalibration and realignment. Variants will be called using an ensembl approach with the GATK UnifiedGenotyper, GATK Haplotype caller and FreeBayes on [[all samples simultaneously | samples grouped by family]]. Combining results in this way from multiple callers can increase overall call confidence (see http://bcbio.wordpress.com/2013/02/06/an-automated-ensemble-method-for-combining-and-evaluating-genomic-variants-from-multiple-callers/). We annotate variant calls using SnpEff (http://snpeff.sourceforge.net/) to predict variant effects, and gemini (http://snpeff.sourceforge.net/) to produce a readily queried database with annotations from external sources such as dbSNP, ENCODE, KEGG and ClinVar. 



Moving to whole genome analysis, from our previous experience with exome and targeted population sequencing, required reworking our analysis approaches. The primary bottleneck is filesystem IO, which requires new scaling approaches. To provide a practical example, writing a single 100Gb aligned genome in BAM format requires approximately 1 day on a reasonably fast filesystem. When writing 10 files concurrently, shared filesystems like NFS slow down considerably
and the time can increase to longer than a week.

Due to these challenges we improved our ability to parallelize these samples by:

- Re-writing our pipeline to avoid file writing steps as much as possible, and running alignment preparation and variant calling steps in smaller chunks of the genome.

- Improving our shared filesystem network throughput by increasing connectivity between cluster nodes.

- Using distributed filesystems to help distribute network and IO traffic. We're currently exploring GlusterFS (http://www.gluster.org/).

The second significant challenge was machine memory, which becomes limiting during variant calling steps with multiple samples. The major challenge is that traditional cluster management systems manage compute core usage but not memory, so computing on shared machines can become problematic when other memory-intensive jobs get scheduled together. There are both technical solutions, like grabbing entire machines for processing to avoid other scheduled jobs, and algorithmic solutions, like Broad's ReduceReads https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_compression_reducereads_ReduceReads.html).



## Cancer variant calling

Sequencing reads will be assed for quality [http://www.bioinformatics.babraham.ac.uk/projects/fastqc/], trimmed and filtered [https://code.google.com/p/cutadapt/], and aligned against the reference genome using Novoalign [http://www.novocraft.com/] before being processed following the GATK's best practice workflow (base recalibration, realignment, and variant calling using the UnifiedGenotyper). In addition, we will apply VarScan2 [Koboldt DC, Zhang Q, Larson DE, Shen D, McLellan MD, Lin L, Miller CA, Mardis ER, Ding L, & Wilson RK (2012). VarScan 2: Somatic mutation and copy number alteration discovery in cancer by exome sequencing. Genome Research PMID: 22300766] to determine somatic status of detected variants. All variants will be annotated using snpEff [ "A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff: SNPs in the genome of Drosophila melanogaster strain w1118; iso-2; iso-3.", Cingolani P, Platts A, Wang le L, Coon M, Nguyen T, Wang L, Land SJ, Lu X, Ruden DM. Fly (Austin). 2012 Apr-Jun;6(2):80-92. PMID: 22728672] before downstream analysis.


## RNA-Seq

HiSeq Illumina sequencing will be performed on our behalf by ###. All samples will be indexed so that pools can be run across ### lanes of an 8-laned Illumina flow cell, providing an estimated 20-30 million single-end reads per sample.

All samples will be processed using an RNA-seq pipeline implemented in the bcbio-nextgen project (https://bcbio-nextgen.readthedocs.org/en/latest/). Raw reads will be examined for quality issues using FastQC (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to ensure library generation and sequencing are suitable for further analysis. Possible ribosomal contamination will be assessed using CollectRNASeqStats (http://picard.sourceforge.net) and cross-sample contamination will be evaluated with ContEST [Cibulskis, Kristian, Aaron McKenna, Tim Fennell, Eric Banks, Mark DePristo, and Gad Getz. “ContEst: Estimating Cross-Contamination of Human Samples in Next-Generation Sequencing Data..” Bioinformatics (Oxford, England) 27, no. 18: 2601–2602. doi:10.1093/bioinformatics/btr446.].

Adapter sequences, other contaminant sequences such as polyA tails and low quality sequences with PHRED quality scores less than five will be trimmed from reads using cutadapt [Martin, M. 2011. “Cutadapt Removes Adapter Sequences From High-Throughput Sequencing Reads.” EMBnet Journal.]. Trimmed reads will be aligned to UCSC build ### of the ## # genome (###), augmented with transcript information from Ensembl release ### using STAR [Dobin, Alexander, Carrie A Davis, Felix Schlesinger, Jorg Drenkow, Chris Zaleski, Sonali Jha, Philippe Batut, Mark Chaisson, and Thomas R Gingeras. 2013. “STAR: Ultrafast Universal RNA-Seq Aligner..” Bioinformatics (Oxford, England) 29 (1). Oxford University Press: 15–21. doi:10.1093/bioinformatics/bts635.].

Alignments will be checked for evenness of coverage, rRNA content, genomic context of alignments (for example, alignments in known transcripts and introns), complexity and other quality checks using a combination of FastQC, RNA-SeQC [DeLuca, David S, Joshua Z Levin, Andrey Sivachenko, Timothy Fennell, Marc-Danie Nazaire, Chris Williams, Michael Reich, Wendy Winckler, and Gad Getz. 2012. “RNA-SeQC: RNA- Seq Metrics for Quality Control and Process Optimization..” Bioinformatics (Oxford, England) 28 (11). Oxford University Press: 1530–32. doi:10.1093/bioinformatics/bts196.] and custom tools. 

Counts of reads aligning to known genes and isoforms will be generated by a combination of featureCounts [Liao, Yang, Gordon K Smyth, and Wei Shi. 2014. “featureCounts: an Efficient General Purpose Program for Assigning Sequence Reads to Genomic Features..” Bioinformatics (Oxford, England) 30 (7). Oxford University Press: 923–30. doi:10.1093/bioinformatics/btt656.], eXpress [Roberts A, Pachter L. (2013). Streaming fragment assignment for real-time analysis of sequencing experiments. Nature Methods, 10(1), 71–3. PMCID: PMC3880119] and DEXSeq [Anders, Simon, Alejandro Reyes, and Wolfgang Huber. 2012. “Detecting Differential Usage of Exons From RNA-Seq Data..” Genome Research 22 (10): 2008–17. doi:10.1101/gr. 133744.111.].

Depending on sequence quality and total counts, novel transcripts will be identified via reference-guided assembly with Cufflinks [Trapnell, Cole, Brian A Williams, Geo Pertea, Ali Mortazavi, Gordon Kwan, Marijke J van Baren, Steven L Salzberg, Barbara J Wold, and Lior Pachter. 2010. “Transcript Assembly and Quantification by RNA-Seq Reveals Unannotated Transcripts and Isoform Switching During Cell Differentiation..” Nature Biotechnology 28 (5). Nature Publishing Group: 511–15. doi: 10.1038/nbt.1621.], with novel transcripts filtered for coding potential agreement [Wang, Liguo, Hyun Jung Park, Surendra Dasari, Shengqin Wang, Jean-Pierre Kocher, and Wei Li. 2013. “CPAT: Coding-Potential Assessment Tool Using an Alignment-Free Logistic Regression Model..” Nucleic Acids Research 41 (6). Oxford University Press: e74–e74. doi: 10.1093/nar/gkt006.] with known genes to reduce false positive assemblies. 

Differential expression at the gene level will be called with DESeq2 [Love, Michael I, Wolfgang Huber, and Simon Anders. 2014. “Moderated Estimation of Fold Change and Dispersion for RNA-Seq Data with DESeq2..” Genome Biology 15 (12): 550. doi:10.1186/PREACCEPT-8897612761307401.], which has been shown to be a robust, conservative differential expression caller. 

Isoform and exon- level calls are much more prone to false positives, and so an ensemble method combining calls that agree from both DEXSeq and EBSeq [Leng, Ning, John A Dawson, James A Thomson, Victor Ruotti, Anna I Rissman, Bart M G Smits, Jill D Haag, Michael N Gould, Ron M Stewart, and Christina Kendziorski. 2013. “EBSeq: an Empirical Bayes Hierarchical Model for Inference in RNA-Seq Experiments..” Bioinformatics (Oxford, England) 29 (8). Oxford University Press: 1035–43. doi:10.1093/bioinformatics/ btt087.] will be used to call differential isoforms.

RNAseq results will be validated by qRT-PCR on a subset of genes. Lists of differentially expressed genes will be examined for gene ontology (GO) term enrichment with g:profiler[Reimand, J., Kull, M., Peterson, H., Hansen, J. & Vilo, J. (2007). g:Profiler--a web-based toolset for functional profiling of gene lists from large-scale experiments. Nucleic Acids Research 35, W193–W200 . doi:10.1093/nar/gkm226]. Functional redundancy in enriched GO terms will be reduced with Revigo [Supek, F., Bošnjak, M., Škunca, N. & Šmuc, T. REVIGO Summarizes and Visualizes Long Lists of Gene Ontology Terms. (2011). PLoS ONE 6, e21800 . doi:10.1371/journal.pone.0021800], which uses GO term semantic similarity and enrichment p-values to group highly similar GO terms under representative terms. In addition, a cut-off-free gene set enrichment analysis will be performed using SeqGSEA [Wang, X., & Cairns, M. J. (2014). SeqGSEA: a Bioconductor package for gene set enrichment analysis of RNA-Seq data integrating differential expression and splicing. Bioinformatics (Oxford, England), btu090. doi:10.1093/bioinformatics/btu090]. SeqGSEA accounts for biological variability and minimizes  biases towards long or highly expressed genes. 


Finally, sample gene expression results will be compared against the Broad’s Connectity Map using gCMAP [http://www.bioconductor.org/packages/devel/bioc/html/gCMAP.html] to identify drug perturbations resulting in highly correlated (or anti-correlated) gene expression changes for subequent validation. We will also explore expanding the drug perturbation database with the LINCS L1000 data [http://lincscloud.org/l1000-data/], subsetting differentially expressed gene lists accordingly.


## Exome Seq (normal or cancer)

To analyze small variants (SNPs) and indels in exome sequencing samples, we
employ a suite of existing software packages combined with methodology for
producing a set of combined ensemble calls from multiple approaches. 

Our exome sample sequencing approach detects small SNPs and indels
using two best practice approaches, via an automated parallelized pipeline
(https://github.com/chapmanb/bcbio-nextgen). We align reads using Novoalign
(http://www.novocraft.com/main/index.php) and post process with de-duplication, 
base quality score recalibration and realignment around indels. We call variants
from these prepared reads using two pipelines: the Broad's GATK tools
(http://gatkforums.broadinstitute.org/discussion/1186/best-practice-variant-detection-with-the-gatk-v4-for-release-2-0)
call variants with the UnifiedGenotyper, following by filtration with the Variant
Quality Score Recalibrator; our second pipeline utilizes FreeBayes from the
Marth Lab (http://gkno.me/pipelines.html) followed by filtration with machine
learning approaches using self-organizing maps. By combining calls from multiple
approaches, we find we can achieve improved sensitivity and specificity compared
to individual methods (http://bcbio.wordpress.com/2013/02/06/an-automated-ensemble-method-for-combining-and-evaluating-genomic-variants-from-multiple-callers/).

We similarly employ two calling methods for variant detection in tumor-normal pairs:
muTect (http://www.broadinstitute.org/cancer/cga/mutect) from the Broad
Institute and VarScan2 (http://varscan.sourceforge.net/) from the Genome
Institute at Washington University. Both provide variant calling and filtering
on matched samples, removing common sequencing artifacts.

To provide biological context for downstream analysis, we annotate variants from
both approaches with predicted effects and associated public data. snpEff 
(http://snpeff.sourceforge.net/) predicts coding effects and Gemini 
(https://github.com/arq5x/gemini) stores variants associated with public
datasets (dbSNP, ENCODE, ClinVar) for query and prioritization.

For copy number detection, we use two freely available packages designed to
identify CNVs from exome data: CoNIFER (http://conifer.sourceforge.net/) and
ExomeCNV (https://secure.genome.ucla.edu/index.php/ExomeCNV_User_Guide). These
provide detection of copy number alterations in single sample and tumor-normal
data, respectively.

Generally, our approach is to utilize multiple best-practice approaches to
identify SNP, indel and copy number variations and use flexible approaches to
combine these into a finalized callset. Our automated installation and pipeline
approaches help ameliorate the difficulties of setting up and running multiple
packages. We make this infrastructure trade off because we value the improvement
in detection sensitivity and specificity, thus passing our most confident variants
into downstream analysis.


[Koboldt DC, Zhang Q, Larson DE, Shen D, McLellan MD, Lin L, Miller CA, Mardis ER, Ding L, & Wilson RK (2012). VarScan 2: Somatic mutation and copy number alteration discovery in cancer by exome sequencing. Genome Research PMID: 22300766]

[ "A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff: SNPs in the genome of Drosophila melanogaster strain w1118; iso-2; iso-3.", Cingolani P, Platts A, Wang le L, Coon M, Nguyen T, Wang L, Land SJ, Lu X, Ruden DM. Fly (Austin). 2012 Apr-Jun;6(2):80-92. PMID: 22728672]



## mRNA

RNA quality will be performed using the Agilent Bioanalyzer (Agilent, Santa Clara, CA). Transcriptional profiling using the ### array will be performed at ###. All data analysis will be performed using the BioConductor framework [1]. Array quality will be assessed using the arrayQualityMetrics package [2], batch adjusted [3], normalized and summarized using Robust Multiarray Averaging (RMA) [4] and tested for genes with significant differential expression using limma [5].

[1]: Reimers, Mark, and Vincent J Carey. “Bioconductor: an Open Source Framework for Bioinformatics and Computational Biology.” Methods in Enzymology 411: 119–134. doi:10.1016/S0076-6879(06)11008-3.
[2]: Kauffmann, Audrey, Robert Gentleman, and Wolfgang Huber. “arrayQualityMetrics--a Bioconductor Package for Quality Assessment of Microarray Data.” Bioinformatics (Oxford, England) 25, no. 3: 415–416. doi:10.1093/bioinformatics/btn647.
[3]: Leek, Jeffrey T, Robert B Scharpf, Héctor Corrada Bravo, David Simcha, Benjamin Langmead, W Evan Johnson, Donald Geman, Keith Baggerly, and Rafael A Irizarry. “Tackling the Widespread and Critical Impact of Batch Effects in High-Throughput Data.” Nature Reviews Genetics 11, no. 10: 733–739. doi:10.1038/nrg2825.
[4]: Wilson, Claire L, and Crispin J Miller. “Simpleaffy: a BioConductor Package for Affymetrix Quality Control and Data Analysis..” Bioinformatics (Oxford, England) 21, no. 18: 3683–3685. doi:10.1093/bioinformatics/bti605.
[5]: Smyth, G. K. (2005). Limma: linear models for microarray data. In: Bioinformatics and Computational Biology Solutions using R and Bioconductor, R. Gentleman, V. Carey, S. Dudoit, R. Irizarry, W. Huber (eds.), Springer, New York, pages 397-420.


## miRNA

Aligned read data will then be postprocessed with the miRDeep2 [Friedländer, Marc R, Sebastian D Mackowiak, Na Li, Wei Chen, and Nikolaus Rajewsky. “miRDeep2 Accurately Identifies Known and Hundreds of Novel microRNA Genes in Seven Animal Clades..” Nucleic Acids Research 40, no. 1 (January 1, 2012): 37–52.], an algorithm that assesses the fit of sequenced RNAs to a biological model of miRNA generation and correct folding. We expect this step to significantly reduce the false-positive rate compared to a simpler analysis of clustered read positions. 

While we expect miRNA to be the most abundant contributor to our sequence data we will also explore other sRNA as potential markers. Aligned reads with a minimum length of 18 nucleotides overlapping genomic loci annotated as snoRNA, tRNA or ncRNA genes in the UCSC Genome Browser [[Kuhn, Robert M, David Haussler, and W James Kent. “The UCSC Genome Browser and Associated Tools..” Briefings in Bioinformatics (August 20, 2012).] will be retained and clustered into distinct groups based on sequence overlap [[Quinlan, Aaron R, and Ira M Hall. “BEDTools: a Flexible Suite of Utilities for Comparing Genomic Features..” Bioinformatics (Oxford, England) 26, no. 6 (March 15, 2010): 841–842.]. Clusters with a minimum coverage of 10 reads in at least 20% of the sample population will be retained and assessed for differential expression 

To correct for underestimation of miRNA expression when using the whole genome background -- due to the higher ambiguity of mapped reads [REF#68] -- we will assess differentially expression of identified miRNA candidates and other sRNA independently from the first step, again using DESeq or DSS as described above.

## miRNA/mRNA integration

To comprehensively integrate the interaction of differentially expressed genes and miRNAs we will use miRTrail [ Laczny, Cedric, Petra Leidinger, Jan Haas, Nicole Ludwig, Christina Backes, Andreas Gerasch, Michael Kaufmann, et al. “miRTrail--a Comprehensive Webserver for Analyzing Gene and miRNA Patterns to Enhance the Understanding of Regulatory Mechanisms in Diseases..” BMC Bioinformatics 13 (2012): 36. doi:10.1186/1471-2105-13-36.] which integrates information on 20.000 genes and almost 1.000 miRNAs. A secondary screen will be performed using MMIA [Nam, Seungyoon, Meng Li, Kwangmin Choi, Curtis Balch, Sun Kim, and Kenneth P Nephew. “MicroRNA and mRNA Integrated Analysis (MMIA): a Web Tool for Examining Biological Functions of microRNA Expression.” Nucleic Acids Research 37, no. Web Server issue (July 1, 2009): W356–62. doi:10.1093/nar/gkp294.], an algorithm incorporating three miRNA prediction algorithms (TargetScan, PITA and PicTar) into the mRNA expression analysis. 


## ChIP-Seq

ChIP-seq data quality will be evaluated using FASTQC [6], and if required filtering and trimming of reads will be performed with Cutadapt [7]. High quality reads will be mapped to the current mouse genome build using Bowtie [8]. We will use MACS2 [9] to call peaks on unique reads only and assess peak quality using phantompeakqualtools, developed as part of the ModENCODE project [10]. Likely peak artifacts will be filtered out using the ModENCODE blacklist. To assess the reproducibility of peaks across replicates, we will run samples through the IDR pipeline [11,12]. For replicates which show high concordance from the IDR analysis, overlapping regions will be used for downstream analysis. Binding sites will be evaluated using the MEME [13] suite of tools for motif discovery and motif enrichment. Functional enrichment of binding site targets will be explored with GREAT [14].

[6]: http://www.bioinformatics.bbsrc.ac.uk/projects/fastqc/
[7]: https://code.google.com/p/cutadapt/
[8]: Langmead, B, C Trapnell, Pop, and S Salzberg. “Ultrafast and Memory-Efficient Alignment of Short DNA Sequences to the Human Genome.” Genome Biology 10, no. 3: R25. doi:10.1186/gb-2009-10-3-r25.
[9]: Zhang, Yong, Tao Liu, Clifford A Meyer, Jérôme Eeckhoute, David S Johnson, Bradley E Bernstein, Chad Nussbaum, et al. “Model-Based Analysis of ChIP-Seq (MACS).” Genome Biology 9, no. 9: R137. doi:10.1186/gb-2008-9-9-r137.
[10]: http://code.google.com/p/phantompeakqualtools/
[11]: https://sites.google.com/site/anshulkundaje/projects/idr
[12]: Landt S.G. et al. "ChIP-seq guidelines and practices of the ENCODE and modENCODE consortia" Genome Res. 2012 Sep; 22(9): 1813–1831
[13]: Bailey, Timothy L, Mikael Bodén, Fabian A Buske, Martin Frith, Charles E Grant, Luca Clementi, Jingyuan Ren, Wilfred W Li, and William S Noble. “MEME SUITE: Tools for Motif Discovery and Searching..” Nucleic Acids Research 37, no. Web Server issue: W202–8. doi:10.1093/nar/gkp335.
[14]: McLean, Cory Y, Dave Bristor, Michael Hiller, Shoa L Clarke, Bruce T Schaar, Craig B Lowe, Aaron M Wenger, and Gill Bejerano. “GREAT Improves Functional Interpretation of Cis-Regulatory Regions..” Nature Biotechnology 28, no. 5: 495–501. doi:10.1038/nbt.1630.

## TF

Likely transcription factors (TF) associated with differentially expressed genes will be identified using oPOSSUM-3[1], a web-accessible software system for identification of over-represented transcription factor binding sites (TFBS) and TFBS families in DNA sequences of co-expressed genes.

[1]: Kwon, Andrew T, David J Arenillas, Rebecca Worsley Hunt, and Wyeth W Wasserman. “oPOSSUM-3: Advanced Analysis of Regulatory Motif Over-Representation Across Genes or ChIP-Seq Datasets..” G3 (Bethesda, Md.) 2, no. 9 (September 2012): 987–1002. doi:10.1534/g3.112.003202.


## RRBS

HiSeq Illumina sequencing will be performed on our behalf by ###. All samples will be indexed so that ### RRBS pools can be run on a single lane of an 8-laned Illumina flow cell, providing an estimated ### million reads per sample.

Illumina sequence quality will be surveyed using FastQC [Andrews, S.: FastQC http://www.bioinformatics.babraham.ac.uk/projects/fastqc/] to ensure library generation and sequencing are suitable for further analysis. Particular attention will be paid to per base sequence quality, overrepresented sequences (an indicator of adapter contamination) and the per base cytosine content (an indicator of conversion efficiency). Using a bioinformatics pipeline constructed in bpipe [Sadedin, S. P., Pope, B., & Oshlack, A. (2012). Bpipe: a tool for running and managing bioinformatics pipelines. Bioinformatics (Oxford, England), 28(11), 1525–1526. doi:10.1093/bioinformatics/bts167], sequences will then be trimmed, aligned and their cytosine methylation status determined. Adapters will be removed and low quality bases (<25–30 Phred quality scores) adaptively trimmed from reads with Trim Galore [Krueger F: Trim Galore!. http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/]. In addition, bases containing a cytosine artificially introduced in the end-repair step during the MspI-RRBS library preparation will be removed. Sequence quality will be re-assessed with FastQC and trimming parameters adjusted if necessary. Trimmed reads will be aligned to the appropriate reference genome with Bismark [Krueger, F., & Andrews, S. R. (2011). Bismark: a flexible aligner and methylation caller for Bisulfite-Seq applications. Bioinformatics (Oxford, England), 27(11), 1571–1572. doi:10.1093/bioinformatics/btr167], a ‘three-letter’ bisulfite aligner based on the gapped aligner Bowtie2 [Langmead, B., & Salzberg, S. L. (2012). Fast gapped-read alignment with Bowtie 2. Nature Methods, 9(4), 357–359. doi:10.1038/nmeth.1923]. Aligned reads will be automatically trimmed with BSeQC to remove nucleotides subject to further technical biases that can result in inaccurate methylation estimation, such as higher bisulfite conversion failure at the 5’ end of reads [Lin, X., Lin, X., Sun, D., Sun, D., Rodriguez, B., Rodriguez, B., et al. (2013). BSeQC: quality control of bisulfite sequencing experiments. Bioinformatics, 29(24), 3227–3229. doi:10.1093/bioinformatics/btt548]. Cytosine methylation will be determined with Bis-SNP [Liu, Y., Siegmund, K. D., Laird, P. W., & Berman, B. P. (2012). Bis-SNP: Combined DNA methylation and SNP calling for Bisulfite-seq data. Genome Biology, 13(7), R61. doi:10.1186/gb–2012–13–7-r61
], a GATK [McKenna, A., Hanna, M., Banks, E., Sivachenko, A., Cibulskis, K., Kernytsky, A., et al. (2010). The Genome Analysis Toolkit: A MapReduce framework for analyzing next-generation DNA sequencing data. Genome Research, 20(9), 1297–1303. doi:10.1101/gr.107524.110] based framework capable of simultaneous genotyping and DNA methylation estimation; methylation at non-CpG cytosines will be used to assess conversion rates. DNA methylation will be validated by bisulfite PCR or pyrosequencing on a subset of CpGs.

