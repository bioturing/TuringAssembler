# TuringAssembler 
The high-throughput barcoding technology, which tags reads that originate from a single DNA molecule,
has been the foundation for generatingSynthetic Long Reads (SLR). 
The accuracy of short-read data combined with the long-range information encoded in the read clouds make protocols 
like **10X GemCode** or **TELL-Seq** the preferred solutions for haplotype genome phasing, and structural variant calling. 
Read clouds technologies also bring new hope to completely solve the long-standing genome and metagenome assembly problems,
at least for simple genome. 

We have developed **TuringAssembler**, the read-cloud assembler which can use only one read-cloud 
data to produce high-quality assembly for a bacterial genome. 

Using a TELL-Seq dataset of the model organism **E.coliK12 MG1655**, TuringAssembler assembled a gapless assembly with **NGA50 4632444**, mis-match rate **5.65/100kbp**, and indel rate **0.47/100kbp**
## Usage
```

=====================================================================================================
TuringAssembler is a program developed by BioTuring for doing genome assembly with read-cloud technology
Please contact support@bioturing.com if you need further support.
=====================================================================================================
TuringAssembler - A genome assembler for read-cloud technology
Version: 0.9-2dec66178be6ed266db4fd4163ac97bda29405c3
Usage:
          assembly3 [options] -1 read_1.fq -2 read_2.fq -l ust/bioturing/sorted
          local_assembly [options] -i graph.bin -1 R1_sorted.fq -2 R2_sorted.fq -l sorted -lk 31 -lc scaffold.full.fasta
Required parameters:
          -1            Forward reads. e.g: R1.fq or R1_lane1.fq R1_lane2.fq ...
          -2            Reverse reads. e.g: R2.fq or R2_lane1.fq R2_lane2.fq ...
          -l            Type of reads library:
                           bioturing: has BX:Z:<barcode> in the comment of the read names
                           sorted   : Sorted reads produced by TuringAssembler. Must be accompanied with -I barcode.idx
                           ust      : reads are generated by TELL-Seq protocols. Must be accompanied with -I I1.fq (I1_lane1.fq I1_lane2.fq)
Optional parameters:
          -t            Number of threads [4]
          -k0           Kmer size for global assembly process [45]
          -lk           Kmer size for local assembly [31]
          -lc           Output file after local assembly step [scaffold.full.fasta]
          -metagenomics Doing assembly for metagenomics dataset [no]
          -o            Output directory [./]
          -i            Input graph binary file. Only use for sub-processes [./graph_k_xx_level_x.bin]
          -sm           Maximum memory size for read sorting (GB) [32]
          -v            Verbose mode. Print log trace and log debug [no]
```

## Authors
-   Hao Tran
-   Huu Che (huu@bioturing.com)
-   Tan Phan (tan@bioturing.com)
-   Thang Tran
-   Son Pham (sonpham@bioturing.com)
