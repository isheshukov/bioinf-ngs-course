# NGS. HW2

## Илья Шешуков

### Часть 1. Анализ ридов Illumina
https://github.com/lh3/seqtk
#### Bwa mem

```shell
$ bwa index data/ref.fasta.gz
[bwa_index] Pack FASTA... 0.05 sec
[bwa_index] Construct BWT for the packed sequence...
[bwa_index] 1.55 seconds elapse.
[bwa_index] Update BWT... 0.03 sec
[bwa_index] Pack forward-only FASTA... 0.05 sec
[bwa_index] Construct SA from BWT and Occ... 0.32 sec
[main] Version: 0.7.17-r1188
[main] CMD: bwa index data/ref.fasta.gz
[main] Real time: 3.029 sec; CPU: 1.990 sec

$ bwa mem data/ref.fasta.gz data/frag.R1.fastq.gz data/frag.R2.fastq.gz > frag.sam
[...]

$ samtools flagstat frag.sam
3432529 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
165 + 0 supplementary
0 + 0 duplicates
3424657 + 0 mapped (99.77% : N/A)
3432364 + 0 paired in sequencing
1716182 + 0 read1
1716182 + 0 read2
3407384 + 0 properly paired (99.27% : N/A)
3416620 + 0 with itself and mate mapped
7872 + 0 singletons (0.23% : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)

$ bwa mem data/ref.fasta.gz data/jump.R1.fastq.gz data/jump.R2.fastq.gz > jump.sam
[...]

$ samtools flagstat jump.sam
4328969 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
13369 + 0 supplementary
0 + 0 duplicates
4279501 + 0 mapped (98.86% : N/A)
4315600 + 0 paired in sequencing
2157800 + 0 read1
2157800 + 0 read2
3764298 + 0 properly paired (87.23% : N/A)
4216854 + 0 with itself and mate mapped
49278 + 0 singletons (1.14% : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)

```
#### Покрытие генома

``` shell
$ samtools sort -o frag.sorted frag.sam
[bam_sort_core] merging from 1 files and 1 in-memory blocks...

$ samtools index frag.sorted

$ samtools depth frag.sorted > frag.depth

```