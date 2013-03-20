Usage

For paired-end Illumina:

./slopmap -1 <PE1 filename> -2 <PE1 filename> -l <Library name> -o <Output prefix> [-k <KMER_SIZE> -d <DISTANCE> -nch <NUM_CONSEQUITIVE_HITS> -nc <NUM_HITS> -t <THRESHOLD>]

For single-end Illumina:

./slopmap -U <SE filename> -l <Library name> -o <Output prefix> [-k <KMER_SIZE> -d <DISTANCE> -nch <NUM_CONSEQUITIVE_HITS> -nc <NUM_HITS> -t <THRESHOLD>]

For Roche 454:

./slopmap -454 <454 filename> -l <Library name> -o <Output prefix> [-k <KMER_SIZE> -d <DISTANCE> -nch <NUM_CONSEQUITIVE_HITS> -nc <NUM_HITS>  -t <THRESHOLD>]


Default parameters:

-k <KMER_SIZE> 15
-d <DISTANCE> 5 (distance between two consequitive kmers)
-nch <NUM_CONSEQUITIVE_HITS> 3
-nc <NUM_HITS> 10
-t <THRESHOLD> 75
Examples

-Paired-end Illumina:
./slopmap -1 ../test_data/SmallTestIllumina_R1.fastq.gz -2 ../test_data/SmallTestIllumina_R2.fastq.gz -o test -l ../test_data/vectors.fasta -t 90

-Single-end Illumina:
./slopmap -U ../test_data/SmallTestIllumina_R1.fastq.gz -o test -l ../test_data/vectors.fasta -t 83

-Roche 454:
./slopmap -454 ../test_data/Small454Test.sff -o test -l ../test_data/vectors.fasta -t 56


