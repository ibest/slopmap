Usage

For paired-end Illumina:

./slopmap -1 <PE1 filename> -2 <PE1 filename> -l <Library name> -o <Output prefix> [-k <KMER_SIZE> -d <DISTANCE> -nch <NUM_CONSEQUITIVE_HITS> -nc <NUM_HITS> --mode2]


Default parameters:

-k <KMER_SIZE> = 15
-d <DISTANCE> = 5 (distance between two consequitive kmers)
-nch <NUM_CONSEQUITIVE_HITS> = 3
-nc <NUM_HITS> = 10
--mode2 = NUM_CONSEQUITIVE_HITS (where --mode2 flag defines two strategies: whether use NUM_CONSEQUITIVE_HITS or NUM_HITS)

Example: 
./slopmap -1 ../test_data/SmallTestIllumina_R1.fastq.gz -2 ../test_data/SmallTestIllumina_R2.fastq.gz -o test -l ../test_data/vectors.fasta