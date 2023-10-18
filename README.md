# Microbiota-and-Antibiotic-Resistance-Genes-at-Human-Pig-Soil-Interface
Our project, 'Microbiota and Antibiotic Resistance Genes at the Human-Pig-Soil Interface,' seeks to investigate the complex interactions and transmission of pathogens among humans, animals, and their environment within two pig farms situated in Chongming District, Shanghai. To achieve this, we have collected a total of 79 samples, comprising employee feces, pig manure, and soil samples. Our primary objective is to gain insights into the microbial composition, diversity, and antibiotic resistance transmission at this critical interface. We use Index 2 × 150 paired-end sequencing through Illumina NovaSeq 6000 for raw data. Our analytical pipeline includes genome assembly programs, taxonomic classifiers, read mappers, and pathway reconstruction tools.

| Step                          | Software             | Software Version | Software Link                                    | Main Parameters                                         |
|-------------------------------|----------------------|-------------------|--------------------------------------------------|----------------------------------------------------------|
| Sequencing Data Quality Control| fastp                | 0.23.2            | [fastp GitHub](https://github.com/OpenGene/fastp) | -l 50 -g -W 5 -5 -q 20 -u 30                           |
| Species Annotations for Unspliced Sequences | Kaiju  | 1.9.0  | [Kaiju GitHub](https://github.com/bioinformatics-centre/kaiju) | -a greedy -x -e 5 -E 0.00001                         |
| Sequence Assembly and Gene Prediction | MEGAHIT | v1.1.2  | [MEGAHIT GitHub](https://github.com/voutcn/megahit) | --k-list 33,55,77,99,127 --min-contig-len 300        |
| Gene Prediction               | Prodigal             | V2.6.3            | [Prodigal GitHub](https://github.com/hyattpd/Prodigal) | -p meta                                                 |
| Gene De-redundancy            | VSEARCH              | v2.7.0_linux_x86_64 | [VSEARCH GitHub](https://github.com/torognes/vsearch) | --fasta_width 0 --notrunclabels                      |
| Protein De-redundancy         | MMseqs easy-linclust | 13.45111          | [MMseqs2 GitHub](https://github.com/soedinglab/MMseqs2) | -c 0.9 --cov-mode 1 --min-seq-id 0.95 --cluster-mode 2 |
| Gene Sequence Species and Function Annotation | MMseqs taxonomy | 13.45111 | [MMseqs2 GitHub](https://github.com/soedinglab/MMseqs2) | --lca-mode 3 -s 2 --max-seq-len 100000 --max-accept 10 --max-seqs 100 |
| CARD Database Comparison      | MMseqs search        | 13.45111          | [MMseqs2 GitHub](https://github.com/soedinglab/MMseqs2) | -s 5.7 --max-seq-len 100000 --max-accept 1 --max-seqs 100 |
