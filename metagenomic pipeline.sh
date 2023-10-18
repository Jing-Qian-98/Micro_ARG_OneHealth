# Activate the environment
conda activate MetaG

threads=2
######################################## Rawdata and metadata prepair ##########################################################


# Create a directory for analysis
mkdir -p  MetaTest
cd  MetaTest

# Prepare raw data and place it in the 1_rawdata directory
## Here we obtain test data from the /share/home/LiuChang/MetaTest_bk/1_rawdata path
## We use the following samples for testing: F01, F02, F03, P01, P02, P03, as the data volume is large, we perform random sampling
mkdir 1_rawdata
for sample in *.fq.gz; do
    #for testing：
    #zcat /share/home/LiuChang/MetaTest_bk/1_rawdata/${sample}_R1.fq.gz | seqkit sample -p 0.05 -s 111 -o 1_rawdata/${sample}_R1.fq.gz &
    #zcat /share/home/LiuChang/MetaTest_bk/1_rawdata/${sample}_R2.fq.gz | seqkit sample -p 0.05 -s 111 -o 1_rawdata/${sample}_R2.fq.gz  

    cp /share/home/LiuChang/MetaTest_bk/1_rawdata/${sample}_R1.fq.gz .  # You can replace 'cp' with 'mv'
    cp /share/home/LiuChang/MetaTest_bk/1_rawdata/${sample}_R2.fq.gz .  # You can replace 'cp' with 'mv'
    # or
    # ln -s /share/home/LiuChang/MetaTest_bk/1_rawdata/${sample}_R1.fq.gz .
    # ln -s /share/home/LiuChang/MetaTest_bk/1_rawdata/${sample}_R2.fq.gz .
    
    echo "$sample" >> name
done

# View how the fastq data looks
less -S 1_rawdata/PLX1_R1.fq.gz

######################################## Reads QualityControl ##########################################################

# Perform fastqc analysis first

mkdir FastQC
fastqc -t ${threads} 1_rawdata/* -o FastQC
multiqc FastQC/* -o MultiQC_fastqc

# Use fastp for quality control
mkdir -p 2_cleandata Fastp
for sample in *.fq.gz; do
    fastp -i 1_rawdata/${sample}_R1.fq.gz -I 1_rawdata/${sample}_R2.fq.gz \
        -o 2_cleandata/${sample}_clean_R1.fq.gz -O 2_cleandata/${sample}_clean_R2.fq.gz \
        -l 50 -g -w ${threads} -W 5 -5 -q 20 -u 30 \
        -j Fastp/${sample}_fastp.json -h Fastp/${sample}_fastp.html 
done

multiqc Fastp/*  -o MultiQC_fastp

# Remove host DNA (optional step)
## We use the bowtie2 mapping method
## First, build the index
cp /share/home/LiuChang/MetaTest_bk/hg38.fa .

bowtie2-build --threads ${threads} hg38.fa hg38 # This step takes some time; you can increase the ${threads} value for actual run

for sample in *.fq.gz; do
    bowtie2 -p ${threads} -x hg38 -1 2_cleandata/${sample}_clean_R1.fq.gz -2 2_cleandata/${sample}_clean_R2.fq.gz -S sample.sam --un-conc-gz 2_cleandata/${sample}_clean_nonhost.fq.gz 
done
rename ".fq.1.gz" "_R1.fq.gz"  2_cleandata/*fq.1.gz
rename ".fq.2.gz" "_R2.fq.gz"  2_cleandata/*fq.2.gz
# Remove ribosomal rRNA (optional step)

rRNA_16SB=/share/home/LiuChang/database/sortmernaDB/silva-bac-16s-id90.fasta 
rRNA_23SB=/share/home/LiuChang/database/sortmernaDB/silva-bac-23s-id98.fasta
rRNA_16SA=/share/home/LiuChang/database/sortmernaDB/silva-arc-16s-id95.fasta
rRNA_23SA=/share/home/LiuChang/database/sortmernaDB/silva-arc-23s-id98.fasta
rRNA_18SE=/share/home/LiuChang/database/sortmernaDB/silva-euk-18s-id95.fasta
rRNA_28SE=/share/home/LiuChang/database/sortmernaDB/silva-euk-28s-id98.fasta
rRNA_5S=/share/home/LiuChang/database/sortmernaDB/rfam-5s-database-id98.fasta
rRNA_5_8S=/share/home/LiuChang/database/sortmernaDB/rfam-5.8s-database-id98.fasta

function rmRNA_fun {
    rm -rf ${1}_sortmerna  
    mkdir -p SortMeRNA
    sortmerna \
        -ref $rRNA_16SB \
        -ref $rRNA_23SB \
        -ref $rRNA_16SA \
        -ref $rRNA_23SA \
        -ref $rRNA_18SE \
        -ref $rRNA_28SE \
        -ref $rRNA_5S \
        -ref $rRNA_5_8S \
        -reads 2_cleandata/${1}_R1.fq.gz \
        -reads 2_cleandata/${1}_R2.fq.gz \
        --fastx \
        --other 2_cleandata/${1}_nonrRNA \
        --aligned SortMeRNA/${1}_rRNA \
        --num_alignments 1 \
        --paired_in \
        -threads ${2} \
        -out2 \
        -workdir ${1}_sortmerna  
    #rm -rf ${1}_sortmerna/idx/
    #rm -rf ${1}_sortmerna/kvdb/
    #rm -f $rmDIR/*_rRNA.log

}

# for sample in *.fq.gz; do
#     rmRNA_fun ${sample}_clean_nonhost 
# done     
 
# Let's try running it, it's quite slow
rmRNA_fun ${sample}_clean_nonhost 20
## For subsequent analysis, we will use 2_cleandata/*_clean_nonhost_R{12}.fq.gz
######################################## Reads based taxonomic analysis ##########################################################

# kraken2
mkdir -p kraken2
cat name | xargs -i -P 1 kraken2 \
    --paired --threads ${threads} --confidence 0.5 --memory-mapping --db ~LiuChang/database/Kraken2db \
    --report kraken2/{}_kraken2.report --output /dev/null \
    2_cleandata/{}_clean_nonhost_R1.fq.gz 2_cleandata/{}_clean_nonhost_R2.fq.gz 

taxonomies="domain kingdom phylum class order family genus species"
for sample in *.fq.gz  do
   for taxo in ${taxonomies};do
       level=$(echo ${taxo:0:1} | tr "a-z" "A-Z")
       echo $level
       bracken \
        -r 150 -l ${level} \
        -d ~LN00/database/Kraken2db \
        -i kraken2/${sample}_kraken2.report \
        -o kraken2/${sample}_${level}.bracken > /dev/null
        # extract taxonID and abundance 
        awk -F "\t" -v sample=${sample} 'NR==1{print "taxonID\t",sample};NR>1{print $2"\t"$6}' kraken2/${sample}_${level}.bracken \
            > kraken2/${sample}_${level}.summary
   done
done


# kaiju 

mkdir -p kaiju
r1=$(cat name | xargs -i echo "2_cleandata/{}_clean_nonhost_R1.fq.gz" | tr "\n" "," | sed 's/,$//')
r2=$(cat name | xargs -i echo "2_cleandata/{}_clean_nonhost_R2.fq.gz" | tr "\n" "," | sed 's/,$//')
out=$(cat name | xargs -i echo "kaiju/kaiju_{}.out" | tr "\n" "," | sed 's/,$//')
nodes=${HOME}/database/kaijudb/nodes.dmp
names=names=${HOME}/database/kaijudb/names.dmp
## To execute exit 1 if kaiju did not succeed and skip the kaiju2summary step

kaiju-multi \
-z ${threads} -a "greedy" -x -e 3 -E 0.00001 \
-i ${r1} -j ${r2} -t ${nodes} \
-f ~/database/kaijudb/kaiju_db_refseq.fmi \
-o ${out}

for sample in F01 F02 F03 P01 P02 P03; do 
    kaiju2table -t ${nodes} -n ${names} -r species -o kaiju/kaiju_${sample}_summary.tsv kaiju/kaiju_${sample}.out -l superkingdom,phylum,class,order,family,genus,species
    sed -i "1s/reads/$sample/" kaiju/kaiju_${sample}_summary.tsv
done

mkdir -p Reads_annotation

conda deactivate
conda activate R4.3
Rscript mergeKaiju.R kaiju/ Reads_annotation # taxa_table_RC.tsv
conda deactivate


# taxonkit
conda activate MetaG
title=$(head -n1 Reads_annotation/taxa_table_RC.tsv)
sed '1d' Reads_annotation/taxa_table_RC.tsv \
    | taxonkit reformat -I 1 --data-dir ${HOME}/database/kaijudb -f "{k};{p};{c};{o};{f};{g};{s}" -r "Unclassified" -P \
    | sed -e"1i ${title}\ttaxonomy" -e "s/k__/d_/g" -e "s/__/_/g" \
    > Reads_annotation/taxa_table_RC.tsv.tmp
mv Reads_annotation/taxa_table_RC.tsv.tmp Reads_annotation/taxa_table_RC.tsv

# modify Unclassified line
sed -i "s/Unclassified;Unclassified;Unclassified;Unclassified;Unclassified;Unclassified;Unclassified/d_Unclassified;p_Unclassified;c_Unclassified;o_Unclassified;f_Unclassified;g_Unclassified;s_Unclassified/" Reads_annotation/taxa_table/taxa_table_RC.tsv


######################################## Reads assembly ##########################################################

# Run MegaHIT assembly

mkdir 3_assembly
# for sample in *.fq.gz; do 
#     megahit -1 2_cleandata/${sample}_clean_nonhost_R1.fq.gz -2 2_cleandata/${sample}_clean_nonhost_R2.fq.gz \
#         --out-prefix ${sample}  \
#         -o ${sample}_megahit_out -t ${threads} # --preset meta-large 
# done  

## --k-min 57 --k-max 137 --k-step 20 
sample=PLX1
megahit -1 2_cleandata/${sample}_clean_nonhost_R1.fq.gz -2 2_cleandata/${sample}_clean_nonhost_R2.fq.gz \
    --out-prefix ${sample}  \
    -o 3_assembly/${sample}_megahit_out -t ${threads} --k-min 57 --k-max 137 --k-step 20 

for sample in F01 F02 F03 P01 P02 P03; do 
    megahit -1 2_cleandata/${sample}_clean_nonhost_R1.fq.gz -2 2_cleandata/${sample}_clean_nonhost_R2.fq.gz \
        --out-prefix ${sample}  \
        -o 3_assembly/${sample}_megahit_out -t 7 &
done

mv 3_assembly/*_megahit_out/*.contigs.fa 3_assembly

# Merging contig files and renaming sequences.
rm -f 3_assembly/all.contigs.fa
cat name | xargs -i sed "s/^>/>{}_/" 3_assembly/{}.contigs.fa >> 3_assembly/all.contigs.fa


# Using mmseqs2 for sequence deduplication with a suggested 0.99 similarity threshold.

~LiuChang/software/mmseqs/bin/mmseqs easy-linclust 3_assembly/all.contigs.fa 3_assembly/all_mmseqs mmseqs_tmp/ \
    -c 0.95 --cov-mode 1 --min-seq-id 0.99 --threads ${threads}  -v 1 --cluster-mode 2 

rm -r mmseqs_tmp/






######################################## gene predict ##########################################################
mkdir -p 4_gene
# prodigal -p meta -f gff -i 3_assembly/all_mmseqs_rep_seq.fasta -o 4_gene/gene.gff -d 4_gene/gene.fa -a 4_gene/pro.fa  -q

# Splitting fasta files
seqkit  split -w 0 -f -s 50000 -O gene_tmp    3_assembly/all_mmseqs_rep_seq.fasta

# Gene prediction
ls gene_tmp/*.part_*.fasta | xargs -i -P ${threads} \
    prodigal -p meta -f gff -i {} -o {}.gff -d {}.gene -a {}.pro  -q



# Processing and reformatting gene prediction files
ls gene_tmp/*.gene | xargs -i -P ${threads} sh -c "
    cat {} | sed 's/ # /;/g' \
    | cut -d';' -f1,2,3,4,6 | sed 's/;/_/g;' \
    | awk  '{if(\$1~/^>/){print \"\n\"\$0}else{printf \$0}}END{print }' > {}_tmp 
    " 

## Intermediate gene prediction files
cat gene_tmp/*.gene_tmp | grep -v '^$' | awk '/^>/{print ">gene_"(NR+1)/2" "$1}!/^>/{print $0}' | sed 's/ >/ /' > 4_gene/genes_predict.fa

# Gene deduplication
vsearch --derep_fulllength 4_gene/genes_predict.fa  --output 4_gene/genes_derep.fa --fasta_width 0 --notrunclabels --uc  4_gene/genes_derep.uc


    
# Processing predicted protein files

ls gene_tmp/*.pro | xargs -i -P4 sh -c "
    cat {} | sed 's/ # /;/g' \
    | cut -d';' -f1,2,3,4,6 | sed 's/;/_/g;' \
    | awk  '{if(\$1~/^>/){print \"\n\"\$0}else{printf \$0}}END{print }' > {}_tmp 
    " 
## Intermediate protein files
cat gene_tmp/*.pro_tmp | grep -v '^$' | awk '/^>/{print ">gene_"(NR+1)/2" "$1}!/^>/{print $0}' | sed 's/ >/ /' > 4_gene/proteins_predict.fa

## Filtering protein sequences based on deduplicated genes
awk 'NR==FNR&&NR%2==1{a[$1]=$1}
    NR!=FNR&&FNR%2==1{if(a[$1]==$1){b=1;print}else{b=0}}
    NR!=FNR&&FNR%2==0{if(b==1)print}' \
    4_gene/genes_derep.fa 4_gene/proteins_predict.fa > 4_gene/proteins_derep.fa 

## Protein deduplication
rm -rf 4_gene/proteins_095 mmseqs_tmp2

~LiuChang/software/mmseqs/bin/mmseqs easy-linclust 4_gene/proteins_derep.fa    4_gene/proteins_095  mmseqs_tmp2 \
-c 0.9 --cov-mode 1 --min-seq-id 0.95 --threads ${threads} --cluster-mode 2  -v 1 

rm -rf mmseqs_tmp2


# Processing gene annotation files gff -> gtf


cat gene_tmp/*.gff | grep -v "^#" | sed '/^$/d' \
    | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$8}' \
    | awk -F "\t" '{ print $0"\tgene_id gene_"++count }' \
    | awk -F" " 'NR==FNR && /^>/ {gsub(">","",$1); $2=gensub(/(.*)_(.*)_(.*)_(.*)_(.*)_(.*)/,"\\1_\\3_\\4","g",$2); a[$2]=$1}
    NR!=FNR {b=$1"_"$4"_"$5; if(b in a){for(i=1;i<9;i++)printf $i"\t";print "gene_id "a[b]} }' \
    4_gene/genes_predict.fa - > 4_gene/genes_predict.gtf


#LC_ALL=C grep "^>" 4_gene/genes_derep.fa | sed "s/>//g" | awk -F" " 'NR==FNR{a[$1]=$1}NR!=FNR{if(a[$10])print}' - 4_gene/predict.gtf_tmp > 4_gene/genes_derep.fa



######################################## Taxonomy assignment of protain_derep ##########################################################

mkdir -p pro2tax

~LiuChang/software/mmseqs/bin/mmseqs createdb 4_gene/proteins_derep.fa 4_gene/proteins_derep.fa_qDB  -v 1
~LiuChang/software/mmseqs/bin/mmseqs createindex       4_gene/proteins_derep.fa_qDB  tmp_QUERY -v 1

~LiuChang/software/mmseqs/bin/mmseqs taxonomy 4_gene/proteins_derep.fa_qDB /PERSONALBIO196/LN00/database/uniref90/uniref90 pro2tax/taxonomyResult tmp_tax --lca-mode 3 \
            --search-type 3 -s 2 --threads ${threads} --db-load-mode 0 --max-seq-len 100000 \
            --lca-ranks superkingdom,phylum,class,order,family,genus,species --tax-lineage 1 \
            --max-accept 10 --max-seqs 100 --split-memory-limit 100G  

~LiuChang/software/mmseqs/bin/mmseqs createtsv 4_gene/proteins_derep.fa_qDB   pro2tax/taxonomyResult  pro2tax/taxonomyResult.tsv --threads ${threads} -v 1 

# the annotation results are filtered to remove unwanted species
#cut -f1,2,3,4,6 pro2tax/taxonomyResult.tsv \
#    | tee pro2tax/pro2taxonomy.tsv \
#    | grep -iv 'Viridiplantae\|Metazoa' \
#    | tee pro2tax/pro2taxonomy_filtered.tsv \
#    | cut -f1  > proteinsID


#seqkit grep -w 0 -f proteinsID  4_gene/proteins_derep.fa > 4_gene/proteins_derep_filtered.fa && rm proteinsID

cut -f1,2,3,4,6 pro2tax/taxonomyResult.tsv \ > pro2tax/pro2taxonomy.tsv 
cut -f1,5 pro2tax/pro2taxonomy.tsv \
    | awk -F"\t" '$2~/d_/{print $0}' \
    | awk -F"[\t;]" '{for(i=1;i<=NF;i++){if(i<NF){printf $i"\t"}else{printf $i"\n"}}}' \
    | sort -k1V,1 | sed '1iQUERY_ID\tsuperkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies' \
    > pro2tax/taxonomy_anno.xls

#cut -f1,5 pro2tax/pro2taxonomy.tsv \
#    | awk -F"[\t;]" '{print $1,$NF}' \
#    | awk -F"[ _]" '{print "gene_"$2,$3}' | awk -v OFS="\t" 'NF==1{print $1,"unclassifued"}NF>1{print $1,$2}' \
#    | sed '1iQUERY_ID\tOTU' > pro2tax/otu.xls


######################################## mapping for Reads Counts of genes ##############################################################

# bowtie2+feature
bowtie2-build --threads ${threads} -q 3_assembly/all_mmseqs_rep_seq.fasta 3_assembly/all_mmseqs_rep_seq

mkdir -p  5_abd
for sample in $(cat name);do
    # add: --no-unal
    bowtie2  --sensitive-local -x 3_assembly/all_mmseqs_rep_seq -1 2_cleandata/${sample}_clean_nonhost_R1.fq.gz -2 2_cleandata/${sample}_clean_nonhost_R2.fq.gz -p 10 | pigz -9 > 5_abd/${sample}.sam.gz
done

for sample in $(cat name);do
    pigz -d -c 5_abd/${sample}.sam.gz | ~LN00/software/subread-2.0.6-Linux-x86_64/bin/featureCounts -p --countReadPairs -t CDS -g gene_id -a 4_gene/genes_predict.gtf -o pre_${sample}.count -T ${threads}
    awk -F "\t" '/^gene_/{print $1"\t"$NF}' pre_${sample}.count > 5_abd/${sample}.count 

done

# Obtain counts for genes dedup
for sample in $(cat name);do
    awk 'NR==FNR&&/^S/{a[$9]=$9}NR==FNR&&/^H/{a[$9]=$11}NR!=FNR{b[a[$1]]+=$2}END{for(i in b)print i"\t"b[i]}' 4_gene/genes_derep.uc 5_abd/${sample}.count > 5_abd/${sample}_genes_derep.count
done

 
4_gene/genes_derep.uc

# coverm # - not used for now
for sample in $(cat name);do
    coverm contig -t 10 -m count --coupled 2_cleandata/${sample}_clean_nonhost_R1.fq.gz  2_cleandata/${sample}_clean_nonhost_R2.fq.gz --reference 3_assembly/all_mmseqs_rep_seq.fasta -o 5_abd/${sample}.coverm
done

######################################## KO anno for  proteins_095  ##############################################################

~LiuChang/software/mmseqs/bin/mmseqs createdb 4_gene/proteins_095_rep_seq.fasta 4_gene/proteins_095_rep_seq.fasta_qDB  -v 1
~LiuChang/software/mmseqs/bin/mmseqs createindex  4_gene/proteins_095_rep_seq.fasta_qDB  tmp_QUERY -v 1

mkdir -p KO_anno
rm -rf  tmp_ko  KO_anno/ko.result.mmseqs*
~LiuChang/software/mmseqs/bin/mmseqs search 4_gene/proteins_095_rep_seq.fasta_qDB ~LN00/database/kobas/ko.mmseqs2_DB  KO_anno/ko.result.mmseqs tmp_ko \
    --threads ${threads} --search-type 3 --split-memory-limit 40G \
    --max-seq-len 100000 --max-accept 1 --max-seqs 100 \
    --start-sens 1 -s 4 --sens-steps 2 -e 1e-5 -v 1

~LiuChang/software/mmseqs/bin/mmseqs convertalis 4_gene/proteins_095_rep_seq.fasta_qDB ~LN00/database/kobas/ko.mmseqs2_DB KO_anno/ko.result.mmseqs KO_anno/ko.mmseqs.m8 --threads 8 -v 1

cut -f-2 KO_anno/ko.mmseqs.m8 | sed '1iQID\tGID' >GID.csv
cp ~LiuChang/database/kobas/sqliteDB/KoGenes.db . ;
kosdb=KoGenes.db

#### query table: GID.csv -> query
#### target table: Kos,KoGenes
python -c \
"import pandas as pd; \
import sqlite3; \
conn= sqlite3.connect('${kosdb}'); \
df = pd.read_csv('GID.csv',sep = '\t'); \
df.to_sql('query', conn, if_exists='replace', index=False); \
cursor = conn.cursor(); \
cursor.execute('SELECT q.QID, ko.koid \
                FROM \
                    (SELECT k.koid, kg.gid \
                    FROM Kos AS k JOIN KoGenes AS kg \
                    ON k.koid = kg.koid) ko, \
                query q \
                WHERE ko.gid=q.GID;'); \
result=cursor.fetchall(); \
cursor.close(); \
conn.close(); \
anno='\n'.join([ '\t'.join(i) for i in result ]) if len(result)>0 else ''; \
print(anno)" > KO_anno/KO_anno.xls 

sed -i '1iQUERY_ID\tKOs'  KO_anno/KO_anno.xls

cut -f2 KO_anno/KO_anno.xls | awk -F"," 'NR>1{for(i=1;i<=NF;i++)print $i}' | sort | sed '1iID' > KO_anno/KO.csv

cp ~LiuChang/database/kobas/sqliteDB/ko.db . 
db=ko.db

#### query table: KO.csv -> query
#### target table: ko_level2
python -c \
"import pandas as pd; \
import sqlite3; \
conn= sqlite3.connect('${db}'); \
df = pd.read_csv('KO_anno/KO.csv'); \
df.to_sql('query', conn, if_exists='replace', index=False); \
cursor = conn.cursor(); \
cursor.execute('SELECT k.pathway_level1, k.pathway_level2 \
                FROM ko_level2 AS k \
                JOIN query AS q \
                ON k.ko_num=q.ID;'); \
result=cursor.fetchall(); \
cursor.close(); \
conn.close(); \
anno='\n'.join([ '\t'.join(i) for i in result ]) if len(result)>0 else ''; \
print(anno)" > result

### ---------------------------------------------------------- ### 

sort result | uniq -c | sed 's/^ {1,}//;s/ /\t/' \
| awk -F "\t" '{print $2"\t"$3"\t"$1}' | sed '1d' > KO_anno/KO_stat.xls



######################################## CAZy anno for  proteins_095  ##############################################################

mkdir -p CAZy_anno
rm -rf  tmp_ko  CAZy_anno/CAZy.result.mmseqs*
~LiuChang/software/mmseqs/bin/mmseqs search 4_gene/proteins_095_rep_seq.fasta_qDB /PERSONALBIO196/LN00/database/CAZyme/CAZyme.mmseqs2_DB  CAZy_anno/CAZy.result.mmseqs tmp_ko \
    --threads ${threads} --search-type 3 --split-memory-limit 40G \
    --max-seq-len 100000 --max-accept 1 --max-seqs 100 \
    --start-sens 1 -s 4 --sens-steps 2 -e 1e-5 -v 1

~LiuChang/software/mmseqs/bin/mmseqs convertalis 4_gene/proteins_095_rep_seq.fasta_qDB /PERSONALBIO196/LN00/database/CAZyme/CAZyme.mmseqs2_DB CAZy_anno/CAZy.result.mmseqs CAZy_anno/CAZy.mmseqs.m8 --threads ${threads} -v 1

cut -f-2 CAZy_anno/CAZy.mmseqs.m8 | sed -re 's/\|[0-9\-]+\.[0-9\-]+\.[0-9\-]+\.[0-9\-]+//g' > CAZy_anno/CAZy.anno
# Converting annotation IDs to GENEs
cat CAZy_anno/CAZy.anno | awk -F "[\t|]" '{{printf $1"\t"}for(i=3;i<=NF;i++){if(i<NF){printf $i","} else{printf $i"\n"}}}' \
| sort -k1V,1 | sed '1iQUERY_ID\tFamilys' > CAZy_anno/Family_anno.xls
