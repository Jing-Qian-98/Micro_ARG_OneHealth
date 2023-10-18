# QIIME2 recommends a minimum of 50 samples for constructing meaningful classifier models, e.g., for medical diagnosis. Therefore, this is not included in standard analysis. When the total number of samples is less than 50, it is suggested to use it only as a means to find marker species. Classifier accuracy can easily overfit and lacks practical reference value. Hence, we perform nested stratified cross-validation analysis to determine the importance scores of all samples' species (ASV/OTU or taxa).

# In this analysis, we use the 'classify-samples-ncv' method from q2-sample-classifier, employing nested cross-validation (NCV) for predictions on all samples. Here, we use 10-fold cross-validation based on the number of samples. We also create confusion matrices and ROC visualizations.

# For user documentation, please refer to https://forum.qiime2.org/t/qiime-2-15-q2-sample-classifier-2019-7/12186

# Create directories and copy necessary files

mkdir -p random_forest && cd random_forest
cp /share/home/LiuChang/MetaTest/data/map.txt ./
cp /share/home/LiuChang/MetaTest/data/1.5.2_TaxonAbdByGenes/taxa_table_genus_RC.tsv ./
## Set variable values
m=`cat map.txt |awk -F '\t'  '{print $4}' |sed -n '1p'`  ## Used as the grouping list for prediction
fm=`cat map.txt|sed 1d |cut -f4 |sort|uniq -c |awk '{print $1}' |sort -n|head -1|awk '{if($1>10){print 50}else{print $1*5}}'`
pms=`cat map.txt|sed 1d |cut -f4 |sort|uniq -c |awk '{print $1}' |sort -n|head -1|awk '{if($1>10){print 10}else{print $1}}'`

## Determine the number of cross-validation iterations
if [ $[`cut -f4 map.txt|sed '1d'|sort|uniq -c|awk '{print $1}'|sort -k 1,1 -n|tail -n 1`-2] -ge 10 ];then 
 pcv=10; 
elif [ $[`cut -f4 map.txt|sed '1d'|sort|uniq -c|awk '{print $1}'|sort -k 1,1 -n|tail -n 1`-2] -ge 5 ];then 
 pcv=5;
else
 pcv=$[`cut -f4 map.txt|sed '1d'|sort|uniq -c|awk '{print $1}'|sort -k 1,1 -n|tail -n 1`-2]
fi
pofs="    "



## Preprocessing abundance data
#cat taxa_table_species_RC.tsv | awk -F "\t" 'BEGIN {OFS="\t"} {$NF=""; print}' | sed 's/.$//g' > l6.txt
#source /PERSONALBIO/work/microbio/m00/software/biosoftware/miniconda3/bin/activate /PERSONALBIO/work/microbio/m13/bin/env/PERSONALBIO/Work/Mb/Mb07/.conda/envs/r-3.3-test
source /share/home/LiuChang/software/biosoftware/miniconda3/bin/deactivate
source /share/home/LiuChang/software/biosoftware/miniconda3/bin/activate /share/home/LiuChang/miniconda3/envs/qiime2-2022.11
qiime dev refresh-cache
biom convert -o otu_table_L6.biom -i taxa_table_genus_RC.tsv --table-type="Taxon table" --process-obs-metadata=taxonomy --to-hdf5
biom convert -i otu_table_L6.biom -o genus.biom --to-hdf5 ##??
qiime tools import --input-path otu_table_L6.biom --type 'FeatureTable[Frequency]'  --output-path genus.qza

#NCVNCV provides results for a mixed nested cross-validation model involving all samples. It outputs feature importance scores and sample predictions.
qiime sample-classifier classify-samples-ncv \
--i-table genus.qza --m-metadata-file map.txt \
--m-metadata-column $m $pofs \
--p-parameter-tuning \
--p-estimator RandomForestClassifier \
--p-cv $pcv \
--p-n-estimators 100 \
--p-random-state 666 \
--output-dir genus_ncv

# Output QZA results as QZV files
qiime sample-classifier confusion-matrix \
--i-predictions genus_ncv/predictions.qza \
--i-probabilities genus_ncv/probabilities.qza \
--m-truth-file map.txt \
--m-truth-column $m \
--o-visualization genus_ncv/accuracy_results.qzv

qiime metadata tabulate \
--m-input-file genus_ncv/feature_importance.qza \
--o-visualization genus_ncv/feature_importance.qzv

qiime metadata tabulate \
--m-input-file genus_ncv/predictions.qza \
--o-visualization genus_ncv/predictions.qzv

# Export QZA results as TSV files
cd genus_ncv
qiime tools export  --input-path probabilities.qza  --output-path probabilities
qiime tools export  --input-path feature_importance.qza  --output-path feature_importance
sed "1s/^/ID/" probabilities/class_probabilities.tsv > metadata.tsv
awk -F "\t" 'NR==FNR{a[$1]=$0}NR>FNR{if($1 in a){print a[$1]}}' ../map.txt metadata.tsv | sed '1i #SampleID\tBarcodeSequence\tLinkerPrimer Sequence\tGroup\tDescription' > map.txt

source /share/home/LiuChang/software/biosoftware/miniconda3/bin/deactivate
source /share/home/LiuChang/software/biosoftware/miniconda3/bin/activate /share/home/LiuChang/project/MetaDenovo_V1/envs/miniconda3/envs/R3.6.1
# Generate feature importance plots
cd feature_importance
Rscript /share/home/LiuChang/script/feature_importance.R importance.tsv 10
cd ..
## ROC
Rscript /share/home/LiuChang/script/ROC_plot.R
gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=pngalpha -r300 -sOutputFile=ROC.png ROC.pdf
# Display accuracy results using QIIME2 view (https://view.qiime2.org/)
mv metadata.tsv probabilities.tsv
mkdir accuracy_results ;cd accuracy_results
unzip ../accuracy_results.qzv 
cd ../
cp accuracy_results/*/data/*png . 
cp accuracy_results/*/data/*pdf . 
cp accuracy_results/*/data/*tsv .
cd ../

source /share/home/LiuChang/software/biosoftware/miniconda3/bin/deactivate
source /share/home/LiuChang/software/biosoftware/miniconda3/bin/activate /share/home/LiuChang/miniconda3/envs/qiime2-2022.11
qiime dev refresh-cache
# Test and Train: 80% of samples are used as the training set to build the model, and the remaining 20% are used as a test set for evaluation.
qiime sample-classifier classify-samples \
  --i-table genus.qza \
  --m-metadata-file map.txt  \
  --m-metadata-column $m \
  --p-parameter-tuning \
  --p-estimator RandomForestClassifier \
  --p-n-estimators 100 \
  --p-random-state 666 \
  --p-cv $pcv \
  --p-test-size 0.2 \
  --output-dir genus

cd genus
qiime tools export  --input-path probabilities.qza  --output-path probabilities
qiime tools export  --input-path feature_importance.qza  --output-path feature_importance
sed "1s/^/ID/" probabilities/class_probabilities.tsv > metadata.tsv
awk -F "\t" 'NR==FNR{a[$1]=$0}NR>FNR{if($1 in a){print a[$1]}}' ../map.txt metadata.tsv \
|sed '1i #SampleID\tBarcodeSequence\tLinkerPrimer Sequence\tGroup\tDescription' > map.txt

source /share/home/LiuChang/software/biosoftware/miniconda3/bin/deactivate
source /share/home/LiuChang/software/biosoftware/miniconda3/bin/activate /share/home/LiuChang/project/MetaDenovo_V1/envs/miniconda3/envs/R3.6.1
## feature_importance
cd feature_importance
Rscript /share/home/LiuChang/script/feature_importance.R importance.tsv 10
cd ..
## ROC
Rscript /share/home/LiuChang/script/ROC_plot.R
gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=pngalpha -r300 -sOutputFile=ROC.png ROC.pdf
# Display accuracy results using QIIME2 view (https://view.qiime2.org/)
mv metadata.tsv probabilities.tsv
mkdir accuracy_results ;cd accuracy_results
unzip ../accuracy_results.qzv 
cd ../
cp accuracy_results/*/data/*png . 
cp accuracy_results/*/data/*pdf . 
cp accuracy_results/*/data/*tsv .
cd ../
