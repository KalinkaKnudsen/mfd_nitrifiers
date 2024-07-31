
# load software modules or environments

eval "$(conda shell.bash hook)"
CONDA_ENV_NAME="graftm_env"
# Activate the Conda environment
conda activate "$CONDA_ENV_NAME"
cd ~/conda/graftM/bin/
export PATH=$PWD:$PATH

module load parallel


WD=path_to_WD
MFD=path_to_samples

### Now to the grafting

# make directories for output

mkdir -p $WD/graftm
mkdir -p $WD/log

# temporary folder for parallel
temp=/home/bio.aau.dk/vj52ou/temp


THREADS=4
GRAFTM_PACKAGE=path_to_graftM_package

cd $WD


ls $MFD | grep "R1" | sed 's/_R1.fastq.gz//' | parallel -j40 --tmpdir $temp graftM graft --forward $MFD/{}_R1.fastq.gz \
 --graftm_package $GRAFTM_PACKAGE --output_directory $WD/graftm/{} --verbosity 5 \
 --threads $THREADS --force --input_sequence_type nucleotide  --search_method hmmsearch+diamond '&>' $WD/log/{}.log


############ Now to extracting the coverage profiles ##########

mkdir $WD/position_files
odir_pos=$WD/position_files


#First, the HMM output file needs to be sorted and the HMM position needs to be extracted 
for line in $(ls $WD/graftm); do
  file_path="$WD/graftm/$line/${line}_R1/${line}_R1.hmmout.txt"

  if [ -e "$file_path" ]; then
    grep -v "#" "$file_path" \
    | tr -s ' ' | awk '{ if ($12 < 1e-10) {print $0}}' \
    | awk -F' ' -v OFS='\t' "{print \"$line\", \$1, \$16, \$17}" \
    >> "$odir_pos/hmm_hit_file.tsv"
    echo "$line"
  else
    echo "File does not exist: $file_path"
  fi
done

## Then sorting this file based on the read_ID
sort -k2 $odir_pos/hmm_hit_file.tsv > $odir_pos/hmm_hit_file_sorted.tsv

#Then I want to concatenate all the read_tax files
for line in $(ls $WD/graftm); do
  file_path="$WD/graftm/$line/${line}_R1/${line}_R1_read_tax.tsv"

  if [ -e "$file_path" ]; then
    cat "$file_path" >> "$odir_pos/tax_file.tsv"
  else
    echo "File does not exist: $file_path"
  fi
done


sort -k1 $odir_pos/tax_file.tsv > $odir_pos/tax_file_sorted.tsv


### Joining and printing header below. Since I join by the tax file, the high E-value reads should be removed, as they are not reported to the read_tax file ###
(join -t $'\t' -1 1 -2 2 -o '2.1 1.1 1.2 2.3 2.4' $odir_pos/tax_file_sorted.tsv $odir_pos/hmm_hit_file_sorted.tsv \
 | awk -v OFS='\t' 'BEGIN{print "Sample\tSequence\tTax\thmm_start\thmm_end"} {print $0}') \
> $odir_pos/position_of_hits.tsv

### Removing all files but the position of hits file ###
rm -r $odir_pos/tax_file_sorted.tsv $odir_pos/hmm_hit_file_sorted.tsv $odir_pos/tax_file.tsv $odir_pos/hmm_hit_file.tsv

echo "position files finished"

##### Then, making quick coverage profiles ###
Rscript $WD/coverage_profiles.R "$WD"




#################### Making e-10 df ###############

for line in $(ls $WD/graftm); do
  output_file="$WD/graftm/$line/combined_count_table_e10.txt"

  # Check if required files exist
  if [ -e "$WD/graftm/$line/"$line"_R1/"$line"_R1.hmmout.txt" ] && [ -e "$WD/graftm/$line/"$line"_R1/"$line"_R1_read_tax.tsv" ]; then
    echo -e "#ID\t$line\tConsensusLineage" > "$output_file"
    count=1
    awk '{ if ($12 < 1e-10) {print $1}}' "$WD/graftm/$line/"$line"_R1/"$line"_R1.hmmout.txt" \
      | grep -f - "$WD/graftm/$line/"$line"_R1/"$line"_R1_read_tax.tsv" \
      | cut -f 2 | sort | uniq -c | sed -e 's/^[ \t]*//' | sed 's/\([[:digit:]]\)\s\+Root/\1\tRoot/' \
      | while read -r line2; do
        echo -e "$count\t$line2" >> "$output_file"
        ((count++))
      done
  else
    echo "Required files do not exist for $line"
  fi
done

### Making combined count table for both e10 and regular
module load SciPy-bundle/2022.05-foss-2020b
python3 join_files.py -e txt -f $WD/combined_count_table_e10 \
-n 0 -s '\t' -p combined_count_table_e10 $WD/graftm ConsensusLineage

module purge

sed -i 's/_R1//g' $WD/combined_count_table_e10.txt

