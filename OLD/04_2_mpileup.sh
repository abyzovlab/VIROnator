#!/usr/bin/bash


#bash ${mpileup_script} ${file} ${mpileup_output} ${sample} ${virus} ${virus_size} ${virus_name} ${dataset} ${reads} ${viral_sizes_bed} ${sample_depth} ${nbins} ${mpileup_cov_plots}

file=$1
mpileup_output=$2
sample=$3
virus=$4
virus_size=$5
virus_name=$6
dataset=$7
reads=$8
viral_sizes_bed=$9
bins=${10}
sample_depth=${11}
nbins=${12}
mpileup_cov_plots=${13}


echo $file
echo $virus

# MPILEUP
# ===================================
# all reads
#samtools mpileup -A $file > $mpileup_output/${dataset}.${sample}.${virus}.mpileup.all.temp
# anomalous reads excluded
#samtools mpileup $file > $mpileup_output/${dataset}.${sample}.${virus}.mpileup.clean.temp



# AVERAGE COVERAGE FOR BINS
# ===================================
# VIRUS DEPTH FOR THE SAMPLE
VIRUS_DEPTH="${mpileup_output}/${dataset}.${sample}.${virus}.depth"
grep ${virus} ${sample_depth} > ${VIRUS_DEPTH}

# CONVERT INTO BED FORMAT
VIRUS_DEPTH_BED="${mpileup_output}/${dataset}.${sample}.${virus}.depth.bed"
awk 'BEGIN{OFS="\t"} {contig=$1; pos=$2+0; depth=$3+0; print contig, pos-1, pos, depth}' $mpileup_output/${dataset}.${sample}.${virus}.depth > ${VIRUS_DEPTH_BED}

# FETCH BINS FOR THE VIRUS
VIRUS_BINS="${mpileup_output}/${dataset}.${sample}.${virus}.nb_${nbins}.bed"
grep ${virus} ${bins} > ${VIRUS_BINS}

# AVERAGE COVERAGE FOR EACH VIRAL BIN
VIRUS_AVG_COV="${mpileup_output}/${dataset}.${sample}.${virus}.nb_${nbins}_avgcov.bed"
bedtools map -a ${VIRUS_BINS} -b ${VIRUS_DEPTH_BED} -c 4 -o mean | awk 'BEGIN{OFS="\t"} {if ($4=="." || $4=="") { $4=0 } printf "%s\t%s\t%s\t%06.3f\n", $1, $2, $3, $4}' > ${VIRUS_AVG_COV}



# BREADTH OF COVERAGE (FLAT COVERAGE)
# ===================================
# FLATTENNING
awk 'BEGIN{OFS="\t"} {if ($4>0) $4=1; else $4=0; print}' ${VIRUS_DEPTH_BED} > ${VIRUS_DEPTH_BED%.bed}.flat.bed
bedtools map \
  -a ${VIRUS_BINS} \
  -b ${VIRUS_DEPTH_BED%.bed}.flat.bed \
  -c 4 -o mean \
| awk 'BEGIN{OFS="\t"} {if ($4=="." || $4=="") $4=0; print $1,$2,$3,$4}' \
> ${VIRUS_BINS%.bed}.flatBinFrac.bed

# THRESHOLDING
t=0.10
tlabel="t$(printf "%.0f" "$(echo "$t*100" | bc)")pct"
VIRUS_AVG_COV_FLAT=${mpileup_output}/${dataset}.${sample}.${virus}.nb_${nbins}_flatBinFrac_${tlabel}.tsv
awk -v dataset="$dataset" -v sample="$sample" -v virus="$virus" \
    -v nbins="$nbins" -v t="$threshold" '
BEGIN {
    OFS="\t"; passing=0; total=0
}
{
    total++
    if ($4 >= t) passing++
}
END {
    printf "%s\t%s\t%s\t%d\t%.2f\t%d\t%.6f\n",
           dataset, sample, virus, nbins, t, passing, passing/total
}' ${VIRUS_BINS%.bed}.flatBinFrac.bed \
> ${VIRUS_AVG_COV_FLAT}



# BREADTH OF COVERAGE (MEAN COVERAGE)
# ===================================
VIRUS_AVG_COV_BREADTH="${mpileup_output}/${dataset}.${sample}.${virus}.nb_${nbins}_avgcov_breadth"
t="0.1"
awk -v t="$t" -v virus="$virus" -v sample="$sample" -v dataset="$dataset" 'BEGIN{ if (t=="") t=0.2; total=0; covered=0 } {cov=$4; total++; if(cov>=t) covered++} END {printf "%s\t%s\t%s\t%.3f\tcovered=%d\ttotal=%d\t%06.3f\n", dataset, sample, virus, t, covered, total, covered/total}' ${VIRUS_AVG_COV} > ${VIRUS_AVG_COV_BREADTH}



# EVENESS OF COVERAGE
# ===================================
awk -v dataset="$dataset" sample="$sample" -v virus="$virus" -v nbins="$nbins" '
BEGIN {
    OFS="\t"
    total=0
    print "sample","virus","nbins","H","J"
}
{
    c[NR]=$4
    total += $4
}
END {
    B = NR
    H = 0
    for (i=1; i<=B; i++) {
        if (c[i] > 0) {
            p = c[i] / total
            H += -p * log(p)
        }
    }
    J = H / log(B)
    printf "%s\t%s\t%s\t%d\t%.6f\t%.6f\n", dataset, sample, virus, nbins, H, J
}
' $mpileup_output/${dataset}.${sample}.${virus}.nb_${nbins}_avgcov.bed > $mpileup_output/${dataset}.${sample}.${virus}.nb_${nbins}_avgcov.evenness


if false; then
awk '
BEGIN {
    OFS="\t"
    total=0;
}
{
    c[NR]=$4;      # store coverage value
    total += $4;   # sum coverage
}
END {
    B = NR
    H = 0
    for(i=1; i<=B; i++){
        if (c[i] > 0){
            p = c[i] / total
            H += -p * log(p)
        }
    }
    J = H / log(B)
    print "Shannon_evenness_J:", J
}
' $mpileup_output/${dataset}.${sample}.${virus}.nb_${nbins}_avgcov.bed > $mpileup_output/${dataset}.${sample}.${virus}.nb_${nbins}_avgcov.evenness
fi # false for turned off


# selected reads
#samtools mpileup $file > $mpileup_output/${dataset}.${sample}.${virus}.pileup.temp


#less $mpileup_output/${dataset}.${sample}.${virus}.pileup.temp | grep "$virus" | awk '$3!="N"' | awk '{print $2"\t"$4}' > $mpileup_output/${dataset}.${sample}.${virus}.pileup
# grep "$virus" $mpileup_output/${dataset}.${sample}.${virus}.mpileup.all.temp | awk '{print $2"\t"$4}' > $mpileup_output/${dataset}.${sample}.${virus}.pileup.all
# grep "$virus" $mpileup_output/${dataset}.${sample}.${virus}.mpileup.clean.temp | awk '{print $2"\t"$4}' > $mpileup_output/${dataset}.${sample}.${virus}.pileup.clean
#less $mpileup_output/${dataset}.${sample}.${virus}.pileup.temp | grep "$virus" | awk '{print $2"\t"$4}' > $mpileup_output/${dataset}.${sample}.${virus}.pileup
#ls $mpileup_output/${sample}.${virus}.pileup



# PLOTTING COVERAGE PLOTS
# ===================================
# python ${mpileup_cov_plots} "$mpileup_output/${dataset}.${sample}.${virus}.pileup.all" "$mpileup_output" "$virus_size" "$sample" "$virus" "$virus_name" "$dataset" "$reads"
# python ${mpileup_cov_plots} "$mpileup_output/${dataset}.${sample}.${virus}.pileup.clean" "$mpileup_output" "$virus_size" "$sample" "$virus" "$virus_name" "$dataset" "$reads"
#python /home/mayo/m277455/h2/genomes/brains/exogene_output/scripts/plots_individual.py "$mpileup_output/${dataset}.${sample}.${virus}.pileup" "$mpileup_output" "$virus_size" "$sample" "$virus" "$virus_name" "$dataset" "$reads"

#find $mpileup_output -type f -name '*pileup*' -empty -delete
#find $mpileup_output -type f -name '*pileup' -empty -delete
