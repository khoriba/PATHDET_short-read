#!/bin/bash

THREADS=20
SPLIT=50000
DEBUG=False
MODE=short

while [ $# -gt 0 ] ; do
  case "$1" in
    -1             ) FQ1=$2
                     shift
                     shift ;;
    -2             ) FQ2=$2
                     shift
                     shift ;;
    -o | --out     ) OUT=$2
                     shift
                     shift ;;
    -t | --threads ) THREADS=$2
                     shift
                     shift ;;
    -s | --split   ) SPLIT=$2
                     shift
                     shift ;;
    --debug        ) DEBUG=True
                     shift ;;
                 * ) echo "ERROR : illegal option"
                     echo "you need to supply path to fastq using -1 and -2 option, and output name using -o/--out option"
                     echo "you could overwrite thread number using -t/--threads option, 64 by default"
                     echo "you could overwrite split number using -s/--split option, 50000 by default"
                     exit ;;
  esac
done

if [ -z "$FQ1" ] || [ -z "FQ2" ] || [ -z "$OUT" ]; then
  echo "you need to supply path to fastq using -1 and -2 option, and output name using -o option"
  exit
elif [ ! -f "$FQ1" ] || [ ! -f "FQ2" ]; then
  echo "fastq file does NOT exists"
  exit
fi

blast_threads=8
if [ "$THREADS" -le "$blast_threads" ]; then
  blast_threads=$THREADS
  parallel_jobs=1
else
  parallel_jobs=`expr $THREADS \/ $blast_threads`
fi



START=`date '+%y%m%d%H%M%S'`
TIMEA=`date +%s`

## Directory PATH
path0=output/$OUT
path1=$path0/${OUT}_qc
path2=$path0/${OUT}_hgs
path3=$path0/${OUT}_blast
path4=$path0/${OUT}_tbl
path5=$path0/${OUT}_map
path_rapid=$path0/${OUT}_rapid

##DataBase
datadir="/path/to/datadir"
# kraken2
PATHDB1=$datadir/kraken2/k2_pluspf_20241228 # all pathogen
GRCH38=$datadir/kraken2/human # GRCh38
# blast
NT=$datadir/blast/nt/nt
HOST2=$datadir/blast/t2t/GCF_009914755.1_T2T-CHM13v2.0_genomic
# bowtie2
HOST=$datadir/bowtie2/GRCh38_noalt_as/GRCh38_noalt_as
# hisat2
GRCH38R=$datadir/hisat2/grch38_tran/genome_tran
# ncbi refseq
path_ng=$datadir/ncbi_genome

export TAXONKIT_DB=$datadir/taxonkit/taxdump

# load spark module
spark=$CONDA_PREFIX/opt/git/SparK/SparK.py



function p01_qc() {
  if [ -f "$path1/.done" ]; then
    echo "QC    process already finished, skipped"
    return
  fi

  mkdir -p $path1

  ##LOG
  START1=`date '+%y%m%d%H%M%S'`
  OUTL1=$path1/${OUT}_log_p1.txt
  echo "${START1} ${OUT} raw FQ stats" > $OUTL1
  seqkit stats -j $THREADS ${FQ1} ${FQ2} >> $OUTL1

  # fastp output
  out1=$path1/${OUT}_trimmed_R1.fastq.gz
  out2=$path1/${OUT}_trimmed_R2.fastq.gz
  # bbduk output
  out3=$path1/${OUT}_trimmed2_R1.fastq.gz
  out4=$path1/${OUT}_trimmed2_R2.fastq.gz

  ##QC
  fastp \
    -i ${FQ1} -I ${FQ2} -o $out1 -O $out2 \
    -h $path1/${OUT}_fp.html -j $path1/${OUT}.json \
    -3 -q 30 -n 20 -l 80 -f 15 -t 1 -T 1 \
    --thread $THREADS 

  bbduk.sh \
    in=$out1 in2=$out2 \
    out=$out3 out2=$out4 \
    entropy=0.5 entropywindow=10 entropyk=5 
  if [ $? -eq 0 ]; then touch $path1/.done.bbduk ; fi

  QC=`date '+%y%m%d%H%M%S'`
  echo -e "\n${QC} ${OUT} trimmed FQ stats" >> $OUTL1
  seqkit stats -j $THREADS $out3 $out4 --quiet >> $OUTL1
  echo "PATHDET_QC_FINISHED" | tee -a $OUTL1

  if [ -f "$path1/.done.bbduk" ]; then
    touch $path1/.done
  else
    echo "an error occurred"
    echo "the analysis was interrupted at the QC process"
    exit
  fi
}

function p02_hgs() {
  if [ -f "$path2/.done" ]; then
    echo "HGS   process already finished, skipped"
    return
  fi

  mkdir -p $path2/HGS

  # input: bbduk output
  in1=$path1/${OUT}_trimmed2_R1.fastq.gz
  in2=$path1/${OUT}_trimmed2_R2.fastq.gz

  ##LOG
  START2=`date '+%y%m%d%H%M%S'`
  OUTL2=$path2/${OUT}_log_p2.txt
  echo "${START2} ${OUT} trimmed FQ stats" > $OUTL2
  cat $path1/${OUT}_log_p1.txt | sed -n 7,9p >> $OUTL2

  # kraken output
  out1=$path2/HGS/${OUT}_hgs_1.fastq.gz
  out2=$path2/HGS/${OUT}_hgs_2.fastq.gz
  # bowtie output
  out3=$path2/HGS/${OUT}_hgs2.1.fastq
  out4=$path2/HGS/${OUT}_hgs2.2.fastq
  # hisat output
  out5=$path2/HGS/${OUT}_hgs.1.fastq
  out6=$path2/HGS/${OUT}_hgs.2.fastq

  ##Human Genome Subtraction
  kraken2 \
    --db ${GRCH38} \
    --paired $in1 $in2 \
    --unclassified-out $path2/HGS/${OUT}_hgs"#".fastq.gz \
    --threads $THREADS \
    > $path2/HGS/${OUT}.kraken.txt 

  bowtie2 \
    -x ${HOST} \
    -1 $out1 -2 $out2 \
    --un-conc $path2/HGS/${OUT}_hgs2.fastq \
    -S $path2/HGS/${OUT}_hgs2.sam \
    -p $THREADS 

  hisat2 \
    -x ${GRCH38R} \
    -1 $out3 -2 $out4 \
    --un-conc $path2/HGS/${OUT}_hgs.fastq \
    -S $path2/HGS/${OUT}_hgs3.sam \
    -p $THREADS 

  #Host Genome Subtraction plus
  #FQ2FA
  seqkit fq2fa $out5 -o ${out5%.fastq}.fasta
  blastn \
    -db ${HOST2} \
    -query ${out5%.fastq}.fasta \
    -evalue 1.0e-10 \
    -max_target_seqs 1 -max_hsps 5 \
    -outfmt "7 std" -num_threads $THREADS > $path2/${OUT}_hgs.txt
  cat $path2/${OUT}_hgs.txt \
    | grep -B 2 "# 0 hits found" \
    | grep "Query" \
    | sed -e 's/# Query: //g' \
    > $path2/${OUT}_nohitID.txt
  seqkit grep -n -w0 -f $path2/${OUT}_nohitID.txt -o $path2/HGS/${OUT}_nohit.fasta ${out5%.fastq}.fasta
  if [ $? -eq 0 ]; then touch $path2/.done.seqkit ; fi  

  HGS=`date '+%y%m%d%H%M%S'`
  echo -e "\n${HGS} ${OUT} subtracted FQ/FA stats" >> $OUTL2
  seqkit stats -j $THREADS \
    $out5 $out6 $path2/HGS/${OUT}_nohit.fasta >> $OUTL2
  #mv $path2/HGS/${OUT}_nohit.fasta $path2/${OUT}_nohit.fasta
  # the one in HGS directory will be refered to in downnstream
  cp $path2/HGS/${OUT}_nohit.fasta $path2/${OUT}_nohit.fasta
  echo "PATHDET_HGS_FINISHED" | tee -a $OUTL2

  if [ -f "$path2/.done.seqkit" ]; then
    touch $path2/.done
  else
    echo "an error occurred"
    echo "the analysis was interrupted at the HGS process"
    exit
  fi
}

function rapid() {
  if [ -f "$path_rapid/.done" ] ; then
    echo "Rapid process already finished, skipped"
    return
  fi

  mkdir -p $path_rapid/KK2

  OUTP=$path_rapid/${OUT}_pre-report.csv
  OUTLr=$path_rapid/${OUT}_log_Rapid.txt

  ##PATHOGEN detection :kraken2
  input=$path2/HGS/${OUT}_nohit.fasta
  kraken2 \
    --report $path_rapid/KK2/${OUT}_p1.kreport \
    --unclassified-out $path_rapid/KK2/${OUT}_unclass2.fastq \
    --classified-out $path_rapid/KK2/${OUT}_class2.fastq \
    --db ${PATHDB1} \
    --confidence 0.5 \
    --output $path_rapid/KK2/${OUT}_kraken_p1.txt \
    --threads $THREADS \
    $input 

  # WCA1: the read number at the beggining
  WCA1=`cat $path1/${OUT}_log_p1.txt | awk 'NR==3 {print $4}' | tr -d ","`
  WCA1=$(( WCA1 * 2 ))
  # RMDA2: the number of reads remained after trimming
  RMDA2=`cat $path1/${OUT}_log_p1.txt | awk 'NR==8 {print $4}' | tr -d ","`
  RMDA2=$(( RMDA2 * 2 ))
  # HG2: the number of reads remained after human genome subtraction
  HG2=`cat $path2/${OUT}_log_p2.txt | awk 'NR==8 {print $4}' | tr -d ","`
  PRstd=`echo "scale=0; ${HG2} * 1000000 / ${RMDA2}" | bc`

  KK2S=`date '+%y%m%d%H%M%S'`
  PKK1=`cat $path_rapid/KK2/${OUT}_class2.fastq | wc -l`
  PKK2=`expr ${PKK1} \/ 4`
  UKK1=`cat $path_rapid/KK2/${OUT}_unclass2.fastq | wc -l`
  UKK2=`expr ${UKK1} \/ 4`
  KKAR=`echo "scale=4; ${PKK2} / ${HG2} * 100" | bc`
  UCR=`echo "scale=4; ${UKK2} / ${WCA1} * 100" | bc`
  echo "${KK2S} kraken2 finished: ${PKK2}reads (${KKAR}% of karaken2 target)" > $OUTLr
  echo -e "\t* un-classified read: ${UKK2}reads (${UCR}%)" >> $OUTLr

  ##Prompt report
  echo -e "Sample_ID\t${OUT}" > $OUTP
  echo -e "Total_read(reads)\t${WCA1}" >> $OUTP
  echo -e "Available_read(reads)\t${RMDA2}" >> $OUTP
  echo -e "Non-human_read(reads)\t${HG2}" >> $OUTP
  echo -e "Non-human_read(RPM)\t${PRstd}" >> $OUTP

  mkdir -p $path_rapid/${OUT}_pRep
  ktImportTaxonomy \
    -q 2 -t 3 -s 4 \
    -o $path_rapid/${OUT}_pRep/${OUT}.plot.html \
    $path_rapid/KK2/${OUT}_kraken_p1.txt 
  if [ $? -eq 0 ]; then touch $path_rapid/.done.ktit ; fi

  grep -v 'unclassified' $path_rapid/KK2/${OUT}_p1.kreport | cut -f 3,5 | awk -F, '$1 > 1{ print $0 }' > $path_rapid/KK2/${OUT}_p2.kreport
  cat $path_rapid/KK2/${OUT}_p2.kreport \
    | awk '{ tax[$2] += $1; } END { for (i in tax) print tax[i] "\t" i; }' \
    | sort -k 1nr,1 \
    | taxonkit lineage -i 2 -o $path_rapid/KK2/${OUT}_tax.kreport
  cat $path_rapid/KK2/${OUT}_tax.kreport \
    | awk -v "r=${RMDA2}" '{print $1/r*1000000"\t"$0}' | cut -f 1,3- > $path_rapid/KK2/${OUT}_tax2.kreport

  # C. acnes is NOT removed

  #Create CSV file for each taxonomic rank
  mkdir -p $path_rapid/CSV
  #Kingdom
  rank=k
  taxonkit reformat -j $THREADS -i 3 -f {$rank} $path_rapid/KK2/${OUT}_tax2.kreport \
    | cut -f 1,4- | sed 's/ /_/g' | awk '{ tax[$2] += $1; } END { for (i in tax) print i "\t" tax[i]; }' \
    | sort -k 2nr,2 \
    | (s=$(cat)&&awk 'NR==FNR{a=a+$2;next}{p=$2/a;print $0"\t"p}' <(echo -e "$s") <(echo -e "$s")) \
    > $path_rapid/CSV/${OUT}_mp${rank^^}.csv
  taxonkit reformat -j $THREADS -i 3 -f {$rank} $path_rapid/KK2/${OUT}_tax.kreport \
    | cut -f 1,4- | sed 's/ /_/g' | awk '{ tax[$2] += $1; } END { for (i in tax) print i "\t" tax[i]; }' \
    | sort -k 2nr,2 \
    > $path_rapid/CSV/${OUT}_rp${rank^^}.csv
  paste $path_rapid/CSV/${OUT}_rp${rank^^}.csv $path_rapid/CSV/${OUT}_mp${rank^^}.csv | cut -f 1,2,4,5 \
    > $path_rapid/${OUT}_pRep/${OUT}_p${rank^^}1.csv

  #Phylum, Family, Genus, Species
  for rank in p f g s ; do
    taxonkit reformat -j $THREADS -i 3 -f {$rank} $path_rapid/KK2/${OUT}_tax2.kreport \
      | cut -f 1,4- | sed 's/ /_/g' | awk '{ tax[$2] += $1; } END { for (i in tax) print i "\t" tax[i]; }' \
      | sort -k 2nr,2 \
      | grep -v -e '^\s' \
      | (s=$(cat)&&awk 'NR==FNR{a=a+$2;next}{p=$2/a;print $0"\t"p}' <(echo -e "$s") <(echo -e "$s")) \
      > $path_rapid/CSV/${OUT}_mp${rank^^}.csv
    taxonkit reformat -j $THREADS -i 3 -f {$rank} $path_rapid/KK2/${OUT}_tax.kreport \
      | cut -f 1,4- | sed 's/ /_/g' | awk '{ tax[$2] += $1; } END { for (i in tax) print i "\t" tax[i]; }' \
      | sort -k 2nr,2 \
      | grep -v -e '^\s' \
      > $path_rapid/CSV/${OUT}_rp${rank^^}.csv
    paste $path_rapid/CSV/${OUT}_rp${rank^^}.csv $path_rapid/CSV/${OUT}_mp${rank^^}.csv \
      | cut -f 1,2,4,5 \
      > $path_rapid/${OUT}_pRep/${OUT}_p${rank^^}1.csv
  done

  #Calculate parameters(RIV,DiversityIndex D,H)
  PF1=`head -n3 $path_rapid/${OUT}_pRep/${OUT}_pF1.csv`
  PG1=`head -n3 $path_rapid/${OUT}_pRep/${OUT}_pG1.csv`
  PS1=`head -n3 $path_rapid/${OUT}_pRep/${OUT}_pS1.csv`

  PDF=`cut -f3 $path_rapid/CSV/${OUT}_mpF.csv | perl -nle '$D += $_*$_; print 1-$D if(eof);'`
  PDG=`cut -f3 $path_rapid/CSV/${OUT}_mpG.csv | perl -nle '$D += $_*$_; print 1-$D if(eof);'`
  PDS=`cut -f3 $path_rapid/CSV/${OUT}_mpS.csv | perl -nle '$D += $_*$_; print 1-$D if(eof);'`

  PHF=`cut -f3 $path_rapid/CSV/${OUT}_mpF.csv | perl -nle '$H += $_*log($_); print $H*(-1) if(eof);'`
  PHG=`cut -f3 $path_rapid/CSV/${OUT}_mpG.csv | perl -nle '$H += $_*log($_); print $H*(-1) if(eof);'`
  PHS=`cut -f3 $path_rapid/CSV/${OUT}_mpS.csv | perl -nle '$H += $_*log($_); print $H*(-1) if(eof);'`

  echo -e "\n[PROMPT_REPORT]" >> $OUTP
  echo -e ">Diversity_Index" >> $OUTP
  echo -e "Simpson_Index(family)\t${PDF}" >> $OUTP
  echo -e "Simpson_Index(genus)\t${PDG}" >> $OUTP
  echo -e "Simpson_Index(species)\t${PDS}" >> $OUTP
  echo -e "\nShannon_Index(family)\t${PHF}" >> $OUTP
  echo -e "Shannon_Index(genus)\t${PHG}" >> $OUTP
  echo -e "Shannon_Index(species)\t${PHS}" >> $OUTP
  echo -e "\n>Top_3hit_pathogen(family) Name,reads,RPM,RA\n${PF1}" >> $OUTP
  echo -e "\n>Top_3hit_pathogen(genus) Name,reads,RPM,RA\n${PG1}" >> $OUTP
  echo -e "\n>Top_3hit_pathogen(species) Name,reads,RPM,RA\n${PS1}" >> $OUTP

  for i in K P F G S ; do
    echo "Microorganisms,reads,RPM,RA" > $path_rapid/${OUT}_pRep/${OUT}_p$i.csv
    cat $path_rapid/${OUT}_pRep/${OUT}_p${i}1.csv >> $path_rapid/${OUT}_pRep/${OUT}_p$i.csv
    rm -f $path_rapid/${OUT}_pRep/${OUT}_p${i}1.csv
  done

  TIMEC=`date +%s`
  PT1=`expr ${TIMEC} \- ${TIMEA}`
  H1=`expr ${PT1} \/ 3600`
  PT1=`expr ${PT1} \% 3600`
  M1=`expr ${PT1} \/ 60`
  S1=`expr ${PT1} \% 60`
  BRKS=`date '+%y%m%d%H%M%S'`
  echo "${BRKS} PROMPT REPORT finished (${H1}:${M1}:${S1})" >> $OUTLr

  echo "PATHDET_Rapid_FINISHED"

  if [ -f "$path_rapid/.done.ktit" ]; then
    touch $path_rapid/.done
  else
    echo "an error occurred"
    echo "the analysis was interrupted at the Rapid process"
    exit
  fi
}

function tidyup() {
  intermediates=(
    $path1/${OUT}_trimmed_R1.fastq.gz
    $path1/${OUT}_trimmed_R2.fastq.gz
    $path2/HGS
    $path3/${OUT}_FA
    $path3/${OUT}_BLAST
    $path3/${OUT}_hgs.fasta
    $path4/${OUT}_CSV
    $path4/${OUT}_CSV2
  )
  if [ -f "$path5/.done" ] && [ "$DEBUG" == "False" ]; then
    rm -fr ${intermediates[@]}
  fi
}



# import common functions
. $CONDA_PREFIX/bin/p03_blast.sh
. $CONDA_PREFIX/bin/p04_krona.sh
. $CONDA_PREFIX/bin/p05_map.sh

# call workflows
functions=(p01_qc p02_hgs rapid p03_blast p04_krona p05_map)
for func in ${functions[@]} ; do $func ; done

tidyup
