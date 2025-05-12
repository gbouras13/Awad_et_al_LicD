# Awad_et_al_LicD
Repository to hold scripts for Awad et al 

This repository holds the code commands to reproduce the bioinformatics analyses in Awad et al. It is not intended to be fully one-click end-to-end reproducible (I lack the time for that) but should be pretty straightforward to follow for users with some knowledge of git, bash and conda. Please make an issue if anything is unclear or if you have questions.

All software was installed with conda using the bioconda Linux installation with version mentioned unless otherwise specificed. 

* To get started with this repository

```bash
git clone https://github.com/gbouras13/Awad_et_al_LicD.git
cd Awad_et_al_LicD
```

# Downloading Raw reads

* To download all raw sequencing reads from the SRA, I use [fastq-dl](https://github.com/rpetit3/fastq-dl)

```bash
mkdir sra_fastqs
cd sra_fastqs
fastq-dl --accession PRJNA1260167 --provider SRA
cd ..
```

* This will download each Run accession (i.e. SRRxxxxxxxx) FASTQ files. See `SRA_metadata.csv` in this repository for more information linking the Run accession to the isolate name and BioSample ID.

# Bacterial Isolate Analyses

* All relevant data can be found in `EF_Bacteria`

## Assembly

* These were done in 2 [hybracter](https://github.com/gbouras13/hybracter) commands with v0.7.3. 
* They were split because we only had long-reads for DFI-250 and DFI-266, whereas we had hybrid for the rest.
* The long-reads were nearly Q20 quality, so I didn't use Medaka (see [this blog](https://rrwick.github.io/2023/10/24/ont-only-accuracy-update.html) for justification)
* This assumes you have installed hybracter and downloaded the databases using `hybracter install`.

#### For all isolates other than DFI-250 and DFI-266

```bash
hybracter hybrid --datadir "../sra_fastqs" -i hybracter_hybrid_sheet.csv  -o hybracter_out_non_dfi -t 32 
```

#### For DFI-250 and DFI-266

```bash
hybracter long --datadir "../sra_fastqs" -i hybracter_long_sheet.csv -t 32 -o hyracter_out_dfi  --no_medaka  --min_quality  15
```

## Annotation

* [Bakta](https://github.com/oschwengers/bakta) v1.9.4 was used to annotate the assemblies
* The hybracter output from above can be found in `ef_isolate_genomes` in this directory, with combined, chromosome and plasmid subdirectories.
* It was not uploaded to this repository as it is too large and not especially relevant other than the genomes
* You can find the annotations in `batka_output`

```bash
# this takes ~70GB and a few hours (in Australia)
bakta_db  download -o bakta_dbs/full_db

mkdir -p bakta_output
for file in ef_isolate_genomes/combined/*.fasta; do
    base_name=$(basename "$file" ".fasta")
    bakta --db bakta_dbs/full_db/db/ --threads 32 --output bakta_output/${base_name}.fasta --complete --prefix ${base_name} genomes/combined/${base_name}.fasta
done
```

## Panaroo pan-genome

* Panaroo v1.5.1 was used to create the pan-genome

```bash
mkdir -p gffs 
cp bakta_output/*/*gff3 gffs
panaroo -i gffs/*  -o ./panaroo/ --clean-mode strict -a core --aligner clustal --core_threshold 0.98 -t 10 --remove-invalid-genes
```

## Mashtree

* Mashtree v1.4.6 was used

```bash
mkdir -p mashtree
mashtree ef_isolate_genomes/combined/*fasta > mashtree/mashtree.dnd
```

## AMRFinderPlus

* AMRFinderPlus v4.0.3 was used

```bash
mkdir -p amrfinderplus
for file in ef_isolate_genomes/combined/*.fasta; do
base_name=$(basename "$file" ".fasta")
amrfinder -p ${base_name}.fasta/${base_name}.faa    -O Enterococcus_faecalis    --plus -o amrfinderplus/${base_name}.txt

done
```

## Mob Suite

* Mob Suite v3.1.9 was used

```bash
mob_typer --multi --infile  ef_isolate_genomes/plasmids/all_plasmids_combined.fasta --out_file mobtyper/sample_mobtyper_results.txt
```

# Phage Isolates

## Sphae 

* To generate the long-read assemblies and Pharokka annotations, we used [sphae](https://github.com/linsalrob/sphae) v1.4.2

```bash
mkdir -p phage_nanopore_fastqs
cp sra_fastqs/SRR33468901.fastq.gz phage_nanopore_fastqs/Ef-P10-208.fastq.gz
cp sra_fastqs/SRR33468899.fastq.gz phage_nanopore_fastqs/Ef-P20-208.fastq.gz
cp sra_fastqs/SRR33468900.fastq.gz phage_nanopore_fastqs/Ef-P16-208.fastq.gz
sphae run --input phage_nanopore_fastqs --sequencing longread --output sphae_output -k --use-conda --conda-frontend mamba
```

* sphae makes a lot of outputs, so I have placed the assembles in `Phages/Fastas` and the Pharokka annotations in `Pharokka`
* Not that APTC-Efa.10 is Ef-P10-208, APTC-Efa.16 is Ef-P16-208 and APTC-Efa.20 is Ef-P20-208

## ColabFold and Protein Structure generation

* This is the most compute intensive and difficult step to reproduce, requiring a large server with more than 1TB of disk space and a GPU (or realistically, a HPC core).
* Protein structure predictions were generated for all predicted CDS from Pharokka (phanotate) with (ColabFold)[https://github.com/sokrypton/ColabFold] via a local server was installed with https://github.com/YoshitakaMo/localcolabfold 
* Specifically, all predictions were run on Pawsey's [Setonix](https://pawsey.org.au/systems/setonix/) using a custom ColabFold v1.5.5 container adapted for Setonix's AMD MI250x GPUs (thanks Sarah Beecroft as always https://quay.io/repository/sarahbeecroft9/colabfold?tab=tags - specifically the rocm6_cpuTF docker file).
* For MSA generation, both the Uniref30 and ColabFold's environment databases were used, without templates.
* For structure prediction, 5 models were generated with 3 recycles. The highest ranking model by pLDDT was chosen as the best model for each protein.
* **All MSAs, plddt and pae json files and top ranking predicted structures are available in `Phages/ColabFold`, along with the pLDDT pTM and length metadata for each protein in `plddt_ptm_len.tsv`** - These are included as tarballs for each phage

* To reproduce for Ef-P10

```bash
# on CPU node to generate MSAs in bulk
THREADS=256
DB_DIR="path_to_where_you_downloaded_the_colabfold_databases"
MSADIR="EF-P10-msa"
PREDIR="EF-P10-predictions"
colabfold_search --threads $THREADS EF-10-phanotate.faa  $DB_DIR  $MSADIR

# on GPU node to generate the structural predictions
colabfold_batch --num-models 5  $MSADIR  $PREDIR

# to select the top ranking model and generate the metadata etc
# needs bioconda to be installed (best to use the same conda env as pharokka)
python process_colabfold_output.py -i EF-10-phanotate.faa -p $PREDIR -o EF-P10-processed -f
```

## Phold

* [Phold](https://github.com/gbouras13/phold) v0.2.0 was used to annotate these phages
* You can find all structures in `Phages/Phold/*_filtered_structures`

```bash
cd Phages
phold compare -t 32 -d phold_db --structure_dir Phold/Ef-P10-208_filtered_structures/ --structures --input Pharokka/Ef-P10-208/pharokka.gbk -o Ef-P10-208
phold compare -t 32 -d phold_db --structure_dir Phold/Ef-P16-208_filtered_structures/ --structures --input Pharokka/Ef-P16-208/pharokka.gbk -o Ef-P16-208
phold compare -t 32 -d phold_db --structure_dir Phold/Ef-P20-208_filtered_structures/ --structures --input Pharokka/Ef-P20-208/pharokka.gbk -o Ef-P20-208
```

## Deposcope and Depolymerase Analysis

* [Deposcope](https://github.com/dimiboeckaerts/DepoScope) is a little more difficult to run and not available on bioconda (but it still nonetheless open source and quite straightforward to install!), so to create a conda environment with the necessary dependencies:

```bash
cd Phages/DepoScope
# download the models
wget 'https://zenodo.org/records/10957073/files/Deposcope.esm2_t12_35M_UR50D.2203.full.model?download=1' -O 'DepoDetection.model'
wget 'https://zenodo.org/records/10957073/files/esm2_t12_finetuned_depolymerases.zip?download=1' -O 'esm2_t12_35M_UR50D-finetuned-depolymerase.labels_4.zip'
unzip esm2_t12_35M_UR50D-finetuned-depolymerase.labels_4.zip
# create conda env - if no NVIDIA gpu, remove '=*=cuda*
conda create -n deposcope transformers biopython pandas tqdm pytorch=*=cuda*
conda activate deposcope
```
* To run deposcope with protein input, I put together a little script `run_deposcope.py` available in this repository in `Phages/Deposcope` modifying the Deposcope-provided [Google Colab notebook](https://colab.research.google.com/drive/1A2XJ_oUtlmIfU3XXmev5dzJUNxqR6VV9?usp=sharing#scrollTo=d8bbf189) to run it only on already gene-called protein sequences:

```bash
for PHAGE in (Ef-P10-208 Ef-P16-208 Ef-20-208); do
python run_deposcope.py -i ../Pharokka/$PHAGE/phanotate.faa -o $PHAGE -e esm2_t12_35M_UR50D__fulltrain__finetuneddepolymerase.2103.4_labels/checkpoint-2255  -d DepoDetection.model  -f
done
```

* To look at the depolymerase domain in more detail see `Depolymerase_Domain`. You can find the the three likely depolymerase containing protein in 1 file (domains.fasta). The MSA was made with MAFFT v7.526, while the MSA was visualised with rich-msa

```bash
cd Phages/Depolymerase_Domain
#MAFFT v7.526 
mafft domains.fasta > domains.msa
mamba create -n msa
conda activate msa
pip install pytrimal
pip install rich-msa
pip install biopython
pip install pandas
python make_alignment_plots.py

```

## Lovis4u

* Phage synteny plots were created with [Lovis4u](https://github.com/art-egorov/lovis4u) v0.0.7 (install from pip not conda)

```bash
cd Phages/lovis4u
# get the initial
lovis4u --gb gbks  -hl -fv-off -o lovis4u_output --set-category-colour 
# change the coordinates manually then
lovis4u --gb gbks  -hl -fv-off -o lovis4u_output_structural_module --set-category-colour --locus-annotation-file locus_annotation_table_structural.tsv
```


# BIM Analysis 

## Assemble the ATCC 700802 reference genome

* Used Trycycler v0.5.4 for this, along with Flye v2.9.3, Raven v1.8 and Canu v2.2 for the consituent assemblies and Filtlong v0.2.1 and porechop ABI v0.5.0 for QC
* The intermediate Raven assemblies were actually discarded as they were of low contiguity/quality

```bash


filtlong  --min_mean_q 10 --min_length 1000 sra_fastqs/SRR33468901.fastq.gz | pigz > temp.fastq.gz 
porechop_abi -i temp.fastq.gz  -o ATCC_700802_filt_trim.fastq.gz -t 32 
trycycler subsample --reads ATCC_700802_filt_trim.fastq.gz --out_dir read_subsets

threads=32
mkdir -p assemblies
genome_size="3000000"

flye --nano-hq read_subsets/sample_01.fastq --threads "$threads" --out-dir assembly_01 && cp assembly_01/assembly.fasta assemblies/assembly_01.fasta && cp assembly_01/assembly_graph.gfa assemblies/assembly_01.gfa && rm -r assembly_01
flye --nano-hq read_subsets/sample_04.fastq --threads "$threads" --out-dir assembly_04 && cp assembly_04/assembly.fasta assemblies/assembly_04.fasta && cp assembly_04/assembly_graph.gfa assemblies/assembly_04.gfa && rm -r assembly_04
flye --nano-hq read_subsets/sample_07.fastq --threads "$threads" --out-dir assembly_07 && cp assembly_07/assembly.fasta assemblies/assembly_07.fasta && cp assembly_07/assembly_graph.gfa assemblies/assembly_07.gfa && rm -r assembly_07
flye --nano-hq read_subsets/sample_10.fastq --threads "$threads" --out-dir assembly_10 && cp assembly_10/assembly.fasta assemblies/assembly_10.fasta && cp assembly_10/assembly_graph.gfa assemblies/assembly_10.gfa && rm -r assembly_10

raven --threads "$threads" --disable-checkpoints --graphical-fragment-assembly assemblies/assembly_03.gfa read_subsets/sample_03.fastq > assemblies/assembly_03.fasta
raven --threads "$threads" --disable-checkpoints --graphical-fragment-assembly assemblies/assembly_06.gfa read_subsets/sample_06.fastq > assemblies/assembly_06.fasta
raven --threads "$threads" --disable-checkpoints --graphical-fragment-assembly assemblies/assembly_09.gfa read_subsets/sample_09.fastq > assemblies/assembly_09.fasta
raven --threads "$threads" --disable-checkpoints --graphical-fragment-assembly assemblies/assembly_12.gfa read_subsets/sample_12.fastq > assemblies/assembly_12.fasta


for i in 02 05 08 11; do
    canu -p canu -d canu_temp_${i} -fast genomeSize="$genome_size" useGrid=false maxThreads="$threads" -nanopore read_subsets/sample_"$i".fastq
    python canu_trim.py canu_temp_${i}/canu.contigs.fasta > assemblies/assembly_"$i".fasta
    #rm -rf canu_temp_${i}
done
```

* cluster

```bash
trycycler cluster --assemblies assemblies/*.fasta --reads ATCC_700802_filt_trim.fastq.gz   --out_dir trycycler

```

* reconcile - manual curation
* There were clearly 4 clusters (corresponding to the expected chromosome + 3 plasmids)

```bash
trycycler reconcile --reads ATCC_700802_filt_trim.fastq.gz  --cluster_dir trycycler/cluster_001 --max_add_seq 20000 
# the flye plasmids weren't great
rm trycycler/cluster_002/1_contigs/*contig*.fasta 
trycycler reconcile --reads ATCC_700802_filt_trim.fastq.gz  --cluster_dir trycycler/cluster_002
trycycler reconcile --reads ATCC_700802_filt_trim.fastq.gz  --cluster_dir trycycler/cluster_003
trycycler reconcile --reads ATCC_700802_filt_trim.fastq.gz  --cluster_dir trycycler/cluster_004

trycycler dotplot --cluster_dir trycycler/cluster_001
trycycler dotplot --cluster_dir trycycler/cluster_002
trycycler dotplot --cluster_dir trycycler/cluster_003
trycycler dotplot --cluster_dir trycycler/cluster_004
```

* MSA

```bash
trycycler msa --cluster_dir trycycler/cluster_001
trycycler msa --cluster_dir trycycler/cluster_002
trycycler msa --cluster_dir trycycler/cluster_003
trycycler msa --cluster_dir trycycler/cluster_004
```

* partition

```bash
trycycler partition --reads ATCC_700802_filt_trim.fastq.gz --cluster_dirs trycycler/cluster_*
```

* consensus 

```bash
trycycler consensus --cluster_dir trycycler/cluster_001
trycycler consensus --cluster_dir trycycler/cluster_002
trycycler consensus --cluster_dir trycycler/cluster_003
trycycler consensus --cluster_dir trycycler/cluster_004
cat trycycler/cluster_*/7_final_consensus.fasta > trycycler/consensus.fasta
cp trycycler/consensus.fasta ATCC_700802_trycycler.fasta
```

## Annotating 

* Bakta v1.9.4

```bash
bakta --db bakta_dbs/full_db/db/ --threads 32 --output ATCC_700802_bakta --complete --prefix ATCC_700802 ATCC_700802_trycycler.fasta
```

## Align BIM data to assembled genome

* Note - need to include the phages in the reference or else some reads (to proteins annotated as LysM and other phage proteins) are shared, and as there are potentially some contaminant phage reads in the BIM FASTQs, these will be misaligned and give false positive variant calls

e.g. this shows lots of hits

```bash
mmseqs easy-search ../Phages/Phold/Ef-P10-208/Ef-P10-208_aa.fasta ATCC_700802_bakta/ATCC_700802.faa phage_vs_bacteria.m8 tmp
```

* minimap2 v2.28-r1209 was used, along with samtools v1.15.1
* All BIMs plus the reference ATCC 700802 readset was used as a baseline for variants

```bash

cat ATCC_700802_trycycler.fasta ../Phages/Fastas/* > ATCC_700802_trycycler_plus_phages_reference.fasta

cd BIMS

cp ../sra_fastqs/SRR33468901.fastq.gz EF_original.fastq.gz
cp ../sra_fastqs/SRR33469033.fastq.gz BIM_10.fastq.gz
cp ../sra_fastqs/SRR33469032.fastq.gz BIM_16.fastq.gz
cp ../sra_fastqs/SRR33469031.fastq.gz BIM_20.fastq.gz
cp ../sra_fastqs/SRR33469030.fastq.gz BIM_1610.fastq.gz
cp ../sra_fastqs/SRR33469029.fastq.gz BIM_2010.fastq.gz
cp ../sra_fastqs/SRR33469028.fastq.gz BIM_2016.fastq.gz

# alignment

BAMS="BAMS"
mkdir -p $BAMS

for prefix in BIM_10 BIM_1610 BIM_16 BIM_2010 BIM_20 BIM_2016 EF_original
do
    minimap2 -ax map-ont -t 31 ATCC_700802_trycycler_plus_phages_reference.fasta ${prefix}.fastq.gz > $BAMS/${prefix}_aln.sam
done

# sorting and conversion to BAM

for prefix in BIM_10 BIM_1610 BIM_16 BIM_2010 BIM_20 BIM_2016  EF_original
do

echo $prefix
samtools view -S -b $BAMS/${prefix}_aln.sam >  $BAMS/${prefix}_aln.bam
samtools sort $BAMS/${prefix}_aln.bam -o $BAMS/${prefix}_aln_sorted.bam
samtools index $BAMS/${prefix}_aln_sorted.bam
done

```
## Variant calling

* Note that the Clair3 output is available in `BIMS/CLAIR3_OUT`, but the aligned BAM files are too large to host on github (and not very hard to reproduce)
* I used Clair3 per [Hall et al](https://doi.org/10.7554/eLife.98300.2) as the best performing tool
* ONT basecaller model was R9.4.1 v3.3 for this data - therefore use that model for Clair3
* To download the model

```bash
mkdir clair3_model
cd clair3_model
wget "https://cdn.oxfordnanoportal.com/software/analysis/models/clair3/r941_prom_sup_g5014.tar.gz"
tar -xzf r941_prom_sup_g5014.tar.gz 
```

* Then run clair3

```bash
samtools faidx ATCC_700802_trycycler_plus_phages_reference.fasta 

for prefix in EF_original BIM_10 BIM_1610 BIM_16 BIM_2010 BIM_20 BIM_2016
do

aln="$BAMS/${prefix}_aln_sorted.bam"
ref="ATCC_700802_trycycler_plus_phages_reference.fasta"
model_path="clair3_model/r941_prom_sup_g5014"
outdir="clair3_${prefix}"
sample="${prefix}"

run_clair3.sh \
    --bam_fn="$aln" \
    --ref_fn="$ref" \
    --threads=32 \
    --platform="ont" \
    --model_path="$model_path" \
    --output="$outdir" \
    --sample_name="$sample" \
    --include_all_ctgs \
    --haploid_precise \
    --no_phasing_for_fa \
    --enable_long_indel
done
```
* To get only the highly confident variant calls (allele frequency > 0.7)

```bash
CLAIR3_OUT="CLAIR3_OUT"
mkdir -p $CLAIR3_OUT
bcftools view -i 'FORMAT/AF > 0.7' clair3_EF_original/merge_output.vcf.gz | bgzip  > $CLAIR3_OUT/clair3_EF_original_high_qual_errors.vcf.gz
tabix -p vcf $CLAIR3_OUT/clair3_EF_original_high_qual_errors.vcf.gz

for prefix in BIM_10 BIM_1610 BIM_16 BIM_2010 BIM_20
do

bcftools view -i 'FORMAT/AF > 0.7' clair3_${prefix}/merge_output.vcf.gz | bgzip  > $CLAIR3_OUT/clair3_${prefix}_high_qual.vcf.gz
tabix  $CLAIR3_OUT/clair3_${prefix}_high_qual.vcf.gz
bcftools isec -C  $CLAIR3_OUT/clair3_${prefix}_high_qual.vcf.gz $CLAIR3_OUT/clair3_EF_original_high_qual_errors.vcf.gz -o $CLAIR3_OUT/clair3_${prefix}_high_qual_filtered.vcf

done
```

* The final variant calls with manually curated notes can be found in `variant_table_summary.xlsx`
* You can also find all the BAM pileups as supplementary figures in this manuscript to verify the existence of the variants
* Code to reproduce these figures can be found in `BIMS/gviz`