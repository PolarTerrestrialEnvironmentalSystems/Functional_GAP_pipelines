<div align="center">
  <img src="https://jobs.awi.de/pubhtml/logo.gif" alt="Your Logo" width="350">
</div>

<h1 align="center">
  <div style="background-color: #008ACC; padding: 10px; border-radius: 10px; color: #FFFFFF;">
    <b>Functional annotation of eukaryotic genes from sedimentary ancient DNA</b>
    <br>
  </div>
  <div style="background-color: #008ACC; padding: 10px; border-radius: 10px; margin-top: 10px;">
  </div>
</h1>

![gap_pipeline_final_v2](https://github.com/ucabuk/functional_eukaryotes_lama/assets/36818676/f1d1e887-cd32-450b-bf43-d1749af70ee0)


This repository includes bioinfornatic scripts for benchmarking gene annoation pipelines for sedimentary ancient DNA shotgun data. Description of the data analysis provided in this repository and the results of benchmarking gene annottaion pipelines are published in Cabuk et al. XXXX. The bioinformatic computations were done on HPC (High Performance Clustering).

A few helping scripts were used to perform data manipulation in between the main bioinformatic scripts, such  helper scripts are provided below.

Here are the full list of the software tools applied for gene annotation benchmarking :


| Tool         | Version   | Link |
|--------------|-----------|------| 
| FastQC       | 0.11.9    | https://github.com/s-andrews/FastQC |
| BBTools      | 38.87     | https://jgi.doe.gov/data-and-tools/software-tools/bbtools |
| fastp        | 0.20.1    | https://github.com/OpenGene/fastp |
| MEGAHIT      | 1.2.9     | https://github.com/voutcn/megahit |
| MetaEuk      | 6.a5d39d9 | https://github.com/soedinglab/metaeuk |
| Pydamage     | 0.72      | https://github.com/maxibor/pydamage |
| Prodigal     | 2.6.3     | https://github.com/hyattpd/Prodigal|
| Tiara        | v1.0.3    | https://github.com/ibe-uw/tiara|
| CD-HIT       | 4.8.1     | https://sites.google.com/view/cd-hit |
| salmon       | v1.10.1   | https://github.com/COMBINE-lab/salmon |
| eggNOG-mapper| 2.1.2     | https://github.com/eggnogdb/eggnog-mapper | 
| MMseqs2      | 14.7e284  | https://github.com/soedinglab/MMseqs2|
| QUAST        | 5.2.0     | https://github.com/ablab/quast| 
| Seqtk        | v1.4      | https://github.com/lh3/seqtk|
| Python        | v3.6.8    | [https://www.r-project.org](https://www.python.org/) |
| R            | 4.2.2     | https://www.r-project.org |
| R studio     | 2022.12.0+353| https://posit.co/download/rstudio-desktop |
| SeqKit       | 2.5.1     |https://github.com/shenwei356/seqkit|

-------------------------------------------

## Databases of Gene Annotation pipeline
### Metaeuk Gene Prediction - Uniref50 database

Description : Uniref50 database is used to predict potential eukaryotic genes in the assembly. To build the database, mmseqs2 was used with the following script. 

Run: ```mmseqs databases UniRef50 uniref50 tmp_uniref50```

### MMseqs2 - NCBI NR

Description : We used NCBI NR database to assign the eukaryotic Orthologs proteins into detailed taxonomic information. To build the database, following script was used.

Run: ```mmseqs databases NR NR_ncbi tmp_nr```

## Data analyses for the Gene Annotation pipelines

### Step1 Pre-processing of raw sequencing data

Description : Initially, the quality of reads before and after processing was evaluated using FastQC (v0.11.9). Subsequently, fastp (v0.20.1) was utilized to remove low-quality reads, adapter sequences, and merge paired-end reads. The removal of PCR duplicates from the shotgun sequence data was performed using clumpify from BBtools (v38.87).

Input File: Raw data from the data source, "Barbara's upload link"

Run   : ```step1_preprocessing.sl```

### Step1.2 Taxonomic assignment with Kraken2 

Description : The paired and merged reads from step1_preprocessing.sl were taxonomically classified using Kraken2 against the nt database at a confidence threshold of 0.8 and a k-mer minimizer of 31nt.

Run   : ```step1.2_kraken_nt22.sl```

### Step2 Metagenome assembly

Description : To improve the assembly quality, processed reads, then, underwent error-correction using tadpole from BBtools (v38.87, (Bushnell, 2014). Following this, quality checked error-corrected merged reads and paired reads from the samples are assembled into contigs using MEGAHIT (v1.2.9).

Run   : ```step2_de_novo_assembly.sl```

Output : ```out.megahit/${sample_name}/final.contigs.fa```

### QUAST Assessment (Basic Assembly Statistics)

Description: Here we assess our assemblies using QUAST. We used ```final.contigs.fa``` as input from MEGAHIT output and it will generate a folder ```quast_out_all_assemblies``` that stores assembly basic statistics.

Input : ```out.megahit/[all_samples_assemblies].final.contigs.fa```

Run : ```quast out.megahit/[all_samples_assemblies]*.fa -o quast_out_all_assemblies```

Output : ```report.tsv ```

Once getting report.tsv file, n50, largest contig, and number of contigs columns can be retrieved and combined with the read count information. To get the publication figure:

Input : ```R_script/assembly/lake_lama_quast_report.csv```
Run : ```R_script/assembly/assembly_stats_plot.R```

### Step3 Gene Prediction

**prokGAP:**

Description : In the prodigal outputs, we have ${sample_name}.fna, ${sample_name}.faa. The ${sample_name}.faa file contains predicted proteins, while ${sample_name}.codon.fna contains the nucleotide sequences of translated amino acids sequences

Input : ```out.megahit/${sample_name}/final.contigs.fa```

Run : ```step3_prodigal_gene_pred.sl```

Output : ```out.prodigal/${sample_name}/${sample_name}.fna, out.prodigal/${sample_name}/${sample_name}.faa, out.prodigal/${sample_name}/${sample_name}.gff ```

Description  : Before running step 4, we create a redundant gene database. To do this, we merge all *.fna and *.faa outputs into one file for salmon input.

bash command : ```cat out.prodigal/*.fna > out.prodigal/redundant_PROKGAP_pCDS.fna```

bash command : ```cat out.prodigal/*.faa > out.prodigal/redundant_PROKGAP_protein.faa```

Here, prodigal IDs should be fixed in these files using ```sed``` script before CD-HIT clustering step.

bash command : ```sed -e "s/ //g" < out.prodigal/redundant_PROKGAP_pCDS.fna ```

bash command : ```sed -e "s/ //g" < out.prodigal/redundant_PROKGAP_protein.faa```

**eukGAP:**

Description : In the metaeuk outputs, we have ${sample_name}.codon.fas, ${sample_name}.fas, ${sample_name}.gff, and ${sample_name}.tsv. The ${sample_name}.fas file contains predicted proteins, while ${sample_name}.codon.fas contains the nucleotide sequences of translated amino acid sequences.

Input : ```out.megahit/${sample_name}/final.contigs.fa```

Run   : ```step3_metaeuk_gene_pred.sl```

Output : ```out.metaeuk/${sample_name}/${sample_name}.codon.fas, out.metaeuk/${sample_name}/${sample_name}.fas, out.metaeuk/${sample_name}/${sample_name}.gff, out.metaeuk/${sample_name}/${sample_name}.tsv```

Description  : Here, redundant catalogs are created using the predicted protein and corresponding nucleotide from all samples.

bash command :

```cat out.metaeuk/*/*.codon.fas > out.metaeuk/redundant_EUKGAP_pCDS.fna```

```cat out.metaeuk/*/*.fas > out.metaeuk/redundant_EUKGAP_protein.faa```


Before running step 4, we create a redundant gene database. To do this, we merge all *.fas outputs into one file for the next step. We perform the same process for predicted protein outputs with .codon.fas extension files.

**Contig pre-classification-based GAP**


Description : We used Tiara deep learning-based contig classification tool to separate eukaryotic origin contigs. Subsequently, we used step3_metaeuk_gene_pred.sl for eukaryotic origin contigs and step3_prodigal_gene_pred.sl for the rest contigs.

Input : ```out.megahit/${sample_name}/final.contigs.fa```

Run   : ```step3_1_tiara_classification.sl```

Output : ```{OUTDIR}/${OUT_TIARA/eukaryota_{ID}_300.fna,  {OUTDIR}/${OUT_TIARA/prokaryotic_{ID}_300.bin.fna```

Description : Running step3_metaeuk_gene_pred.sl for eukaryotic bins (eukGAP)

Input : ```${sample_name}.eukaryotic_300.fna```

Run   : ```step3_metaeuk_gene_pred.sl```

Output : ```out.tiara/out.metaeuk/${sample_name}/${sample_name}.codon.fas, out.tiara/out.metaeuk/${sample_name}/${sample_name}.fas, out.tiara/out.metaeuk/${sample_name}/${sample_name}.gff, out.tiara/out.metaeuk/${sample_name}/${sample_name}.tsv```

bash command : 

```cat out.tiara/eukaryotes/out.metaeuk/*/*.codon.fas > out.tiara/out.metaeuk/redundant_preclass_EUKGAP_pCDS.fna```
   
```cat out.tiara/out.metaeuk/*/*.fas > out.tiara/out.metaeuk/redundant_preclass_EUKGAP_protein.faa```


Description : Running step3_prodigal_gene_pred.sl for prokaryotic bins (prokGAP)

Input : ```${sample_name}.prokaryotic_300.fna```

Run   : ```step3_prodigal_gene_pred.sl```

Output : ```out.tiara/out.prodigal/${sample_name}/${sample_name}.fna, out.tiara/out.prodigal/${sample_name}/${sample_name}.faa, out.tiara/out.prodigal/${sample_name}/${sample_name}.gff ```
 

bash command : 

```cat out.tiara/out.prodigal/*/*.fna > out.tiara/out.prodigal/redundant_preclass_prokGAP_pCDS.fna```
        
```cat out.tiara/out.prodigal/*/*.faa > out.tiara/out.prodigal/redundant_preclass_prokGAP_protein.faa```

Here, prodigal IDs should be fixed in these files using ```sed``` script before CD-HIT clustering step.

bash command : 

```sed -i "s/ //g" out.tiara/out.prodigal/redundant_preclass_prokGAP_pCDS.fna```

```sed -i "s/ //g" out.tiara/out.prodigal/redundant_preclass_prokGAP_protein.faa```

### Step4 Functional Gene/Protein Catalog 

Description : Functional redundant catalogs were generated from the prediction pipelines and present the functional diversity on gene and protein level aggregated from all sediment core samples. 

INPUTS (respectively): 

```out.metaeuk/redundant_EUKGAP_protein.faa, out.prodigal/redundant_PROKGAP_protein.faa, out.tiara/out.prodigal/redundant_preclass_prokGAP_protein.fna, out.tiara/out.metaeuk/redundant_preclass_EUKGAP_protein.faa```

Run : ```step4_cd_hit.sl```

Description: In the CD-HIT output, we extract gene IDs from the non-redundant gene catalog. These IDs are used to extract proteins with the same IDs from redundant_protein.fas.

OUTPUT (respectivel):  

```out.metaeuk/non_redundant_EUKGAP_protein.faa, out.prodigal/non_redundant_PROKGAP_protein.faa, out.tiara/out.prodigal/non_redundant_preclass_prokGAP_protein.faa, out.tiara/out.metaeuk/non_redundant_preclass_eukGAP_protein.faa```

Description: Here we extract their respective nucleotides for the gene abundance estimation with salmon. Example for 'redundant_PROKGAP_protein.faa' and should be implemented it on all redundant protein dataset from the pipelines

Bash command: ```grep '>' out.prodigal/redundant_PROKGAP_protein.faa | cut -c 2- > out.prodigal/non_redundant_PROKGAP_protein_to_pCDS.id```

Bash command: ```seqtk subseq out.prodigal/redundant_PROGAP_pCDS.fna non_redundant_EUKGAP_protein_to_pCDS.id > out.prodigal/non_redundant_PROKGAP_protein_to_pCDS.fna ```

### Step5 Gene Abundance Estimation

Description : Salmon was use for the quality-checked merged and paired reads of the samples to quantify the abundance of each gene in the non-redundant gene catalog 


INDEX FILE INPUTS: (respectively): 

```out.prodigal/non_redundant_PROGAP_protein_to_pCDS.fna, out.metaeuk/non_redundant_EUKGAP_protein_to_pCDS.fna, out.tiara/out.prodigal/non_redundant_preclass_prokGAP_protein_to_pCDS.fna, out.tiara/out.metaeuk/non_redundant_preclass_prokGAP_protein_to_pCDS.fna```

Run : ```step5_abundance_salmon.sl```

Description : The paired and merged gene counts were aggregated for each sample. 

Description : First, sum up all merged and read counts using the python script for each pipeline.

Run : ```python3 python_script/sum_up_qc_merge_paired.py```

After running the batch and python script, we merge the output from the samples, considering only the 'numreads'  and 'len' column in the outputs instead of TPM (transcripts per million).

Run : ```salmon quantmerge --quants directory_of_salmon_output -c numreads -o prokGAP_all_lake_lama_gene_quant.raw.count.qf```

Run : ```salmon quantmerge --quants directory_of_salmon_output -c len -o prokGAP_all_lake_lama_gene_quant.raw.count.len```

Description : Fix the 'all_lake_lama_gene_quant.raw.count.qf' and 'prokGAP_all_lake_lama_gene_quant.raw.count.len'

Run : ```python_scripts/fix_merge_len.py```

Description : To calculate the normalized gene count (NGC) over samples using raw counts, we used custom python script.

Input : ```prokGAP_all_lake_lama_gene_quant.raw.count.len```

Run: ```python3 python_script/ngc_calculation.py```

Output : ```prokGAP_all_lake_lama_gene_quant.CPM.fixed.tsv```

To get number table for the gene abundance, you will have lots of table in here for non-redundant, eukaryotes/prokaryotes gene subset, KOs subset. You need to merge the table as I already prepared in the folder; gene_catalog_visualization.r.:

Run: ```gene_catalog/gene_catalog_abundance_calc.R```

To visualization the gene abundance figure, see R_script/Gene_catalog.

Run: ```gene_catalog/gene_catalog_visualization.r```

### Step6 Functional annotation with EggNOG

Description : EggNOG is used to get the functional annotation of the proteins in the non-redundant protein datasets from the pipelines.

INPUTS (Respectively) :  

```out.metaeuk/non_redundant_EUKGAP_protein.faa,out.prodigal/non_redundant_PROKGAP_protein.faa,out.tiara/out.prodigal/non_redundant_preclass_prokGAP_protein.faa,out.tiara/out.metaeuk/non_redundant_preclass_eukGAP_protein.faa```

Run   : ```step6_eggnog_annnotation.sl```

OUTPUTS (Respectively) :  

```out.eggnog/non_redundant_EUKGAP_protein_eggNOG.emapper.annotations, out.eggnog/non_redundant_PROKGAP_protein_eggNOG.emapper.annotations, out.eggnog/non_redundant_preclass_prokGAP_protein_eggNOG.emapper.annotations, out.eggnog/non_redundant_preclass_eukGAP_protein_eggNOG.emapper.annotations```

Description : From the output of the eggNOG result, we select eukaryotic origin protein and extract them from the non-redundant protein file using their gene IDs. Example from prokGAP input. In that step there will be 8 datasets from all pipelines.

#for ProkGAP example

Bash command: ```grep "Eukaryota" out.eggnog/non_redundant_PROKGAP_protein_eggNOG.emapper.annotations | awk '{print $1} > out.eggnog/eggNOG_eukaryotes_OG_prokGAP.id```

Bash command: ```seqtk subseq out.eggnog/non_redundant_PROKGAP_protein.faa eggNOG_eukaryotes_OG_prokGAP.id > out.eggnog/out.eggnog/eggnog_eukaryote_OG_proteins_prokGAP.faa```

### Diversity and abundance of KEGG Orthologs 

Description : Protein sequences with KEGG Orthologs (KO) annotations were extracted from the EggNOG outputs and proteins associated with multiple KO identifiers to reduce ambiguous annotations were discarded using an R script. 

eukgap_input: ```out.eggnog/non_redundant_EUKGAP_protein_eggNOG.emapper.annotations, eukGAP_all_lake_lama_gene_quant.CPM.fixed.tsv```

prokgap_input: ```out.eggnog/non_redundant_PROKGAP_protein_eggNOG.emapper.annotations,prokGAP_all_lake_lama_gene_quant.CPM.fixed.tsv ```

preclass_eukgap_input: ```out.eggnog/non_redundant_preclass_prokGAP_protein_eggNOG.emapper.annotations,preclass_eukGAP_all_lake_lama_gene_quant.CPM.fixed.tsv```

preclass_prokgap_Input: ```out.eggnog/non_redundant_preclass_eukGAP_protein_eggNOG.emapper.annotations,preclass_prokGAP_all_lake_lama_gene_quant.CPM.fixed.tsv```

Run: ```R_scripts/kegg_diversity.R```

Output : ```prokGAP_final_df_row_counts_lama.csv,eukGAP_final_df_row_KO_counts_lama.csv ,preclass_eukGAP_final_df_row_KO_counts_lama.csv, preclass_prokGAP_final_df_row_KO_counts_lama.csv, observed_count.png```

Description: For KEGG pathway analysis, there is another script and files are needed. Inputs are the same.

KEGG reference files: ```kegg_ko_to_map.txt, kegg_pathway_mod_draft.tsv, ko_numbers_reference_list.tsv```

Run: ```R_scripts/KOs_to_pathway.R```

KEGG Photosynthesis pathwway visualization:

Run: ```R_scripts/photosynthesis_pathway/KEGG_pathway_visualization.R```

Description : To perform resampling and calculating KEGG diversity, run the Resampling_KOs.R

Resampling analysis:  ```R_scripts/Resampling_KOs.R```

Visualization : ```R_scripts/KOs_richness.R```

Correlation analysis between KOs and Assembly: ```R_scripts/corr_analysis_KOs/corr_analysis_KOs```

KOs Venn Diagra: ```R_scripts/KO_Venn_diagram/KO_venn.R```

### Step7 Taxonomic assignments of eukaryotic proteins 

Description : The taxonomic origin of eukaryotic proteins from the four datasets (based on the eggNOG result) were cross-checked against the NCBI protein database. This script should be iterated across pipelines.

prokgap_input : ```eggnog_eukaryote_OG_proteins_prokGAP.faa```

eukgap_input : ```eggnog_eukaryote_OG_proteins_eukGAP.faa```

preclass_eukgap_input: ```eggnog_eukaryote_OG_proteins_preclass_eukGAP.faa```

preclass_prokgap_Input: ```eggnog_eukaryote_OG_proteins_preclass_prokGAP.faa```

Run   : ```step7_mmseqs2_taxonomy.sl```

Output : ```eggnog_eukaryote_OG_proteins_prokGAP_taxonomy.tsv, eggnog_eukaryote_OG_proteins_eukGAP_taxonomy.tsv, eggnog_eukaryote_OG_proteins_preclass_eukGAP_taxonomy.tsv, eggnog_eukaryote_OG_proteins_preclass_prokGAP_taxonomy.tsv```

Description : Taxonomy abundance.R should be run to get the protein taxonomy abundance. 

eukgap_input: ```out.eggnog/non_redundant_EUKGAP_protein_eggNOG.emapper.annotations, eggnog_eukaryote_OG_proteins_eukGAP_taxonomy.tsv,eukGAP_all_lake_lama_gene_quant.CPM.fixed.tsv ```

prokgap_input: ```out.eggnog/non_redundant_PROKGAP_protein_eggNOG.emapper.annotations,eggnog_eukaryote_OG_proteins_prokGAP_taxonomy.tsv,prokGAP_all_lake_lama_gene_quant.CPM.fixed.tsv ```

preclass_eukgap_input: ```out.eggnog/non_redundant_preclass_prokGAP_protein_eggNOG.emapper.annotations, eggnog_eukaryote_OG_proteins_preclass_eukGAP_taxonomy.tsv,preclass_eukGAP_all_lake_lama_gene_quant.CPM.fixed.tsv```

preclass_prokgap_Input: ```out.eggnog/non_redundant_preclass_eukGAP_protein_eggNOG.emapper.annotations, eggnog_eukaryote_OG_proteins_preclass_prokGAP_taxonomy.tsv,preclass_prokGAP_all_lake_lama_gene_quant.CPM.fixed.tsv```

Kraken input (To compare with kraken result) : ```lama_shotgun_join_0.8.csv```

Run   : ```taxonomy_abund_mmseqs2.R```

Correlation analysis: ```spearman_corr_test.R```

### step8 aDNA Assessment (PyDamage)

Description: Here we assess aDNA contigs in our assemblies using PyDamage. This process include several subprocesses, bwa, pydamage and kraken2. 

Input : ```out.megahit/[sample_name].final.contigs.fa```

Run : ```aDNA_pydamage_assessment.sl```

Output : ```$sample_name.pydamage_result.csv ```

To bring it all together, ```cat *_pydamage_result.csv > all_samples_pydamage_result.csv ```

Visualization : ```R_script/pydamage_workflow_v1.R```

