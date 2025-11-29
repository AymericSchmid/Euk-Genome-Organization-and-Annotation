# Eukaryotic Genome Organization and Annotation

## Project Overview

This repository contains all scripts and results from the **Eukaryotic Genome Organization and Annotation** course project, focused on the *Arabidopsis thaliana* accession **Ms-0**. The project encompasses comprehensive genome annotation, including transposable element (TE) annotation, gene prediction, functional annotation, and comparative genomics analyses.

## Biological Context

- **Organism**: *Arabidopsis thaliana*
- **Accession**: Ms-0
- **Genome Assembly**: PacBio HiFi-based assembly (produced in a previous course using hifiasm)
- **Assembly Quality**: High completeness (BUSCO ~97.6%)
- **Comparative Genomes**: TAIR10 (reference), Altai-5, Are-6, Ice-1

## 1. Transposable Element (TE) Annotation

### Overview
Comprehensive annotation and analysis of transposable elements in the Ms-0 genome, including classification, age estimation, and spatial distribution.

### Methods

#### 1.1 TE Annotation with EDTA
- **Tool**: EDTA (Extensive de-novo TE Annotator) v2.2
- **Parameters**: 
  - Species: others (generic parameters)
  - Sensitive mode enabled
  - CDS filtering using TAIR10 representative gene models

#### 1.2 TE Classification Refinement
- **TEsorter**: Classification of LTR retrotransposons into clades
- **Database**: rexdb-plant (plant-specific TE database)
- **Focus**: Copia and Gypsy superfamilies, with clade-level classification

#### 1.3 Divergence and Age Analysis
- **RepeatMasker**: Genome masking and divergence calculation
- **parseRM.pl**: Parsing RepeatMasker output for landscape generation
- **Divergence bins**: 1% bins up to 50% divergence

### Results

#### TE Content Summary
- Complete TE annotation across the genome
- Classification by order, superfamily, and clade
- Divergence profiles for different TE families

#### Key Findings
- **Centromere-enriched families**: Athila and CRM clades (known centromeric TEs in Brassicaceae)
- **LTR-RT distribution**: Copia and Gypsy superfamilies show distinct spatial patterns
- **Divergence landscape**: Shows TE insertion history and activity periods

### Scripts
- `01_EDTA_TE_annotation.slurm`: Main TE annotation with EDTA
- `02_LTRRT_identity_clade_analysis.slurm`: LTR-RT classification and identity analysis
- `03_class_refinment_TEsorter.slurm`: Refinement of TE classification with clade information
- `04_visualize_TE_distribution.slurm`: Circular plot generation for TE distribution
- `05_a_TE_divergence.slurm`: Divergence landscape data generation
- `05_b_plot_landscape.slurm`: Divergence landscape plot generation
- `scripts/perlScripts/parseRM.pl`: RepeatMasker output parser
- `scripts/RAnalysis/02_full_length_LTRs_identity.R`: LTR-RT identity distribution plots
- `scripts/RAnalysis/04_annotation_circlize.R`: Circular genome plots
- `scripts/RAnalysis/05_plot_div.R`: TE divergence landscape plots

### Output Figures
- TE insertion age distributions
- TE divergence landscape (sequence amount in Mbp by divergence percentage)
- Circos-style TE distribution plot (showing Gypsy, Copia, top superfamilies, Athila, CRM clades)

---

## 2. Gene Annotation (MAKER)

### Overview
*Ab initio* and evidence-based gene prediction using the MAKER pipeline, followed by quality filtering and annotation refinement.

### Methods

#### 2.1 MAKER Pipeline
- **Tool**: MAKER v3.01.03
- **Evidence sources**:
  - Protein homology (UniProt Viridiplantae)
  - EST evidence (Trinity transcriptome assembly)
  - *Ab initio* predictions (AUGUSTUS)
- **Repeat masking**: EDTA-annotated TEs

#### 2.2 Quality Filtering
- **AED filtering**: Annotation Edit Distance-based quality filtering
- **InterProScan integration**: Functional domain information used in filtering
- **Feature extraction**: Genes, mRNAs, CDS, exons, UTRs

#### 2.3 Quality Assessment
- **AGAT**: Gene feature statistics
- **BUSCO**: Completeness assessment on protein and transcript sets

### Results

#### Gene Statistics
- **Total predicted genes**: 35,718
- **After quality filtering**: 35,718 (all passed filters)
- **Mono-exonic genes**: 11,405 (31.9%)
- **Mean gene length**: 1,999 bp
- **Longest gene**: 194,792 bp

#### BUSCO Results
- **Genome-level BUSCO**: ~97.6% (from assembly QC)
- **Protein-level BUSCO**: Results available in output files
- **Transcript-level BUSCO**: Results available in output files


### Scripts
- `06_a_maker_prepare.slurm`: Generate MAKER control files
- `06_b_gene_anno_maker.slurm`: Run MAKER gene annotation pipeline
- `06_c_maker_output_prep.slurm`: Merge MAKER output files
- `07_a_rename_genes_transcr.slurm`: Rename gene and transcript IDs
- `07_c_maker_aed_filter_and_subset.slurm`: AED-based quality filtering
- `07_e_plot_AED_distribution.slurm`: AED distribution visualization
- `08_a_busco_QC.slurm`: BUSCO quality assessment
- `08_b_AGAT_annotation_stats.slurm`: Gene annotation statistics
- `scripts/RAnalysis/07_genes_AED_distribution.R`: AED distribution plots

### Output Figures
- AED cumulative distribution (histogram and cumulative curve)
- BUSCO barplots (genome, protein, transcript levels)

---

## 3. Functional Annotation

### Overview
Comprehensive functional annotation of predicted proteins using domain databases and homology searches.

### Methods

#### 3.1 InterProScan
- **Tool**: InterProScan
- **Output format**: TSV and GFF integration

#### 3.2 BLAST Searches
- **Tool**: BLAST+ v2.15.0
- **Databases**:
  1. **TAIR10**: *A. thaliana* representative proteome
  2. **UniProt Viridiplantae**: Reviewed plant proteins
- **Parameters**: 
  - E-value threshold: 1e-5
  - Max target sequences: 10
  - Best hit per query retained

### Results

#### InterProScan Summary
- **Total proteins**: 44,184
- **InterProScan annotated**: 29,131 (66%)
- **Pfam domains**: 28,829 (65%)
- **GO terms**: 18,771 (42%)

#### BLAST Results

**TAIR10 (same species)**
- **Genes with hit**: 35,213
- **Genes without hit**: 505
- **Coverage**: ~98.6%

**UniProt Viridiplantae (curated database)**
- **Genes with hit**: 28,290
- **Genes without hit**: 7,428
- **Coverage**: ~79.2%

### Scripts
- `07_b_iprscan_annotation.slurm`: InterProScan functional annotation
- `09_prot_sequ_homology.slurm`: BLAST searches against TAIR10 and UniProt
- `scripts/RAnalysis/07_iprscan_summary.R`: InterProScan summary statistics

### Output Files
- InterProScan summary table
- BLAST summary table (TAIR10 and UniProt)
- Functional annotation GFF files

---

## 4. Pangenome Analysis

### Overview
Comparative genomics analysis across five *A. thaliana* accessions to identify core, accessory, and species-specific genes.

### Methods

#### 4.1 Data Preparation
- **Accessions**: Altai-5, Are-6, Ice-1, Ms-0, TAIR10
- **Input files**: BED files (gene positions) and FASTA files (protein sequences)
- **Preparation**: Gene ID standardization and format conversion

#### 4.2 Orthogroup Inference
- **Tool**: GENESPACE
- **Synteny detection**: MCScanX
- **Reference genome**: TAIR10

#### 4.3 Pangenome Classification
- **Core orthogroups**: Present in all genomes
- **Accessory orthogroups**: Present in some but not all genomes
- **Species-specific orthogroups**: Present in only one genome

### Results

#### Orthogroup Statistics
- **Total orthogroups**: 33,550
- **Core orthogroups**: 22,061 (65.8%)
- **Species-specific orthogroups**: 5,453 (16.3%)
- **Accessory orthogroups**: 6,036 (18.0%)

#### Gene Counts
- **Genes in core genome**: 111,562
- **Species-specific genes**: 7,005
- **Total genes across all accessions**: ~150,000+

#### Per-Accession Breakdown
Gene counts by category (core/accessory/specific) are available in the output files.

### Scripts
- `10_a_genespace_files_preparation.slurm`: Prepare BED and FASTA files for GENESPACE
- `10_b_genespace_run.slurm`: Run GENESPACE analysis
- `10_c_process_pangenome.slurm`: Process pangenome results and generate plots
- `scripts/RAnalysis/10_b_run_genespace.R`: GENESPACE R interface
- `scripts/RAnalysis/10_c_process_pangenome.R`: Pangenome analysis and visualization

### Output Figures
- Pangenome distribution barplot (orthogroups and genes vs. number of genomes)
- TAIR10 conservation plot (genes shared vs. not shared with TAIR10 per genome)
- Per-accession gene counts table (core/accessory/specific)

---

## 5. Synteny Analysis

### Overview
Synteny analysis between Ms-0 contigs and TAIR10 chromosomes to assess genome organization and identify potential misassemblies.

### Methods
- **Tool**: GENESPACE synteny visualization
- **Reference**: TAIR10 chromosomes
- **Query**: Ms-0 contigs

### Results and Limitations

#### Observations
- Alignment of Ms-0 contigs to TAIR10 chromosomes shows general synteny patterns
- Some contigs align to multiple chromosomes, suggesting potential misassemblies
- High assembly fragmentation limits detailed synteny interpretation

#### Important Limitations

**Conclusions are limited due to**:
1. **High assembly fragmentation**: Many small contigs make comprehensive synteny analysis challenging
2. **Potential misassemblies**: Some contigs include segments from multiple chromosomes, indicating assembly artifacts
3. **Contig-level analysis**: Without chromosome-level assembly, detailed synteny interpretation is restricted

### Scripts
- Synteny visualization is part of the GENESPACE output (generated in `10_b_genespace_run.slurm`)

### Output Figures
- Synteny plot (Ms-0 contigs aligned to TAIR10 chromosomes)

---

## Scripts and Workflow

### Directory Structure
```
scripts/
├── 00_dir_set_up.slurm                    # Directory setup and file linking
├── 01_EDTA_TE_annotation.slurm            # TE annotation with EDTA
├── 02_LTRRT_identity_clade_analysis.slurm  # LTR-RT classification
├── 03_class_refinment_TEsorter.slurm      # TE classification refinement
├── 04_visualize_TE_distribution.slurm     # TE distribution visualization
├── 05_a_TE_divergence.slurm               # TE divergence analysis
├── 05_b_plot_landscape.slurm              # Divergence landscape plots
├── 06_a_maker_prepare.slurm               # MAKER control file generation
├── 06_b_gene_anno_maker.slurm             # MAKER gene annotation
├── 06_c_maker_output_prep.slurm           # MAKER output merging
├── 07_a_rename_genes_transcr.slurm        # Gene ID renaming
├── 07_b_iprscan_annotation.slurm          # InterProScan annotation
├── 07_c_maker_aed_filter_and_subset.slurm # AED filtering
├── 07_d_visualize_TE_and_genes_distribution.slurm # Combined visualization
├── 07_e_plot_AED_distribution.slurm       # AED distribution plots
├── 08_a_busco_QC.slurm                    # BUSCO quality assessment
├── 08_b_AGAT_annotation_stats.slurm       # Gene annotation statistics
├── 09_prot_sequ_homology.slurm            # BLAST homology searches
├── 10_a_genespace_files_preparation.slurm # GENESPACE file preparation
├── 10_b_genespace_run.slurm               # GENESPACE analysis
├── 10_c_process_pangenome.slurm          # Pangenome processing
├── perlScripts/
│   └── parseRM.pl                         # RepeatMasker parser
└── RAnalysis/
    ├── 02_full_length_LTRs_identity.R     # LTR-RT identity plots
    ├── 04_annotation_circlize.R           # Circular genome plots
    ├── 05_plot_div.R                      # Divergence landscape plots
    ├── 07_genes_AED_distribution.R        # AED distribution plots
    ├── 07_iprscan_summary.R               # InterProScan summaries
    ├── 10_b_run_genespace.R               # GENESPACE R interface
    └── 10_c_process_pangenome.R           # Pangenome analysis
```

### Workflow Overview

#### Phase 1: Setup and TE Annotation
1. `00_dir_set_up.slurm`: Set up directory structure and link input files
2. `01_EDTA_TE_annotation.slurm`: Annotate TEs with EDTA
3. `02_LTRRT_identity_clade_analysis.slurm`: Classify LTR-RTs into clades
4. `03_class_refinment_TEsorter.slurm`: Refine TE classification
5. `04_visualize_TE_distribution.slurm`: Visualize TE distribution
6. `05_a_TE_divergence.slurm`: Calculate TE divergence
7. `05_b_plot_landscape.slurm`: Plot divergence landscape

#### Phase 2: Gene Annotation
8. `06_a_maker_prepare.slurm`: Prepare MAKER control files
9. `06_b_gene_anno_maker.slurm`: Run MAKER annotation pipeline
10. `06_c_maker_output_prep.slurm`: Merge MAKER outputs
11. `07_a_rename_genes_transcr.slurm`: Rename gene IDs
12. `07_c_maker_aed_filter_and_subset.slurm`: Filter by AED
13. `07_e_plot_AED_distribution.slurm`: Plot AED distribution

#### Phase 3: Functional Annotation
14. `07_b_iprscan_annotation.slurm`: InterProScan annotation
15. `09_prot_sequ_homology.slurm`: BLAST searches

#### Phase 4: Quality Assessment
16. `08_a_busco_QC.slurm`: BUSCO analysis
17. `08_b_AGAT_annotation_stats.slurm`: Gene statistics

#### Phase 5: Comparative Genomics
18. `10_a_genespace_files_preparation.slurm`: Prepare comparative files
19. `10_b_genespace_run.slurm`: Run GENESPACE
20. `10_c_process_pangenome.slurm`: Process pangenome results

### Reproducibility

All scripts are designed to run on a SLURM cluster with the following requirements:
- **Container system**: Apptainer/Singularity
- **Required modules**: R, Perl, BLAST+, BUSCO, SeqKit, and others (see individual scripts)
- **Data dependencies**: Course data directory with reference genomes, databases, and software containers

To reproduce the analyses:
1. Ensure all input files are in the correct locations (as specified in `00_dir_set_up.slurm`)
2. Run scripts in numerical order (00 → 10_c)
3. Adjust paths in scripts to match your environment
4. Ensure all required software modules and containers are available

---

## Output Files

### Key Figures Generated

1. **TE Analysis**
   - `02_LTR_Copia_Gypsy_cladelevel.png/pdf`: LTR-RT identity distribution by clade
   - `04_TE_distribution.pdf`: Circular plot of TE distribution
   - `05_TE_divergence_landscape.pdf`: TE divergence landscape
   - `07_TE_and_genes_distribution.pdf`: Combined TE and gene distribution

2. **Gene Annotation**
   - `07_AED_distribution.pdf`: AED score distribution (histogram and cumulative)
   - BUSCO plots: Genome, protein, and transcript completeness

3. **Functional Annotation**
   - InterProScan summary tables
   - BLAST summary tables

4. **Comparative Genomics**
   - `10_c_pangenome_frequency_plot.pdf`: Pangenome distribution plot
   - TAIR10 conservation plot
   - Synteny visualization

### Data Files

All intermediate and final annotation files are stored in organized directories:
- `TE_annotation/`: TE annotation results
- `gene_annotation/`: Gene annotation results
- `annotation_qc/`: Quality assessment results
- `comparative_genomics/`: Pangenome and synteny results
- `plots/`: All generated figures

---

## References

### Databases and Resources
- **TAIR10**: [The Arabidopsis Information Resource](https://www.arabidopsis.org/)
- **UniProt**: [UniProt Knowledgebase](https://www.uniprot.org/)
- **InterPro**: [InterPro Protein Classification](https://www.ebi.ac.uk/interpro/)
- **Pfam**: [Pfam Protein Families Database](https://pfam.xfam.org/)
- **Gene Ontology**: [GO Consortium](http://geneontology.org/)

### Software Tools
- **EDTA**: [EDTA GitHub Repository](https://github.com/oushujun/EDTA)
- **MAKER**: [MAKER Website](http://www.yandell-lab.org/software/maker.html)
- **BUSCO**: [BUSCO Documentation](https://busco.ezlab.org/)
- **InterProScan**: [InterProScan Documentation](https://www.ebi.ac.uk/interpro/download/)
- **GENESPACE**: [GENESPACE GitHub Repository](https://github.com/jtlovell/GENESPACE)
- **TEsorter**: [TEsorter GitHub Repository](https://github.com/zhangrengang/TEsorter)
- **RepeatMasker**: [RepeatMasker Website](http://www.repeatmasker.org/)

### Key Publications
- **Lian et al. 2024**: Pangenome analysis of *Arabidopsis thaliana*

### Additional Resources
- **AUGUSTUS**: [AUGUSTUS Gene Prediction](http://augustus.gobics.de/)
- **AGAT**: [AGAT GitHub Repository](https://github.com/NBISweden/AGAT)
- **MCScanX**: [MCScanX GitHub Repository](https://github.com/wyp1125/MCScanX)