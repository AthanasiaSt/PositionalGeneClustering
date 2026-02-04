# LinearGenomeClustering

## Overview

**LinearGenomeClustering** is a Python-based computational pipeline for detecting non-random positional clustering of genes along linear chromosomes. The method implements a permutation-based statistical framework to assess whether a given set of genes is organized more closely (or more sparsely) along the genome than expected by chance.

The pipeline is designed to be general and reusable across organisms and genomic categorizations. It enables the systematic study of how linear genome organization relates to functional, regulatory, evolutionary, and expression-related gene properties.

Code is provided as used in the analysis and require path and input adaptation.
---

## Key Concepts

The pipeline evaluates two complementary aspects of genomic organization:

1. **Chromosomal enrichment / avoidance**  
   Tests whether a given gene set is preferentially enriched or avoided on specific chromosomes compared to random expectation.

2. **Linear intergenic distance structure**  
   Tests whether consecutive genes in the input set are positioned closer together (or farther apart) than expected under a chromosome-aware random model.

Together, these analyses allow detection of both global chromosomal biases and fine-scale positional clustering along chromosomes.

---

## Method Summary

Starting from a user-defined gene set with genomic coordinates, the pipeline:

1. Computes observed intergenic distances between consecutive genes per chromosome.
2. Generates chromosome-aware null distributions using random permutation of gene positions.
3. Computes empirical z-scores and p-values for positional clustering per chromosome.
4. Identifies significantly clustered (or dispersed) chromosomes.
5. Segments clustered regions into coordinate-level **sub-clusters** based on intra-cluster distance structure.
6. Computes density and structural properties of detected sub-clusters.
7. Optionally compares and overlaps sub-clusters across multiple genomic categories.

---

## Outputs

The pipeline returns multiple result tables, including:

- Chromosome-level z-scores and empirical p-values
- Chromosome enrichment / avoidance statistics
- Significant positionally clustered chromosomes
- Coordinate-level sub-clusters of clustered genes
- Sub-cluster density and structural annotations

These outputs can be used for downstream statistical analysis, visualization, and genome browser tracks.

---

## Applications

This framework was applied to a wide range of genomic annotations in *Saccharomyces cerevisiae*, including:

- Transcription factor target sets
- Gene Ontology categories
- Evolutionary age partitions
- Conservation level categories
- Transcriptional plasticity groups

These analyses revealed pervasive and functionally informative patterns of linear genome organization, linking positional clustering to transcriptional regulation, expression variability, and evolutionary constraints.

---

## Input Requirements

The main input is a gene coordinate table with the following required columns:

- `chr` — chromosome
- `start` — genomic start coordinate
- `end` — genomic end coordinate
- `gene_name` — gene identifier

Additional inputs may include:

- Precomputed intergenic distance dictionaries
- Gene category annotations (e.g. TF targets, GO terms)

---

## Main Components

Key functions in the pipeline include:

- `linearly_clust()` — main positional clustering analysis
- `compute_z_distances()` — intra-cluster distance normalization
- `subclustering()` — segmentation of clustered regions into sub-clusters
- `sub_clusters_elaborate()` — expansion of sub-clusters to gene-level annotations
- `subclusters_plot_across_categories()` — visualization of sub-clusters across chromosomes

---

## Citation

If you use this pipeline in your work, please cite or acknowledge this repository. The method was originally developed as part of a Master’s thesis, under the supervision of Christoforos Nikolaou, focused on linear genome organization and gene clustering in *Saccharomyces cerevisiae*.

---

## Author

Developed by Athanasia Stavropoulou / Supervisor: Christoforos Nikolaou, BSRC 'Alexander Fleming'

