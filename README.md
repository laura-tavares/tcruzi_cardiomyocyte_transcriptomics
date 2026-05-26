# RNA-seq Analysis of Human Cardiomyocytes Infected with Trypanosoma cruzi

Differential gene expression analysis of human cardiomyocytes infected with *Trypanosoma cruzi* using bulk RNA-seq and DESeq2.

---

## Dataset

- GEO accession: GSE223600
- Model: Human cardiomyocytes
- Conditions:
  - 0 hpi (control)
  - 24 hpi infected

---

## Workflow

```text
GEO download
→ count matrix construction
→ DESeq2 normalization
→ differential expression analysis
→ PCA and clustering
→ volcano plot visualization
→ GO enrichment analysis
→ targeted hypoxia-related heatmap
