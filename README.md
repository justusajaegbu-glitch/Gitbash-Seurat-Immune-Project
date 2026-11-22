# Gitbash-Seurat-Immune-Project
Aligning control and interferon-beta-stimulated PBMCs to find which cell types respond, and which genes drive the response.


## 1. Objective

**Question:** Once control and stimulated PBMCs are aligned, do interferon-response genes (ISG15, IFI6, MX1, IFIT3, CXCL10) rise in every immune cell type, or are some responses specific to particular cell types?

**What I tested:**
(i) Can the two conditions be integrated so that cells group by type and not by treatment?
(ii) Which genes mark each cell type in both conditions (conserved markers)?
(iii) Which genes change with stimulation within a cell type (CD4 naive T, CD14 monocytes, B cells)?

## 2. Dataset

PBMCs from [Kang et al.](https://www.nature.com/articles/nbt.4042), split into a control group and a group stimulated with interferon-beta. Two count matrices were used: [`immune_control_expression_matrix.txt.gz` and `immune_stimulated_expression_matrix.txt.gz`](https://www.dropbox.com/s/79q6dttg8yl20zg/immune_alignment_expression_matrices.zip?dl=1).

## 3. Tools

R, Seurat (log-normalisation, anchor-based integration, PCA, UMAP, clustering, differential expression), ggplot2 and patchwork for figures. Run from the command line in Git Bash with `Rscript`.

## 4. Workflow

1. **Load and filter:** build one Seurat object per condition (genes in at least 5 cells), add a `stim` label (CTRL or STIM), keep cells with more than 500 detected genes.
2. **Normalise and select features:** log-normalisation, then the top 2000 variable genes (vst) per condition.
3. **Integrate:** find integration anchors between the two conditions (dims 1:20) and integrate the data.
4. **Reduce and cluster:** scale, PCA (30 components), UMAP and nearest-neighbour graph on dims 1:20, clustering at resolution 0.5. This gave 14 clusters (0-13).
5. **Annotate:** assign cell types from canonical markers.
6. **Conserved markers:** find genes that mark a cluster in both conditions.
7. **Stimulation response:** compare STIM vs CTRL within cell types.
   <img width="1190" height="586" alt="image" src="https://github.com/user-attachments/assets/fc1c1263-f4ac-4dd9-83c6-f00f8968f307" />


## 5. QC and integration results

- Cells with 500 or fewer detected genes were removed in each condition.
- After integration, control and stimulated cells overlap within the same clusters, so clusters reflect cell type and not treatment (`umap_split`, `immune_alignment_umap`).
- Cell counts before and after filtering are printed by the script and can be added here.

## 6. Cluster annotation

| Cluster | Cell type | Cluster | Cell type |
|---|---|---|---|
| 0 | CD14 monocytes | 7 | NK cells |
| 1 | CD4 naive T | 8 | Dendritic cells |
| 2 | CD4 memory T | 9 | Activated B |
| 3 | CD16 monocytes | 10 | Megakaryocytes |
| 4 | B cells | 11 | Plasmacytoid DCs |
| 5 | CD8 T | 12 | Erythrocytes |
| 6 | Activated T | 13 | Monocyte/megakaryocyte doublets |

Annotation used canonical markers (CD3D, SELL, CREM, CD8A, GNLY, CD79A, FCGR3A, CCL2, PPBP) shown in `canonical_markers_featureplot` and `celltype_markers_dotplot`.

## 7. Results

- **NK cluster (7):** conserved markers across both conditions include GNLY, NKG7 and GZMB (`cluster7_nk_markers`).
- **Interferon response is broad:** ISG15, IFI6 and MX1 are higher in stimulated CD4 naive T cells and CD14 monocytes (`avg_expression_scatter`).
- **B cells:** the strongest stimulation changes (B STIM vs B CTRL) are ISG15, IFIT3, IFI6, ISG20 and CXCL10 (`b_interferon_response`).
- **Cell-type markers are stable:** CD3D and GNLY keep their cell-type pattern in both conditions, while IFI6 rises across clusters (`interferon_response_featureplot`).
- **CXCL10** responds mainly in monocytes and B cells, and CD14 falls after stimulation in monocytes (`interferon_response_violin`).

## 8. Biological interpretation

Interferon-beta triggers a shared antiviral programme (ISG15, IFI6, MX1, IFIT3) in almost every PBMC type, and only some genes are cell-type specific: CXCL10 is strongest in monocytes and B cells. Integration first removes the treatment effect, so cell types can be named once and then compared between conditions. This pattern (a common interferon core with selective extras) is what makes cell-type-aware analysis more informative than pooling all cells.

## 9. Limitations

- Annotation is based on a small set of canonical markers and one clustering resolution (0.5).
- Cluster 13 (monocyte/megakaryocyte doublets) is a technical group and is not interpreted biologically.
- Pooled comparisons of cells are not replicate-aware, so p-values from cell-level tests overstate confidence. A pseudobulk analysis with replicates would be the next step.
- Only a single stimulation time point and dose are available.

## 10. Reproducibility

```bash
Rscript scripts/run_alignment.R
# or run the full pipeline:
bash scripts/run_full_pipeline.sh

Place the two input matrices in data/. Outputs are written to output/, and sessionInfo() is saved with the results to record package versions.

```
Gitbash-Seurat-Immune-Project/
├── README.md
├── .gitignore                        # *.rds, data/*.gz, .RData, .Rhistory
├── scripts/
│   ├── run_alignment.R
│   ├── run_full_pipeline.sh
│   └── report.html
├── data/                             # input matrices (not tracked)
└── output/
    ├── umap_celltypes, umap_split, immune_alignment_umap
    ├── canonical_markers_featureplot, celltype_markers_dotplot
    ├── cluster7_nk_markers, b_interferon_response
    ├── avg_expression_scatter
    ├── interferon_response_featureplot, interferon_response_violin
    └── immune_combined_seurat.rds, immune_annotated_seurat.rds   (not tracked)
```

## 11. References

- Kang et al., multiplexed droplet single-cell RNA-sequencing using natural genetic variation, Nature Biotechnology (2018). Source of the PBMC data.
- Stuart et al., Comprehensive integration of single-cell data, Cell (2019). Anchor-based integration method.
- Seurat integration vignette, Satija Lab: https://satijalab.org/seurat/archive/v3.0/immune_alignment.html. The analysis steps follow this workflow, and I ran it end to end in Git Bash.
