
library(Seurat)
library(tidyverse)
library(edgeR)        # bulk RNA-seq DE, used on the pseudobulk
library(patchwork)    # put two ggplots side by side
library(ggrepel)      # gene labels that don't overlap
# DESeq2 must be installed for FindMarkers(test.use = "DESeq2")
# fgsea + msigdbr are needed for the GSEA part (section 6)

set.seed(1)
dir.create("plots", showWarnings = FALSE)
dir.create("res",   showWarnings = FALSE)

theme_set(theme_classic(base_size = 16))               # one style for all ggplots
col_inj <- c(host = "grey55", injected = "#d6604d")    # one colour per condition

# helper: save a ggplot as a PNG in plots/
save_plot <- function(p, file, w = 10, h = 5.5) {
  ggsave(file.path("plots", file), p, width = w, height = h, dpi = 200, bg = "white")
}

## ---- 1. Load the Day 2 object and define the condition ----------------------
seu <- readRDS("/varidata/researchtemp/hpctmp/single_cell_workshop_2026/day2/data/day2_session2_seu.rds")

DefaultAssay(seu) <- "RNA"     # DE uses the RNA assay, never the integrated values
seu <- JoinLayers(seu)         # Day 2 split the layers by sample; join them again

# remove cells the published annotation flags as doublets or "stripped" nuclei
seu <- subset(seu, subset = !celltype.mapped %in% c("Doublet", "Stripped"))

# tomato = TRUE  -> cell comes from the injected td-Tomato+ stem cells
# tomato = FALSE -> cell comes from the host embryo
seu$injection <- ifelse(seu$tomato, "injected", "host")
seu$injection <- factor(seu$injection, levels = c("host", "injected"))   # host = reference

meta <- seu[[]]                         # cell metadata table (same as seu@meta.data)
numbers <- list(n_cells = ncol(seu))    # numbers quoted on the slides

## ---- Figure: cells per sample, grouped by pool --------------------------------
cells_per_sample <- meta |> dplyr::count(sample, pool, injection)
cells_per_sample$sample <- factor(cells_per_sample$sample)

p <- ggplot(cells_per_sample, aes(x = sample, y = n, fill = injection)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n), vjust = -0.3, size = 5) +
  facet_grid(~ pool, scales = "free_x", space = "free_x",
             labeller = labeller(pool = function(x) paste("pool", x))) +
  scale_fill_manual(values = col_inj) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "sample", y = "cells after QC", fill = NULL,
       title = "Each pool contributes one host and one injected sample")
save_plot(p, "design_cells_per_sample.png", 10, 5)

## ---- Figure: host and injected cells on the UMAP -------------------------------
p <- DimPlot(seu, reduction = "umap", group.by = "injection",
             shuffle = TRUE,                    # random plotting order: no group hides the other
             cols = col_inj, pt.size = 0.2) +
  labs(title = "Host and injected cells share one map", x = "UMAP 1", y = "UMAP 2")
save_plot(p, "umap_injection.png", 8, 6)

## ---- Figure: where are the Mesenchyme cells? ----------------------------------

mes_cells <- WhichCells(seu, expression = celltype.mapped == "Mesenchyme")
p_umap <- DimPlot(seu, reduction = "umap", cells.highlight = mes_cells,
                  cols.highlight = "#2166ac", sizes.highlight = 0.6, pt.size = 0.2) +
  NoLegend() +
  labs(title = "Mesenchyme (published label)", x = "UMAP 1", y = "UMAP 2")

if ("RNA_snn_res.0.3" %in% colnames(meta)) {
  clust_df <- meta |>
    group_by(cluster = RNA_snn_res.0.3) |>
    summarise(mes = sum(celltype.mapped == "Mesenchyme"),   # Mesenchyme cells in the cluster
              total = n()) |>                               # all cells in the cluster
    filter(mes > 0) |>
    mutate(share = round(100 * mes / total)) |>
    slice_max(mes, n = 6)
  p_bar <- ggplot(clust_df, aes(x = reorder(cluster, -mes), y = mes)) +
    geom_col(fill = "#2166ac", width = 0.7) +
    geom_text(aes(label = paste0(share, "% of cluster")), vjust = -0.4, size = 4.5) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(x = "Day 2 cluster (resolution 0.3)", y = "Mesenchyme cells",
         title = "Which clusters contain them?")
  save_plot(p_umap | p_bar, "mesenchyme_location.png", 13, 5.5)
} else {
  save_plot(p_umap, "mesenchyme_location.png", 8, 6)
}

## ---- 2. Per-cell test (Seurat FindMarkers), Mesenchyme -----------------------
mes <- subset(seu, subset = celltype.mapped == "Mesenchyme")
numbers$mes_cells <- table(mes$sample)

# injected vs host; every cell is treated as an independent observation
Idents(mes) <- "injection"
sc_de <- FindMarkers(mes, ident.1 = "injected", ident.2 = "host",
                     logfc.threshold = 0,   # test all genes, no fold-change cut-off
                     min.pct = 0.1)         # gene detected in at least 10 % of cells in one group
numbers$sc_degs <- sum(sc_de$p_val_adj < 0.05)

# control comparison: two HOST samples (6 and 8), no condition difference
host_only <- subset(mes, subset = injection == "host" & sample %in% c(6, 8))
Idents(host_only) <- "sample"
fake_de <- FindMarkers(host_only, ident.1 = "6", ident.2 = "8",
                       logfc.threshold = 0, min.pct = 0.1)
numbers$fake_degs <- sum(fake_de$p_val_adj < 0.05)

## ---- 3. Pseudobulk with Seurat: AggregateExpression ---------------------------
# Sum the RAW counts of all cells that share a cell type AND a sample.
# pool and injection are fixed within a sample, so listing them does not create
# extra groups: it only keeps them as metadata columns of the result.
pseudo <- AggregateExpression(seu,
                              assays = "RNA",
                              return.seurat = TRUE,
                              group.by = c("celltype.mapped", "sample", "pool", "injection"))
pseudo$injection <- factor(pseudo$injection, levels = c("host", "injected"))
head(pseudo[[]])                          # one row per pseudobulk sample

# How many cells were summed into each pseudobulk sample?
n_cells <- table(seu$celltype.mapped, seu$sample)                  # cell type x sample
pseudo$ncells <- as.vector(n_cells[cbind(pseudo$celltype.mapped,   # look up row = cell type
                                         pseudo$sample)])          #        column = sample
write.csv(pseudo[[]], "res/pseudobulk_ncells.csv", row.names = FALSE)

all_celltypes <- sort(unique(pseudo$celltype.mapped))

# Drop pseudobulk samples built from fewer than 10 cells (too noisy)
pseudo <- subset(pseudo, subset = ncells >= 10)

## ---- 4a. Seurat vignette way: FindMarkers with DESeq2 on the pseudobulk -------
mes_pseudo <- subset(pseudo, subset = celltype.mapped == "Mesenchyme")
Idents(mes_pseudo) <- "injection"
bulk_de <- FindMarkers(mes_pseudo, ident.1 = "injected", ident.2 = "host",
                       test.use = "DESeq2")     # 3 vs 3 samples, pairing by pool NOT used
numbers$pb_deseq2_unpaired <- sum(bulk_de$p_val_adj < 0.05, na.rm = TRUE)
write.csv(bulk_de, "res/DE_Mesenchyme_pseudobulk_FindMarkers_DESeq2.csv")

## ---- 4b. Paired model with edgeR: ~ pool + injection --------------------------
y <- DGEList(counts  = as.matrix(LayerData(mes_pseudo, assay = "RNA", layer = "counts")),
             samples = mes_pseudo[[]])

# keep genes with enough counts to be tested
keep <- filterByExpr(y, group = y$samples$injection)
numbers$genes_kept  <- sum(keep)
numbers$genes_total <- length(keep)
y <- y[keep, , keep.lib.sizes = FALSE]

# TMM normalization, exactly as for bulk RNA-seq
y <- calcNormFactors(y)

## ---- Figure: MDS ------------------------------------------------------------------
mds <- plotMDS(y, plot = FALSE)        
mds_df <- data.frame(
  dim1      = mds$x,                   
  dim2      = mds$y,
  injection = y$samples$injection,
  pool      = factor(y$samples$pool),
  label     = paste0("s", y$samples$sample, " (", y$samples$ncells, " cells)")
)
p <- ggplot(mds_df, aes(x = dim1, y = dim2, colour = injection, shape = pool)) +
  geom_point(size = 5) +
  geom_text_repel(aes(label = label), size = 5, show.legend = FALSE, box.padding = 0.8) +
  scale_colour_manual(values = col_inj) +
  labs(x = sprintf("Leading logFC dim 1 (%.0f%%)", 100 * mds$var.explained[1]),
       y = sprintf("Leading logFC dim 2 (%.0f%%)", 100 * mds$var.explained[2]),
       title = "Mesenchyme pseudobulk samples (MDS)")
save_plot(p, "mds_mesenchyme.png", 9, 6)

## ---- Model and test ---------------------------------------------------------------
# expression ~ pool + injection: compare injected vs host WITHIN each pool (paired)
design <- model.matrix(~ factor(pool) + injection, data = y$samples)

y   <- estimateDisp(y, design)                    # biological variation between replicates
fit <- glmQLFit(y, design, robust = TRUE)         # quasi-likelihood negative binomial model
res <- glmQLFTest(fit, coef = "injectioninjected")   # test the injection effect

mes_table <- topTags(res, n = Inf)$table          # all genes, sorted by p-value
mes_table$gene <- rownames(mes_table)
write.csv(mes_table, "res/DE_Mesenchyme_pseudobulk_edgeR.csv", row.names = FALSE)

numbers$pb_up      <- sum(mes_table$FDR < 0.05 & mes_table$logFC > 0)
numbers$pb_down    <- sum(mes_table$FDR < 0.05 & mes_table$logFC < 0)
numbers$common_bcv <- sqrt(y$common.dispersion)

## ---- Figure: BCV -------------------------------------------------------------------
png("plots/bcv_mesenchyme.png", width = 1800, height = 1200, res = 200)
par(cex = 1.3, mar = c(4.5, 4.5, 3, 1))
plotBCV(y, main = sprintf("Mesenchyme: common BCV = %.2f", numbers$common_bcv))
dev.off()

## ---- Figure: one point per sample, lines join the two samples of a pool -----------
logcpm <- cpm(y, log = TRUE)                     # normalized log2 counts per million
genes  <- intersect(c("tomato-td", "Xist", "Cdkn1c"), rownames(logcpm))

paired_df <- data.frame(
  gene      = rep(genes, each = ncol(logcpm)),
  logCPM    = as.vector(t(logcpm[genes, , drop = FALSE])),
  injection = rep(y$samples$injection, times = length(genes)),
  pool      = factor(rep(y$samples$pool, times = length(genes)))
)
paired_df$gene <- factor(paired_df$gene, levels = genes)

p <- ggplot(paired_df, aes(x = injection, y = logCPM, group = pool, colour = pool)) +
  geom_line(linewidth = 1) +
  geom_point(size = 4) +
  facet_wrap(~ gene, scales = "free_y") +
  labs(x = NULL, y = "log2 CPM",
       title = "Mesenchyme pseudobulk: one point per sample, lines join each pool") +
  theme(strip.text = element_text(face = "bold", size = 16))
save_plot(p, "paired_mesenchyme.png", 12, 5)

## ---- Figure: how many significant genes does each approach find? ------------------
deg_counts <- data.frame(
  test = c("Per cell (Wilcoxon)\ninjected vs host",
           "Per cell (Wilcoxon)\nhost vs host, s6 vs s8",
           "Pseudobulk, FindMarkers\nDESeq2 (unpaired)",
           "Pseudobulk, edgeR\n~ pool + injection"),
  n = c(numbers$sc_degs, numbers$fake_degs,
        numbers$pb_deseq2_unpaired, numbers$pb_up + numbers$pb_down)
)
deg_counts$test <- factor(deg_counts$test, levels = deg_counts$test)

p <- ggplot(deg_counts, aes(x = test, y = n, fill = test)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = n), vjust = -0.3, size = 7) +
  scale_fill_manual(values = c("#d6604d", "#f4a582", "#92c5de", "#2166ac")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "significant genes (5 %)", title = "Mesenchyme: per-cell vs pseudobulk tests")
save_plot(p, "deg_counts_comparison.png", 12, 5.5)

## ---- Figure: same genes, very different p-values --------------------------------------
pb_degs  <- mes_table$gene[mes_table$FDR < 0.05]
cell_top <- setdiff(rownames(sc_de), pb_degs)[1:8]          # top per-cell genes not DE in pseudobulk
genes_show <- c(pb_degs, cell_top)

pv <- data.frame(
  gene       = genes_show,
  per_cell   = -log10(sc_de[genes_show, "p_val"] + 1e-300),
  pseudobulk = -log10(mes_table[match(genes_show, mes_table$gene), "PValue"])
)
pv <- pv[complete.cases(pv), ]
pv$label <- ifelse(pv$gene %in% pb_degs, paste0(pv$gene, " *"), pv$gene)

pv_long <- pivot_longer(pv, c(per_cell, pseudobulk), names_to = "test", values_to = "mlog10p")
p <- ggplot(pv_long, aes(x = mlog10p, y = reorder(label, mlog10p, FUN = max))) +
  geom_line(aes(group = label), colour = "grey70") +
  geom_point(aes(colour = test), size = 3.5) +
  geom_vline(xintercept = -log10(0.05), lty = 2) +
  scale_x_log10() +
  scale_colour_manual(values = c(per_cell = "#d6604d", pseudobulk = "#2166ac"),
                      labels = c(per_cell = "per cell (n = cells)", pseudobulk = "pseudobulk (n = 3 vs 3)")) +
  labs(x = "-log10(p)   (log scale; dashed = p 0.05)", y = NULL, colour = NULL,
       title = "Same genes, very different p-values",
       subtitle = "* = significant in the pseudobulk analysis (FDR < 0.05)") +
  theme(legend.position = "top", axis.text.y = element_text(size = 12))
save_plot(p, "pvalues_per_gene.png", 10, 7)

## ---- 5. All cell types: the same edgeR steps in a function -----------------------
# pb = pseudobulk Seurat object for ONE cell type. Returns NULL if it can't be tested.
run_edger <- function(pb) {
  # testable only with >= 2 pseudobulk samples per condition ...
  if (min(table(pb$injection)) < 2) return(NULL)
  y <- DGEList(counts  = as.matrix(LayerData(pb, assay = "RNA", layer = "counts")),
               samples = pb[[]])
  design <- model.matrix(~ factor(pool) + injection, data = y$samples)
  # ... and more samples than model terms (replicates left to estimate variation)
  if (nrow(design) <= ncol(design)) return(NULL)
  
  keep <- filterByExpr(y, group = y$samples$injection)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design, robust = TRUE)
  res <- glmQLFTest(fit, coef = "injectioninjected")
  
  tab <- topTags(res, n = Inf)$table
  tab$gene <- rownames(tab)
  tab
}

de_list <- list()
for (ct in unique(pseudo$celltype.mapped)) {
  pb_ct  <- subset(pseudo, subset = celltype.mapped == ct)
  res_ct <- run_edger(pb_ct)
  if (!is.null(res_ct)) de_list[[ct]] <- res_ct
}
numbers$untested <- setdiff(all_celltypes, names(de_list))   # cell types that could not be tested
numbers$n_tested <- length(de_list)

de_all <- bind_rows(de_list, .id = "celltype")               # one long table, all genes
write.csv(de_all, "res/DE_injected_vs_host_pseudobulk_edgeR.csv", row.names = FALSE)

## ---- Figure: DEGs per cell type ---------------------------------------------------------
per_celltype <- de_all |>
  group_by(celltype) |>
  summarise(up = sum(FDR < 0.05 & logFC > 0), down = sum(FDR < 0.05 & logFC < 0)) |>
  pivot_longer(c(up, down), names_to = "direction", values_to = "n")
per_celltype$n_signed <- ifelse(per_celltype$direction == "down", -per_celltype$n, per_celltype$n)

p <- ggplot(per_celltype, aes(x = n_signed, y = reorder(celltype, n, FUN = sum), fill = direction)) +
  geom_col() +
  geom_vline(xintercept = 0) +
  scale_fill_manual(values = c(down = "#4393c3", up = "#d6604d")) +
  labs(x = "DEGs at FDR < 0.05 (down | up)", y = NULL, fill = NULL,
       title = "Injected vs host, per cell type (pseudobulk edgeR)",
       caption = paste("Not testable:", paste(numbers$untested, collapse = ", "))) +
  theme(axis.text.y = element_text(size = 12), plot.caption = element_text(size = 11, hjust = 1))
save_plot(p, "degs_per_celltype.png", 11, 7)

## ---- Figure: genes that change in many cell types ----------------------------------------
gene_class <- c(
  "tomato-td" = "sorting marker (positive control)",
  Xist = "sex chromosomes", Erdr1 = "sex chromosomes", Mid1 = "sex chromosomes",
  Cdkn1c = "imprinted", Phlda2 = "imprinted", Kcnq1ot1 = "imprinted", H13 = "imprinted",
  Mcts2 = "imprinted", Nnat = "imprinted", Impact = "imprinted", Slc38a4 = "imprinted",
  Zdbf2 = "imprinted", Peg3 = "imprinted")

recurring <- de_all |>
  filter(FDR < 0.05) |>
  mutate(direction = ifelse(logFC > 0, "up", "down")) |>
  dplyr::count(gene, direction, sort = TRUE) |>           # in how many cell types?
  head(15)
recurring$class <- ifelse(recurring$gene %in% names(gene_class), gene_class[recurring$gene], "other")
recurring$gene  <- paste0(recurring$gene, " (", recurring$direction, ")")

p <- ggplot(recurring, aes(x = n, y = reorder(gene, n), fill = class)) +
  geom_col() +
  geom_vline(xintercept = numbers$n_tested, lty = 2) +    # dashed line = all tested cell types
  scale_fill_manual(values = c("sorting marker (positive control)" = "#1b9e77",
                               "sex chromosomes" = "#7570b3", imprinted = "#d95f02",
                               other = "grey60")) +
  labs(x = "number of cell types where the gene is DE (dashed = all tested)", y = NULL, fill = NULL,
       title = "Genes that change in (almost) every cell type") +
  theme(legend.position = "bottom", axis.text.y = element_text(size = 12))
save_plot(p, "recurring_genes.png", 11, 7)

## ---- 6. GSEA on the Mesenchyme pseudobulk results ------------------------------------------
# Rank ALL tested genes by a signed statistic (direction x strength of evidence),
# then ask whether the genes of a pathway sit at the top or bottom of the ranking.
tryCatch({                       # if anything in the GSEA part fails, the script carries on
  if (requireNamespace("fgsea", quietly = TRUE) && requireNamespace("msigdbr", quietly = TRUE)) {
    library(fgsea)
    library(msigdbr)
    
    # mouse versions of the 50 MSigDB Hallmark gene sets
    hallmark <- tryCatch(msigdbr(species = "Mus musculus", collection = "H"),            # msigdbr >= 10
                         error = function(e) msigdbr(species = "Mus musculus", category = "H"))  # older
    pathways <- split(hallmark$gene_symbol, hallmark$gs_name)    # list: pathway -> genes
    
    ranks <- sign(mes_table$logFC) * -log10(mes_table$PValue)    # + up in injected, - down
    names(ranks) <- mes_table$gene
    ranks <- sort(ranks, decreasing = TRUE)
    
    gsea_res <- fgsea(pathways = pathways, stats = ranks, minSize = 15, maxSize = 500)
    gsea_res <- gsea_res[order(gsea_res$pval), ]
    numbers$gsea_sig <- sum(gsea_res$padj < 0.05)
    
    gsea_out <- as.data.frame(gsea_res)
    gsea_out$leadingEdge <- sapply(gsea_out$leadingEdge, paste, collapse = ";")
    write.csv(gsea_out, "res/GSEA_Mesenchyme_hallmark.csv", row.names = FALSE)
    
    top_gsea <- head(gsea_out, 12)
    top_gsea$pathway <- gsub("HALLMARK_", "", top_gsea$pathway)
    p <- ggplot(top_gsea, aes(x = NES, y = reorder(pathway, NES), fill = padj < 0.05)) +
      geom_col() +
      geom_vline(xintercept = 0) +
      scale_fill_manual(values = c(`FALSE` = "grey70", `TRUE` = "#b2182b"),
                        labels = c(`FALSE` = "padj >= 0.05", `TRUE` = "padj < 0.05")) +
      labs(x = "normalized enrichment score (NES)\n< 0: lower in injected   |   > 0: higher in injected",
           y = NULL, fill = NULL, title = "Mesenchyme: Hallmark GSEA (12 smallest p-values)") +
      theme(axis.text.y = element_text(size = 11), legend.position = "top")
    save_plot(p, "gsea_mesenchyme.png", 11, 6.5)
  } else {
    message("fgsea and/or msigdbr not installed: GSEA figure skipped")
  }
}, error = function(e) message("GSEA skipped: ", conditionMessage(e)))

## ---- 7. Differential abundance: do cell-type proportions change? ---------------------------
# Count cells per cell type (rows) and sample (columns), then use the same
# edgeR model on the cell numbers.
abundances <- unclass(table(seu$celltype.mapped, seu$sample))

sample_info <- meta |> distinct(sample, pool, injection) |> arrange(sample)
sample_info <- sample_info[match(colnames(abundances), sample_info$sample), ]

y.ab <- DGEList(abundances, samples = sample_info)
keep <- filterByExpr(y.ab, group = y.ab$samples$injection)
y.ab <- y.ab[keep, ]
design.ab <- model.matrix(~ factor(pool) + injection, y.ab$samples)
y.ab   <- estimateDisp(y.ab, design.ab, trend = "none")      # few "genes" (cell types): no trend
fit.ab <- glmQLFit(y.ab, design.ab, robust = TRUE, abundance.trend = FALSE)
da <- topTags(glmQLFTest(fit.ab, coef = "injectioninjected"), n = Inf)$table
da$celltype <- rownames(da)
write.csv(da, "res/DA_injected_vs_host.csv", row.names = FALSE)
numbers$da_sig <- da$celltype[da$FDR < 0.05]

## ---- Figure: proportions per sample for the top DA cell types ---------------------------
show_ct <- unique(c(head(da$celltype, 4), "Mesenchyme"))
prop_df <- as.data.frame(as.table(prop.table(table(seu$celltype.mapped, seu$sample), margin = 2)))
colnames(prop_df) <- c("celltype", "sample", "proportion")
prop_df$celltype <- as.character(prop_df$celltype)
prop_df$sample   <- as.character(prop_df$sample)
prop_df <- prop_df |>
  filter(celltype %in% show_ct) |>
  left_join(sample_info |> mutate(sample = as.character(sample)), by = "sample") |>
  mutate(pct = 100 * proportion,
         celltype = factor(paste0(celltype, "\nFDR = ", signif(da$FDR[match(celltype, da$celltype)], 2)),
                           levels = paste0(show_ct, "\nFDR = ", signif(da$FDR[match(show_ct, da$celltype)], 2))))

p <- ggplot(prop_df, aes(x = injection, y = pct, group = factor(pool), colour = factor(pool))) +
  geom_line(linewidth = 1) +
  geom_point(size = 3.5) +
  facet_wrap(~ celltype, scales = "free_y", nrow = 1) +
  labs(x = NULL, y = "% of the sample's cells", colour = "pool",
       title = "Differential abundance: cell-type proportions per sample") +
  theme(strip.text = element_text(size = 12))
save_plot(p, "da_proportions.png", 14, 5)

## ---- 8. Part B figure: CV^2 = 1/mu + BCV^2 ------------------------------------------------
mu <- 10^seq(0, 3, length.out = 200)                     # pseudobulk count of a gene
cv_df <- rbind(
  data.frame(mu = mu, organism = "Mouse (BCV 0.1)", bcv = 0.1),
  data.frame(mu = mu, organism = "Human (BCV 0.4)", bcv = 0.4)
)
cv_df$total <- 1 / cv_df$mu + cv_df$bcv^2                # technical + biological variation

p <- ggplot(cv_df, aes(x = mu, y = total, colour = organism)) +
  geom_line(linewidth = 1.4) +
  geom_hline(aes(yintercept = bcv^2, colour = organism), lty = 2, show.legend = FALSE) +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values = c("Mouse (BCV 0.1)" = "#1b9e77", "Human (BCV 0.4)" = "#d95f02")) +
  labs(x = "pseudobulk count per gene  (more cells x more depth  →)",
       y = "CV²  (lower = more precise)", colour = NULL,
       title = "More counts help until biological variation (dashed) dominates",
       subtitle = sprintf("This dataset (Mesenchyme): common BCV = %.2f", numbers$common_bcv)) +
  theme(legend.position = "top")
save_plot(p, "power_cv2_mouse_vs_human.png", 10, 6)

## ---- 9. Save the numbers quoted on the slides + session info --------------------------------
saveRDS(numbers, "res/day3_session2_numbers.rds")
writeLines(capture.output(sessionInfo()), "res/sessionInfo_day3_session2.txt")
str(numbers)