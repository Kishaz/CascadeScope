# Provenance / build script for the shipped CascadeScope assets. NOT required to
# run the app (the app already ships the frozen reference and example below).
# Rebuilding needs the reference-cohort expression matrix, the 98-gene panel
# definition, and an example-source cohort, which are NOT distributed, plus
# Bioconductor's org.Hs.eg.db. Point the three paths below at your local copies,
# then run from the repo root:  Rscript make_app_assets.R
#   -> references/ICI_v1.rds   (frozen locked reference + panel + metadata)
#   -> example_cohort.csv      (example upload: a public IMvigor210 subset)
suppressPackageStartupMessages({ library(data.table) })
source("ici_locked_reference.R"); source("score_cohort.R")

## inputs (not distributed; set these to your local copies) -------------------
panel_csv   <- "ICI_98gene_panel.csv"       # panel definition (columns: gene, layer, direction)
ref_expr    <- "reference_baseline.txt.gz"  # reference cohort: genes x samples, HGNC symbols, *_Baseline columns
example_rds <- "example_source.rds"         # a public cohort for the example (list with $E = genes x samples)

PANEL <- load_panel(panel_csv)
PG <- unique(toupper(unlist(PANEL)))

## 1. freeze the locked reference (a validated external baseline cohort) ----
nt <- as.data.frame(fread(ref_expr))
rownames(nt) <- toupper(nt$HGNC_approved_symbol)
refmat <- as.matrix(nt[, grep("_Baseline$", colnames(nt)), drop=FALSE])
REF <- fit_locked_reference(refmat, PANEL, "ICI_locked_v1")

## alias / Ensembl -> canonical-symbol map for the panel genes (org.Hs.eg.db) ----
suppressMessages(library(org.Hs.eg.db))
eg <- AnnotationDbi::mapIds(org.Hs.eg.db, PG, "ENTREZID", "SYMBOL"); eg <- eg[!is.na(eg)]
al <- AnnotationDbi::mapIds(org.Hs.eg.db, eg, "ALIAS",   "ENTREZID", multiVals="list")
en <- AnnotationDbi::mapIds(org.Hs.eg.db, eg, "ENSEMBL", "ENTREZID", multiVals="list")
mp <- do.call(rbind, lapply(names(eg), function(sym){ e<-eg[[sym]]
  syns <- unique(toupper(c(sym, unlist(al[[e]]), unlist(en[[e]])))); syns <- syns[!is.na(syns) & syns!=""]
  data.frame(syn=syns, canon=sym, row.names=NULL) }))
mp <- mp[!(mp$syn %in% names(which(table(mp$syn) > 1))), ]   # drop ambiguous synonyms
alias_map <- setNames(mp$canon, mp$syn)
cat(sprintf("alias/Ensembl map: %d synonyms -> %d panel genes\n", length(alias_map), length(unique(alias_map))))

## reference sub-module score distribution (for the per-patient drill-down) -----
SMref <- submod_scores(refmat, REF)
submod_stats <- data.frame(submodule=colnames(SMref), mean=colMeans(SMref, na.rm=TRUE),
                           sd=apply(SMref, 2, sd, na.rm=TRUE), row.names=NULL)
cat("submodule stats:", paste(submod_stats$submodule, collapse=", "), "\n")

dir.create("references", showWarnings=FALSE, recursive=TRUE)
saveRDS(list(ref=REF, panel=PANEL, panel_genes=PG, alias_map=alias_map, submod_stats=submod_stats,
             version="ICI panel v1", panel_version="v1", n_genes=length(PG),
             layers=names(PANEL), n_ref=ncol(refmat),
             ref_label="NeoTRIP breast baseline (n=241)", built=Sys.Date()),
        "references/ICI_v1.rds")
cat(sprintf("froze reference 'ICI panel v1': %d panel genes, %d reference samples -> references/ICI_v1.rds\n",
            length(PG), ncol(refmat)))

## 2. example cohort (a public IMvigor210 subset: panel genes + random extras) ----
obj <- readRDS(example_rds); E <- obj$E
rownames(E) <- toupper(rownames(E))
set.seed(1)
samp <- sample(colnames(E), 60)
extra <- sample(setdiff(rownames(E), PG), 400)
genes <- intersect(c(PG, extra), rownames(E))
ex <- E[genes, samp, drop=FALSE]
out <- data.frame(gene=rownames(ex), round(ex, 3), check.names=FALSE)
fwrite(out, "example_cohort.csv")
cat(sprintf("example cohort: %d genes x %d samples -> example_cohort.csv\n", nrow(out), ncol(ex)))
