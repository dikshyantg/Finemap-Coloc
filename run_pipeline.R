### 📄 File 2: `run_pipeline.R`

```R
#!/usr/bin/env Rscript
# Batch Fine-mapping and Colocalization Pipeline

# Date: 2026-03-1

# --- 1. Global Configuration ---
options(stringsAsFactors = FALSE)
options(warn = -1)

# --- 2. Load Required Packages ---
libs <- c("argparse", "data.table", "dplyr", "ggplot2", "plyr", "devtools", 
          "stringr", "tibble", "coloc", "writexl")

# Install missing packages automatically
for (lib in libs) {
  if (!requireNamespace(lib, quietly = TRUE)) {
    install.packages(lib, repos = "http://cran.us.r-project.org")
  }
  suppressMessages(library(lib, character.only = TRUE))
}

# --- 3. Argument Parsing ---
parser <- ArgumentParser(description="Run Batch Fine-mapping and Coloc")
parser$add_argument("--z_file", required = TRUE, help = "Path to the GWAS Z-file")
parser$add_argument("--qtl_file", default = NULL, help = "Optional: Path to QTL data for coloc")
parser$add_argument("--ld_dir", default = "./ld_matrices/", help = "Directory containing LD matrices")
parser$add_argument("--finemap_exe", default = "finemap", help = "Path to FINEMAP executable")
parser$add_argument("--study_id", default = "study_result", help = "Prefix for output files")
parser$add_argument("--sample_n", type = "integer", default = 482730, help = "GWAS Sample size")
parser$add_argument("--pip_threshold", type = "double", default = 0.6, help = "PIP threshold")
parser$add_argument("--h4_threshold", type = "double", default = 0.6, help = "Coloc H4 threshold")
parser$add_argument("--rsid_col", default = "ID", help = "Column name for RSIDs")
parser$add_argument("--freq_col", default = "no_freq", help = "Column name for allele frequency")
parser$add_argument("--locus_window", type = "integer", default = 1500000, help = "Window size (+/-)")

args <- parser$parse_args()

# --- 4. Data Processing & Standardization ---
cat("Loading GWAS data from:", args$z_file, "\n")
gwas.raw <- fread(args$z_file)
col_names <- colnames(gwas.raw)

if (args$freq_col == "no_freq" || !(args$freq_col %in% col_names)) {
    gwas.raw$freq <- NA
    cat("Note: Frequency column missing, setting to NA\n")
}

gwas.dt <- gwas.raw %>% 
    rename(rsid = !!sym(args$rsid_col)) %>%
    mutate(
        sample_n = args$sample_n,
        allele1 = toupper(A1),
        allele2 = toupper(A2)
    ) %>%
    select(rsid, chr, pos, allele1, allele2, beta, se, p, freq, sample_n)

# QC: Remove missing values
i1 <- nrow(gwas.dt)
gwas.dt <- gwas.dt[complete.cases(gwas.dt[, c("chr", "pos", "beta", "se", "p")]), ]
cat((i1 - nrow(gwas.dt)), "SNPs removed due to missing values\n")

# Load QTL data if provided
qtl.dt <- NULL
if (!is.null(args$qtl_file)) {
    cat("Loading QTL data from:", args$qtl_file, "\n")
    qtl.dt <- fread(args$qtl_file)
}

# --- 5. Locus Partitioning Logic (Greedy Clumping) ---
get_independent_loci <- function(gwas_data, window, p_thresh = 5e-8) {
    sig_snps <- gwas_data %>% filter(p < p_thresh) %>% arrange(p)
    loci <- list()
    
    while(nrow(sig_snps) > 0) {
        lead_snp <- sig_snps[1, ]
        locus_name <- paste0("locus_", lead_snp$chr, "_", lead_snp$pos)
        
        # Define window boundaries
        start_pos <- lead_snp$pos - window
        end_pos <- lead_snp$pos + window
        
        # Extract locus data
        locus_data <- gwas_data %>% 
            filter(chr == lead_snp$chr, pos >= start_pos, pos <= end_pos)
        
        loci[[locus_name]] <- list(name = locus_name, lead_rsid = lead_snp$rsid, data = locus_data)
        
        # Remove these SNPs from the pool to find the next independent locus
        sig_snps <- sig_snps %>% filter(!(rsid %in% locus_data$rsid))
    }
    return(loci)
}

# --- 6. Finemapping & Coloc Execution ---
run_batch_locus <- function(locus_info, qtl_data, args) {
    locus_name <- locus_info$name
    locus_data <- locus_info$data
    cat("\nProcessing:", locus_name, "| SNPs:", nrow(locus_data), "\n")
    
    # 6a. Generate FINEMAP inputs
    z_file <- paste0(locus_name, ".z")
    master_file <- paste0(locus_name, ".master")
    ld_file <- file.path(args$ld_dir, paste0(locus_info$lead_rsid, "_ld.txt"))
    
    # Format Z file for FINEMAP (rsid, chromosome, position, allele1, allele2, maf, beta, se)
    finemap_in <- locus_data %>% 
        mutate(maf = ifelse(is.na(freq), 0.05, freq)) %>% # Dummy MAF if missing, FINEMAP needs it
        select(rsid, chromosome = chr, position = pos, allele1, allele2, maf, beta, se)
    
    fwrite(finemap_in, z_file, sep = " ")
    
    # Create master file
    master_content <- paste("z;ld;snp;config;cred;log;n_samples",
                            paste(z_file, ld_file, paste0(locus_name, ".snp"), 
                                  paste0(locus_name, ".config"), paste0(locus_name, ".cred"), 
                                  paste0(locus_name, ".log"), args$sample_n, sep=";"), 
                            sep="\n")
    writeLines(master_content, master_file)
    
    # 6b. Run FINEMAP
    finemap_res <- data.frame()
    if (file.exists(ld_file)) {
        cmd <- paste(args$finemap_exe, "--sss --in-files", master_file, "--n-causal-snps 5")
        system(cmd, ignore.stdout = TRUE)
        
        if (file.exists(paste0(locus_name, ".snp"))) {
            finemap_res <- fread(paste0(locus_name, ".snp")) %>% filter(prob > args$pip_threshold)
        }
    } else {
        cat("  -> Warning: LD file not found (", ld_file, "), skipping FINEMAP.\n")
    }
    
    # 6c. Run Coloc (if QTL data available)
    coloc_res <- NULL
    if (!is.null(qtl_data)) {
        qtl_locus <- qtl_data %>% filter(rsid %in% locus_data$rsid)
        
        if (nrow(qtl_locus) > 50) { # Ensure enough overlapping SNPs
            dataset1 <- list(beta = locus_data$beta, varbeta = locus_data$se^2, 
                             type = "cc", snp = locus_data$rsid, N = args$sample_n)
            
            dataset2 <- list(beta = qtl_locus$beta, varbeta = qtl_locus$se^2, 
                             type = "quant", snp = qtl_locus$rsid, N = max(qtl_locus$N, na.rm=TRUE))
            
            res <- coloc.abf(dataset1, dataset2)
            if (res$summary["PP.H4.abf"] > args$h4_threshold) {
                coloc_res <- data.frame(locus = locus_name, H4 = res$summary["PP.H4.abf"])
            }
        }
    }
    
    # Cleanup temp files
    unlink(c(z_file, master_file, paste0(locus_name, ".*")))
    
    return(list(finemap = finemap_res, coloc = coloc_res))
}

# --- 7. Execution Loop & Output ---
loci_list <- get_independent_loci(gwas.dt, args$locus_window)
cat("Identified", length(loci_list), "independent loci.\n")

all_finemap <- list()
all_coloc <- list()

for (locus in loci_list) {
    results <- run_batch_locus(locus, qtl.dt, args)
    if (nrow(results$finemap) > 0) all_finemap[[locus$name]] <- results$finemap
    if (!is.null(results$coloc)) all_coloc[[locus$name]] <- results$coloc
}

# Combine results
df_finemap <- bind_rows(all_finemap, .id = "Locus")
df_coloc <- bind_rows(all_coloc)

# Write outputs
cat("\nWriting results...\n")
if(nrow(df_finemap) > 0) {
    fwrite(df_finemap, paste0(args$study_id, "_credible_set_combine.txt"), sep="\t")
}

if(nrow(df_coloc) > 0) {
    fwrite(df_coloc, paste0(args$study_id, "_finemap_coloc_combine_summary.csv"))
    write_xlsx(df_coloc, paste0(args$study_id, "_causal_gene_summary.xlsx"))
} else {
    cat("No significant colocalization results to write.\n")
}

cat("Pipeline completed successfully!\n")