#setwd('~/Projects/09_whiffin_rotation/whiffin-rotation/')
devtools::load_all()
library(data.table)
library(reticulate)
library(argparse)

# load LOEUF deciles
loeuf_path <- "/well/lindgren/flassen/ressources/genesets/genesets/data/gnomad/karczewski2020/supplementary_dataset_11_full_constraint_metrics.tsv"
loeuf <- fread(loeuf_path)
loeuf <- loeuf[,c("gene","transcript","gene_id", "oe_lof_upper_bin")]
loeuf <- loeuf[!is.na(loeuf$oe_lof_upper_bin),]
loeuf$decile <- unlist(lapply(loeuf$oe_lof_upper_bin, function(l) paste0(l,"-",l+1)))
loeuf$decile <- factor(loeuf$decile, levels = paste0(0:9,'-',1:10))
loeuf_mapping <- loeuf$decile 
names(loeuf_mapping) <- loeuf$transcript


# ensure that correct env is used
binary <- "/gpfs3/well/lindgren/users/mmq446/conda/skylake/envs/reticulate/bin/python3"
Sys.setenv(RETICULATE_PYTHON = binary)
print(reticulate::py_config())

# ensure that packages can be loaded
ushuffle <- reticulate::import('ushuffle')
source_python('python/shuffle_utrs.py')

main <- function(args) {

    path_features = args$path_features
    path_mane = args$path_mane
    out_prefix = args$out_prefix
    iterations = as.numeric(args$iterations)
    subset = args$subset
    interval = TRUE # only test 1:10, for example, TRUE=test all
 
    # import 5' or 3' UTR data
    #d <- fread('../../210629_MANE.GRCh38.v0.95.combined-table.txt', sep = '\t')
    d <- fread(path_mane)
    stopifnot(subset %in% d$type)
    stopifnot("enstid_version" %in% colnames(d))
    d <- d[d$type == subset,]
    stopifnot(nrow(d) > 10)
    d$transcript <- gsub("\\.[0-9]+$","",d$enstid_version)
    features <- fread(path_features)
    #features <- fread('derived/tables/210629_MANE.v0.95.UTR_features.txt', sep = '\t')
    #ensgids <- features$ensgid #[features$u5_AUG > 0]

    # codons to evalaute    
    all_codons <- generate_codons()
    codons <- all_codons

    # count how many observed codons are in each sequence ----
    outfile_obs <- paste0(out_prefix, "_observed_codons.csv")
    if (!file.exists(outfile_obs)){
        print('counting observed codons..')
        res_obs <- lapply(d$seq, function(seq){
          mat <- do.call(cbind, lapply(codons, function(codon) count_codon(seq, codon)))
          colnames(mat) <- codons
          return(mat)
        })
        # prettify and export
        res_mat_obs <- as.data.frame(do.call(rbind, res_obs))
        colnames(res_mat_obs) <- paste0('obs.',codons)
        res_mat_obs$ensgid_version <- d$ensgid_version
        res_mat_obs$enstid_version <- d$enstid_version
        fwrite(res_mat_obs, outfile_obs, sep = ',')
    }

    # simulate expected codons given sequence context
    write('simulating expected codons..',stdout())
    i <- as.numeric(args$array_replicate)
    write(paste0(get_time(), ' - Replicate ',i),stdout())
    set.seed(i)
    sim_seq(d$seq[interval], f = function(x) find_codon(x))
    res_expt <- sim_expected_codons(d$seq[interval], k = 2, iter = iterations, codons = codons, parallel = T, seed = i)
  
    # save confidence intervals
    outfile_ci <- paste0(out_prefix, '_codons_expt_ci_rep',i,'.csv')
    res_mat_expt <- as.data.frame(do.call(rbind, res_expt))
    colnames(res_mat_expt) <- paste0('expt.',codons)
    res_mat_expt$ensgid_version <- d$ensgid_version[interval]
    res_mat_expt$enstid_version <- d$enstid_version[interval]
    res_mat_expt$transcript <- d$transcript[interval]
    res_mat_expt$loeuf <- loeuf_mapping[res_mat_expt$transcript]
    fwrite(res_mat_expt, outfile_ci, sep = ',')
  
    # save only estimates
    outfile_est <- paste0(out_prefix, 'codons_expt_rep',i,'.csv')
    mat_split_expt <- as.data.frame(matrixsplit(res_mat_expt, ';', as.numeric, 3))
    mat_split_expt$ensgid_version <- d$ensgid_version[interval]
    mat_split_expt$enstid_version <- d$enstid_version[interval]
    mat_split_expt$transcript <- d$transcript[interval]
    mat_split_expt$loeuf <- loeuf_mapping[mat_split_expt$transcript]
    fwrite(mat_split_expt, outfile_est, sep = ',')
    #})

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--path_mane", default=NULL, help = "?")
parser$add_argument("--path_features", default=NULL, help = "?")
parser$add_argument("--subset", default="five_prime_UTR", help = "?")
parser$add_argument("--iterations", default=5, help = "?")
parser$add_argument("--array_replicate", default=1, help = "?")
parser$add_argument("--out_prefix", default=NULL, help = "?")
args <- parser$parse_args()

main(args)

