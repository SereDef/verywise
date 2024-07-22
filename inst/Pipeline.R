library(verywise)
library(bigreadr)

hemi = "lh"
measure = "thickness"
fwhm = 10
target = "fsaverage"
n_cores = 1
subj_dir = "/Users/Serena/Desktop/Infant2adult/Package/try"

# Simulate FreSurfer and phenotype dataset
# simulate_data(subj_dir)

# Read phenotype data
if (!exists("pheno")){
  pheno <- read.csv(file.path(subj_dir, "phenotype.csv"))
}

# Read vertex data
ss_backfile <- file.path(subj_dir,
                         paste0(hemi, ".", measure, ".supersubject"))

if (!file.exists(paste0(ss_backfile,'.bk'))) {
  ss <- build_supersubject(subj_dir,
                           folder_id = pheno$folder_id,
                           files_list = list.dirs.till(subj_dir, n = 2),
                           measure = measure,
                           hemi = hemi,
                           fwhmc = paste0("fwhm",fwhm),
                           target = target,
                           backing = ss_backfile,
                           n_cores = n_cores)
} else {
  ss <- bigstatsr::big_attach(bigstatsr::sub_bk(ss_backfile, ".rds"))
}

# Mask vertex data
masked <- mask_cortex(ss, hemi = hemi, target = target, n_cores = n_cores)

# Prepare chunk sequence
iv <- which(!masked)
chunk_seq <- make_chunk_sequence(iv)
