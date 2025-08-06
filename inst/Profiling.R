library(verywise)
library(foreach)

subj_dir <- test_path("fixtures", "fs7")
pheno <- read.csv(file.path(subj_dir, "phenotype.csv"))

n_obs <- nrow(pheno)
test_formula <- vw_area ~ sex + age + (1 | id)
fs_home = "/Applications/freesurfer/7.4.1"

outp_dir <- file.path(subj_dir,'results')

pheno_list <- list(pheno, pheno, pheno)

system.time({
  result <- run_vw_lmm(
    formula = test_formula,
    pheno = pheno_list,
    folder_id = 'folder_id',
    subj_dir = subj_dir,
    outp_dir = outp_dir,
    hemi = "lh",
    seed = 123,
    n_cores = 4,
    FS_HOME = fs_home,
    # prioritize speed over accuracy
    lmm_control = lme4::lmerControl(calc.derivs = FALSE,
                                    use.last.params = TRUE,
                                    check.rankX = "ignore",
                                    check.nobs.vs.rankZ = "ignore",
                                    check.nobs.vs.nlev = "ignore",
                                    check.nlev.gtreq.5 = "ignore",
                                    check.nlev.gtr.1 = "ignore",
                                    check.nobs.vs.nRE = "ignore",
                                    check.formula.LHS = "ignore",
                                    check.scaleX = "ignore",
                                    check.conv.grad = "ignore",
                                    check.conv.singular = "ignore",
                                    check.conv.hess = "ignore"),
    use_model_template = TRUE,
    verbose = TRUE
  )
})
