#!/usr/bin/env Rscript
#-------------------------------------------------------------------------------
#-This R script implements the ProMix tool for improving design and normalization of multiplex proteomics study.
#-Contact: Hua Tang (huatang@stanford.edu), Huaying Fang (hyfang@cnu.edu.cn)
#-Date: 20141214
#-Import: R package lme4
#-------------------------------------------------------------------------------
require(lme4);
#---------------------------------------
promix_wrap = function(file_in = "input_promix.txt", file_out = "output_promix.txt") {
  data_in = read.table(file_in, header = T, sep = " ");
  #-Input file should have no less than 4 columns and the first three columns are Samp, Subj and Run.
  if(!all(c("Samp", "Subj", "Run") == names(data_in)[1:3])) {
	  print("The first three columns of input file should be Samp, Subj and Run!");
	  quit("no");
  };
  design_dat = data_in[, c("Samp", "Subj", "Run")];
  prt_dat = t(data.matrix(data_in[, -(1:3)]));
  colnames(prt_dat) = design_dat$Samp;
  medincl_prt = rownames(prt_dat);
  lmmincl_prt = rownames(prt_dat);
  mega_id = names(which.max(table(design_dat$Subj)));
  output = promix_func(prt_dat, design_dat, medincl_prt, lmmincl_prt, mega_id = mega_id);
  IDsubj = setdiff(sort(unique(design_dat$Subj)), mega_id);
  output = data.frame(Subj = IDsubj, t(output[, IDsubj]));
  write.table(output, file = file_out, quote = F, sep = " ", row.names = F, col.names = T);
};
#---------------------------------------
#-ProMix: Improving design and normalization of multiplex proteomics study
#-Input:
#-  prt_dat: a Nprt X Nsamp data matrix with rownames and column names; rownames is protein names; colnames is sample names. 
#-  design_dat: a Nsamp X 3 data frame; each row represents a sample; names are Samp for sample id (SA1-SA_Nsamp), Subj for subject id (SU0 for mega-pool, SU1-SU_Nsubj), Run for run id (R1-R_Nrun).
#-  Nsamp = Nrun + Nsubj X 2 (each subject for SU1-SU_Nsubj has 2 replicates while mega-pool has Nrun replicates)
#-  medincl_prt: subset of proteins for calculating the median of each sample
#-  lmmincl_prt: subset of proteins for fitting linear mixed effect models
#-  mega_id: Subject id for mega-pool
#-  niter: iteration time for ProMix
promix_func = function(prt_dat, design_dat, medincl_prt, lmmincl_prt, mega_id = "SU0", niter = 5) {
  prt_dat = prt_dat[, design_dat$Samp];
  Nprt = nrow(prt_dat);
  Nsamp = ncol(prt_dat);
  all_prt = rownames(prt_dat);
  IDsamp = design_dat$Samp;
  IDsubj1 = unique(design_dat$Subj);
  IDrun = unique(design_dat$Run);
  Nsubj1 = length(IDsubj1); #-Nsubj1 = Nsubj + 1
  Nrun = length(IDrun);
  samp2run = with(design_dat, setNames(Run, Samp));
  samp2subj =  with(design_dat, setNames(Subj, Samp));
  runid_samp = samp2run[IDsamp];
  subjid_samp = samp2subj[IDsamp];
  sampid_mega = with(design_dat, Samp[Subj == mega_id]);
  #-Run effect matrix: Nprt X Nrun
  b_gr = matrix(0, nrow = Nprt, ncol = Nrun);
  rownames(b_gr) = all_prt;
  colnames(b_gr) = IDrun;
  #-Subject effect matrix: Nprt X Nsubj1
  #-mu_g is implicitly included in gamma_gl (gamma_gl[, mega_id]) (intercept term for fixed effect in linear mixed effect models)
  gamma_gl = matrix(0, nrow = Nprt, ncol = Nsubj1);
  rownames(gamma_gl) = all_prt;
  colnames(gamma_gl) = IDsubj1;
  #-Iteration for ProMix
  lmm_dat = cbind(design_dat[, c("Samp", "Subj", "Run")], prt = NA);
  lmm_dat$Subj = with(lmm_dat, relevel(factor(Subj), ref = mega_id));
  for(iter in 1:niter) {
    mu_g = gamma_gl[, mega_id];
    gamma_gl[, mega_id] = 0;
    prt_dat1 = prt_dat - (mu_g + b_gr[, runid_samp] + gamma_gl[, subjid_samp]);
    #-1. Calculate alpha_k
    alpha_k = apply(prt_dat1[medincl_prt, ], 2, median, na.rm = T);
	alpha_k = alpha_k - mean(alpha_k);
    zgk = prt_dat - rep(alpha_k, each = Nprt);
    #-Remove outliers in mega-pool: threshold of z-score is 3 for outliers.
    thresh = 3;
    zgk_mega = zgk[, sampid_mega];
    rowmu = apply(zgk_mega, 1, mean, na.rm = T);
    rowsd = apply(zgk_mega, 1, sd, na.rm = T);
    rowz = (zgk_mega - rowmu) / (rowsd * 0.8 + quantile(rowsd, 0.5, na.rm = T) * 0.2);
    zgk_mega[abs(rowz) > thresh] = NA;
    zgk[, sampid_mega] = zgk_mega;
    #-2. Fit linear mixed effect models
    for(prti in lmmincl_prt) {
      lmm_dat$prt = zgk[prti, ];
      lmm_fit = suppressMessages(lmer(prt ~ (1|Run) + Subj, data = lmm_dat, control = lmerControl(optimizer = "Nelder_Mead")));
      #-Random effects for runs
      eff_run = ranef(lmm_fit)$Run;
      b_gr[prti, rownames(eff_run)] = eff_run[[1]];
      #-Fixed effects for subjects
      eff_subj = fixef(lmm_fit);
      subj_name = gsub("Subj", "", names(eff_subj));
      gamma_gl[prti, subj_name[-1]] = eff_subj[-1];
	  gamma_gl[prti, mega_id] = eff_subj[1];
    };
  };
  lmmexcl_prt = setdiff(all_prt, lmmincl_prt);
  if(length(lmmexcl_prt) > 0) {
    #-PQN for proteins not used in linear mixed model
    gamma_gl[lmmexcl_prt, ] = pqn_func(zgk[lmmexcl_prt, ], design_dat, medincl_prt = NULL, mega_id = mega_id);
  };
  return(gamma_gl[, setdiff(IDsubj1, mega_id)]);
};
#---------------------------------------
average2mats_func = function(mat1, mat2, name_rc = NULL) {
  if(length(name_rc) == 0) {
    name_rc = list(rownames(mat1), colnames(mat1));
  };
  outdat = (mat1 + mat2)/2;
  naind1 = is.na(mat1);
  naind2 = is.na(mat2);
  outdat[!naind1 & naind2] = mat1[!naind1 & naind2];
  outdat[naind1 & !naind2] = mat2[naind1 & !naind2];
  rownames(outdat) = name_rc[[1]];
  colnames(outdat) = name_rc[[2]];
  return(outdat);
};
#---------------------------------------
#-Probability quotient normalization (PQN)
pqn_func = function(prt_dat, design_dat, medincl_prt, mega_id = "SU0") {
    if(length(medincl_prt) == 0) {
      subset_median = 0;
    } else {
      subset_median = apply(prt_dat[medincl_prt, ], 2, median, na.rm = T);
    };
    IDsubj1 = unique(design_dat$Subj);
    Nprt = nrow(prt_dat);
    prt_dat_ = prt_dat - rep(subset_median, each = Nprt);
    sampid_mega = with(design_dat, Samp[Subj == mega_id]);
    prt_mega = rowMeans(prt_dat_[, sampid_mega], na.rm = T);
    design_nomega = design_dat[design_dat$Subj != mega_id, ];
    design_nomega = design_nomega[order(design_nomega$Subj), ];
    id1 = 2*(1:(nrow(design_nomega)/2));
    id2 = id1 - 1;
    prt_nomega = average2mats_func(mat1 = prt_dat_[, design_nomega$Samp[id1]], mat2 = prt_dat_[, design_nomega$Samp[id2]], name_rc = list(rownames(prt_dat), design_nomega$Subj[id1]));
    outdat = cbind(prt_mega, prt_nomega);
    colnames(outdat)[1] = mega_id;
    return(outdat[, IDsubj1]);
};
#---------------------------------------#---------------------------------------
#-Arguments from command line: file_in file_out
args_cmd = commandArgs(trailingOnly = T);
if(length(args_cmd) == 0) {
	print("Example is used!");
	file_in = "input_promix.txt";
	file_out = "output_promix.txt"
};
if(length(args_cmd) == 1) {
	file_in = args_cmd[1];
	file_out = "output_promix.txt"	
};
if(length(args_cmd) >= 2) {
	file_in = args_cmd[1];
	file_out = args_cmd[2];
};
promix_wrap(file_in = file_in, file_out = file_out); 
quit("no");
#---------------------------------------#---------------------------------------
#---------------------------------------