library(tidyverse)
library(fgem)
library(ggrepel)
library(patchwork)
library(tidyverse)
library(fgem)
library(firatheme)

file_n <- fs::dir_ls("data/new_cv_gobp_models_conservative_relaxed/")
safe_r <- safely(qs::qread)
fit_l <- map(file_n, safe_r)
fit_df <- map(fit_l, "result") %>%
  compact() %>%
  bind_rows() %>%
  group_by(cancer)

fit_df <- map_df(file_n, qs::qread) %>%
  group_by(cancer)

sfit_df <- fgem:::summarise_cv_lik(fit_df)
sfit_df <- filter(sfit_df, l2 == 0) %>%
    transmute(cancer = cancer,
              cv_sum0 = cv_sum,
              l0mean0 = l0_mean,
              l1mean0 = l1_mean,
              l2mean0 = l2_mean) %>%
  inner_join(sfit_df)
sfit_df <- sfit_df %>% mutate(is_max_cv = cv_sum == max(cv_sum), max_l2 = l2[is_max_cv])

GO_def <- readRDS("data/GO_Definitions.RDS") %>%
  select(feature_name, description = Term)

## filter(sfit_df,l2>0) %>%
## ggplot(aes(x=l2_mean,y=cv_sum)) +
##   geom_point(aes(col=is_max_cv)) +
##   facet_wrap(~cancer, scales = "free") +
##   scale_x_log10()

## filter(sfit_df) %>%
##   ggplot(aes(x=l2,y=cv_sum)) +
##   geom_point(aes(col=is_max_cv)) +
##   facet_wrap(~cancer, scales = "free") +
##   geom_hline(aes(yintercept=cv_sum0)) +
##   scale_x_log10()

## unnest(sfit_df,Beta) %>%
##   group_by(cancer, l2, is_max_cv) %>%
##   filter(feature_name!="Intercept") %>%
##   summarise(max_abs_Beta=max(abs(Beta))) %>% filter(l2>0) %>%
##   ggplot(aes(x=l2,y=max_abs_Beta))+geom_point(aes(col=is_max_cv))+facet_wrap(~cancer,scales="free")+scale_x_log10()

## unnest(sfit_df,Beta) %>%
##   filter(feature_name!="Intercept")  %>%
##   filter(l2>0 | is_max_cv) %>%
##   ggplot(aes(x=l2,y=Beta,col=feature_name))+geom_line()+scale_x_log10()+facet_wrap(~cancer,scales="free")+theme(legend.position = "none")+geom_vline(aes(xintercept=max_l2))

## sfit_df %>%
##   mutate(is_max_cv=cv_sum == max(cv_sum)) %>%
##   filter(is_max_cv) %>%
##   unnest(Beta)  %>%
##   filter(Beta != 0) %>%
##   filter(feature_name != "Intercept") %>%
##   filter(abs(Beta) == max(abs(Beta)))  %>%
##   arrange(desc(abs(Beta)))  %>%
##   inner_join(GO_def) %>%
##   select(-is_max_cv) %>%
##   as.data.frame()


ocancer_df <- qs::qread("data/cleanest_driverMAPS_results_20TumorTypes.qs")
tcancer_df <- ocancer_df %>%
  select(Gene, cancer, BF) %>%
  filter(cancer == cancer[1])
fgnf <- fgem:::fgem_null_fit(tcancer_df$BF, log_BF = TRUE)
tpr <- fgem:::predict_fgem(fgnf, fgem:::empty_x(tcancer_df$BF, tcancer_df$Gene), log = TRUE)
xb <- stats::qlogis(tpr, log.p = TRUE)
tbf <- fgem:::log_BF2_logpost(tcancer_df$BF, xb)
cancer_df <- ocancer_df %>%
  select(Gene, cancer, BF) %>%
  group_by(cancer) %>%
  mutate(uniform_posterior = c(BF2posterior(BF, Gene, prior = fgem:::predict_uniform_prior(BF, Gene, log_BF = TRUE), log_BF = TRUE)))

all_genes <- unique(cancer_df$Gene)

cancer_dfl <- group_nest(cancer_df, .key = "data")



go_df <- qs::qread("data/GO_df.qs")
long_df <- bind_rows(go_df)
blong_df <- bind_rows(
  long_df,
  tibble(Gene = all_genes, feature_name = "Intercept")
)
all_def <- bind_rows(mutate(GO_def, category = "GO", sub_category = "BP"))
xMS <- fgem:::df2sparse(long_df, total_rownames = all_genes)

comp_df <- read_tsv("data/2020-02-02_IntOGen-Drivers-20200213/Compendium_Cancer_Genes.tsv")
cohort_df <- read_tsv("data/2020-02-02_IntOGen-Cohorts-20200213/cohorts.tsv") %>%
    mutate(is_TCGA = SOURCE == "TCGA")

is_tcga <- select(cohort_df, cohort = COHORT, is_TCGA) %>%
    distinct()

cohort_tcga <- filter(cohort_df, is_TCGA) %>%
  rename(cohort = COHORT) %>%
  mutate(cancer = str_remove(cohort, "TCGA_W.S_")) %>%
  select(intogen_cancer = CANCER_TYPE, cancer)
# cohort_notcga <- filter(cohort_df,SOURCE!="TCGA") %>% rename(cohort=COHORT)

sub_comp_df <- select(comp_df,
                      Gene = SYMBOL,
                      cohort = COHORT,
                      intogen_cancer = CANCER_TYPE,
                      qval_cons = QVALUE_COMBINATION) %>%
  inner_join(is_tcga) %>%
  left_join(cohort_tcga) %>%
  mutate(cancer = if_else(is.na(cancer), intogen_cancer, cancer)) %>%
  select(-intogen_cancer) %>%
  distinct(cancer, Gene, is_TCGA)


intogen_l <- nest(sub_comp_df,
    intogen_l = -cancer
)
qs::qsave(intogen_l, "data/intogen_df.qs")

assign_intogen <- function(data, cancer, include_tcga = TRUE) {
  if ("intogen_gene" %in% colnames(data)) {
    return(data)
  }
  if (!include_tcga) {
    scdf <- filter(sub_comp_df, !is_TCGA)
  } else {
    scdf <- sub_comp_df
  }
  intogen_l <- mutate(scdf, intogen_gene = if_else(cancer == {{ cancer }},
    if_else(is_TCGA, 2L, 1L),
    if_else(is_TCGA, 4L, 3L)
  )) %>%
    group_by(Gene) %>%
    summarise(intogen_gene = min(intogen_gene), .groups = "drop")
  rdf <- left_join(data, intogen_l, by = "Gene") %>%
    replace_na(list(intogen_gene = 5L)) %>%
    group_by(Gene) %>%
    filter(intogen_gene == min(intogen_gene)) %>%
    slice(1) %>%
    ungroup()
  return(rdf)
}

kt <- c(
  "Known Type-Specific Cancer Gene (Non-TCGA)",
  "Known Type-Specific Cancer Gene (TCGA)",
  "Known Cancer Gene (Non-TCGA)",
  "Known Cancer Gene (TCGA)",
  "Unknown Cancer Gene"
)

skt <- c(
  "Type-Specific(Non-TCGA)",
  "Type-Specific(TCGA)",
  "Generic(Non-TCGA)",
  "Generic(TCGA)",
  "Unknown"
)
igdf <- tibble(
  intogen_gene = 1:5,
  IntogenClass = skt
) %>% mutate(
  TypeSpecific = str_detect(IntogenClass, "Type-Specific"),
  isCancer = intogen_gene < 5,
  is_non_tcga = intogen_gene %in% c(1L, 3L)
)


group_mapn <- function(.x, .f, ...) {
  gv <- group_vars(.x)
  .x <- ungroup(.x)
  split(.x, .x[[gv]]) %>% map(.f = .f, ... = ...)
}

cancer_dfl <- nest(cancer_df, data = -cancer) %>% mutate(data = map2(data, cancer, assign_intogen, include_tcga = TRUE))

add_pred <- function(beta_df, cancer_df, log = FALSE) {
  ret_df <- mutate(cancer_df,
    functional_prior = c(fgem:::predict_long_df(beta_df, xMS[cancer_df$Gene, ], log = log)),
    functional_posterior = c(fgem:::BF2posterior(BF, Gene, functional_prior, log_BF = log))
  )
}

log_pbinom <- function(lp, x) {
  lp <- pmin(lp, -.Machine$double.eps)
  sum(x * lp + (1 - x) * log(expm1(-lp)))
}


plot_post <- function(data,
                      log = TRUE,
                      func_post_cutoff = 0.5,
                      post_cutoff = 0.99,
                      strict_func_post_cutoff = 0.9) {
  pdata <- filter(data) %>%
    arrange(desc(functional_posterior - uniform_posterior)) %>%
    mutate(
      idn = 1:n(),
      `intOGen Validation Status` = factor(kt[intogen_gene], levels = kt)
    )
  if (!log) {
      pdata <- mutate(pdata,
                      functional_posterior = exp(functional_posterior),
                      uniform_posterior = exp(uniform_posterior)
                      )
  } else {
      func_post_cutoff <- log(func_post_cutoff)
      post_cutoff <- log(func_post_cutoff)
      strict_func_post_cutoff <- log(strict_func_post_cutoff)
  }
  top_diff <- filter(
      pdata,
      functional_posterior > func_post_cutoff
  ) %>%
      slice(1:25)
  top_g <- filter(pdata,
                  idn >= 15,
                  functional_posterior > post_cutoff | uniform_posterior > post_cutoff) %>%
      group_by(`intOGen Validation Status`) %>%
      arrange(desc(pmax(functional_posterior, uniform_posterior))) %>%
      slice(1:15) %>%
      ungroup()
  ttg <- filter(pdata, intogen_gene == 5, functional_posterior > strict_func_post_cutoff)
  top_g <- bind_rows(top_g, ttg) %>% distinct()

  rp <- pdata %>%
    ggplot(aes(
      x = uniform_posterior,
      y = functional_posterior,
      col = `intOGen Validation Status`,
      label = Gene
    )) +
      geom_label_repel(data = top_diff, max.overlaps = 10) +
      geom_label_repel(data = top_g, max.overlaps = 10) +
      geom_point() +
      geom_abline(slope = 1)
  return(rp)
}


plot_beta <- function(Beta, cancer,width=25,k=10,description=FALSE) {
  tBeta <- filter(Beta, feature_name != "Intercept", Beta != 0)
  nr <- ceiling(nrow(tBeta) / k)
  BR <- range(tBeta$Beta)
  l <- nrow(tBeta)
  if (description) {
      Beta_p <- tBeta %>%
          arrange(abs(Beta)) %>%
          mutate(
              label = str_wrap(paste0(feature_name, ":\n", description), width = width),
              label = factor(label, levels = label), Beta_ct = gl(n = nr, k = k, length = l),
              feature_name = factor(feature_name, levels = feature_name),
              idx = 1:n()
          )
  } else {
      Beta_p <- tBeta %>%
          arrange(abs(Beta)) %>%
          mutate(
              label = feature_name,
              label = factor(label, levels = label), Beta_ct = gl(n = nr, k = k, length = l),
              feature_name = factor(feature_name, levels = feature_name),
              idx = 1:n()
          )
  }
  
  ggplot(
    tail(Beta_p, k),
    aes(y = label, x = Beta)
  ) +
    geom_point() +
    ggtitle(cancer)
}


plt_fun <- function(Beta, data, cancer, l0_mean, ...) {
  l0n <- l0_mean
  pta <- proc.time()
  out_png_e <- fs::path(paste0("org/conservative_", "gobp", "_enrichment"))
  fs::dir_create(out_png_e)

  out_png_p <- fs::path(paste0("org/conservative_", "gobp", "_posterior"))
  fs::dir_create(out_png_p)
  lout_png_p <- fs::path(out_png_p, paste0(cancer, "_log"), ext = "png")
  out_png_p <- fs::path(out_png_p, cancer, ext = "png")

  out_png_e <- fs::path(out_png_e, cancer, ext = "png")

  pa <- plot_beta(Beta, cancer) + theme_fira() + plot_annotation(title = cancer, subtitle = paste0("L0: ", as.integer(l0n)))
  ggsave(out_png_e, pa, dpi = 320, width = 13.33, height = 7.5)

  pb <- plot_post(data, log = FALSE) + theme_fira() + plot_annotation(title = cancer, subtitle = paste0("L0: ", as.integer(l0n)))
  ggsave(out_png_p, pb, dpi = 320, width = 13.33, height = 7.5)
  lpb <- plot_post(data, log = TRUE) + theme_fira() + plot_annotation(title = cancer, subtitle = paste0("L0: ", as.integer(l0n)))
  ggsave(lout_png_p, lpb, dpi = 320, width = 13.33, height = 7.5)
  return(proc.time() - pta)
}

bf <- filter(sfit_df,
             is_max_cv) %>%
  inner_join(cancer_dfl) %>%
  mutate(
    data = map2(Beta, data, add_pred, log = TRUE),
    Beta = map(Beta, inner_join, y = all_def, by = "feature_name")
  ) %>%
  dplyr::mutate(data = map2(data, cancer, assign_intogen, include_tcga = TRUE))
qs::qsave(bf,"data/plot_data.qs")
                                        # library(furrr)
                                        # plan(multisession)
                                        
ucec_bf_df <- filter(bf, cancer == "UCEC") 
ucec_l0n <- ucec_bf_df$l0_mean
plot_ucec <- plot_post(ucec_bf_df$data[[1]], log = FALSE) +
    theme_fira() + ggtitle("UCEC")




brca_bf_df <- filter(bf, cancer == "BRCA")
brca_l0n <- brca_bf_df$l0_mean
plot_brca <- plot_post(brca_bf_df$data[[1]], log = FALSE)  +
    theme_fira() +
    ggtitle("BRCA") 

(plot_brca + plot_ucec) + plot_layout(guides = "collect")

pb_brca <- plot_beta(brca_bf_df$Beta[[1]], "BRCA",20,25)
pb_ucec <- plot_beta(ucec_bf_df$Beta[[1]],"UCEC",20,25)
## (plot_beta(brca_bf_df$Beta[[1]], "BRCA",25) +
##  plot_beta(ucec_bf_df$Beta[[1]],"UCEC",25))

plots_post <- ((plot_brca + plot_ucec) + plot_layout(guides = "collect"))

plots_enr <- (pb_brca + pb_ucec) & theme_fira()

ggsave("org/fgem_posterior_plot.png",plot=plots_post,width=2250/300,height=(2625/300)/2,units="in")
ggsave("org/fgem_enrichment_plot.png",plot=plots_enr,width=2250/300,height=(2625/300)/2,units="in")
                                        #/ plots_enr



    
    ## plot_annotation(title = cancer)

pwalk(bf, plt_fun)

cancer_null <- group_by(cancer_df, cancer) %>%
    summarise(null_lik = fgem:::fgem_null_lik(BF, log_BF = TRUE),
              Intercept = fgem:::fgem_null(BF, log_BF = TRUE))

feat_df <- qs::qread("data/single_df.qs") %>%
    select(feature_name,
           cancer,
           univariate_p = pval_fgem,
           fisher_p = pval_fisher_0_1,
           univariate_q = qval_fgem
           )

bfdf <- select(feat_df, cancer, feature_name, contains("fisher")) %>% gather(key = "cutoff", value = "pval", -cancer, -feature_name)

group_by(bfdf, cancer, cutoff) %>%
  mutate(pa = p.adjust(pval, "fdr")) %>%
  group_by(cutoff) %>%
  summarise(n_sig = mean(pa < 0.05)) %>%
  arrange(desc(n_sig))

## feat_df <- map_df(fs::dir_ls("data/fgem_m_results"), qs::qread) %>%
##   inner_join(select(cancer_null,-Intercept)) %>%
##   mutate(pval=stats::pchisq(-2 * (null_lik - (lik)), df=1, lower.tail=F)) %>%
##   group_by(cancer) %>%
##   mutate(qval = p.adjust(dplyr::if_else(convergence == 0, pval, 1), method = "fdr")) %>%
##   unnest(Beta) %>%
##   filter(feature_name != "Intercept") %>%
##

head(feature_odf)

post_df <- unnest(select(bf, cancer, data), data) %>%
  mutate(
    functional_posterior = exp(functional_posterior),
    uniform_posterior = exp(uniform_posterior),
    post_ratio = functional_posterior / uniform_posterior,
    post_diff = functional_posterior - uniform_posterior,
    rank_diff = rank(functional_posterior) - rank(uniform_posterior)
  ) %>%
  inner_join(igdf)
pda <- post_df %>%
  group_by(cancer) %>%
  summarise(ave_diff_overall = mean(functional_posterior - uniform_posterior))
ggplot(post_df, aes(x = cancer, y = rank_diff, col = cancer, fill = cancer)) +
  facet_wrap(~IntogenClass, scales = "free_y") +
  geom_violin()
filter(post_df, is_non_tcga | !isCancer) %>%
  mutate(is_ooc = IntogenClass == "Type-Specific(Non-TCGA)") %>%
  mutate(
    rank_unif = rank(uniform_posterior),
    rank_func = rank(functional_posterior),
    rank_diff = rank_func - rank_unif
  ) %>%
  group_by(cancer, is_ooc) %>%
  summarise(ave_diff = mean()) %>%
  spread(key = "is_ooc", value = "ave_diff") %>%
  arrange(desc(`TRUE` - `FALSE`))
spread(key = "intogen_gene", value = "ave_diff") %>%
  write_tsv("data/posterior_diff.tsv")

huge_bf <- unnest(select(bf, -data), Beta) %>%
  filter(Beta != 0) %>%
  inner_join(blong_df) %>%
  inner_join(unnest(select(filter(bf), -Beta), data)) %>%
  mutate(intogen_gene = kt[intogen_gene])

feature_odf <- unnest(select(bf, -data), Beta) %>%
  filter(Beta != 0) %>%
  select(cancer, feature_name, Beta, category, sub_category, description) %>%
  inner_join(feat_df) %>%
  arrange(cancer, desc(abs(Beta)))
BPGOdf <- select(BPGOdf, Gene, feature_name = feature)


feature_odf <- inner_join(select(feature_odf, cancer, feature_name), BPGOdf) %>%
  inner_join(cancer_df) %>%
  mutate(uniform_posterior = exp(uniform_posterior)) %>%
  group_by(cancer, feature_name) %>%
  summarise(
    num_genes = n_distinct(Gene),
    n_sig = sum(uniform_posterior > 0.9)
  ) %>%
  inner_join(feature_odf)

feature_odf %>%
  write_tsv("data/all_conservative_features.tsv")

filter(huge_bf, Beta != 0) %>%
  filter(Beta < 0)

w_bf <- distinct(huge_bf, Gene, cancer, BF, functional_prior, functional_posterior, uniform_posterior, intogen_gene) %>%
  mutate(
    uniform_posterior = exp(uniform_posterior),
    functional_posterior = exp(functional_posterior),
    functional_prior = exp(functional_prior)
  ) %>%
  inner_join(
    filter(huge_bf, Beta != 0) %>%
      group_by(cancer, Gene, intogen_gene) %>%
      summarise(features = paste0(feature_name, collapse = ";"))
  ) %>%
  arrange(cancer, desc(functional_posterior)) %>%
  write_tsv("data/all_conservative_genes.tsv")

comp_df <- read_tsv("data/2020-02-02_IntOGen-Drivers-20200213/Compendium_Cancer_Genes.tsv")
cohort_df <- read_tsv("data/2020-02-02_IntOGen-Cohorts-20200213/cohorts.tsv")

cohort_notcga <- filter(cohort_df, SOURCE != "TCGA") %>% rename(cohort = COHORT)

sub_comp_df <- select(comp_df, Gene = SYMBOL, cohort = COHORT, cancer = CANCER_TYPE, qval_cons = QVALUE_COMBINATION) %>%
  semi_join(cohort_notcga)

tcomp_df <- select(comp_df, Gene = SYMBOL, cohort = COHORT, cancer = CANCER_TYPE, qval_cons = QVALUE_COMBINATION) %>%
  semi_join(cohort_notcga) %>%
  distinct(cancer, Gene) %>%
  mutate(value = 1L)
