# ══════════════════════════════════════════════════════════════════════════════
# BACHELOR THESIS — FULL STATISTICAL ANALYSIS (revised)
# ══════════════════════════════════════════════════════════════════════════════

library(ez)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(moments)
library(psych)
library(janitor)

# ── PATHS ─────────────────────────────────────────────────────────────────────
bids_root <- "C:/Users/alisa/OneDrive/Desktop/HRV_Biotrace/bids_data"
hrv_path  <- file.path(bids_root, "derivatives", "hrv-analysis")
eda_path  <- file.path(bids_root, "derivatives", "eda-analysis")
bf_path   <- file.path(bids_root, "derivatives", "bf-analysis")
fig_dir   <- file.path(bids_root, "derivatives", "group_figures_revised")
dir.create(fig_dir, recursive=TRUE, showWarnings=FALSE)

# ── EXCLUSION SETS ────────────────────────────────────────────────────────────
excl_ecg  <- c("sub-021","sub-010","sub-002","sub-007","sub-020","sub-023")
excl_eda  <- c("sub-021","sub-003","sub-005","sub-011")
excl_temp <- c("sub-021","sub-008","sub-012","sub-022","sub-025","sub-001")

S2_TO_MS2    <- 1e6
LN_SHIFT_MS2 <- log(1e6)

colors       <- c("control"="#2166AC", "experimental"="#D6604D")
phase_labels <- c("pre"="PRE","early"="EARLY","late"="LATE","post"="POST")

MyTheme <- theme(
  axis.title.x     = element_text(size=13),
  axis.text.x      = element_text(size=11),
  axis.title.y     = element_text(size=13, vjust=3),
  axis.text.y      = element_text(size=10),
  plot.title       = element_text(size=13, face="bold"),
  plot.subtitle    = element_text(size=10, color="grey40"),
  panel.background = element_rect(fill="white", colour="grey80"),
  panel.grid.major = element_line(linewidth=0.4, colour="grey92"),
  panel.grid.minor = element_blank(),
  plot.margin      = margin(0.5,0.5,0.5,0.8,"cm"),
  legend.position  = "top"
)

make_long_complete <- function(data, metric_col, excl_set) {
  df <- data %>%
    dplyr::filter(!VP %in% excl_set) %>%
    dplyr::select(VP, ses_id, condition, phase, dplyr::all_of(metric_col)) %>%
    dplyr::filter(!is.na(.data[[metric_col]])) %>%
    dplyr::rename(value = dplyr::all_of(metric_col))
  complete_VPs <- df %>%
    dplyr::group_by(VP) %>%
    dplyr::summarise(n=dplyr::n(), .groups="drop") %>%
    dplyr::filter(n==8) %>% dplyr::pull(VP)
  cat(sprintf("  %s: n=%d complete\n", metric_col, length(complete_VPs)))
  df %>% dplyr::filter(VP %in% complete_VPs) %>%
    dplyr::mutate(VP=droplevels(as.factor(VP)),
                  condition=as.factor(condition),
                  phase=factor(phase, levels=c("pre","early","late","post")))
}

check_normality <- function(data, col="value", label="") {
  if (nchar(label)) cat(sprintf("\nNormalverteilung %s:\n", label))
  for (ph in c("pre","early","late","post")) {
    vals <- data %>% dplyr::filter(phase==ph) %>% dplyr::pull(all_of(col))
    sw   <- shapiro.test(vals)
    cat(sprintf("  %-6s: W=%.4f, p=%.4f%s\n",
                toupper(ph), sw$statistic, sw$p.value,
                ifelse(sw$p.value < 0.05, " (!)", "")))
  }
}

check_diff_normality <- function(wide_df, label="") {
  if (nchar(label)) cat(sprintf("\nNormalverteilung Differenzvariablen — %s:\n", label))
  diff_vars <- names(wide_df)[grepl("^diff_", names(wide_df))]
  for (dv in diff_vars) {
    vals <- wide_df[[dv]][!is.na(wide_df[[dv]])]
    sw   <- shapiro.test(vals)
    cat(sprintf("  %-20s: p=%.4f%s\n", dv, sw$p.value,
                ifelse(sw$p.value < 0.05, " (!)", "")))
  }
  for (dv in diff_vars) {
    vals <- wide_df[[dv]][!is.na(wide_df[[dv]])]
    sw   <- shapiro.test(vals)
    if (sw$p.value < 0.05) {
      sk <- skewness(vals); ku <- kurtosis(vals)
      cat(sprintf("  %-20s skew=%.3f (z=%.2f) kurt=%.3f (z=%.2f)\n",
                  dv, sk, sk/0.455, ku, ku/0.91))
    }
  }
}

run_anova_block <- function(df_long, dv_label, h1_direction="greater") {
  cat(sprintf("\n%s\n%s\n", strrep("─",65), dv_label))
  anova_ez <- ezANOVA(data=df_long, wid=VP, dv=value,
                      within=.(phase, condition),
                      detailed=TRUE, return_aov=TRUE)
  cat("\nMauchly:\n");     print(anova_ez$Mauchly)
  cat("\nSpherizität:\n"); print(anova_ez$Sphericity)
  cat("\nANOVA (ez):\n");  print(anova_ez$ANOVA)
  anova_rs <- anova_test(data=df_long, dv=value, wid=VP,
                         within=c(phase, condition))
  cat("\nrstatix (HF-Korrektur):\n")
  print(get_anova_table(anova_rs, correction="HF"))
  cat("\nPost-hoc Paarvergleiche Phase (Bonferroni, zweiseitig):\n")
  bon <- df_long %>%
    rstatix::pairwise_t_test(value ~ phase, pool.sd=FALSE, paired=TRUE,
                             p.adjust.method="bonferroni", alternative="two.sided")
  print(bon)
  pairwise.t.test(df_long$value, df_long$phase, p.adj="bonferroni", paired=TRUE)
  dir_str <- ifelse(h1_direction=="greater","exp > ctrl",
                    ifelse(h1_direction=="less","exp < ctrl","zweiseitig"))
  cat(sprintf("\nHaupthypothese (%s) pro Phase:\n", dir_str))
  for (ph in c("pre","early","late","post")) {
    ctrl <- df_long %>% dplyr::filter(phase==ph, condition=="control") %>% dplyr::pull(value)
    exp  <- df_long %>% dplyr::filter(phase==ph, condition=="experimental") %>% dplyr::pull(value)
    tt   <- t.test(exp, ctrl, paired=TRUE, alternative=h1_direction)
    cat(sprintf("  %-6s: t(%.0f)=%.3f, p=%.4f%s\n",
                toupper(ph), tt$parameter, tt$statistic, tt$p.value,
                ifelse(tt$p.value < 0.05, " *", "")))
  }
  invisible(anova_ez)
}

plot_signal_split <- function(df_long, signal_label, ylab, fname_base) {
  means <- ezStats(data=df_long, wid=VP, dv=value, within=.(phase, condition))
  means$phase <- factor(means$phase, levels=c("pre","early","late","post"))
  p1 <- ggplot(means, aes(x=phase, y=Mean, color=condition, group=condition)) +
    geom_line(linewidth=1.2) +
    geom_point(size=5, shape=15) +
    geom_errorbar(aes(ymin=Mean-0.5*FLSD, ymax=Mean+0.5*FLSD), width=0.1, linewidth=1) +
    scale_color_manual(values=colors, labels=c("Control","Experimental")) +
    scale_x_discrete(labels=phase_labels) +
    labs(x="Phase", y=ylab,
         title=paste(signal_label, "— Group Means ± FLSD"),
         caption="Error bars: FLSD", color="Condition") +
    MyTheme
  print(p1)
  ggsave(file.path(fig_dir, paste0(fname_base,"_means.png")), p1,
         width=8, height=6, dpi=300, bg="white")
  p2 <- ggplot() +
    geom_line(data=df_long %>% dplyr::mutate(
      phase=factor(phase,levels=c("pre","early","late","post"))),
      aes(x=phase, y=value, group=interaction(VP, condition), color=condition),
      alpha=0.25, linewidth=0.6) +
    scale_color_manual(values=colors, labels=c("Control","Experimental")) +
    scale_x_discrete(labels=phase_labels) +
    labs(x="Phase", y=ylab,
         title=paste(signal_label, "— Individual Trajectories"),
         color="Condition") +
    MyTheme
  print(p2)
  ggsave(file.path(fig_dir, paste0(fname_base,"_indiv.png")), p2,
         width=8, height=6, dpi=300, bg="white")
  cat(sprintf("Saved: %s_means.png + %s_indiv.png\n", fname_base, fname_base))
}


# ══════════════════════════════════════════════════════════════════════════════
# 1.  HRV
# ══════════════════════════════════════════════════════════════════════════════
cat("\n", strrep("=",65), "\n")
cat("HRV ANALYSIS (ECG group, n=15)\n")

hrv_time <- read_tsv(file.path(hrv_path, "task-BBSIG_per_phase_hrv-time.tsv"),
                     show_col_types=FALSE) %>%
  dplyr::rename(VP=subj_id) %>%
  dplyr::mutate(VP=as.factor(VP), condition=as.factor(condition),
                phase=factor(phase,levels=c("pre","early","late","post")))

hrv_freq <- read_tsv(file.path(hrv_path, "task-BBSIG_per_phase_hrv-freq.tsv"),
                     show_col_types=FALSE) %>%
  dplyr::rename(VP=subj_id) %>%
  dplyr::mutate(VP=as.factor(VP), condition=as.factor(condition),
                phase=factor(phase,levels=c("pre","early","late","post")))

hrv_all <- hrv_time %>%
  dplyr::left_join(hrv_freq, by=c("VP","ses_id","condition","phase"))

hrv_all <- hrv_all %>% dplyr::mutate(
  RMSSD_log  = log(RMSSD + 1),
  HFn_final  = HFn,
  HF_ms2     = HF * S2_TO_MS2,
  HF_ms2_log = log(HF_ms2 + 1),
  LnHF_ms2   = LnHF + LN_SHIFT_MS2
)

cat("\nSanity check — HF units:\n")
cat(sprintf("  HF raw:       %.6f – %.6f s²\n",
            min(hrv_all$HF,na.rm=TRUE), max(hrv_all$HF,na.rm=TRUE)))
cat(sprintf("  HF converted: %.1f – %.1f ms²\n",
            min(hrv_all$HF_ms2,na.rm=TRUE), max(hrv_all$HF_ms2,na.rm=TRUE)))
cat(sprintf("  RMSSD:        %.1f – %.1f ms\n",
            min(hrv_all$RMSSD,na.rm=TRUE), max(hrv_all$RMSSD,na.rm=TRUE)))

hrv_ecg <- hrv_all %>% dplyr::filter(!VP %in% excl_ecg) %>%
  dplyr::mutate(VP=droplevels(VP))

# ── RMSSD ─────────────────────────────────────────────────────────────────────
cat("\n── RMSSD ──────────────────────────────────────────────────\n")

# >>> ADDED: print range and summary to check for implausible values
cat("\nRMSSD raw summary (ms):\n")
print(summary(hrv_ecg$RMSSD))

ggdensity(hrv_ecg$RMSSD,
          title="RMSSD raw",
          xlab="RMSSD (ms)")
ggqqplot(hrv_ecg$RMSSD, title="RMSSD raw Q-Q")

check_normality(hrv_ecg %>% dplyr::rename(value=RMSSD), label="RMSSD (roh)")

ggdensity(hrv_ecg$RMSSD_log,
          title="RMSSD log(x+1) — expect bell shape",
          xlab="log RMSSD")
ggqqplot(hrv_ecg$RMSSD_log, title="RMSSD log Q-Q")

check_normality(hrv_ecg %>% dplyr::rename(value=RMSSD_log), label="RMSSD (log)")

rmssd_wide <- hrv_ecg %>%
  dplyr::select(VP,ses_id,condition,phase,RMSSD_log) %>%
  tidyr::pivot_wider(names_from=phase, values_from=RMSSD_log, names_prefix="RMSSD_") %>%
  dplyr::mutate(diff_pre_early=RMSSD_pre-RMSSD_early, diff_pre_late=RMSSD_pre-RMSSD_late,
                diff_pre_post=RMSSD_pre-RMSSD_post, diff_early_late=RMSSD_early-RMSSD_late,
                diff_early_post=RMSSD_early-RMSSD_post, diff_late_post=RMSSD_late-RMSSD_post)
check_diff_normality(rmssd_wide, "RMSSD")
saveRDS(rmssd_wide, file.path(hrv_path, "data_rmssd_wide_n15.Rds"))

df_RMSSD <- make_long_complete(hrv_all, "RMSSD_log", excl_ecg)

cat("\nDeskriptivstatistik RMSSD:\n")
ezStats(data=df_RMSSD, wid=VP, dv=value, within=.(phase, condition))

run_anova_block(df_RMSSD, "RMSSD  [log(ms+1)]  ECG group  n=15",
                h1_direction="greater")
saveRDS(df_RMSSD, file.path(hrv_path, "data_df_RMSSD_n15.Rds"))
plot_signal_split(df_RMSSD, "RMSSD (n=15)", "log RMSSD (ms)", "fig_RMSSD")


# ── HFn ───────────────────────────────────────────────────────────────────────
cat("\n── HFn ────────────────────────────────────────────────────\n")

# >>> ADDED: summary to check boundary values (should stay well between 0 and 1)
cat("\nHFn raw summary:\n")
print(summary(hrv_ecg$HFn))
cat(sprintf("  Values < 0.05: %d\n", sum(hrv_ecg$HFn < 0.05, na.rm=TRUE)))
cat(sprintf("  Values > 0.95: %d\n", sum(hrv_ecg$HFn > 0.95, na.rm=TRUE)))

ggdensity(hrv_ecg$HFn_final,
          title="HFn raw",
          xlab="HFn")
ggqqplot(hrv_ecg$HFn_final, title="HFn raw Q-Q")

check_normality(hrv_ecg %>% dplyr::rename(value=HFn),       label="HFn (roh)")
check_normality(hrv_ecg %>% dplyr::mutate(value=sqrt(HFn)), label="HFn (sqrt)")
# sqrt introduces violation in EARLY — retain raw HFn

hrv_ecg <- hrv_ecg %>% dplyr::mutate(HFn_final=HFn)

hfn_wide <- hrv_ecg %>%
  dplyr::select(VP,ses_id,condition,phase,HFn_final) %>%
  tidyr::pivot_wider(names_from=phase, values_from=HFn_final, names_prefix="HFn_") %>%
  dplyr::mutate(diff_pre_early=HFn_pre-HFn_early, diff_pre_late=HFn_pre-HFn_late,
                diff_pre_post=HFn_pre-HFn_post, diff_early_late=HFn_early-HFn_late,
                diff_early_post=HFn_early-HFn_post, diff_late_post=HFn_late-HFn_post)
check_diff_normality(hfn_wide, "HFn")
saveRDS(hfn_wide, file.path(hrv_path, "data_hfn_wide_n15.Rds"))

df_HFn <- make_long_complete(hrv_all, "HFn_final", excl_ecg)
cat("\nDeskriptivstatistik HFn:\n")
ezStats(data=df_HFn, wid=VP, dv=value, within=.(phase, condition))
run_anova_block(df_HFn, "HFn  [raw 0–1]  ECG group  n=15", h1_direction="greater")
saveRDS(df_HFn, file.path(hrv_path, "data_df_HFn_n15.Rds"))
plot_signal_split(df_HFn, "HFn (n=15)", "Normalised HF Power (HFn)", "fig_HFn")


# ── HF absolute ───────────────────────────────────────────────────────────────
cat("\n── HF absolute (ms²) ──────────────────────────────────────\n")

check_normality(hrv_ecg %>% dplyr::mutate(value=HF*S2_TO_MS2,
                                          phase=factor(phase,levels=c("pre","early","late","post"))),
                label="HF_ms2 (roh)")
check_normality(hrv_ecg %>% dplyr::mutate(value=log(HF*S2_TO_MS2+1),
                                          phase=factor(phase,levels=c("pre","early","late","post"))),
                label="HF_ms2 (log+1)")

df_HF <- make_long_complete(hrv_all, "HF_ms2_log", excl_ecg)
cat("\nDeskriptivstatistik HF absolute:\n")
ezStats(data=df_HF, wid=VP, dv=value, within=.(phase, condition))
run_anova_block(df_HF, "HF absolute  [log(ms²+1)]  ECG group  n=15", h1_direction="greater")
saveRDS(df_HF, file.path(hrv_path, "data_df_HF_abs_n15.Rds"))
plot_signal_split(df_HF, "HF absolute (n=15)", "log HF Power (ms²)", "fig_HF_abs")


# ── LnHF ──────────────────────────────────────────────────────────────────────
cat("\n── LnHF (ln ms²) ──────────────────────────────────────────\n")

check_normality(hrv_ecg %>% dplyr::mutate(value=LnHF+LN_SHIFT_MS2,
                                          phase=factor(phase,levels=c("pre","early","late","post"))),
                label="LnHF_ms2")

df_LnHF <- make_long_complete(hrv_all, "LnHF_ms2", excl_ecg)
run_anova_block(df_LnHF, "LnHF  [ln(ms²)]  ECG group  n=15", h1_direction="greater")
saveRDS(df_LnHF, file.path(hrv_path, "data_df_LnHF_n15.Rds"))


# ── SENSITIVITY: HFn n=15 vs n=20 ────────────────────────────────────────────
cat("\n", strrep("=",65), "\n")
cat("SENSITIVITY ANALYSIS — HFn n=15 vs n=20\n")

hrv_n20 <- hrv_all %>%
  dplyr::filter(!VP %in% c("sub-021")) %>%
  dplyr::mutate(VP=droplevels(as.factor(VP)),
                condition=as.factor(condition),
                phase=factor(phase,levels=c("pre","early","late","post")))

df_HFn_n20 <- make_long_complete(hrv_n20, "HFn_final", c())
run_anova_block(df_HFn_n20, "HFn n=20 sensitivity", h1_direction="greater")


# ══════════════════════════════════════════════════════════════════════════════
# 2.  HEART RATE
# ══════════════════════════════════════════════════════════════════════════════
cat("\n", strrep("=",65), "\n")
cat("HEART RATE ANALYSIS (ECG group, n=15)\n")

data_HR <- read_tsv(file.path(hrv_path, "task-BBSIG_per_phase_hr_nk.tsv"),
                    show_col_types=FALSE) %>%
  dplyr::rename(VP=subj_id) %>%
  dplyr::mutate(VP=as.factor(VP), condition=as.factor(condition),
                phase=factor(phase,levels=c("pre","early","late","post")),
                HR_log=log(HR_mean+2))

# >>> ADDED: identify who is causing the second bump around 90 bpm
# The density showed a main peak and a second bump around 90 — likely 1-2
# participants with consistently high resting HR. Check here:
cat("\nHR raw summary (bpm):\n")
print(summary(data_HR$HR_mean))

cat("\nMean HR per participant (check for outliers above ~90 bpm):\n")
data_HR %>%
  dplyr::filter(!VP %in% excl_ecg) %>%
  dplyr::group_by(VP) %>%
  dplyr::summarise(mean_HR = round(mean(HR_mean, na.rm=TRUE), 1),
                   max_HR  = round(max(HR_mean,  na.rm=TRUE), 1),
                   .groups="drop") %>%
  dplyr::arrange(desc(mean_HR)) %>%
  print()

ggdensity(data_HR$HR_mean,
          title="HR raw — second bump ~90 bpm may indicate high-HR participants",
          xlab="HR (bpm)")
ggqqplot(data_HR$HR_mean, title="HR raw Q-Q")

check_normality(data_HR %>% dplyr::rename(value=HR_mean), label="HR (roh)")

ggdensity(data_HR$HR_log,
          title="HR log(x+2) — second bump should be reduced",
          xlab="log HR")
ggqqplot(data_HR$HR_log, title="HR log Q-Q")

check_normality(data_HR %>% dplyr::rename(value=HR_log), label="HR (log)")

# >>> ADDED: per-phase density for log HR to check phase-level distributions
data_HR %>%
  dplyr::filter(!VP %in% excl_ecg) %>%
  ggplot(aes(x=HR_log)) +
  geom_density(fill="steelblue", alpha=0.4) +
  facet_wrap(~phase, nrow=1) +
  labs(title="HR log(x+2) — density per phase",
       x="log HR", y="density") +
  theme_minimal()

data_HR_log_wide <- data_HR %>%
  dplyr::filter(!VP %in% excl_ecg) %>%
  dplyr::select(VP,ses_id,condition,phase,HR_log) %>%
  tidyr::pivot_wider(names_from=phase, values_from=HR_log, names_prefix="HR_") %>%
  dplyr::mutate(diff_pre_early=HR_pre-HR_early, diff_pre_late=HR_pre-HR_late,
                diff_pre_post=HR_pre-HR_post, diff_early_late=HR_early-HR_late,
                diff_early_post=HR_early-HR_post, diff_late_post=HR_late-HR_post)
check_diff_normality(data_HR_log_wide, "HR")
saveRDS(data_HR_log_wide, file.path(hrv_path, "data_HR_log_wide_n15.Rds"))

df_HR <- make_long_complete(data_HR, "HR_log", excl_ecg)
cat("\nDeskriptivstatistik HR (log):\n")
ezStats(data=df_HR, wid=VP, dv=value, within=.(phase, condition))

cat("\nRohe Mittelwerte HR (bpm):\n")
data_HR %>%
  dplyr::filter(VP %in% unique(df_HR$VP)) %>%
  dplyr::group_by(phase, condition) %>%
  dplyr::summarise(M=round(mean(HR_mean,na.rm=TRUE),2),
                   SD=round(sd(HR_mean,na.rm=TRUE),2), .groups="drop") %>% print()

run_anova_block(df_HR, "Heart Rate  [log(bpm+2)]  ECG group  n=15",
                h1_direction="less")
saveRDS(df_HR, file.path(hrv_path, "data_df_HR_n15.Rds"))

data_HR_raw_complete <- data_HR %>%
  dplyr::filter(VP %in% unique(df_HR$VP)) %>%
  dplyr::select(VP,ses_id,condition,phase,HR_mean) %>%
  dplyr::rename(value=HR_mean) %>%
  dplyr::mutate(VP=droplevels(as.factor(VP)), condition=as.factor(condition),
                phase=factor(phase,levels=c("pre","early","late","post")))
plot_signal_split(data_HR_raw_complete, "Heart Rate (n=15)", "Heart Rate (bpm)", "fig_HR")
cat("Heart Rate analysis complete\n")


# ══════════════════════════════════════════════════════════════════════════════
# 3.  BREATHING FREQUENCY
# ══════════════════════════════════════════════════════════════════════════════
cat("\n", strrep("=",65), "\n")
cat("BREATHING FREQUENCY ANALYSIS (ECG group, n=15)\n")

data_BF <- read_tsv(file.path(bf_path, "task-BBSIG_per_phase_bf.tsv"),
                    show_col_types=FALSE) %>%
  dplyr::rename(VP=subj_id) %>%
  dplyr::mutate(VP=as.factor(VP), condition=as.factor(condition),
                phase=factor(phase,levels=c("pre","early","late","post")))

# >>> ADDED: range check — confirm no values outside exclusion bounds remain
cat("\nBF raw range check (should be 9–24 br/min after exclusions):\n")
cat(sprintf("  Min: %.2f  Max: %.2f\n",
            min(data_BF$BF_bpm, na.rm=TRUE),
            max(data_BF$BF_bpm, na.rm=TRUE)))
cat(sprintf("  Values below 9 br/min:  %d\n", sum(data_BF$BF_bpm < 9,  na.rm=TRUE)))
cat(sprintf("  Values above 24 br/min: %d\n", sum(data_BF$BF_bpm > 24, na.rm=TRUE)))


# >>> CHANGED: density on long-format data so the full distribution is visible
# (previously this was on wide format which had no single plottable column)
ggdensity(data_BF$BF_bpm,
          title="BF raw",
          xlab="BF (br/min)")
ggqqplot(data_BF$BF_bpm, title="BF raw Q-Q")

# >>> ADDED: per-phase density raw
data_BF %>%
  dplyr::filter(!VP %in% excl_ecg) %>%
  ggplot(aes(x=BF_bpm)) +
  geom_density(fill="steelblue", alpha=0.4) +
  facet_wrap(~phase, nrow=1) +
  labs(title="BF raw — density per phase", x="BF (br/min)", y="density") +
  theme_minimal()

data_BF_wide <- data_BF %>%
  dplyr::filter(!VP %in% excl_ecg, !is.na(BF_bpm)) %>%
  dplyr::select(VP,ses_id,condition,phase,BF_bpm) %>%
  tidyr::pivot_wider(names_from=phase, values_from=BF_bpm, names_prefix="BF_")

cat("\nNormalverteilung BF (roh):\n")
for (ph in c("pre","early","late","post")) {
  vals <- data_BF_wide[[paste0("BF_",ph)]]
  sw <- shapiro.test(vals[!is.na(vals)])
  cat(sprintf("  %-6s: W=%.4f, p=%.4f%s\n", toupper(ph), sw$statistic, sw$p.value,
              ifelse(sw$p.value<0.05," (!)","")))}

data_BF_log <- data_BF_wide %>% dplyr::mutate(
  BF_pre=log(BF_pre+2), BF_early=log(BF_early+2),
  BF_late=log(BF_late+2), BF_post=log(BF_post+2))

cat("\nNormalverteilung BF (log):\n")
for (ph in c("pre","early","late","post")) {
  vals <- data_BF_log[[paste0("BF_",ph)]]
  sw <- shapiro.test(vals[!is.na(vals)])
  cat(sprintf("  %-6s: W=%.4f, p=%.4f%s\n", toupper(ph), sw$statistic, sw$p.value,
              ifelse(sw$p.value<0.05," (!)","")))}

# >>> ADDED: per-phase density after log transformation using long format
data_BF_log_long_check <- data_BF_log %>%
  tidyr::pivot_longer(cols=c(BF_pre,BF_early,BF_late,BF_post),
                      names_to="phase", values_to="BF_log") %>%
  dplyr::mutate(phase=factor(gsub("BF_","",phase),
                             levels=c("pre","early","late","post")))

ggdensity(data_BF_log_long_check$BF_log,
          title="BF log(x+2) — full distribution, expect bell shape",
          xlab="log BF")
ggqqplot(data_BF_log_long_check$BF_log, title="BF log Q-Q overall")

data_BF_log_long_check %>%
  ggplot(aes(x=BF_log)) +
  geom_density(fill="steelblue", alpha=0.4) +
  facet_wrap(~phase, nrow=1) +
  labs(title="BF log(x+2) — density per phase", x="log BF", y="density") +
  theme_minimal()

data_BF_log <- data_BF_log %>% dplyr::mutate(
  diff_pre_early=BF_pre-BF_early, diff_pre_late=BF_pre-BF_late,
  diff_pre_post=BF_pre-BF_post, diff_early_late=BF_early-BF_late,
  diff_early_post=BF_early-BF_post, diff_late_post=BF_late-BF_post)
check_diff_normality(data_BF_log, "BF")
saveRDS(data_BF_log, file.path(bf_path, "data_BF_log_n15.Rds"))

# >>> FIXED: removed typo (data_BF_lodata_BF_lodata_BF_log_long)
data_BF_log_long <- data_BF_log %>%
  tidyr::pivot_longer(cols=c(BF_pre,BF_early,BF_late,BF_post),
                      names_to="phase", values_to="bf_values") %>%
  dplyr::mutate(phase=factor(gsub("BF_","",phase),levels=c("pre","early","late","post")))

complete_VPs_BF <- data_BF_log_long %>%
  dplyr::filter(!is.na(bf_values)) %>%
  dplyr::group_by(VP) %>%
  dplyr::summarise(n=dplyr::n(), .groups="drop") %>%
  dplyr::filter(n==8) %>% dplyr::pull(VP)
cat(sprintf("\nComplete BF participants: %d\n", length(complete_VPs_BF)))

df_BF <- data_BF_log_long %>%
  dplyr::filter(VP %in% complete_VPs_BF) %>%
  dplyr::mutate(VP=droplevels(VP)) %>%
  dplyr::rename(value=bf_values)

cat("\nMittelwerte BF (log):\n")
ezStats(data=df_BF, wid=VP, dv=value, within=.(phase, condition))

cat("\nRohe Mittelwerte BF (br/min):\n")
data_BF %>%
  dplyr::filter(VP %in% complete_VPs_BF, !is.na(BF_bpm)) %>%
  dplyr::group_by(phase, condition) %>%
  dplyr::summarise(M=round(mean(BF_bpm,na.rm=TRUE),2),
                   SD=round(sd(BF_bpm,na.rm=TRUE),2), .groups="drop") %>% print()

run_anova_block(df_BF, "Breathing Frequency  [log(br/min+2)]  ECG group  n=15",
                h1_direction="less")
saveRDS(df_BF, file.path(bf_path, "data_df_BF_n15.Rds"))

data_BF_raw_long <- data_BF %>%
  dplyr::filter(VP %in% complete_VPs_BF, !is.na(BF_bpm)) %>%
  dplyr::select(VP,ses_id,condition,phase,BF_bpm) %>%
  dplyr::rename(value=BF_bpm) %>%
  dplyr::mutate(VP=droplevels(as.factor(VP)), condition=as.factor(condition),
                phase=factor(phase,levels=c("pre","early","late","post")))
plot_signal_split(data_BF_raw_long, "Breathing Frequency (n=15)",
                  "Breathing Frequency (br/min)", "fig_BF")
cat("BF analysis complete\n")


# ══════════════════════════════════════════════════════════════════════════════
# 4.  EDA & TEMPERATURE
# ══════════════════════════════════════════════════════════════════════════════
cat("\n", strrep("=",65), "\n")
cat("EDA & TEMPERATURE ANALYSIS\n")

data_raw <- read_tsv(
  file.path(eda_path, "task-BBSIG_eda-temp-bf_per_phase_final.tsv"),
  show_col_types=FALSE) %>%
  dplyr::rename(VP=subj_id) %>%
  dplyr::mutate(VP=as.factor(VP), condition=as.factor(condition),
                phase=factor(phase,levels=c("pre","early","late","post")))

janitor::get_dupes(data_raw)

# ── EDA ───────────────────────────────────────────────────────────────────────
data_EDA_wide <- data_raw %>%
  dplyr::filter(!VP %in% excl_eda, !EDA_excluded) %>%
  dplyr::select(VP,ses_id,condition,phase,EDA_mean) %>%
  tidyr::pivot_wider(names_from=phase, values_from=EDA_mean, names_prefix="EDA_")

# >>> ADDED: full distribution density before transformation using long format
data_EDA_raw_long <- data_raw %>%
  dplyr::filter(!VP %in% excl_eda, !EDA_excluded) %>%
  dplyr::select(VP,ses_id,condition,phase,EDA_mean)

cat("\nEDA raw summary (µS):\n")
print(summary(data_EDA_raw_long$EDA_mean))

ggdensity(data_EDA_raw_long$EDA_mean,
          title="EDA raw",
          xlab="EDA (µS)")
ggqqplot(data_EDA_raw_long$EDA_mean, title="EDA raw Q-Q")

# >>> ADDED: per-phase density raw
data_EDA_raw_long %>%
  ggplot(aes(x=EDA_mean)) +
  geom_density(fill="steelblue", alpha=0.4) +
  facet_wrap(~phase, nrow=1) +
  labs(title="EDA raw — density per phase", x="EDA (µS)", y="density") +
  theme_minimal()

cat("\nNormalverteilung EDA (roh):\n")
for (ph in c("pre","early","late","post")) {
  vals <- data_EDA_wide[[paste0("EDA_",ph)]]
  sw <- shapiro.test(vals[!is.na(vals)])
  cat(sprintf("  %-6s: W=%.4f, p=%.4f%s\n", toupper(ph), sw$statistic, sw$p.value,
              ifelse(sw$p.value<0.05," (!)","")))}

data_EDA_log <- data_EDA_wide %>% dplyr::mutate(
  EDA_pre=log(EDA_pre+2), EDA_early=log(EDA_early+2),
  EDA_late=log(EDA_late+2), EDA_post=log(EDA_post+2))

# >>> ADDED: full distribution density after transformation using long format
data_EDA_log_long_check <- data_EDA_log %>%
  tidyr::pivot_longer(cols=c(EDA_pre,EDA_early,EDA_late,EDA_post),
                      names_to="phase", values_to="EDA_log") %>%
  dplyr::mutate(phase=factor(gsub("EDA_","",phase),
                             levels=c("pre","early","late","post")))

ggdensity(data_EDA_log_long_check$EDA_log,
          title="EDA log(x+2) — full distribution, expect bell shape",
          xlab="log EDA")
ggqqplot(data_EDA_log_long_check$EDA_log, title="EDA log Q-Q overall")

data_EDA_log_long_check %>%
  ggplot(aes(x=EDA_log)) +
  geom_density(fill="steelblue", alpha=0.4) +
  facet_wrap(~phase, nrow=1) +
  labs(title="EDA log(x+2) — density per phase", x="log EDA", y="density") +
  theme_minimal()

cat("\nNormalverteilung EDA (log):\n")
for (ph in c("pre","early","late","post")) {
  vals <- data_EDA_log[[paste0("EDA_",ph)]]
  sw <- shapiro.test(vals[!is.na(vals)])
  cat(sprintf("  %-6s: W=%.4f, p=%.4f%s\n", toupper(ph), sw$statistic, sw$p.value,
              ifelse(sw$p.value<0.05," (!)","")))}

data_EDA_log <- data_EDA_log %>% dplyr::mutate(
  diff_pre_early=EDA_pre-EDA_early, diff_pre_late=EDA_pre-EDA_late,
  diff_pre_post=EDA_pre-EDA_post, diff_early_late=EDA_early-EDA_late,
  diff_early_post=EDA_early-EDA_post, diff_late_post=EDA_late-EDA_post)
check_diff_normality(data_EDA_log, "EDA")
saveRDS(data_EDA_log, file.path(eda_path, "data_EDA_log_n17.Rds"))

data_EDA_log_long <- data_EDA_log %>%
  tidyr::pivot_longer(cols=c(EDA_pre,EDA_early,EDA_late,EDA_post),
                      names_to="phase", values_to="eda_values") %>%
  dplyr::mutate(phase=factor(gsub("EDA_","",phase),levels=c("pre","early","late","post")))

complete_VPs_EDA <- data_EDA_log_long %>%
  dplyr::filter(!is.na(eda_values)) %>%
  dplyr::group_by(VP) %>%
  dplyr::summarise(n=dplyr::n(), .groups="drop") %>%
  dplyr::filter(n==8) %>% dplyr::pull(VP)
cat(sprintf("\nComplete EDA participants: %d\n", length(complete_VPs_EDA)))

data_EDA_log_long1 <- data_EDA_log_long %>%
  dplyr::filter(VP %in% complete_VPs_EDA) %>%
  dplyr::mutate(VP=droplevels(VP))

df_EDA <- data_EDA_log_long1 %>% dplyr::rename(value=eda_values)

cat("\nMittelwerte EDA (log):\n")
ezStats(data=data_EDA_log_long1, wid=VP, dv=eda_values, within=.(phase, condition))

anova_EDA <- ezANOVA(data=data_EDA_log_long1, wid=VP, dv=eda_values,
                     within=.(phase, condition), detailed=TRUE, return_aov=TRUE)
cat("\nMauchly EDA:\n");     print(anova_EDA$Mauchly)
cat("\nSpherizität EDA:\n"); print(anova_EDA$Sphericity)
cat("\nANOVA EDA:\n");       print(anova_EDA$ANOVA)

rep_anova_EDA <- anova_test(data=data_EDA_log_long1, dv=eda_values,
                            wid=VP, within=c(phase, condition))
cat("\nrstatix EDA (HF-Korrektur):\n")
print(get_anova_table(rep_anova_EDA, correction="HF"))

pairwise.t.test(data_EDA_log_long1$eda_values, data_EDA_log_long1$phase,
                p.adj="bonferroni", paired=TRUE)
bon_EDA <- data_EDA_log_long1 %>%
  rstatix::pairwise_t_test(eda_values ~ phase, pool.sd=FALSE, paired=TRUE,
                           p.adjust.method="bonferroni", alternative="two.sided")
print(bon_EDA)

cat("\nHaupthypothese EDA — exp < ctrl pro Phase:\n")
for (ph in c("pre","early","late","post")) {
  ctrl <- data_EDA_log_long1 %>% dplyr::filter(phase==ph, condition=="control") %>%
    dplyr::pull(eda_values)
  exp  <- data_EDA_log_long1 %>% dplyr::filter(phase==ph, condition=="experimental") %>%
    dplyr::pull(eda_values)
  tt <- t.test(exp, ctrl, paired=TRUE, alternative="less")
  cat(sprintf("  %-6s: t(%.0f)=%.3f, p=%.4f%s\n",
              toupper(ph), tt$parameter, tt$statistic, tt$p.value,
              ifelse(tt$p.value<0.05," *","")))}

saveRDS(data_EDA_log_long1, file.path(eda_path, "data_df_EDA_n17.Rds"))

data_EDA_long <- data_raw %>%
  dplyr::filter(!VP %in% excl_eda, !EDA_excluded) %>%
  dplyr::select(VP,ses_id,condition,phase,EDA_mean) %>%
  dplyr::rename(value=EDA_mean) %>%
  dplyr::filter(VP %in% complete_VPs_EDA) %>%
  dplyr::mutate(VP=droplevels(as.factor(VP)), condition=as.factor(condition),
                phase=factor(phase,levels=c("pre","early","late","post")))
plot_signal_split(data_EDA_long, "EDA (n=17)", "EDA (µS)", "fig_EDA")


# ── TEMPERATURE ───────────────────────────────────────────────────────────────
cat("\n── TEMPERATURE ─────────────────────────────────────────────\n")

data_TEMP_wide <- data_raw %>%
  dplyr::filter(!VP %in% excl_temp, !temp_excluded) %>%
  dplyr::select(VP,ses_id,condition,phase,temp_mean) %>%
  tidyr::pivot_wider(names_from=phase, values_from=temp_mean, names_prefix="TEMP_")

# >>> ADDED: full raw distribution and per-phase summary
data_TEMP_raw_long <- data_raw %>%
  dplyr::filter(!VP %in% excl_temp, !temp_excluded) %>%
  dplyr::select(VP,ses_id,condition,phase,temp_mean)

cat("\nTEMP raw summary (°C):\n")
print(summary(data_TEMP_raw_long$temp_mean))

# >>> ADDED: per-participant mean to spot residual low-value outliers

cat("\nMean temperature per participant (check for values near 28°C boundary):\n")
data_TEMP_raw_long %>%
  dplyr::group_by(VP) %>%
  dplyr::summarise(mean_T = round(mean(temp_mean, na.rm=TRUE), 2),
                   min_T  = round(min(temp_mean,  na.rm=TRUE), 2),
                   .groups="drop") %>%
  dplyr::arrange(min_T) %>%
  print()

ggdensity(data_TEMP_raw_long$temp_mean,
          title="TEMP raw",
          xlab="Skin Temperature (°C)")
ggqqplot(data_TEMP_raw_long$temp_mean, title="TEMP raw Q-Q overall")

# >>> ADDED: per-phase density raw — the PRE non-normality should be visible here
data_TEMP_raw_long %>%
  ggplot(aes(x=temp_mean)) +
  geom_density(fill="steelblue", alpha=0.4) +
  facet_wrap(~phase, nrow=1) +
  labs(title="TEMP raw — density per phase (PRE non-normality visible?)",
       x="Temperature (°C)", y="density") +
  theme_minimal()

cat("\nNormalverteilung TEMP (roh):\n")
for (ph in c("pre","early","late","post")) {
  vals <- data_TEMP_wide[[paste0("TEMP_",ph)]]
  sw <- shapiro.test(vals[!is.na(vals)])
  cat(sprintf("  %-6s: W=%.4f, p=%.4f%s\n", toupper(ph), sw$statistic, sw$p.value,
              ifelse(sw$p.value<0.05," (!)","")))}

data_TEMP_log <- data_TEMP_wide %>% dplyr::mutate(
  TEMP_pre=log(TEMP_pre+2), TEMP_early=log(TEMP_early+2),
  TEMP_late=log(TEMP_late+2), TEMP_post=log(TEMP_post+2))

# >>> ADDED: full distribution and per-phase density after log transformation
data_TEMP_log_long_check <- data_TEMP_log %>%
  tidyr::pivot_longer(cols=c(TEMP_pre,TEMP_early,TEMP_late,TEMP_post),
                      names_to="phase", values_to="TEMP_log") %>%
  dplyr::mutate(phase=factor(gsub("TEMP_","",phase),
                             levels=c("pre","early","late","post")))

ggdensity(data_TEMP_log_long_check$TEMP_log,
          title="TEMP log(x+2) — full distribution",
          xlab="log Temperature")
ggqqplot(data_TEMP_log_long_check$TEMP_log, title="TEMP log Q-Q overall")

# >>> ADDED: this is the critical plot — PRE Q-Q after transformation
# Points deviating bottom-left indicate low-value outliers not caught by exclusion
data_TEMP_log_long_check %>%
  ggplot(aes(x=TEMP_log)) +
  geom_density(fill="steelblue", alpha=0.4) +
  facet_wrap(~phase, nrow=1) +
  labs(title="TEMP log(x+2) — density per phase (check PRE for outlier influence)",
       x="log Temperature", y="density") +
  theme_minimal()

# >>> ADDED: Q-Q per phase separately for TEMP — most important diagnostic
data_TEMP_log_long_check %>%
  dplyr::group_by(phase) %>%
  dplyr::group_map(~ {
    ggqqplot(.x$TEMP_log,
             title=paste("TEMP log Q-Q —", .y$phase,
                         "(outlier points bottom-left = residual sensor issue)"))
  }, .keep=TRUE) %>%
  invisible()

cat("\nNormalverteilung TEMP (log):\n")
for (ph in c("pre","early","late","post")) {
  vals <- data_TEMP_log[[paste0("TEMP_",ph)]]
  sw <- shapiro.test(vals[!is.na(vals)])
  cat(sprintf("  %-6s: W=%.4f, p=%.4f%s\n", toupper(ph), sw$statistic, sw$p.value,
              ifelse(sw$p.value<0.05," (!)","")))}

data_TEMP_log <- data_TEMP_log %>% dplyr::mutate(
  diff_pre_early=TEMP_pre-TEMP_early, diff_pre_late=TEMP_pre-TEMP_late,
  diff_pre_post=TEMP_pre-TEMP_post, diff_early_late=TEMP_early-TEMP_late,
  diff_early_post=TEMP_early-TEMP_post, diff_late_post=TEMP_late-TEMP_post)
check_diff_normality(data_TEMP_log, "TEMP")
saveRDS(data_TEMP_log, file.path(eda_path, "data_TEMP_log_n15.Rds"))

data_TEMP_log_long <- data_TEMP_log %>%
  tidyr::pivot_longer(cols=c(TEMP_pre,TEMP_early,TEMP_late,TEMP_post),
                      names_to="phase", values_to="temp_values") %>%
  dplyr::mutate(phase=factor(gsub("TEMP_","",phase),levels=c("pre","early","late","post")))

complete_VPs_TEMP <- data_TEMP_log_long %>%
  dplyr::filter(!is.na(temp_values)) %>%
  dplyr::group_by(VP) %>%
  dplyr::summarise(n=dplyr::n(), .groups="drop") %>%
  dplyr::filter(n==8) %>% dplyr::pull(VP)
cat(sprintf("\nComplete TEMP participants: %d\n", length(complete_VPs_TEMP)))

data_TEMP_log_long1 <- data_TEMP_log_long %>%
  dplyr::filter(VP %in% complete_VPs_TEMP) %>%
  dplyr::mutate(VP=droplevels(VP))

cat("\nMittelwerte TEMP (log):\n")
ezStats(data=data_TEMP_log_long1, wid=VP, dv=temp_values, within=.(phase, condition))

anova_TEMP <- ezANOVA(data=data_TEMP_log_long1, wid=VP, dv=temp_values,
                      within=.(phase, condition), detailed=TRUE, return_aov=TRUE)
cat("\nMauchly TEMP:\n");     print(anova_TEMP$Mauchly)
cat("\nSpherizität TEMP:\n"); print(anova_TEMP$Sphericity)
cat("\nANOVA TEMP:\n");       print(anova_TEMP$ANOVA)

rep_anova_TEMP <- anova_test(data=data_TEMP_log_long1, dv=temp_values,
                             wid=VP, within=c(phase, condition))
cat("\nrstatix TEMP (HF-Korrektur):\n")
print(get_anova_table(rep_anova_TEMP, correction="HF"))

pairwise.t.test(data_TEMP_log_long1$temp_values, data_TEMP_log_long1$phase,
                p.adj="bonferroni", paired=TRUE)
bon_TEMP <- data_TEMP_log_long1 %>%
  rstatix::pairwise_t_test(temp_values ~ phase, pool.sd=FALSE, paired=TRUE,
                           p.adjust.method="bonferroni", alternative="two.sided")
print(bon_TEMP)

cat("\nBedingungsvergleich TEMP — zweiseitig pro Phase:\n")
for (ph in c("pre","early","late","post")) {
  ctrl <- data_TEMP_log_long1 %>% dplyr::filter(phase==ph, condition=="control") %>%
    dplyr::pull(temp_values)
  exp  <- data_TEMP_log_long1 %>% dplyr::filter(phase==ph, condition=="experimental") %>%
    dplyr::pull(temp_values)
  tt <- t.test(exp, ctrl, paired=TRUE, alternative="two.sided")
  cat(sprintf("  %-6s: t(%.0f)=%.3f, p=%.4f%s\n",
              toupper(ph), tt$parameter, tt$statistic, tt$p.value,
              ifelse(tt$p.value<0.05," *","")))}

saveRDS(data_TEMP_log_long1, file.path(eda_path, "data_df_TEMP_n15.Rds"))

data_TEMP_long <- data_raw %>%
  dplyr::filter(!VP %in% excl_temp, !temp_excluded) %>%
  dplyr::select(VP,ses_id,condition,phase,temp_mean) %>%
  dplyr::rename(value=temp_mean) %>%
  dplyr::filter(VP %in% complete_VPs_TEMP) %>%
  dplyr::mutate(VP=droplevels(as.factor(VP)), condition=as.factor(condition),
                phase=factor(phase,levels=c("pre","early","late","post")))

# Identify who is in the lower temperature cluster
cat("\nTEMP per-participant minimum and mean — identify lower cluster:\n")
data_raw %>%
  dplyr::filter(!VP %in% excl_temp, !temp_excluded) %>%
  dplyr::group_by(VP) %>%
  dplyr::summarise(
    mean_T  = round(mean(temp_mean,  na.rm=TRUE), 2),
    min_T   = round(min(temp_mean,   na.rm=TRUE), 2),
    max_T   = round(max(temp_mean,   na.rm=TRUE), 2),
    .groups = "drop"
  ) %>%
  dplyr::arrange(mean_T) %>%
  print(n=Inf)

plot_signal_split(data_TEMP_long, "Skin Temperature (n=13)",
                  "Skin Temperature (°C)", "fig_TEMP")

cat("\n", strrep("=",65), "\n")
cat("ALL ANALYSES COMPLETE\n\n")
cat("n per measure:\n")
cat("  RMSSD, HFn, HF_abs, LnHF, HR, BF  → n=15\n")
cat("  EDA                                → n=16\n")
cat("  Skin Temperature                   → n=13\n")

