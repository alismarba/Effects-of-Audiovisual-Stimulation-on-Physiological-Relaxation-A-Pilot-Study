# ==============================================================================
# Bacheloarbeit
# ==============================================================================


### Load packages ##############################################################

library(ez)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(moments)
library(psych)
library(janitor)


### Paths ######################################################################

bids_root <- "C:/Users/alisa/OneDrive/Desktop/HRV_Biotrace/bids_data"
hrv_path  <- file.path(bids_root, "derivatives", "hrv-analysis")
eda_path  <- file.path(bids_root, "derivatives", "eda-analysis")
bf_path   <- file.path(bids_root, "derivatives", "bf-analysis")
fig_dir   <- file.path(bids_root, "derivatives", "group_figures_revised")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)


### Exclusion sets #############################################################

excl_ecg  <- c("sub-021","sub-010","sub-002","sub-007","sub-020","sub-023")
excl_eda  <- c("sub-021","sub-003","sub-005","sub-011")
excl_temp <- c("sub-021","sub-001","sub-008","sub-012","sub-022","sub-025")
# sub-001 added manually: PRE temp = 26.3°C in ses-02 below 28°C threshold,
# not flagged by automated procedure; both sessions required -> excluded.


### Constants and style ########################################################

S2_TO_MS2    <- 1e6
LN_SHIFT_MS2 <- log(1e6)

colors <- c("control" = "#2166AC", "experimental" = "#D6604D")

phase_labels <- c("pre" = "PRE", "early" = "EARLY",
                  "late" = "LATE", "post" = "POST")

MyTheme <- theme(
  axis.title.x     = element_text(size = 13),
  axis.text.x      = element_text(size = 11),
  axis.title.y     = element_text(size = 13, vjust = 3),
  axis.text.y      = element_text(size = 10),
  plot.title       = element_text(size = 13, face = "bold"),
  plot.subtitle    = element_text(size = 10, color = "grey40"),
  panel.background = element_rect(fill = "white", colour = "grey80"),
  panel.grid.major = element_line(linewidth = 0.4, colour = "grey92"),
  panel.grid.minor = element_blank(),
  plot.margin      = margin(0.5, 0.5, 0.5, 0.8, "cm"),
  legend.position  = "top"
)


### Helper functions ###########################################################

make_long_complete <- function(data, metric_col, excl_set) {
  df <- data %>%
    dplyr::filter(!VP %in% excl_set) %>%
    dplyr::select(VP, ses_id, condition, phase, dplyr::all_of(metric_col)) %>%
    dplyr::filter(!is.na(.data[[metric_col]])) %>%
    dplyr::rename(value = dplyr::all_of(metric_col))
  complete_VPs <- df %>%
    dplyr::group_by(VP) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::filter(n == 8) %>%
    dplyr::pull(VP)
  df %>%
    dplyr::filter(VP %in% complete_VPs) %>%
    dplyr::mutate(
      VP        = droplevels(as.factor(VP)),
      condition = as.factor(condition),
      phase     = factor(phase, levels = c("pre","early","late","post"))
    )
}

plot_signal_split <- function(df_long, signal_label, ylab, fname_base) {
  means       <- ezStats(data = df_long, wid = VP, dv = value,
                         within = .(phase, condition))
  means$phase <- factor(means$phase, levels = c("pre","early","late","post"))

  p1 <- ggplot(means, aes(x = phase, y = Mean,
                           color = condition, group = condition)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 5, shape = 15) +
    geom_errorbar(aes(ymin = Mean - 0.5 * FLSD, ymax = Mean + 0.5 * FLSD),
                  width = 0.1, linewidth = 1) +
    scale_color_manual(values = colors, labels = c("Control","Experimental")) +
    scale_x_discrete(labels = phase_labels) +
    labs(x = "Phase", y = ylab,
         title   = paste(signal_label, "— Group Means ± FLSD"),
         caption = "Error bars: FLSD", color = "Condition") +
    MyTheme
  ggsave(file.path(fig_dir, paste0(fname_base, "_means.png")), p1,
         width = 8, height = 6, dpi = 300, bg = "white")

  p2 <- ggplot() +
    geom_line(
      data = df_long %>%
        dplyr::mutate(phase = factor(phase,
                                     levels = c("pre","early","late","post"))),
      aes(x = phase, y = value,
          group = interaction(VP, condition), color = condition),
      alpha = 0.25, linewidth = 0.6) +
    scale_color_manual(values = colors, labels = c("Control","Experimental")) +
    scale_x_discrete(labels = phase_labels) +
    labs(x = "Phase", y = ylab,
         title = paste(signal_label, "— Individual Trajectories"),
         color = "Condition") +
    MyTheme
  ggsave(file.path(fig_dir, paste0(fname_base, "_indiv.png")), p2,
         width = 8, height = 6, dpi = 300, bg = "white")
}


# ==============================================================================
# 1.  HRV ANALYSIS  (ECG group, n = 15)
# ==============================================================================

hrv_time <- read_tsv(
  file.path(hrv_path, "task-BBSIG_per_phase_hrv-time.tsv"),
  show_col_types = FALSE) %>%
  dplyr::rename(VP = subj_id) %>%
  dplyr::mutate(VP = as.factor(VP), condition = as.factor(condition),
                phase = factor(phase, levels = c("pre","early","late","post")))

hrv_freq <- read_tsv(
  file.path(hrv_path, "task-BBSIG_per_phase_hrv-freq.tsv"),
  show_col_types = FALSE) %>%
  dplyr::rename(VP = subj_id) %>%
  dplyr::mutate(VP = as.factor(VP), condition = as.factor(condition),
                phase = factor(phase, levels = c("pre","early","late","post")))

hrv_all <- hrv_time %>%
  dplyr::left_join(hrv_freq, by = c("VP","ses_id","condition","phase"))

# NOTE: HR_log is NOT added here — HR is read from a separate file in Section 2
hrv_all <- hrv_all %>% dplyr::mutate(
  RMSSD_log  = log(RMSSD + 1),
  HF_ms2     = HF * S2_TO_MS2,
  HF_ms2_log = log(HF_ms2 + 1),
  LnHF_ms2   = LnHF + LN_SHIFT_MS2
)

hrv_ecg <- hrv_all %>%
  dplyr::filter(!VP %in% excl_ecg) %>%
  dplyr::mutate(VP = droplevels(VP))


### 1a. RMSSD ##################################################################
# Non-normally distributed in raw form -> log(x + 1) applied

shapiro.test(hrv_ecg$RMSSD[hrv_ecg$phase == "pre"])    # p = 0.0027 (!)
shapiro.test(hrv_ecg$RMSSD[hrv_ecg$phase == "early"])  # p = 0.0231 (!)
shapiro.test(hrv_ecg$RMSSD[hrv_ecg$phase == "late"])   # p = 0.0050 (!)
shapiro.test(hrv_ecg$RMSSD[hrv_ecg$phase == "post"])   # p = 0.0148 (!)

shapiro.test(hrv_ecg$RMSSD_log[hrv_ecg$phase == "pre"])    # p = 0.2533
shapiro.test(hrv_ecg$RMSSD_log[hrv_ecg$phase == "early"])  # p = 0.9455
shapiro.test(hrv_ecg$RMSSD_log[hrv_ecg$phase == "late"])   # p = 0.6984
shapiro.test(hrv_ecg$RMSSD_log[hrv_ecg$phase == "post"])   # p = 0.1383

# Normality of difference variables (log-transformed)
rmssd_wide <- hrv_ecg %>%
  dplyr::select(VP, ses_id, condition, phase, RMSSD_log) %>%
  tidyr::pivot_wider(names_from = phase, values_from = RMSSD_log,
                     names_prefix = "RMSSD_") %>%
  dplyr::mutate(
    diff_pre_early  = RMSSD_pre   - RMSSD_early,
    diff_pre_late   = RMSSD_pre   - RMSSD_late,
    diff_pre_post   = RMSSD_pre   - RMSSD_post,
    diff_early_late = RMSSD_early - RMSSD_late,
    diff_early_post = RMSSD_early - RMSSD_post,
    diff_late_post  = RMSSD_late  - RMSSD_post
  )

shapiro.test(rmssd_wide$diff_pre_early)   # p = 0.0318 (!) skew = -0.753, kurt = 4.360
shapiro.test(rmssd_wide$diff_pre_late)    # p = 0.9984
shapiro.test(rmssd_wide$diff_pre_post)    # p = 0.9956
shapiro.test(rmssd_wide$diff_early_late)  # p = 0.0845
shapiro.test(rmssd_wide$diff_early_post)  # p = 0.1242
shapiro.test(rmssd_wide$diff_late_post)   # p = 0.3291
# diff_pre_early: skew = -0.753 (z = -1.65), kurt = 4.360 (z = 4.79) — retained

skewness(rmssd_wide$diff_pre_early);  skewness(rmssd_wide$diff_pre_early) / 0.455
kurtosis(rmssd_wide$diff_pre_early); kurtosis(rmssd_wide$diff_pre_early) / 0.91

saveRDS(rmssd_wide, file.path(hrv_path, "data_rmssd_wide_n15.Rds"))

df_RMSSD <- make_long_complete(hrv_all, "RMSSD_log", excl_ecg)  # n = 15

ezStats(data = df_RMSSD, wid = VP, dv = value, within = .(phase, condition))

anova_RMSSD <- ezANOVA(data = df_RMSSD, wid = VP, dv = value,
                        within = .(phase, condition),
                        detailed = TRUE, return_aov = TRUE)
anova_RMSSD$Mauchly
anova_RMSSD$Sphericity
anova_RMSSD$ANOVA

get_anova_table(anova_test(df_RMSSD, dv = value, wid = VP,
                            within = c(phase, condition)), correction = "HF")

df_RMSSD %>%
  rstatix::pairwise_t_test(value ~ phase, pool.sd = FALSE, paired = TRUE,
                           p.adjust.method = "bonferroni", alternative = "two.sided")

# Primary hypothesis (exp > ctrl) per phase
t.test(df_RMSSD$value[df_RMSSD$phase=="pre"   & df_RMSSD$condition=="experimental"],
       df_RMSSD$value[df_RMSSD$phase=="pre"   & df_RMSSD$condition=="control"],
       paired = TRUE, alternative = "greater")  # p = 0.7059
t.test(df_RMSSD$value[df_RMSSD$phase=="early" & df_RMSSD$condition=="experimental"],
       df_RMSSD$value[df_RMSSD$phase=="early" & df_RMSSD$condition=="control"],
       paired = TRUE, alternative = "greater")  # p = 0.7155
t.test(df_RMSSD$value[df_RMSSD$phase=="late"  & df_RMSSD$condition=="experimental"],
       df_RMSSD$value[df_RMSSD$phase=="late"  & df_RMSSD$condition=="control"],
       paired = TRUE, alternative = "greater")  # p = 0.7474
t.test(df_RMSSD$value[df_RMSSD$phase=="post"  & df_RMSSD$condition=="experimental"],
       df_RMSSD$value[df_RMSSD$phase=="post"  & df_RMSSD$condition=="control"],
       paired = TRUE, alternative = "greater")  # p = 0.6699

saveRDS(df_RMSSD, file.path(hrv_path, "data_df_RMSSD_n15.Rds"))
plot_signal_split(df_RMSSD, "RMSSD (n = 15)", "log RMSSD (ms)", "fig_RMSSD")


### 1b. HFn ####################################################################
# Approximately normal in raw form -> entered untransformed

shapiro.test(hrv_ecg$HFn[hrv_ecg$phase == "pre"])    # p = 0.3632
shapiro.test(hrv_ecg$HFn[hrv_ecg$phase == "early"])  # p = 0.0963
shapiro.test(hrv_ecg$HFn[hrv_ecg$phase == "late"])   # p = 0.8724
shapiro.test(hrv_ecg$HFn[hrv_ecg$phase == "post"])   # p = 0.0376 (!) borderline

# Normality of difference variables (raw)
hfn_wide <- hrv_ecg %>%
  dplyr::select(VP, ses_id, condition, phase, HFn) %>%
  tidyr::pivot_wider(names_from = phase, values_from = HFn,
                     names_prefix = "HFn_") %>%
  dplyr::mutate(
    diff_pre_early  = HFn_pre   - HFn_early,
    diff_pre_late   = HFn_pre   - HFn_late,
    diff_pre_post   = HFn_pre   - HFn_post,
    diff_early_late = HFn_early - HFn_late,
    diff_early_post = HFn_early - HFn_post,
    diff_late_post  = HFn_late  - HFn_post
  )

shapiro.test(hfn_wide$diff_pre_early)   # p = 0.4223
shapiro.test(hfn_wide$diff_pre_late)    # p = 0.6907
shapiro.test(hfn_wide$diff_pre_post)    # p = 0.5932
shapiro.test(hfn_wide$diff_early_late)  # p = 0.2995
shapiro.test(hfn_wide$diff_early_post)  # p = 0.1215
shapiro.test(hfn_wide$diff_late_post)   # p = 0.2008

saveRDS(hfn_wide, file.path(hrv_path, "data_hfn_wide_n15.Rds"))

df_HFn <- make_long_complete(hrv_all, "HFn", excl_ecg)  # n = 15

ezStats(data = df_HFn, wid = VP, dv = value, within = .(phase, condition))

anova_HFn <- ezANOVA(data = df_HFn, wid = VP, dv = value,
                      within = .(phase, condition),
                      detailed = TRUE, return_aov = TRUE)
anova_HFn$Mauchly
anova_HFn$Sphericity
anova_HFn$ANOVA

get_anova_table(anova_test(df_HFn, dv = value, wid = VP,
                            within = c(phase, condition)), correction = "HF")

df_HFn %>%
  rstatix::pairwise_t_test(value ~ phase, pool.sd = FALSE, paired = TRUE,
                           p.adjust.method = "bonferroni", alternative = "two.sided")

# Primary hypothesis (exp > ctrl) per phase
t.test(df_HFn$value[df_HFn$phase=="pre"   & df_HFn$condition=="experimental"],
       df_HFn$value[df_HFn$phase=="pre"   & df_HFn$condition=="control"],
       paired = TRUE, alternative = "greater")  # p = 0.6052
t.test(df_HFn$value[df_HFn$phase=="early" & df_HFn$condition=="experimental"],
       df_HFn$value[df_HFn$phase=="early" & df_HFn$condition=="control"],
       paired = TRUE, alternative = "greater")  # p = 0.9165
t.test(df_HFn$value[df_HFn$phase=="late"  & df_HFn$condition=="experimental"],
       df_HFn$value[df_HFn$phase=="late"  & df_HFn$condition=="control"],
       paired = TRUE, alternative = "greater")  # p = 0.1332
t.test(df_HFn$value[df_HFn$phase=="post"  & df_HFn$condition=="experimental"],
       df_HFn$value[df_HFn$phase=="post"  & df_HFn$condition=="control"],
       paired = TRUE, alternative = "greater")  # p = 0.7865

saveRDS(df_HFn, file.path(hrv_path, "data_df_HFn_n15.Rds"))
plot_signal_split(df_HFn, "HFn (n = 15)", "Normalised HF Power (HFn)", "fig_HFn")


### 1c. HF absolute (ms²) ######################################################

shapiro.test(log(hrv_ecg$HF[hrv_ecg$phase=="pre"]   * S2_TO_MS2 + 1))  # p = 0.6608
shapiro.test(log(hrv_ecg$HF[hrv_ecg$phase=="early"] * S2_TO_MS2 + 1))  # p = 0.1579
shapiro.test(log(hrv_ecg$HF[hrv_ecg$phase=="late"]  * S2_TO_MS2 + 1))  # p = 0.0113 (!)
shapiro.test(log(hrv_ecg$HF[hrv_ecg$phase=="post"]  * S2_TO_MS2 + 1))  # p = 0.2634

df_HF <- make_long_complete(hrv_all, "HF_ms2_log", excl_ecg)  # n = 15

ezStats(data = df_HF, wid = VP, dv = value, within = .(phase, condition))

anova_HF <- ezANOVA(data = df_HF, wid = VP, dv = value,
                     within = .(phase, condition),
                     detailed = TRUE, return_aov = TRUE)
anova_HF$Mauchly
anova_HF$Sphericity
anova_HF$ANOVA

get_anova_table(anova_test(df_HF, dv = value, wid = VP,
                            within = c(phase, condition)), correction = "HF")

df_HF %>%
  rstatix::pairwise_t_test(value ~ phase, pool.sd = FALSE, paired = TRUE,
                           p.adjust.method = "bonferroni", alternative = "two.sided")

saveRDS(df_HF, file.path(hrv_path, "data_df_HF_abs_n15.Rds"))
plot_signal_split(df_HF, "HF absolute (n = 15)", "log HF Power (ms²)", "fig_HF_abs")


### 1d. Sensitivity: HFn n = 15 vs n = 20 #####################################

hrv_n20 <- hrv_all %>%
  dplyr::filter(!VP %in% c("sub-021")) %>%
  dplyr::mutate(VP = droplevels(as.factor(VP)), condition = as.factor(condition),
                phase = factor(phase, levels = c("pre","early","late","post")))

df_HFn_n20 <- make_long_complete(hrv_n20, "HFn", c())  # n = 20

anova_HFn_n20 <- ezANOVA(data = df_HFn_n20, wid = VP, dv = value,
                           within = .(phase, condition),
                           detailed = TRUE, return_aov = TRUE)
anova_HFn_n20$ANOVA
get_anova_table(anova_test(df_HFn_n20, dv = value, wid = VP,
                            within = c(phase, condition)), correction = "HF")


# ==============================================================================
# 2.  HEART RATE ANALYSIS  (ECG group, n = 15)
# ==============================================================================
# HR is read from a separate file — HR_log added here, not in hrv_all

data_HR <- read_tsv(
  file.path(hrv_path, "task-BBSIG_per_phase_hr_nk.tsv"),
  show_col_types = FALSE) %>%
  dplyr::rename(VP = subj_id) %>%
  dplyr::mutate(VP = as.factor(VP), condition = as.factor(condition),
                phase = factor(phase, levels = c("pre","early","late","post")),
                HR_log = log(HR_mean + 2))

### 2a. Normality ##############################################################
# Raw: non-normal at PRE -> log(x + 2) applied

shapiro.test(data_HR$HR_mean[data_HR$phase == "pre"])    # p = 0.0148 (!)
shapiro.test(data_HR$HR_mean[data_HR$phase == "early"])  # p = 0.0604
shapiro.test(data_HR$HR_mean[data_HR$phase == "late"])   # p = 0.1080
shapiro.test(data_HR$HR_mean[data_HR$phase == "post"])   # p = 0.0752

shapiro.test(data_HR$HR_log[data_HR$phase == "pre"])    # p = 0.1857
shapiro.test(data_HR$HR_log[data_HR$phase == "early"])  # p = 0.4801
shapiro.test(data_HR$HR_log[data_HR$phase == "late"])   # p = 0.3336
shapiro.test(data_HR$HR_log[data_HR$phase == "post"])   # p = 0.2081

data_HR_log_wide <- data_HR %>%
  dplyr::filter(!VP %in% excl_ecg) %>%
  dplyr::select(VP, ses_id, condition, phase, HR_log) %>%
  tidyr::pivot_wider(names_from = phase, values_from = HR_log,
                     names_prefix = "HR_") %>%
  dplyr::mutate(
    diff_pre_early  = HR_pre   - HR_early,
    diff_pre_late   = HR_pre   - HR_late,
    diff_pre_post   = HR_pre   - HR_post,
    diff_early_late = HR_early - HR_late,
    diff_early_post = HR_early - HR_post,
    diff_late_post  = HR_late  - HR_post
  )

shapiro.test(data_HR_log_wide$diff_pre_early)   # p = 0.1624
shapiro.test(data_HR_log_wide$diff_pre_late)    # p = 0.6213
shapiro.test(data_HR_log_wide$diff_pre_post)    # p = 0.0294 (!) skew = -0.806, kurt = 3.590
shapiro.test(data_HR_log_wide$diff_early_late)  # p = 0.6454
shapiro.test(data_HR_log_wide$diff_early_post)  # p = 0.1532
shapiro.test(data_HR_log_wide$diff_late_post)   # p = 0.0053 (!) skew = -1.290, kurt = 5.106

skewness(data_HR_log_wide$diff_pre_post);  skewness(data_HR_log_wide$diff_pre_post)  / 0.455
kurtosis(data_HR_log_wide$diff_pre_post);  kurtosis(data_HR_log_wide$diff_pre_post)  / 0.91
skewness(data_HR_log_wide$diff_late_post); skewness(data_HR_log_wide$diff_late_post) / 0.455
kurtosis(data_HR_log_wide$diff_late_post); kurtosis(data_HR_log_wide$diff_late_post) / 0.91

saveRDS(data_HR_log_wide, file.path(hrv_path, "data_HR_log_wide_n15.Rds"))

### 2b. ANOVA ##################################################################

df_HR <- make_long_complete(data_HR, "HR_log", excl_ecg)  # n = 15

ezStats(data = df_HR, wid = VP, dv = value, within = .(phase, condition))

# Raw means for reporting
data_HR %>%
  dplyr::filter(VP %in% unique(df_HR$VP)) %>%
  dplyr::group_by(phase, condition) %>%
  dplyr::summarise(M  = round(mean(HR_mean, na.rm = TRUE), 2),
                   SD = round(sd(HR_mean,   na.rm = TRUE), 2), .groups = "drop")

anova_HR <- ezANOVA(data = df_HR, wid = VP, dv = value,
                     within = .(phase, condition),
                     detailed = TRUE, return_aov = TRUE)
anova_HR$Mauchly
anova_HR$Sphericity
anova_HR$ANOVA

get_anova_table(anova_test(df_HR, dv = value, wid = VP,
                            within = c(phase, condition)), correction = "HF")

df_HR %>%
  rstatix::pairwise_t_test(value ~ phase, pool.sd = FALSE, paired = TRUE,
                           p.adjust.method = "bonferroni", alternative = "two.sided")

# Primary hypothesis (exp < ctrl) per phase
t.test(df_HR$value[df_HR$phase=="pre"   & df_HR$condition=="experimental"],
       df_HR$value[df_HR$phase=="pre"   & df_HR$condition=="control"],
       paired = TRUE, alternative = "less")  # p = 0.9036
t.test(df_HR$value[df_HR$phase=="early" & df_HR$condition=="experimental"],
       df_HR$value[df_HR$phase=="early" & df_HR$condition=="control"],
       paired = TRUE, alternative = "less")  # p = 0.9345
t.test(df_HR$value[df_HR$phase=="late"  & df_HR$condition=="experimental"],
       df_HR$value[df_HR$phase=="late"  & df_HR$condition=="control"],
       paired = TRUE, alternative = "less")  # p = 0.9689
t.test(df_HR$value[df_HR$phase=="post"  & df_HR$condition=="experimental"],
       df_HR$value[df_HR$phase=="post"  & df_HR$condition=="control"],
       paired = TRUE, alternative = "less")  # p = 0.7465

saveRDS(df_HR, file.path(hrv_path, "data_df_HR_n15.Rds"))

data_HR_raw_complete <- data_HR %>%
  dplyr::filter(VP %in% unique(df_HR$VP)) %>%
  dplyr::select(VP, ses_id, condition, phase, HR_mean) %>%
  dplyr::rename(value = HR_mean) %>%
  dplyr::mutate(VP = droplevels(as.factor(VP)), condition = as.factor(condition),
                phase = factor(phase, levels = c("pre","early","late","post")))
plot_signal_split(data_HR_raw_complete, "Heart Rate (n = 15)",
                  "Heart Rate (bpm)", "fig_HR")


# ==============================================================================
# 3.  BREATHING FREQUENCY  (ECG group, n = 15)
# ==============================================================================
# Normally distributed in raw form (all p > .08) -> no transformation applied

data_BF <- read_tsv(
  file.path(bf_path, "task-BBSIG_per_phase_bf.tsv"),
  show_col_types = FALSE) %>%
  dplyr::rename(VP = subj_id) %>%
  dplyr::mutate(VP = as.factor(VP), condition = as.factor(condition),
                phase = factor(phase, levels = c("pre","early","late","post")))

data_BF_wide <- data_BF %>%
  dplyr::filter(!VP %in% excl_ecg, !is.na(BF_bpm)) %>%
  dplyr::select(VP, ses_id, condition, phase, BF_bpm) %>%
  tidyr::pivot_wider(names_from = phase, values_from = BF_bpm,
                     names_prefix = "BF_")

### 3a. Normality ##############################################################

shapiro.test(data_BF_wide$BF_pre)    # p = 0.3960
shapiro.test(data_BF_wide$BF_early)  # p = 0.0843
shapiro.test(data_BF_wide$BF_late)   # p = 0.1389
shapiro.test(data_BF_wide$BF_post)   # p = 0.9060

data_BF_wide <- data_BF_wide %>% dplyr::mutate(
  diff_pre_early  = BF_pre   - BF_early,
  diff_pre_late   = BF_pre   - BF_late,
  diff_pre_post   = BF_pre   - BF_post,
  diff_early_late = BF_early - BF_late,
  diff_early_post = BF_early - BF_post,
  diff_late_post  = BF_late  - BF_post
)

shapiro.test(data_BF_wide$diff_pre_early)   # p = 0.7751
shapiro.test(data_BF_wide$diff_pre_late)    # p = 0.1651
shapiro.test(data_BF_wide$diff_pre_post)    # p = 0.6401
shapiro.test(data_BF_wide$diff_early_late)  # p = 0.1184
shapiro.test(data_BF_wide$diff_early_post)  # p = 0.2378
shapiro.test(data_BF_wide$diff_late_post)   # p = 0.1493

saveRDS(data_BF_wide, file.path(bf_path, "data_BF_raw_wide_n15.Rds"))

### 3b. ANOVA ##################################################################

data_BF_long <- data_BF_wide %>%
  tidyr::pivot_longer(cols = c(BF_pre, BF_early, BF_late, BF_post),
                      names_to = "phase", values_to = "bf_values") %>%
  dplyr::mutate(phase = factor(gsub("BF_", "", phase),
                               levels = c("pre","early","late","post")))

complete_VPs_BF <- data_BF_long %>%
  dplyr::filter(!is.na(bf_values)) %>%
  dplyr::group_by(VP) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n == 8) %>%
  dplyr::pull(VP)

df_BF <- data_BF_long %>%
  dplyr::filter(VP %in% complete_VPs_BF) %>%
  dplyr::mutate(VP = droplevels(VP)) %>%
  dplyr::rename(value = bf_values)  # n = 15

ezStats(data = df_BF, wid = VP, dv = value, within = .(phase, condition))

anova_BF <- ezANOVA(data = df_BF, wid = VP, dv = value,
                     within = .(phase, condition),
                     detailed = TRUE, return_aov = TRUE)
anova_BF$Mauchly
anova_BF$Sphericity
anova_BF$ANOVA

get_anova_table(anova_test(df_BF, dv = value, wid = VP,
                            within = c(phase, condition)), correction = "HF")

df_BF %>%
  rstatix::pairwise_t_test(value ~ phase, pool.sd = FALSE, paired = TRUE,
                           p.adjust.method = "bonferroni", alternative = "two.sided")

# Primary hypothesis (exp < ctrl) per phase
t.test(df_BF$value[df_BF$phase=="pre"   & df_BF$condition=="experimental"],
       df_BF$value[df_BF$phase=="pre"   & df_BF$condition=="control"],
       paired = TRUE, alternative = "less")  # p = 0.3226
t.test(df_BF$value[df_BF$phase=="early" & df_BF$condition=="experimental"],
       df_BF$value[df_BF$phase=="early" & df_BF$condition=="control"],
       paired = TRUE, alternative = "less")  # p = 0.9422
t.test(df_BF$value[df_BF$phase=="late"  & df_BF$condition=="experimental"],
       df_BF$value[df_BF$phase=="late"  & df_BF$condition=="control"],
       paired = TRUE, alternative = "less")  # p = 0.8842
t.test(df_BF$value[df_BF$phase=="post"  & df_BF$condition=="experimental"],
       df_BF$value[df_BF$phase=="post"  & df_BF$condition=="control"],
       paired = TRUE, alternative = "less")  # p = 0.4778

saveRDS(df_BF, file.path(bf_path, "data_df_BF_n15.Rds"))

data_BF_raw_long <- data_BF %>%
  dplyr::filter(VP %in% complete_VPs_BF, !is.na(BF_bpm)) %>%
  dplyr::select(VP, ses_id, condition, phase, BF_bpm) %>%
  dplyr::rename(value = BF_bpm) %>%
  dplyr::mutate(VP = droplevels(as.factor(VP)), condition = as.factor(condition),
                phase = factor(phase, levels = c("pre","early","late","post")))
plot_signal_split(data_BF_raw_long, "Breathing Frequency (n = 15)",
                  "Breathing Frequency (br/min)", "fig_BF")


# ==============================================================================
# 4.  EDA & SKIN TEMPERATURE
# ==============================================================================

data_raw <- read_tsv(
  file.path(eda_path, "task-BBSIG_eda-temp-bf_per_phase_final.tsv"),
  show_col_types = FALSE) %>%
  dplyr::rename(VP = subj_id) %>%
  dplyr::mutate(VP = as.factor(VP), condition = as.factor(condition),
                phase = factor(phase, levels = c("pre","early","late","post")))

janitor::get_dupes(data_raw)  # no duplicates expected


### 4a. EDA ####################################################################
# Severely non-normal in raw form -> log(x + 2) applied

data_EDA_wide <- data_raw %>%
  dplyr::filter(!VP %in% excl_eda, !EDA_excluded) %>%
  dplyr::select(VP, ses_id, condition, phase, EDA_mean) %>%
  tidyr::pivot_wider(names_from = phase, values_from = EDA_mean,
                     names_prefix = "EDA_")

shapiro.test(data_EDA_wide$EDA_pre)    # p = 0.0000 (!)
shapiro.test(data_EDA_wide$EDA_early)  # p = 0.0002 (!)
shapiro.test(data_EDA_wide$EDA_late)   # p = 0.0028 (!)
shapiro.test(data_EDA_wide$EDA_post)   # p = 0.0212 (!)

data_EDA_log <- data_EDA_wide %>% dplyr::mutate(
  EDA_pre   = log(EDA_pre   + 2),
  EDA_early = log(EDA_early + 2),
  EDA_late  = log(EDA_late  + 2),
  EDA_post  = log(EDA_post  + 2)
)

shapiro.test(data_EDA_log$EDA_pre)    # p = 0.4325
shapiro.test(data_EDA_log$EDA_early)  # p = 0.7623
shapiro.test(data_EDA_log$EDA_late)   # p = 0.0702
shapiro.test(data_EDA_log$EDA_post)   # p = 0.3857

data_EDA_log <- data_EDA_log %>% dplyr::mutate(
  diff_pre_early  = EDA_pre   - EDA_early,
  diff_pre_late   = EDA_pre   - EDA_late,
  diff_pre_post   = EDA_pre   - EDA_post,
  diff_early_late = EDA_early - EDA_late,
  diff_early_post = EDA_early - EDA_post,
  diff_late_post  = EDA_late  - EDA_post
)

shapiro.test(data_EDA_log$diff_pre_early)   # p = 0.0668
shapiro.test(data_EDA_log$diff_pre_late)    # p = 0.1053
shapiro.test(data_EDA_log$diff_pre_post)    # p = 0.1992
shapiro.test(data_EDA_log$diff_early_late)  # p = 0.0233 (!) skew = 0.829, kurt = 2.814
shapiro.test(data_EDA_log$diff_early_post)  # p = 0.0049 (!) skew = -1.352, kurt = 6.135
shapiro.test(data_EDA_log$diff_late_post)   # p = 0.1831

skewness(data_EDA_log$diff_early_late);  skewness(data_EDA_log$diff_early_late) / 0.455
kurtosis(data_EDA_log$diff_early_late);  kurtosis(data_EDA_log$diff_early_late) / 0.91
skewness(data_EDA_log$diff_early_post);  skewness(data_EDA_log$diff_early_post) / 0.455
kurtosis(data_EDA_log$diff_early_post);  kurtosis(data_EDA_log$diff_early_post) / 0.91

saveRDS(data_EDA_log, file.path(eda_path, "data_EDA_log_n16.Rds"))

data_EDA_log_long <- data_EDA_log %>%
  tidyr::pivot_longer(cols = c(EDA_pre, EDA_early, EDA_late, EDA_post),
                      names_to = "phase", values_to = "eda_values") %>%
  dplyr::mutate(phase = factor(gsub("EDA_", "", phase),
                               levels = c("pre","early","late","post")))

complete_VPs_EDA <- data_EDA_log_long %>%
  dplyr::filter(!is.na(eda_values)) %>%
  dplyr::group_by(VP) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n == 8) %>%
  dplyr::pull(VP)

data_EDA_log_long1 <- data_EDA_log_long %>%
  dplyr::filter(VP %in% complete_VPs_EDA) %>%
  dplyr::mutate(VP = droplevels(VP))  # n = 16

ezStats(data = data_EDA_log_long1, wid = VP, dv = eda_values,
        within = .(phase, condition))

anova_EDA <- ezANOVA(data = data_EDA_log_long1, wid = VP, dv = eda_values,
                     within = .(phase, condition),
                     detailed = TRUE, return_aov = TRUE)
anova_EDA$Mauchly
anova_EDA$Sphericity
anova_EDA$ANOVA

get_anova_table(anova_test(data_EDA_log_long1, dv = eda_values, wid = VP,
                            within = c(phase, condition)), correction = "HF")

data_EDA_log_long1 %>%
  rstatix::pairwise_t_test(eda_values ~ phase, pool.sd = FALSE, paired = TRUE,
                           p.adjust.method = "bonferroni", alternative = "two.sided")

# Primary hypothesis (exp < ctrl) per phase
t.test(data_EDA_log_long1$eda_values[data_EDA_log_long1$phase=="pre"   & data_EDA_log_long1$condition=="experimental"],
       data_EDA_log_long1$eda_values[data_EDA_log_long1$phase=="pre"   & data_EDA_log_long1$condition=="control"],
       paired = TRUE, alternative = "less")  # p = 0.0345 *
t.test(data_EDA_log_long1$eda_values[data_EDA_log_long1$phase=="early" & data_EDA_log_long1$condition=="experimental"],
       data_EDA_log_long1$eda_values[data_EDA_log_long1$phase=="early" & data_EDA_log_long1$condition=="control"],
       paired = TRUE, alternative = "less")  # p = 0.1487
t.test(data_EDA_log_long1$eda_values[data_EDA_log_long1$phase=="late"  & data_EDA_log_long1$condition=="experimental"],
       data_EDA_log_long1$eda_values[data_EDA_log_long1$phase=="late"  & data_EDA_log_long1$condition=="control"],
       paired = TRUE, alternative = "less")  # p = 0.3027
t.test(data_EDA_log_long1$eda_values[data_EDA_log_long1$phase=="post"  & data_EDA_log_long1$condition=="experimental"],
       data_EDA_log_long1$eda_values[data_EDA_log_long1$phase=="post"  & data_EDA_log_long1$condition=="control"],
       paired = TRUE, alternative = "less")  # p = 0.1323

saveRDS(data_EDA_log_long1, file.path(eda_path, "data_df_EDA_n16.Rds"))

data_EDA_long <- data_raw %>%
  dplyr::filter(!VP %in% excl_eda, !EDA_excluded) %>%
  dplyr::select(VP, ses_id, condition, phase, EDA_mean) %>%
  dplyr::rename(value = EDA_mean) %>%
  dplyr::filter(VP %in% complete_VPs_EDA) %>%
  dplyr::mutate(VP = droplevels(as.factor(VP)), condition = as.factor(condition),
                phase = factor(phase, levels = c("pre","early","late","post")))
plot_signal_split(data_EDA_long, "EDA (n = 16)", "EDA (µS)", "fig_EDA")


### 4b. Skin Temperature #######################################################
# Normally distributed in raw form after exclusions (all p > .14)
# -> no transformation applied

data_TEMP_wide <- data_raw %>%
  dplyr::filter(!VP %in% excl_temp, !temp_excluded) %>%
  dplyr::select(VP, ses_id, condition, phase, temp_mean) %>%
  tidyr::pivot_wider(names_from = phase, values_from = temp_mean,
                     names_prefix = "TEMP_")

shapiro.test(data_TEMP_wide$TEMP_pre)    # p = 0.1769
shapiro.test(data_TEMP_wide$TEMP_early)  # p = 0.5785
shapiro.test(data_TEMP_wide$TEMP_late)   # p = 0.3748
shapiro.test(data_TEMP_wide$TEMP_post)   # p = 0.4845

data_TEMP_wide <- data_TEMP_wide %>% dplyr::mutate(
  diff_pre_early  = TEMP_pre   - TEMP_early,
  diff_pre_late   = TEMP_pre   - TEMP_late,
  diff_pre_post   = TEMP_pre   - TEMP_post,
  diff_early_late = TEMP_early - TEMP_late,
  diff_early_post = TEMP_early - TEMP_post,
  diff_late_post  = TEMP_late  - TEMP_post
)

shapiro.test(data_TEMP_wide$diff_pre_early)   # p = 0.2490
shapiro.test(data_TEMP_wide$diff_pre_late)    # p = 0.6323
shapiro.test(data_TEMP_wide$diff_pre_post)    # p = 0.3727
shapiro.test(data_TEMP_wide$diff_early_late)  # p = 0.8310
shapiro.test(data_TEMP_wide$diff_early_post)  # p = 0.4197
shapiro.test(data_TEMP_wide$diff_late_post)   # p = 0.0145 (!) skew = 1.162, kurt = 5.601

skewness(data_TEMP_wide$diff_late_post);  skewness(data_TEMP_wide$diff_late_post) / 0.455
kurtosis(data_TEMP_wide$diff_late_post);  kurtosis(data_TEMP_wide$diff_late_post) / 0.91

saveRDS(data_TEMP_wide, file.path(eda_path, "data_TEMP_raw_wide_n12.Rds"))

data_TEMP_long_analysis <- data_TEMP_wide %>%
  tidyr::pivot_longer(cols = c(TEMP_pre, TEMP_early, TEMP_late, TEMP_post),
                      names_to = "phase", values_to = "temp_values") %>%
  dplyr::mutate(phase = factor(gsub("TEMP_", "", phase),
                               levels = c("pre","early","late","post")))

complete_VPs_TEMP <- data_TEMP_long_analysis %>%
  dplyr::filter(!is.na(temp_values)) %>%
  dplyr::group_by(VP) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n == 8) %>%
  dplyr::pull(VP)

data_TEMP_long1 <- data_TEMP_long_analysis %>%
  dplyr::filter(VP %in% complete_VPs_TEMP) %>%
  dplyr::mutate(VP = droplevels(VP))  # n = 12

ezStats(data = data_TEMP_long1, wid = VP, dv = temp_values,
        within = .(phase, condition))

anova_TEMP <- ezANOVA(data = data_TEMP_long1, wid = VP, dv = temp_values,
                      within = .(phase, condition),
                      detailed = TRUE, return_aov = TRUE)
anova_TEMP$Mauchly
anova_TEMP$Sphericity
anova_TEMP$ANOVA

get_anova_table(anova_test(data_TEMP_long1, dv = temp_values, wid = VP,
                            within = c(phase, condition)), correction = "HF")

pairwise.t.test(data_TEMP_long1$temp_values, data_TEMP_long1$phase,
                p.adj = "bonferroni", paired = TRUE)

data_TEMP_long1 %>%
  rstatix::pairwise_t_test(temp_values ~ phase, pool.sd = FALSE, paired = TRUE,
                           p.adjust.method = "bonferroni", alternative = "two.sided")

# Condition comparison (two-sided) per phase
t.test(data_TEMP_long1$temp_values[data_TEMP_long1$phase=="pre"   & data_TEMP_long1$condition=="experimental"],
       data_TEMP_long1$temp_values[data_TEMP_long1$phase=="pre"   & data_TEMP_long1$condition=="control"],
       paired = TRUE, alternative = "two.sided")  # p = 0.8652
t.test(data_TEMP_long1$temp_values[data_TEMP_long1$phase=="early" & data_TEMP_long1$condition=="experimental"],
       data_TEMP_long1$temp_values[data_TEMP_long1$phase=="early" & data_TEMP_long1$condition=="control"],
       paired = TRUE, alternative = "two.sided")  # p = 0.9931
t.test(data_TEMP_long1$temp_values[data_TEMP_long1$phase=="late"  & data_TEMP_long1$condition=="experimental"],
       data_TEMP_long1$temp_values[data_TEMP_long1$phase=="late"  & data_TEMP_long1$condition=="control"],
       paired = TRUE, alternative = "two.sided")  # p = 0.6577
t.test(data_TEMP_long1$temp_values[data_TEMP_long1$phase=="post"  & data_TEMP_long1$condition=="experimental"],
       data_TEMP_long1$temp_values[data_TEMP_long1$phase=="post"  & data_TEMP_long1$condition=="control"],
       paired = TRUE, alternative = "two.sided")  # p = 0.5885

saveRDS(data_TEMP_long1, file.path(eda_path, "data_df_TEMP_n12.Rds"))

data_TEMP_long <- data_raw %>%
  dplyr::filter(!VP %in% excl_temp, !temp_excluded) %>%
  dplyr::select(VP, ses_id, condition, phase, temp_mean) %>%
  dplyr::rename(value = temp_mean) %>%
  dplyr::filter(VP %in% complete_VPs_TEMP) %>%
  dplyr::mutate(VP = droplevels(as.factor(VP)), condition = as.factor(condition),
                phase = factor(phase, levels = c("pre","early","late","post")))
plot_signal_split(data_TEMP_long, "Skin Temperature (n = 12)",
                  "Skin Temperature (°C)", "fig_TEMP")


# ==============================================================================
# Final sample sizes
#   RMSSD, HFn, HF absolute, HR, BF  -> n = 15
#   EDA                               -> n = 16
#   Skin Temperature                  -> n = 12
# ==============================================================================
