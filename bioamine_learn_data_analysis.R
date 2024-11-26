###############################################################################
# Importing libraries
###############################################################################

library(DHARMa) # for diagnostics of LMMs
library(lmerTest) # for building LMMs using overloaded lmer from lme4 
library(broom.mixed) # for extracting estimates, conf. intervals, and p-values from LMMs
library(ggpubr) # for creating multi-panel figures
library(ggplot2) # for plots
library(ggsci) # for npg palette
library(dplyr) # for wrangling datasets
library(tidyr) # for pivoting dataframes
library(patchwork)
library(stringr)

###############################################################################
# Data manipulation
###############################################################################

# import data
als <- read.csv("BA_Data.csv")

# re-assign factors 
als$hive <- as.factor(als$hive)
als$trial <- as.factor(als$trial)

###############################################################################
# Distribution of max. learning accuracy in individual and social contexts
###############################################################################

# reshape dataframe
als_speedAccuracy <- als %>%
  pivot_longer(c("socM", "indM", "socS", "indS"), names_to = "var", values_to = "score") %>%
  # pivot_longer(c("socS", "indS"), names_to = "speed", values_to = "speed_score") %>%
  select(c("var", "score")) %>% # select learning accuracy scores
  group_by(var) %>%
  mutate(score_mean = mean(score)) %>% # calculate averages
  ungroup %>%
  group_by(speed) %>%
  mutate(speed_mean = mean(speed_score)) %>% # calculate averages
  ungroup

# plot max. learning accuracy in individual and social contexts
acc_dist <- als_speedAccuracy %>%
  filter(str_detect(var, "M")) %>%
  ggplot(aes(x = score, fill = var)) +
  geom_histogram(position = "dodge", bins = 15, color = "black") +
  labs(x = "Accuracy", y = "Number of Bees") +
  scale_fill_manual(limits = c("indM", "socM"), 
                    labels = c("Individual", "Social"), 
                    values = c("white", "black")) +
  scale_y_continuous(expand = c(0,0), 
                     # limits = c(0, 14),
                     breaks = seq(0, 14, 2)) +
  scale_x_continuous(breaks = seq(0, 1, 0.1),
                     labels = as.character(seq(0, 1, 0.1))) +
  geom_vline(aes(xintercept = score_mean, linetype=var)) + 
  geom_vline(aes(xintercept = 0.5), linetype = "dotted") + 
  scale_linetype_manual(limits=c("indM", "socM"),
                        labels=c("Individual", "Social"),
                        values=c("dashed", "solid")) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(size=12),
        axis.text = element_text(size=10, color = "black"))

speed_dist <- als_speedAccuracy %>%
  filter(str_detect(var, "S")) %>%
  ggplot( aes(x = score, fill = var)) +
  geom_histogram(position = "dodge", bins = 30, color = "black") +
  labs(x = "Speed", y = "Number of Bees") +
  scale_fill_manual(limits = c("indS", "socS"), 
                    labels = c("Individual", "Social"), 
                    values = c("white", "black")) +
  scale_y_continuous(expand = c(0,0), 
                     limits = c(0, 9),
                     breaks = seq(0, 9, 1)) +
  scale_x_continuous(breaks = seq(-0.05, 0.13, 0.02)) +
  geom_vline(aes(xintercept = score_mean, linetype = var)) + 
  scale_linetype_manual(limits=c("indS", "socS"),
                        labels=c("Individual", "Social"),
                        values=c("dashed", "solid")) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(size=12),
        axis.text = element_text(size=10, color = "black"))

png("fig2_accuracy_speed_dist.png", width=1500, height=2000, res=300)
speed_dist/ acc_dist
dev.off()

# statistical tests
## differences in distribution
### two-sided KS test
ks.test(als$socM, als$indM)
ks.test(als$socS, als$indS)
## mean max. learning accuracy in individual and social contexts > 0.5
### one-sample t-test
t.test(als$socM, mu = 0.5, alternative = "greater")
t.test(als$indM, mu = 0.5, alternative = "greater")

###############################################################################
# Fast vs. accurate learners in individual and social contexts
###############################################################################

learnSpdAcc_tbl <- data.frame(learnSpd = c("Individual", "Social", "Individual", "Social"),
                              learnAcc = c("Individual", "Individual", "Social", "Social"),
                              count = c(round(prop.table(table(als$learnSpd, als$learnAcc), 1)*100, 2)))

png("fig1_cont_tbl.png", width=1500, height=1000, res=300)
ggplot(learnSpdAcc_tbl, aes(x = learnSpd, y = learnAcc, fill = as.factor(count))) +
  geom_tile() +
  scale_fill_manual(values = c("#4DBBD5FF", "#4DBBD5FF", "#E64B35FF", "#E64B35FF")) +
  geom_text(aes(label = paste0(count, "%")), vjust = 1) +
  xlab("Fast learning speed") +
  ylab("High learning accuracy") +
  ggtitle("") +
  theme_bw() + 
  theme(legend.position = "none", panel.grid = element_blank(),
        # axis.text = element_text(color = "black", size = 12),
        # axis.title = element_text(size = 14)
        )
dev.off()

fisher.test(als$learnSpd, als$learnAcc)

###############################################################################
# Building LMMs with interactions for biogenic amines
###############################################################################

# standardize data for LMMs
scaled_als <- data.frame("DA" = scale(als$DA, TRUE, TRUE),
                         "Ser" = scale(als$Ser, TRUE, TRUE),
                         "OA" = scale(als$OA, TRUE, TRUE),
                         "TA" = scale(als$TA, TRUE, TRUE),
                         "Glu" = scale(als$Glu, TRUE, TRUE),
                         "socM" = scale(als$socM, TRUE, TRUE),
                         "indM" = scale(als$indM, TRUE, TRUE),
                         "socS" = scale(als$socS, TRUE, TRUE),
                         "indS" = scale(als$indS, TRUE, TRUE),
                         "age" = scale(als$age, TRUE, TRUE),
                         "hive" = als$hive,
                         "trial" = als$trial)

set.seed(2024) # for reproducibility

socS_mod <- lmer(socS ~ (DA + Ser + OA + TA + Glu + age)^2
                 + (DA + Ser + OA + TA + Glu |hive)
                 + (DA + Ser + OA + TA + Glu |trial),
                 data = scaled_als)

socS_mod@optinfo$conv # check for optimizer convergence

socS_resid <- simulateResiduals(socS_mod, n=1000) # simulate residuals
testDispersion(socS_resid) # check for over-/underdispersion
testOutliers(socS_resid) # check outliers
plot(socS_resid) # plot diagnostics

indS_mod <- lmer(indS ~ (DA + Ser + OA + TA + Glu + age)^2
                 + (DA + Ser + OA + TA + Glu|hive)
                 + (DA + Ser + OA + TA + Glu|trial),
                 data = scaled_als)

indS_mod@optinfo$conv # check for optimizer convergence

indS_resid <- simulateResiduals(indS_mod, n=1000)
testDispersion(indS_resid)
testOutliers(indS_resid)
plot(indS_resid)

socM_mod <- lmer(socM ~ (DA + Ser + OA + TA + Glu + age)^2
                 + (DA + Ser + OA + TA + Glu |hive)
                 + (DA + Ser + OA + TA + Glu |trial),
                 data = scaled_als)

socM_mod@optinfo$conv # check for optimizer convergence

socM_resid <- simulateResiduals(socM_mod, n=1000)
testDispersion(socM_resid)
testOutliers(socM_resid)
plot(socM_resid)

indM_mod <- lmer(indM ~ (DA + Ser + OA + TA + Glu + age)^2
                 + (DA + Ser + OA + TA + Glu|hive)
                 + (DA + Ser + OA + TA + Glu|trial),
                 data = scaled_als)

indM_mod@optinfo$conv # check for optimizer convergence

indM_resid <- simulateResiduals(indM_mod, n=1000)
testDispersion(indM_resid)
testOutliers(indM_resid)
plot(indM_resid)

###############################################################################
# Building LMMs without interactions
###############################################################################

socS_mod2 <- lmer(socS ~ DA + Ser + OA + TA + Glu + age
                  + (DA + Ser + OA + TA + Glu|hive)
                  + (DA + Ser + OA + TA + Glu|trial),
                  data = scaled_als)

socS_mod2@optinfo$conv # check for optimizer convergence

indS_mod2 <- lmer(indS ~ DA + Ser + OA + TA + Glu + age
                  + (DA + Ser + OA + TA + Glu|hive)
                  + (DA + Ser + OA + TA + Glu|trial),
                  data = scaled_als)

indS_mod2@optinfo$conv # check for optimizer convergence

socM_mod2 <- lmer(socM ~ DA + Ser + OA + TA + Glu + age
                  + (DA + Ser + OA + TA + Glu|hive)
                  + (DA + Ser + OA + TA + Glu|trial),
                  data = scaled_als)

socM_mod2@optinfo$conv # check for optimizer convergence

indM_mod2 <- lmer(indM ~ DA + Ser + OA + TA + Glu + age
                  + (DA + Ser + OA + TA + Glu|hive)
                  + (DA + Ser + OA + TA + Glu|trial),
                  data = scaled_als)

indM_mod2@optinfo$conv # check for optimizer convergence

################################################################################
# Comparing LMMs with interactions with the ones without interactions using LRTs
################################################################################

anova(update(socS_mod, REML = FALSE), update(socS_mod2, REML = FALSE))
anova(update(indS_mod, REML = FALSE), update(indS_mod2, REML = FALSE))
anova(update(socM_mod, REML = FALSE), update(socM_mod2, REML = FALSE))
anova(update(indM_mod, REML = FALSE), update(indM_mod2, REML = FALSE))

################################################################################
# Extracting and plotting fixed effects from LMMs
################################################################################

add_sig_sym <- function(pvalues){
  ifelse(pvalues < 0.001, "***", 
         ifelse(pvalues < 0.01, "**",
                ifelse(pvalues < 0.05, "*", NA)))
}

# extracting estimates, confidence intervals, and p-values from LMMs
socS_est <- tidy(socS_mod, conf.int = TRUE, effects = "fixed") %>%
  mutate(p.stars = add_sig_sym(p.value)) %>% 
  filter(term != "(Intercept)") %>%
  mutate(group = ifelse(estimate > 0, "pos", "neg"))
indS_est <- tidy(indS_mod, conf.int = TRUE, effects = "fixed") %>%
  mutate(p.stars = add_sig_sym(p.value)) %>% 
  filter(term != "(Intercept)") %>%
  mutate(group = ifelse(estimate > 0, "pos", "neg"))
socM_est <- tidy(socM_mod, conf.int = TRUE, effects = "fixed") %>%
  mutate(p.stars = add_sig_sym(p.value)) %>% 
  filter(term != "(Intercept)") %>%
  mutate(group = ifelse(estimate > 0, "pos", "neg"))
indM_est <- tidy(indM_mod, conf.int = TRUE, effects = "fixed") %>%
  mutate(p.stars = add_sig_sym(p.value)) %>% 
  filter(term != "(Intercept)") %>%
  mutate(group = ifelse(estimate > 0, "pos", "neg"))

# term labels for y-axis in plots
y_label <- rev(c("DA" = "DA", 
                 "Ser" = "5-HT", 
                 "OA" = "OA",
                 "TA" = "TA",
                 "Glu" = "Glu", 
                 "DA:Ser" = "DA \u00D7 5-HT",
                 "DA:OA" = "DA \u00D7 OA",
                 "DA:TA" = "DA \u00D7 TA",
                 "DA:Glu" = "DA \u00D7 Glu", 
                 "Ser:OA" = "5-HT \u00D7 OA",
                 "Ser:TA" = "5-HT \u00D7 TA",
                 "Ser:Glu" = "5-HT \u00D7 Glu", 
                 "OA:TA" = "OA \u00D7 TA",
                 "OA:Glu" = "OA \u00D7 Glu",
                 "TA:Glu" = "TA \u00D7 Glu"))

# custom ggplot theme
lmm_plot_theme <- theme(panel.background=element_blank(), 
                        panel.border=element_blank(),
                        axis.line=element_line(color=1, size=0.1, linetype=1),
                        panel.grid.major=element_line(color="gray80", size=0.1, linetype=1), 
                        panel.grid.minor=element_line(color="gray80", size=0.1, linetype=1),
                        strip.background=element_blank(),
                        axis.text.x=element_text(color="black", size=8),
                        axis.text.y=element_text(color="black", size=8),
                        axis.title.x=element_text(size=8),
                        axis.title.y=element_text(size=8),
                        axis.ticks=element_blank(),
                        plot.title=element_text(hjust=0.5, size=10),
                        legend.position="top",
                        legend.text=element_text(size=8))

# socS LMM term estimates and significance
socS_plot <- ggplot(socS_est, aes(x=estimate, y=term, color=group, label=p.stars)) +
  geom_point(size=3) + ylab("") + xlab("Estimate") + 
  geom_text(vjust=0, nudge_y=-0.05, size=6, color=1) +
  geom_errorbar(aes(xmin=conf.low, xmax=conf.high, width=0.3)) +
  geom_vline(xintercept=0, linetype="dashed", color = "black", size=0.3) + 
  scale_color_manual(values=c("neg"="black", pos="grey"), guide="none") +
  scale_y_discrete(limits = names(y_label), labels = y_label) +
  ggtitle("Social Learning Speed") +
  lmm_plot_theme +
  theme(axis.text.y = element_text(face = c("plain", "bold", rep("plain", 3), "bold", rep('plain', 9))))

# indS LMM term estimates and significance
indS_plot <- ggplot(indS_est, aes(x=estimate, y=term, color=group, label=p.stars)) +
  geom_point(size=3) + ylab("") + xlab("Estimate") + 
  geom_text(vjust=0, nudge_y=-0.05, size=6, color=1) +
  geom_errorbar(aes(xmin=conf.low, xmax=conf.high, width=0.3)) +
  geom_vline(xintercept=0, linetype="dashed", color = "black", size=0.3) + 
  scale_color_manual(values=c("neg"="black", pos="grey"), guide="none") +
  scale_y_discrete(limits = names(y_label), labels = y_label) +
  ggtitle("Individual Learning Speed") +
  lmm_plot_theme +
  theme(axis.text.y = element_text(face = c("plain", "bold", rep("plain", 3), "bold", rep('plain', 9))))

# socM LMM term estimates and significance
socM_plot <- ggplot(socM_est, aes(x=estimate, y=term, color=group, label=p.stars)) +
  geom_point(size=3) + ylab("") + xlab("Estimate") + 
  geom_text(vjust=0, nudge_y=-0.05, size=6, color=1) +
  geom_errorbar(aes(xmin=conf.low, xmax=conf.high, width=0.3)) +
  geom_vline(xintercept=0, linetype="dashed", color = "black", size=0.3) + 
  scale_color_manual(values=c("neg"="black", pos="grey"), guide="none") +
  scale_y_discrete(limits = names(y_label), labels = y_label) +
  ggtitle("Social Learning Accuracy") +
  lmm_plot_theme +
  theme(axis.text.y = element_text(face = c("plain", "bold", rep("plain", 3), "bold", rep('plain', 9))))

# indM LMM term estimates and significance
indM_plot <- ggplot(indM_est, aes(x=estimate, y=term, color=group, label=p.stars)) +
  geom_point(size=3) + ylab("") + xlab("Estimate") + 
  geom_text(vjust=0, nudge_y=-0.05, size=6, color=1) +
  geom_errorbar(aes(xmin=conf.low, xmax=conf.high, width=0.3)) +
  geom_vline(xintercept=0, linetype="dashed", color = "black", size=0.3) + 
  scale_color_manual(values=c("neg"="black", pos="grey"), guide="none") +
  scale_y_discrete(limits = names(y_label), labels = y_label) +
  ggtitle("Individual Learning Accuracy") +
  lmm_plot_theme +
  theme(axis.text.y = element_text(face = c("plain", "bold", rep("plain", 3), "bold", rep('plain', 9))))

png("lmm_estimates.png", width=3000, height=3500, res=400)
# arranging aa plots in one figure and removing x-axis title
comb_plots <- ggarrange(indS_plot + rremove("ylab") + rremove("xlab"), 
                        indM_plot + rremove("y.text") + rremove("xlab"), 
                        socS_plot  + rremove("ylab") + rremove("xlab"), 
                        socM_plot  + rremove("y.text") + rremove("xlab"), 
                        labels="AUTO", ncol=2, nrow=2, font.label=list(size=10),
                        widths = c(1.1, 1, 1, 1))
# adding a common x-axis title
annotate_figure(comb_plots, bottom=text_grob("Estimate", hjust=0.5, size=10))
dev.off()
