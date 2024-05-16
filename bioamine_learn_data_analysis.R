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

###############################################################################
# Data manipulation
###############################################################################

# import data
als <- read.csv("bioamine_learn_data.csv")

# re-assign factors 
als$hive <- as.factor(als$hive)
als$trial <- as.factor(als$trial)

###############################################################################
# Behavioral data analysis
###############################################################################

# reshape dataframe
als_accuracy <- als %>%
  pivot_longer(c("socM", "indM"),
               names_to = "accuracy",
               values_to = "value") %>%
  select(c("accuracy", "value")) %>% # select learning accuracy scores
  group_by(accuracy) %>%
  mutate(mean = mean(value)) %>% # calculate averages
  ungroup

# plot max. learning accuracy in individual and social contexts
ggplot(als_accuracy, aes(x = value, fill = accuracy)) +
  geom_histogram(position = "identity", bins = 40, color = "white") +
  labs(title = "Two-sided Kolmogorov-Smirnov test: D = 0.32, p-value = 0.03", 
       x = "Accuracy", 
       y = "Number of bees", 
       fill = "") +
  scale_fill_manual(limits = c("indM", "socM"), 
                    labels = c("Individual", "Social"), 
                    values = c("#E64B35FF", "#4DBBD5FF")) +
  scale_x_continuous(breaks = seq(0, 1, 0.1), 
                     labels = as.character(seq(0, 1, 0.1))) +
  geom_vline(aes(xintercept = mean), 
             linetype = "dashed") + 
  geom_vline(aes(xintercept = 0.5), 
             linetype = "dotted") + 
  facet_wrap(factor(accuracy, 
                    levels = c("indM", "socM"), 
                    labels = c("Individual", "Social"))~., 
             nrow = 2) +
  theme_bw() +
  theme(legend.position = "none", plot.title = element_text(hjust=0.5))

###############################################################################
# Building LMMs with interactions
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

set.seed(1) # for reproducibility

socS_mod <- lmer(socS ~ (DA + Ser + OA + TA + Glu + age)^2
                 + (DA + Ser + OA + TA + Glu |hive)
                 + (DA + Ser + OA + TA + Glu |trial),
                 data = scaled_als, 
                 control = lmerControl(optimizer = "bobyqa",
                             optCtrl = list(maxfun = 1e6)))

socS_mod@optinfo$conv # check for optimizer convergence

socS_resid <- simulateResiduals(socS_mod, n=1000) # simulate residuals
testDispersion(socS_resid) # check for over-/underdispersion
testOutliers(socS_resid) # check outliers
plot(socS_resid) # plot diagnostics

indS_mod <- lmer(indS ~ (DA + Ser + OA + TA + Glu + age)^2
                 + (DA + Ser + OA + TA + Glu|hive)
                 + (DA + Ser + OA + TA + Glu|trial),
                 data = scaled_als,
                 control = lmerControl(optimizer = "bobyqa",
                                       optCtrl = list(maxfun = 1e6)))

indS_mod@optinfo$conv # check for optimizer convergence

indS_resid <- simulateResiduals(indS_mod, n=1000)
testDispersion(indS_resid)
testOutliers(indS_resid)
plot(indS_resid)

socM_mod <- lmer(socM ~ (DA + Ser + OA + TA + Glu + age)^2
                 + (DA + Ser + OA + TA + Glu |hive)
                 + (DA + Ser + OA + TA + Glu |trial),
                 data = scaled_als,
                 control = lmerControl(optimizer = "bobyqa",
                                       optCtrl = list(maxfun = 1e6)))

socM_mod@optinfo$conv # check for optimizer convergence

socM_resid <- simulateResiduals(socM_mod, n=1000)
testDispersion(socM_resid)
testOutliers(socM_resid)
plot(socM_resid)

indM_mod <- lmer(indM ~ (DA + Ser + OA + TA + Glu + age)^2
                 + (DA + Ser + OA + TA + Glu|hive)
                 + (DA + Ser + OA + TA + Glu|trial),
                 data = scaled_als,
                 control = lmerControl(optimizer = "bobyqa",
                                       optCtrl = list(maxfun = 1e6)))

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
                  data = scaled_als,
                  control = lmerControl(optimizer = "bobyqa",
                                        optCtrl = list(maxfun = 1e6)))

socS_mod2@optinfo$conv # check for optimizer convergence

indS_mod2 <- lmer(indS ~ DA + Ser + OA + TA + Glu + age
                  + (DA + Ser + OA + TA + Glu|hive)
                  + (DA + Ser + OA + TA + Glu|trial),
                  data = scaled_als,
                  control = lmerControl(optimizer = "bobyqa",
                                        optCtrl = list(maxfun = 1e6)))

indS_mod2@optinfo$conv # check for optimizer convergence

socM_mod2 <- lmer(socM ~ DA + Ser + OA + TA + Glu + age
                  + (DA + Ser + OA + TA + Glu|hive)
                  + (DA + Ser + OA + TA + Glu|trial),
                  data = scaled_als,
                  control = lmerControl(optimizer = "bobyqa",
                                       optCtrl = list(maxfun = 1e6)))

socM_mod2@optinfo$conv # check for optimizer convergence

indM_mod2 <- lmer(indM ~ DA + Ser + OA + TA + Glu + age
                  + (DA + Ser + OA + TA + Glu|hive)
                  + (DA + Ser + OA + TA + Glu|trial),
                  data = scaled_als,
                  control = lmerControl(optimizer = "bobyqa",
                                        optCtrl = list(maxfun = 1e6)))

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
                 "Glu" = "Glu", 
                 "DA:Ser" = "DA \u00D7 5-HT",
                 "DA:OA" = "DA \u00D7 OA", 
                 "DA:Glu" = "DA \u00D7 Glu", 
                 "Ser:OA" = "5-HT \u00D7 OA", 
                 "Ser:Glu" = "5-HT \u00D7 Glu", 
                 "OA:Glu" = "OA \u00D7 Glu"))

# custom ggplot theme
lmm_plot_theme <- theme(panel.background=element_blank(), 
                        panel.border=element_blank(),
                        axis.line=element_line(color=1, size=0.1, linetype=1),
                        panel.grid.major=element_line(color="gray80", size=0.1, linetype=1), 
                        panel.grid.minor=element_line(color="gray80", size=0.1, linetype=1),
                        strip.background=element_blank(),
                        axis.text.x=element_text(color="black", size=6),
                        axis.text.y=element_text(color="black", size=6),
                        axis.title.x=element_text(size=8),
                        axis.title.y=element_text(size=8),
                        axis.ticks=element_blank(),
                        plot.title=element_text(hjust=0.5, size=8),
                        legend.position="top",
                        legend.text=element_text(size=8))

# socS LMM term estimates and significance
socS_plot <- ggplot(socS_est, aes(x=estimate, y=term, color=group, label=p.stars)) +
  geom_point(size=1.25) + ylab("") + xlab("Estimate") + 
  geom_text(vjust=0, nudge_y=-0.05, size=4, color=1) +
  geom_errorbar(aes(xmin=conf.low, xmax=conf.high, width=0.3)) +
  geom_vline(xintercept=0, linetype="dashed", color = "black", size=0.3) + 
  scale_color_npg(guide="none") +
  scale_y_discrete(limits = names(y_label), labels = y_label) +
  ggtitle("Social Learning Speed") +
  lmm_plot_theme +
  theme(axis.text.y = element_text(face = c("bold", "plain", "bold", rep('plain', 7)),
                       size = c(8, 6, 8, rep(6, 7))))

# indS LMM term estimates and significance
indS_plot <- ggplot(indS_est, aes(x=estimate, y=term, color=group, label=p.stars)) +
  geom_point(size=1.25) + ylab("") + xlab("Estimate") + 
  geom_text(vjust=0, nudge_y=-0.05, size=4, color=1) +
  geom_errorbar(aes(xmin=conf.low, xmax=conf.high, width=0.3)) +
  geom_vline(xintercept=0, linetype="dashed", color = "black", size=0.3) + 
  scale_color_npg(guide="none") +
  scale_y_discrete(limits = names(y_label), labels = y_label) +
  ggtitle("Individual Learning Speed") +
  lmm_plot_theme +
  theme(axis.text.y = element_text(face = c("bold", "plain", "bold", rep('plain', 7)),
                                   size = c(8, 6, 8, rep(6, 7))))

# socM LMM term estimates and significance
socM_plot <- ggplot(socM_est, aes(x=estimate, y=term, color=group, label=p.stars)) +
  geom_point(size=1.25) + ylab("") + xlab("Estimate") + 
  geom_text(vjust=0, nudge_y=-0.05, size=4, color=1) +
  geom_errorbar(aes(xmin=conf.low, xmax=conf.high, width=0.3)) +
  geom_vline(xintercept=0, linetype="dashed", color = "black", size=0.3) + 
  scale_color_npg(guide="none") +
  scale_y_discrete(limits = names(y_label), labels = y_label) +
  ggtitle("Social Learning Accuracy") +
  lmm_plot_theme +
  theme(axis.text.y = element_text(face = c("bold", "plain", "bold", rep('plain', 7)),
                                   size = c(8, 6, 8, rep(6, 7))))

# indM LMM term estimates and significance
indM_plot <- ggplot(indM_est, aes(x=estimate, y=term, color=group, label=p.stars)) +
  geom_point(size=1.25) + ylab("") + xlab("Estimate") + 
  geom_text(vjust=0, nudge_y=-0.05, size=4, color=1) +
  geom_errorbar(aes(xmin=conf.low, xmax=conf.high, width=0.3)) +
  geom_vline(xintercept=0, linetype="dashed", color = "black", size=0.3) + 
  scale_color_npg(guide="none") +
  scale_y_discrete(limits = names(y_label), labels = y_label) +
  ggtitle("Individual Learning Accuracy") +
  lmm_plot_theme +
  theme(axis.text.y = element_text(face = c("bold", "plain", "bold", rep('plain', 7)),
                                   size = c(8, 6, 8, rep(6, 7))))

png("lmm_estimates.png", width=2500, height=2000, res=300)
# arranging aa plots in one figure and removing x-axis title
comb_plots <- ggarrange(indS_plot + rremove("ylab") + rremove("xlab"), 
                        indM_plot + rremove("y.text") + rremove("xlab"), 
                        socS_plot  + rremove("ylab") + rremove("xlab"), 
                        socM_plot  + rremove("y.text") + rremove("xlab"), 
                        labels="AUTO", ncol=2, nrow=2, font.label=list(size=10),
                        widths = c(1.1, 1, 1, 1))
# adding a common x-axis title
annotate_figure(comb_plots, bottom=text_grob("Estimate", hjust=0.5, size=8))
dev.off()
