###############################################################################
# Importing libraries
###############################################################################

library(DHARMa) # for diagnostics of LMMs
library(lme4) # for building LMMs
library(sjPlot) # for extracting estimates, conf. intervals, and p-values from LMMs
library(ggpubr) # for creating multi-panel figures
library(ggplot2) # for plots
library(ggsci) # for npg palette
library(car) # for Type II Wald chi-square tests

###############################################################################
# Data manipulation
###############################################################################

# import data
als <- read.csv("bioamine_learn_data.csv")

# re-assign factors 
als$hive <- as.factor(als$hive)
als$trial <- as.factor(als$trial)
als$learnSpd <- as.factor(als$learnSpd)
als$learnAcc <- as.factor(als$learnAcc)

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
                         "trial" = als$trial,
                         "learnSpd" = als$learnSpd,
                         "learnAcc" = als$learnAcc)

###############################################################################
# Assessing learnSpd-learnAcc trade-off with Fisher's and plotting proportions
###############################################################################

# proportions
prop.table(table(als$learnAcc, als$learnSpd), margin=2)

# Fisher's exact test
fisher.test(als$learnAcc, als$learnSpd)

# custom annotations for ggplot
## significance of Fisher's exact test  
sig <- grobTree(textGrob("**", x=0,  y=0.95, hjust=-8.8,
                         gp=gpar(col="black", fontsize=14)))
# proportion of accurate individual learners that are also fast individual learners
indv_per <- grobTree(textGrob("33.33%", x=0.28,  y=0.77,
                              gp=gpar(col="black", fontsize=12)))
# proportion of accurate social learners that are also fast social learners
soc_per <- grobTree(textGrob("15.38%", x=0.73,  y=0.11,
                             gp=gpar(col="black", fontsize=12)))

# plotting and saving the barplot
png("manuscript/Prop_learnSpd_indivSocLearnSpeed.png", width=1500, height=1000, res=300)
ggplot(als, aes(x = as.factor(learnSpd), fill = as.factor(learnAcc))) + 
  labs(y="Proportions", x="Learning speed", fill="Learning accuracy") +
  ylim(c(0, 1.05)) +
  scale_x_discrete(labels=c("Individual", "Social")) + 
  scale_fill_npg(labels=c("Individual", "Social")) +
  geom_bar(position = "fill") + 
  annotation_custom(sig) +
  annotation_custom(indv_per) +
  annotation_custom(soc_per) +
  theme()
dev.off()

###############################################################################
# Building LMMs with interactions
###############################################################################

socS_mod <- lmer(socS ~ (DA + Ser + OA + TA + Glu + age)^2
                 + (DA + Ser + OA + TA + Glu |hive)
                 + (DA + Ser + OA + TA + Glu |trial),
                 data = scaled_als)

socS_resid <- simulateResiduals(socS_mod, n=1000) # simulate residuals
testDispersion(socS_resid) # check for over-/underdispersion
testOutliers(socS_resid) # check outliers
plot(socS_resid) # plot diagnostics

indS_mod <- lmer(indS ~ (DA + Ser + OA + TA + Glu + age)^2
                 + (DA + Ser + OA + TA + Glu|hive)
                 + (DA + Ser + OA + TA + Glu|trial),
                 data = scaled_als)

indS_resid <- simulateResiduals(indS_mod, n=1000)
testDispersion(indS_resid)
testOutliers(indS_resid)
plot(indS_resid)

socM_mod <- lmer(socM ~ (DA + Ser + OA + TA + Glu + age)^2
                 + (DA + Ser + OA + TA + Glu |hive)
                 + (DA + Ser + OA + TA + Glu |trial),
                 data = scaled_als)

socM_resid <- simulateResiduals(socM_mod, n=1000)
testDispersion(socM_resid)
testOutliers(socM_resid)
plot(socM_resid)

indM_mod <- lmer(indM ~ (DA + Ser + OA + TA + Glu + age)^2
                 + (DA + Ser + OA + TA + Glu|hive)
                 + (DA + Ser + OA + TA + Glu|trial),
                 data = scaled_als)

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

indS_mod2 <- lmer(indS ~ DA + Ser + OA + TA + Glu + age
                  + (DA + Ser + OA + TA + Glu|hive)
                  + (DA + Ser + OA + TA + Glu|trial),
                  data = scaled_als)

socM_mod2 <- lmer(socM ~ DA + Ser + OA + TA + Glu + age
                  + (DA + Ser + OA + TA + Glu|hive)
                  + (DA + Ser + OA + TA + Glu|trial),
                  data = scaled_als)

indM_mod2 <- lmer(indM ~ DA + Ser + OA + TA + Glu + age
                  + (DA + Ser + OA + TA + Glu|hive)
                  + (DA + Ser + OA + TA + Glu|trial),
                  data = scaled_als)

################################################################################
# Comparing LMMs with interactions with the ones without interactions using LRTs
################################################################################

anova(socS_mod, socS_mod2, test = "LRT")
anova(indS_mod, indS_mod2, test = "LRT")
anova(socM_mod, socM_mod2, test = "LRT")
anova(indM_mod, indM_mod2, test = "LRT")

################################################################################
# LMMs with interactions are sig. better -- getting estimates and p-values
################################################################################

# summary of models
summary(socS_mod) 
summary(indS_mod)
summary(socM_mod)
summary(indM_mod)

# Type II Wald test of significance of estimates
Anova(socS_mod) 
Anova(indS_mod)
Anova(socM_mod)
Anova(indM_mod)

################################################################################
# Plotting estimates from LMMs
################################################################################

# extracting estimates, confidence intervals, and p-values from LMMs
socS_est <- data.frame(plot_model(socS_mod)$data)
indS_est <- data.frame(plot_model(indS_mod)$data)
socM_est <- data.frame(plot_model(socM_mod)$data)
indM_est <- data.frame(plot_model(indM_mod)$data)

# term labels for y-axis in plots
est_labe <- rev(c("DA", "5-HT", "OA", "TA", "Glu", "age", "DA \u00D7 5-HT",
                  "DA \u00D7 OA", "DA \u00D7 TA", "DA \u00D7 Glu", "DA \u00D7 age", 
                  "5-HT \u00D7 OA", "5-HT \u00D7 TA", "5-HT \u00D7 Glu", 
                  "5-HT \u00D7 age", "OA \u00D7 TA", "OA \u00D7 Glu", "OA \u00D7 age",
                  "TA \u00D7 Glu", "TA \u00D7 age", "Glu \u00D7 age"))

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
  geom_point(size=0.7) + ylab("") + xlab("Estimate") + 
  geom_text(vjust=0, nudge_y=-0.1, size=2.5, color=1) +
  geom_errorbar(aes(xmin=conf.low, xmax=conf.high, width=0.3)) +
  geom_vline(xintercept=0, linetype="dashed", color = "black", size=0.3) + 
  scale_color_npg(guide="none") +
  scale_y_discrete(labels=est_labe) +
  xlim(-25, 30) +
  lmm_plot_theme +
  theme(axis.text.y = 
          element_text(face = c(rep('plain', 4), 'bold', rep('plain', 4), 'bold', rep('plain', 11)),
                       size = c(rep(6, 4), 8, rep(6, 4), 8, rep(6, 11))))

# indS LMM term estimates and significance
indS_plot <- ggplot(indS_est, aes(x=estimate, y=term, color=group, label=p.stars)) +
  geom_point(size=0.7) + ylab("") + xlab("Estimate") + 
  geom_text(vjust=0, nudge_y=-0.1, size=2.5, color=1) +
  geom_errorbar(aes(xmin=conf.low, xmax=conf.high, width=0.3)) +
  geom_vline(xintercept=0, linetype="dashed", color = "black", size=0.3) + 
  scale_color_npg(guide="none") +
  scale_y_discrete(labels=est_labe) +
  xlim(-25, 30) +
  lmm_plot_theme +
  theme(axis.text.y = 
          element_text(face = c(rep('plain', 4), 'bold', rep('plain', 4), 'bold', rep('plain', 11)),
                       size = c(rep(6, 4), 8, rep(6, 4), 8, rep(6, 11))))

# socM LMM term estimates and significance
socM_plot <- ggplot(socM_est, aes(x=estimate, y=term, color=group, label=p.stars)) +
  geom_point(size=0.7) + ylab("") + xlab("Estimate") + 
  geom_text(vjust=0, nudge_y=-0.1, size=2.5, color=1) +
  geom_errorbar(aes(xmin=conf.low, xmax=conf.high, width=0.3)) +
  geom_vline(xintercept=0, linetype="dashed", color = "black", size=0.3) + 
  scale_color_npg(guide="none") +
  scale_y_discrete(labels=est_labe) +
  xlim(-25, 30) +
  lmm_plot_theme +
  theme(axis.text.y = 
          element_text(face = c(rep('plain', 4), 'bold', rep('plain', 4), 'bold', rep('plain', 11)),
                       size = c(rep(6, 4), 8, rep(6, 4), 8, rep(6, 11))))

# indM LMM term estimates and significance
indM_plot <- ggplot(indM_est, aes(x=estimate, y=term, color=group, label=p.stars)) +
  geom_point(size=0.7) + ylab("") + xlab("Estimate") + 
  geom_text(vjust=0, nudge_y=-0.1, size=2.5, color=1) +
  geom_errorbar(aes(xmin=conf.low, xmax=conf.high, width=0.3)) +
  geom_vline(xintercept=0, linetype="dashed", color = "black", size=0.3) + 
  scale_color_npg(guide="none") +
  scale_y_discrete(labels=est_labe) +
  xlim(-25, 30) +
  lmm_plot_theme +
  theme(axis.text.y = 
          element_text(face = c(rep('plain', 4), 'bold', rep('plain', 4), 'bold', rep('plain', 11)),
                       size = c(rep(6, 4), 8, rep(6, 4), 8, rep(6, 11))))

# saving figure in png format
png("lmm_estimates.png", width=2000, height=2000, res=300)
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