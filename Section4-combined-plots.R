#
#  Supporting material for "Replication in random translation designs for model-robust prediction"  
# by Waite, T.W. (2023+)

# Script to produce plots in the paper from saved results 

library(tidyverse)
results_df <- read_csv("heuristics-results.csv")

# make names etc nicer (for plotting with ggplot)
results_df <- tibble(results_df)
names(results_df ) <- c("tau2", "Model", "Heuristic", "Optimum delta", "minPsi", "MIV", "Max MISB", "Relative sd of discrepancy")
results_df %>% 
  mutate(biasratio = `Max MISB`/MIV) %>% 
  ggplot(aes(x=tau2, y=biasratio, colour=Model, linetype=Heuristic)) + 
  geom_line() +
  xlab(expression(tau^2/sigma^2)) +
  ylab("Ratio of maximum MISB to MIV")

# tau2=0.2 seems like a sensible rough cut-off point below which bias and variance are of a comparable magnitude, above that bias begins to dominate the IMSPE

#  make Heuristic and Model factors, order levels to improve plots 
results_df  %>% mutate(Heuristic= factor(Heuristic, levels=c("Split", "Half split", "Retain")),
                       Model = factor(Model, 
                                      levels=c("Linear (1 factor)", 
                                               "Quadratic (1 factor)", 
                                               "Linear (2 factor)", 
                                               "Quadratic (3 factor)")))  -> results_df

# compute efficiencies 
results_df %>% 
  dplyr::select(tau2, Model, Heuristic, minPsi) %>%
  pivot_wider(names_from=Heuristic, values_from = minPsi) %>% # put different heuristics in different columns 
  transmute(tau2=tau2,
            Model=Model,
            `Half split`  = Retain/`Half split`,
            Split=Retain/Split) %>% # take ratio of Psi of retain vs splitting (or half splitting)
  pivot_longer(cols=c(`Half split`,
                      Split), names_to = "Heuristic", values_to = "Efficiency") %>% # turn data back into long format 
  mutate(Heuristic=factor(Heuristic, levels=c("Split", "Half split"))) %>% 
  filter(!((Heuristic=="Half split")&(Model!="Linear (1 factor)"))) -> efficiency_results


# Headline results 
# summarize efficiency results for different models 
efficiency_results %>% filter(tau2<=0.2) %>% 
  group_by(Model, Heuristic) %>% summarize( EfficiencyMin = min(Efficiency), EfficiencyMax = max(Efficiency))

# efficiency plots - split vs retain
pdf("../figures/efficiencies.pdf", width=5 ,height=3)
efficiency_results %>% 
  filter(Heuristic!="Half split"&tau2<=0.2) %>%
  ggplot(aes(x=tau2, y=Efficiency, colour=Model, shape=Model)) + 
  geom_line() + 
  geom_point() +
  theme_bw() + 
  xlab(expression(tau^2)) +
  scale_colour_grey()
dev.off()

# Bias-variance tradeoffs

base_graph <- results_df %>% 
  filter(Model!="Quadratic (3 factor)"&Heuristic!="Half split"&tau2<=0.2) %>% 
  ggplot(aes(x=tau2, linetype=Heuristic)) + facet_wrap(~Model) + theme_bw() + scale_colour_grey()  + xlab(expression(tau^2))

pdf("../figures/problems1to3-optdelta.pdf", width=8, height=3)
base_graph + geom_line(aes(y=`Optimum delta`))
dev.off()

pdf("../figures/problems1to3-MIV.pdf", width=8, height=3)
base_graph + geom_line(aes(y=MIV))   + ylim(c(0,1.5))
dev.off()

pdf("../figures/problems1to3-MISB.pdf", width=8, height=3)
base_graph + geom_line(aes(y=`Max MISB`))    + ylim(c(0,1.5))
dev.off()


# 3 factor quadratic - plot of optimum deltas, MIV, max MISB
base_graph <- results_df %>% filter(Model=="Quadratic (3 factor)"&tau2<=0.05) %>%
  ggplot(aes(x=tau2, linetype=Heuristic, shape=Heuristic)) + theme_bw() + scale_colour_grey() + xlab(expression(tau^2))

pdf("../figures/problem4-optdelta.pdf", width=5, height=3)
base_graph + geom_line(aes(y=`Optimum delta`)) + geom_point(aes(y=`Optimum delta`))
dev.off()

pdf("../figures/problem4-MIV.pdf", width=5, height=3)
base_graph + geom_line(aes(y=MIV)) 
dev.off()

pdf("../figures/problem4-MISB.pdf", width=5, height=3)
base_graph + geom_line(aes(y=`Max MISB`)) 
dev.off()

