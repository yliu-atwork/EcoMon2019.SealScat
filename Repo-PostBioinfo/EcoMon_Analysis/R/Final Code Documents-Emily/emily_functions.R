# functions for eDNA volume analysis
library(ggplot2)
library(dplyr)
library(reshape2)
library(doBy)
library(ggpubr)

# Create a box plot for any two variables within a data set. 
box_plot <- function(data, x, y) { 
  p <- ggplot(data, aes({{x}}, {{y}})) +
    geom_boxplot() +
    geom_jitter(shape=16, width = 0.1, height = 0.1) + theme_bw() +
    coord_flip()
  
  return(p)
}

# make wide data frame that calculates differences between desired variable between filtration volumes
paired_dif <- function(df, variable){
  dat.123.wide <- dcast(df, Cast + Sampling.Depth.Type + Sampling.Depth.Meter ~ Filtration.Volume, value.var = variable, fun.aggregate=mean)
  dat.123.wide$diff1.2 <- dat.123.wide$`2` - dat.123.wide$`1`
  dat.123.wide$diff2.3 <- dat.123.wide$`3` - dat.123.wide$`2`
  
  return(dat.123.wide)
  
}

# create scatter plot
scatter_plot <- function(df, x, y){
  p <- ggplot(df, aes({{x}}, {{y}})) +
    geom_point() +
    geom_smooth(method = lm) +
    stat_cor(method = "pearson")
  return(p)
}

# create violin plot for species distribution in water column
species_plot <- function(data, x, y) { 
  p <- ggplot(data, aes({{x}}, {{y}})) +
    geom_violin() +
    geom_jitter(shape=16, width = 0.1, height = 0.1) + 
    theme_bw() +
    scale_y_reverse()
  
  return(p)
}

# Create a summary table of count, mean, standard deviation, and 95% confidence intervals based on data
summary_table <- function(data, group, value) {
  data %>%
  group_by({{group}}) %>%
  dplyr::summarize(
    count = n(),
    mean = mean({{value}}, na.rm = TRUE),
    sd = sd({{value}}, na.rm = TRUE),
    min = min({{value}}, na.rm = TRUE),
    max = max({{value}}, na.rm = TRUE)) %>%
  mutate(se = sd / sqrt(count),
         lower.ci = mean - qt(1 - (0.05 / 2), count - 1) * se,
         upper.ci = mean + qt(1 - (0.05 / 2), count - 1) * se)

}

# Create a violin plot for any two variables within a data set. 
violin_plot <- function(data, x, y) { 
  p <- ggplot(data, aes({{x}}, {{y}})) +
    geom_violin() +
    geom_jitter(shape=16, width = 0.1, height = 0.1) + theme_bw() +
    coord_flip()
  
  return(p)
}





