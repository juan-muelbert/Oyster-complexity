library(tidyverse)
library(visreg)
library(gridExtra)
library(lme4)
library(lmerTest)
library(car)
library(scales)

## Three-dimensional niche construction maximizes offspring survival. Esquivel-Muelbert et al. 2025 ##

complex_data <- read.csv("data_complexity_21oct2023.csv")
metrics <- read.csv("tile_metrics.csv")
data <- complex_data %>% left_join(metrics) %>% filter(treatment != "control") %>% na.omit()
data$density <- data$oysters/data$surface_area
data$site<-as.factor(data$site)
data$treatment <- as.factor(data$treatment)
data$sqrt.den <- sqrt(data$density)
reefs <- read.csv("reef_metrics.csv")
reefs$logHeight <- log10(reefs$height)


caged <- data %>% filter(treatment == "caged")
uncaged <- data %>% filter(treatment == "uncaged")

### Complexity mediates predation on oyster recruits
# Models for single variables (logRG, fd and logHeight) comparing treatment effects, using site as random effects

## Rugosity
rg.mod <- lmer(sqrt(density) ~ poly(logRG, 2) * treatment + (1|site),
               data = data)
summary(rg.mod)
summary.rg.mod <- summary(rg.mod)
write.csv(summary.rg.mod$coefficients, "RG_summary.csv")
hist(residuals(rg.mod))
ranef(rg.mod)
r.squaredGLMM(rg.mod)

# Fig2a

rugosity.new <- seq(from = min(data$logRG),
                    to = max(data$logRG),
                    length.out = 200)
predic.rug <- tibble(logRG = rugosity.new) %>%
  expand_grid(tibble(treatment = c("caged", "uncaged"))) %>%
  expand_grid(tibble(site = c("brisbane_waters", "port_hacking", "hawkesbury")))

rg.predict <- predict(object = rg.mod,
                      newdata = predic.rug,
                      se.fit = TRUE,
                      type = "response") %>%
  as_tibble()

plot_data <- bind_cols(predic.rug, rg.predict) %>%
  mutate(confidence95 = se.fit*1.96)

plot_data.summary <- plot_data %>% group_by(logRG, treatment) %>%
  summarise(fit = mean(fit), confidence95 = mean(confidence95))

rg.plot <- ggplot(data = plot_data.summary, aes(x = logRG, y = fit, colour = treatment)) +
  geom_ribbon(aes(ymin = fit - confidence95,
                  ymax = fit + confidence95), alpha = 0.3, linetype = 0) +
  geom_point(data = data, aes(x = logRG, y = sqrt(density)), alpha = 0.3, size = 2) +
  geom_line(linewidth = 1) + scale_colour_manual(values = c("#c54c00ff", "#5F5F5F")) +
  scale_fill_manual(values = alpha(c("#c54c00ff", "#5F5F5F"), 0.5)) + theme_classic() + ylab("Oyster density (per cm2)") +
  xlab("log Rugosity")

rg.plot

## Fractal dimension

fd.mod <- lmer(sqrt(density) ~ poly(fd, 2) * treatment + (1|site),
               data = data)
summary(fd.mod)
summary.fd.mod <- summary(fd.mod)
write.csv(summary.fd.mod$coefficients, "FD_summary.csv")
hist(residuals(fd.mod))
ranef(fd.mod)
r.squaredGLMM(fd.mod)

# Fig2b

fd.new <- seq(from = min(data$fd),
              to = max(data$fd),
              length.out = 200)
predic.fd <- tibble(fd = fd.new) %>%
  expand_grid(tibble(treatment = c("caged", "uncaged"))) %>%
  expand_grid(tibble(site = c("brisbane_waters", "port_hacking", "hawkesbury")))

fd.predict <- predict(object = fd.mod,
                      newdata = predic.fd,
                      se.fit = TRUE,
                      type = "response") %>%
  as_tibble()

fd.plot_data <- bind_cols(predic.fd, fd.predict) %>%
  mutate(confidence95 = se.fit*1.96)

fd.plot_data.summary <- fd.plot_data %>% group_by(fd, treatment) %>%
  summarise(fit = mean(fit), confidence95 = mean(confidence95))

fd.plot <- ggplot(data = fd.plot_data.summary, aes(x = fd, y = fit, colour = treatment)) +
  geom_ribbon(aes(ymin = fit - confidence95,
                  ymax = fit + confidence95), alpha = 0.3, linetype = 0) +
  geom_point(data = data, aes(x = fd, y = sqrt(density)), alpha = 0.3, size = 2) +
  geom_line(linewidth = 1) +
  theme_classic() + scale_colour_manual(values = c("#c54c00ff", "#5F5F5F")) +
  scale_fill_manual(values = alpha(c("#c54c00ff", "#5F5F5F"), 0.5)) + ylab("Oyster density (per cm2)") + xlab("Fractal dimension")

fd.plot

## Height
height.mod <- lmer(sqrt(density) ~ poly(logHeight, 2) * treatment + (1|site),
                   data = data)
summary(height.mod)
summary.height.mod <- summary(height.mod)
write.csv(summary.height.mod$coefficients, "Height_summary.csv")
hist(residuals(height.mod))
r.squaredGLMM(height.mod)


# Height plot

Height.new <- seq(from = min(data$logHeight),
                  to = max(data$logHeight),
                  length.out = 200)
predic.Height <- data_frame(logHeight = Height.new) %>%
  expand_grid(tibble(treatment = c("caged", "uncaged"))) %>%
  expand_grid(tibble(site = c("brisbane_waters", "port_hacking", "hawkesbury")))

Height.predict <- predict(object = height.mod,
                          newdata = predic.Height,
                          se.fit = TRUE,
                          type = "response") %>%
  as_tibble()

Height.plot_data <- bind_cols(predic.Height, Height.predict) %>%
  mutate(confidence95 = se.fit*1.96)

Height.plot_data.summary <- Height.plot_data %>% group_by(logHeight, treatment) %>%
  summarise(fit = mean(fit), confidence95 = mean(confidence95))

Height.plot <- ggplot(data = Height.plot_data.summary, aes(x = logHeight, y = fit, colour = treatment)) +
  geom_ribbon(aes(ymin = fit - confidence95,
                  ymax = fit + confidence95), linetype = 0, alpha = 0.3) +
  geom_point(data = data, aes(x = logHeight, y = sqrt(density)), alpha = 0.3, size = 2) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = c("#c54c00ff", "#5F5F5F")) +
  scale_fill_manual(values = alpha(c("#c54c00ff", "#5F5F5F"), 0.5)) + theme_classic() + scale_x_continuous(breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1)) +
  ylab("Oyster density (per cm2)") + xlab("log Height range (cm)")
Height.plot

# Plot Fig 2 Raw

grid.arrange(rg.plot, fd.plot, Height.plot, ncol=3, nrow=1)

## Three-dimensional niche construction on oyster reefs

# Modelling the interactive effects of Fractal dimension and Height range on oyster density

fd.height_int <- lmer(sqrt(density) ~ poly(fd, 2) * poly(logHeight, 2) + (1|site), data = uncaged)
vif(fd.height_int)
summary.fd_height <- summary(fd.height_int)
summary.fd_height
hist(residuals(fd.height_int))
plot(fd.height_int)
write.csv(summary.fd_height$coefficients, "fd_height_full_mod.csv")


# Plot model (heatmap) and overlay data from natural oyster reefs
# summarizing oyster reef data
oyster_reefs <- reefs %>% group_by(patch) %>% # summarising reef points by patch
  summarise(meanFD = mean(fd),
            seFD   = sd(fd)/sqrt(n()),
            meanlogHeight = mean(logHeight),
            seHeight = sd(logHeight)/sqrt(n()),
            sdFD = sd(fd), sdHeight = sd(logHeight))

# Fig3
heatmap <- visreg2d(fd.height_int, xvar = "fd", yvar = "logHeight", plot.type ="gg", scale = "response",
                 data = uncaged)

heatmap.df <- heatmap$data
max(heatmap.df$z) ## max predicted sqrt(density) = 0.1369969

heatmap.df$scaled_z <- rescale(heatmap.df$z, 
                               to = c(0, 0.1369969), 
                               from = range(heatmap.df$z,
                               na.rm = TRUE, finite = TRUE)) # scaling model from 0 to max predicted density(sqrt) to remove negative values 

## heatmap

ggplot(data = heatmap.df, aes(x = x, y = y, fill = scaled_z)) +
  geom_raster() +
  geom_contour(data = heatmap.df, aes(z=scaled_z), color="grey10", linewidth = 0.3, binwidth = 0.025) +
  geom_point(data = oyster_reefs, aes (x = meanFD, y = meanlogHeight), size = 2, inherit.aes = FALSE) + 
  geom_errorbar(data = oyster_reefs, aes(x = meanFD, y = meanlogHeight,
                                         ymin = meanlogHeight - sdHeight,
                                         ymax = meanlogHeight + sdHeight), width=0.01, inherit.aes = FALSE ) + 
  geom_errorbar(data = oyster_reefs, aes(x = meanFD, y = meanlogHeight,
                                         xmin = meanFD - sdFD,
                                         xmax = meanFD + sdFD), width=0.01, inherit.aes = FALSE) +
  scale_fill_viridis_c(option = "B") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "grey50")) +
  xlab("Fractal dimension") + ylab("log Height range") + labs(fill = "Oysters per cm2 (sqrt)") +
  scale_y_continuous(breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, 1))

### For suplementary materials ##

# Model with and without surface area

### Effects of logRG on oyster counts and density for uncaged and caged treatments

## Count model (effects of surface area included)
# Uncaged
count.mod_uncaged <- lmer(sqrt(oysters) ~ poly(logRG, 2) + (1|site), data = uncaged)
sum_count.mod_uncaged <- summary(count.mod_uncaged)
sum_count.mod_uncaged
write.csv(sum_count.mod_uncaged$coefficients, "area_count_uncaged.csv")
hist(residuals(area.count))
r.squaredGLMM(count.mod_uncaged)

# Caged

count.mod_caged <- lmer(sqrt(oysters) ~ poly(logRG, 2) + (1|site), data = caged)
sum_count.mod_caged <- summary(count.mod_caged)
sum_count.mod_caged
write.csv(sum_count.mod_caged$coefficients, "area_count_caged.csv")
hist(residuals(area.count))
r.squaredGLMM(count.mod_caged)

# Density Model (effects of area excluded)
# Uncaged
density.mod_uncaged <- lmer(sqrt(density) ~ poly(logRG, 2) + (1|site), data = uncaged)
sum_density.mod_uncaged <- summary(density.mod_uncaged)
sum_density.mod_uncaged
write.csv(sum_density.mod_uncaged$coefficients, "area_density_uncaged.csv")
hist(residuals(area.density))
r.squaredGLMM(density.mod_uncaged)

# Caged

density.mod_caged <- lmer(sqrt(density) ~ poly(logRG, 2) + (1|site), data = caged)
sum_density.mod_caged <- summary(density.mod_caged)
sum_density.mod_caged
write.csv(sum_density.mod_caged$coefficients, "area_density_caged.csv")
hist(residuals(area.density))
r.squaredGLMM(density.mod_caged)

# Figure S1

rugosity.new <- seq(from = min(data$logRG),
                    to = max(data$logRG),
                    length.out = 200)
predic.rug_uncaged <- tibble(logRG = rugosity.new) %>%
  expand_grid(tibble(treatment = c("uncaged"))) %>%
  expand_grid(tibble(site = c("brisbane_waters", "port_hacking", "hawkesbury")))

predic.rug_caged <- tibble(logRG = rugosity.new) %>%
  expand_grid(tibble(treatment = c("caged"))) %>%
  expand_grid(tibble(site = c("brisbane_waters", "port_hacking", "hawkesbury")))

# Fig S1a - Count data uncaged

count.predict_uncaged <- predict(object = count.mod_uncaged,
                           newdata = predic.rug_uncaged,
                           se.fit = TRUE,
                           type = "response") %>%
  as_tibble()

plot_count_data_uncaged <- bind_cols(predic.rug_uncaged, count.predict_uncaged) %>%
  mutate(confidence95 = se.fit*1.96)

plot_count_data.summary_uncaged <- plot_count_data_uncaged %>% group_by(logRG, treatment) %>%
  summarise(fit = mean(fit), confidence95 = mean(confidence95))

count.plot_uncaged <- ggplot(data = plot_count_data.summary_uncaged, aes(x = logRG, y = fit, colour = treatment)) +
  geom_ribbon(aes(ymin = fit - confidence95,
                  ymax = fit + confidence95), linetype = 0, alpha = 0.3) +
  geom_point(data = uncaged, aes(x = logRG, y = sqrt(oysters)), alpha = 0.3, size = 1.5) +
  geom_line(linewidth = 0.8) + scale_colour_manual(values = c("#5F5F5F")) +
  scale_fill_manual(values = alpha(c("#5F5F5F"), 0.5)) + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylab("Oyster count") + theme(legend.position = "none") + ylim(0, 14)

count.plot_uncaged

# FigS1b
count.predict_caged <- predict(object = count.mod_caged,
                                 newdata = predic.rug_caged,
                                 se.fit = TRUE,
                                 type = "response") %>%
  as_tibble()

plot_count_data_caged <- bind_cols(predic.rug_caged, count.predict_caged) %>%
  mutate(confidence95 = se.fit*1.96)

plot_count_data.summary_caged <- plot_count_data_caged %>% group_by(logRG, treatment) %>%
  summarise(fit = mean(fit), confidence95 = mean(confidence95))

count.plot_caged <- ggplot(data = plot_count_data.summary_caged, aes(x = logRG, y = fit, colour = treatment)) +
  geom_ribbon(aes(ymin = fit - confidence95,
                  ymax = fit + confidence95), linetype = 0, alpha = 0.3) +
  geom_point(data = caged, aes(x = logRG, y = sqrt(oysters)), alpha = 0.3, size = 1.5) +
  geom_line(linewidth = 0.8) + scale_colour_manual(values = c("#c54c00ff")) +
  scale_fill_manual(values = alpha(c("#c54c00ff"), 0.5)) + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylab("Oyster count") + theme(legend.position = "none") + ylim(0, 14)

count.plot_caged

### Fig S1c
density.predict_uncaged <- predict(object = density.mod_uncaged,
                                 newdata = predic.rug_uncaged,
                                 se.fit = TRUE,
                                 type = "response") %>%
  as_tibble()

plot_density_data_uncaged <- bind_cols(predic.rug_uncaged, density.predict_uncaged) %>%
  mutate(confidence95 = se.fit*1.96)

plot_density_data.summary_uncaged <- plot_density_data_uncaged %>% group_by(logRG, treatment) %>%
  summarise(fit = mean(fit), confidence95 = mean(confidence95))

density.plot_uncaged <- ggplot(data = plot_density_data.summary_uncaged, aes(x = logRG, y = fit, colour = treatment)) +
  geom_ribbon(aes(ymin = fit - confidence95,
                  ymax = fit + confidence95), linetype = 0, alpha = 0.3) +
  geom_point(data = uncaged, aes(x = logRG, y = sqrt(density)), alpha = 0.3, size = 1.5) +
  geom_line(linewidth = 0.8) + scale_colour_manual(values = c("#5F5F5F")) +
  scale_fill_manual(values = alpha(c("#5F5F5F"), 0.5)) + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylab("Oyster density") + theme(legend.position = "none") + ylim(0, 0.55)

density.plot_uncaged

# FigS1d
density.predict_caged <- predict(object = density.mod_caged,
                               newdata = predic.rug_caged,
                               se.fit = TRUE,
                               type = "response") %>%
  as_tibble()

plot_density_data_caged <- bind_cols(predic.rug_caged, density.predict_caged) %>%
  mutate(confidence95 = se.fit*1.96)

plot_density_data.summary_caged <- plot_density_data_caged %>% group_by(logRG, treatment) %>%
  summarise(fit = mean(fit), confidence95 = mean(confidence95))

density.plot_caged <- ggplot(data = plot_density_data.summary_caged, aes(x = logRG, y = fit, colour = treatment)) +
  geom_ribbon(aes(ymin = fit - confidence95,
                  ymax = fit + confidence95), linetype = 0, alpha = 0.3) +
  geom_point(data = caged, aes(x = logRG, y = sqrt(density)), alpha = 0.3, size = 1.5) +
  geom_line(linewidth = 0.8) + scale_colour_manual(values = c("#c54c00ff")) +
  scale_fill_manual(values = alpha(c("#c54c00ff"), 0.5)) + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylab("Oyster density") + theme(legend.position = "none") + ylim(0, 0.55)

density.plot_caged

# Create Fig S1. Plot panels a, b, c and d together

grid.arrange(count.plot_uncaged, count.plot_caged, density.plot_uncaged, density.plot_caged, ncol=2, nrow=2)


# Relationship between Rugosity, Height range and Fractal dimension

mod_metrics <- lm(logRG ~ fd * logHeight, data = metrics)
summary_mod_metrics <- summary(mod_metrics)
summary_mod_metrics
write.csv(summary_mod_metrics$coefficients, "mod_metrics.csv")
visreg2d(mod_metrics, "fd", "logHeight", plot.type = "persp")

# Model testing for caging artifacts
library(mvabund)
library(pscl)

control_data <- read.csv("control.csv")
control.mod <- zeroinfl(oysters ~ treatment, data = control_data, dist = "negbin")
summary_control.mod <- summary(control.mod) # no difference between cage control
summary_control
Anova(control.mod)
hist(control.mod$residuals)
write.csv(summary_control.mod$coefficients, "control_mod.csv")


control_mod <- manyglm(oysters ~ treatment, data = control_data, family = "poisson")
summary_control_mod <- summary(control_mod)
hist(control_mod$residuals)
plot.manyglm(control_mod)
write.csv(summary_control_mod$coefficients, "summary_control.csv")
