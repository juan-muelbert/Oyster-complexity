library(tidyverse)
library(visreg)
library(lme4)
library(lmerTest)
library(car)
library(scales)
library(glmmTMB)
library(MuMIn)
library(pscl)

## The natural architecture of oyster reefs maximises recruit survival. Esquivel-Muelbert et al. 2025 ##

# Data loading and setup ====

complex_data <- read.csv("data_complexity_21oct2023.csv")
metrics <- read.csv("tile_metrics.csv")
metrics$logSurface_area <- log10(metrics$surface_area_cm2)
data <- complex_data %>% left_join(metrics) %>% filter(treatment != "control") %>% na.omit()
data$density <- data$oysters/data$surface_area_cm2
data$site<-as.factor(data$site)
data$treatment <- as.factor(data$treatment)
reefs <- read.csv("reef_metrics.csv")
reefs$logHeight <- log10(reefs$height)
reefs$logRG <- log10(reefs$rugosity)

caged <- data %>% filter(treatment == "caged")
uncaged <- data %>% filter(treatment == "uncaged")

# Main ====

### Complexity mediates predation on oyster recruits
## Surface area
# Modelling the effect of surface area on oyster counts, comparing treatment effects, using site as random effects

sa.count.mod <- glmmTMB(oysters ~ poly(logSurface_area, 2) * treatment + (1|site),
                        family = nbinom2(link='log'),
                        data = data)
summary(sa.count.mod)
summary.sa.count.mod <- summary(sa.count.mod)
write.csv(summary.sa.count.mod$coefficients$cond, "SA_count_summary.csv")
hist(residuals(sa.count.mod))
ranef(sa.count.mod)
r.squaredGLMM(sa.count.mod)
Anova(sa.count.mod)

# Fig. 2a

surface_area.new <- seq(from = min(data$logSurface_area),
                        to = max(data$logSurface_area),
                        length.out = 200)

predic.sa <- tibble(logSurface_area = surface_area.new) %>%
  expand_grid(tibble(treatment = c("caged", "uncaged"))) %>%
  expand_grid(tibble(site = c("brisbane_waters", "port_hacking", "hawkesbury")))

sa.count.predict <- predict(object = sa.count.mod,
                            newdata = predic.sa,
                            se.fit = TRUE,
                            type = "response") %>%
  as_tibble()

plot_data <- bind_cols(predic.sa, sa.count.predict) %>%
  mutate(confidence95 = se.fit*1.96)

plot_data.summary <- plot_data %>% group_by(logSurface_area, treatment) %>%
  summarise(fit = mean(fit), confidence95 = mean(confidence95))

sa.plot <- ggplot(data = plot_data.summary, aes(x = logSurface_area, y = fit, colour = treatment)) +
  geom_ribbon(aes(ymin = fit - confidence95,
                  ymax = fit + confidence95), alpha = 0.3, linetype = 0) +
  geom_point(data = data, aes(x = logSurface_area, y = oysters), alpha = 0.3, size = 2) +
  geom_line(linewidth = 1) + scale_colour_manual(values = c("#c54c00ff", "#5F5F5F")) +
  scale_fill_manual(values = alpha(c("#c54c00ff", "#5F5F5F"), 0.5)) + theme_classic() + ylab("Oyster count") +
  xlab("log10 surface area")

sa.plot

## Modelling the effect of Fractal dimension and Height Range on oyster oyster densities, comparing treatment effects, using site as random effects
## Fractal dimension

fd.mod <- lmer(sqrt(density) ~ poly(fd, 2) * treatment + (1|site),
               data = data)
summary(fd.mod)
summary.fd.mod <- summary(fd.mod)
Anova(fd.mod, test = "F")
hist(residuals(fd.mod))

# Fig 2b

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
  scale_fill_manual(values = alpha(c("#c54c00ff", "#5F5F5F"), 0.5)) + 
  ylab("Oyster density (per cm2)") + xlab("Fractal dimension") + theme(legend.position = "none")

fd.plot

## Height range
height.mod <- lmer(sqrt(density) ~ poly(logHeight, 2) * treatment + (1|site),
                   data = data)
summary(height.mod)
summary.height.mod <- summary(height.mod)
Anova(height.mod, test = "F")
hist(residuals(height.mod))

# Fig 2c

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
  ylab("Oyster density (per cm2)") + xlab("log Height range (cm)") + theme(legend.position = "none")

Height.plot

### Ecosystem engineering in oyster reefs

# Modelling the interactive effects of Fractal dimension and Height range on oyster density in uncaged artifical habitat designs

fd.height_int <- lmer(sqrt(density) ~ poly(fd, 2) * poly(logHeight, 2) + (1|site), data = uncaged)
vif(fd.height_int)
summary.fd_height <- summary(fd.height_int)
summary.fd_height
Anova(fd.height_int, test = "F")
hist(residuals(fd.height_int))

# Plot model (heatmap) and overlay data from natural oyster reefs
# summarizing oyster reef data
oyster_reefs <- reefs %>% group_by(patch) %>% # summarising reef points by patch
  summarise(meanFD = mean(fd),
            seFD   = sd(fd)/sqrt(n()),
            meanlogHeight = mean(logHeight),
            seHeight = sd(logHeight)/sqrt(n()),
            sdFD = sd(fd), sdHeight = sd(logHeight,),
            meanlogRG = mean(logRG),
            seRG = sd(logRG)/sqrt(n()))

# Fig3
# Generate raw fig with visreg2d
heatmap <- visreg2d(fd.height_int, xvar = "fd", yvar = "logHeight", plot.type ="gg", scale = "response",
                    data = uncaged)
heatmap
heatmap.df <- heatmap$data # extract heatmap data for plotting using ggplot 
heatmap.df$z[heatmap.df$z < 0] <- NA

# Plot heatmap in ggplot

inc <- seq(0, max(sqrt(heatmap.df$z), na.rm=TRUE), length=20)^2 # create ribbons of 5% increments.

ggplot(data = heatmap.df, aes(x = x, y = y, fill = z)) +
  geom_raster() +
  geom_contour(data = heatmap.df, aes(z=z), color="black", linewidth = 0.3, breaks = inc) +
  geom_point(data = oyster_reefs, aes (x = meanFD, y = meanlogHeight), size = 1.5, inherit.aes = FALSE) + 
  geom_errorbar(data = oyster_reefs, aes(x = meanFD, y = meanlogHeight,
                                         ymin = meanlogHeight - sdHeight,
                                         ymax = meanlogHeight + sdHeight), width=0.01, inherit.aes = FALSE) + 
  geom_errorbar(data = oyster_reefs, aes(x = meanFD, y = meanlogHeight,
                                         xmin = meanFD - sdFD,
                                         xmax = meanFD + sdFD), width=0.01, inherit.aes = FALSE) +
  scale_fill_viridis_c(option = "B", na.value="white", begin = 0.05) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "grey50")) +
  xlab("Fractal dimension") + ylab("log Height range") + labs(fill = "Oysters per cm2 (sqrt)") +
  scale_y_continuous(breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, 1)) + 
  theme(text = element_text(size = 7))

# Suplementary materials ====
### Model testing for caging artifacts

control_data <- read.csv("control.csv")
control.mod <- zeroinfl(oysters ~ treatment, data = control_data, dist = "negbin")
summary_control.mod <- summary(control.mod) # no difference between cage control and uncaged.
summary_control.mod ## Supplementary Table 6
