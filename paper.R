## Complexity script ## Juan Muelbert et al. 2024

# Set up ====
# ..... Libraries ====

library(tidyverse)
library(visreg)
library(gridExtra)

# Main ====

# ..... Load and modify data ====
complex_data <- read.csv("data_complexity_21oct2023.csv")
metrics <- read.csv("tile_metrics.csv")
data <- complex_data %>% left_join(metrics) %>% filter(treatment != "control") %>% na.omit()
data$density <- data$oysters/data$surface_area
data$logHeight <- log1p(data$height)
data$site<-as.factor(data$site)
data$sqrt.den <- sqrt(data$density)
data$treatment <- as.factor(data$treatment)
reefs <- read.csv("reef_metrics.csv")
reefs$logHeight <- log1p(reefs$height)

data.bw <- data %>% filter(site == "brisbane_waters")
data.hr <- data %>% filter(site == "hawkesbury")
data.ph <- data %>% filter(site == "port_hacking")
caged <- data %>% filter(treatment == "caged")
uncaged <- data %>% filter(treatment == "uncaged")


# ..... Modelling ====

sites <- uncaged$site %>% unique()

# .......... Effects of logRG on oyster counts and density ==== 
# caged treatment only

area.count <- lm(sqrt(oysters) ~ poly(logRG, 2) * site,
                 data = uncaged)
sum_area.count <- summary(area.count)
sum_area.count
anova(area.count)
write.csv(sum_area.count$coefficients, "area_count.csv")
visreg::visreg(area.count, "logRG", by = "site", overlay=T)
plot(area.count$residuals, area.count$fitted.values)
hist(area.count$residuals)



area.count.sites <- lapply(1:length(sites), function(x) {
  site_x <- sites[x]
  
  lm(sqrt(oysters) ~ poly(logRG, 2),
     data = uncaged %>% filter(site == site_x))
  
})

names(area.count.sites) <- sites
                 

lapply(area.count.sites, summary)




area.density <- lm(sqrt(density) ~ poly(logRG, 2) * site,
                   data = uncaged)
sum_area.density <- summary(area.density)
sum_area.density
anova(area.density)
write.csv(sum_area.count$coefficients, "area_density.csv")
visreg::visreg(area.density, "logRG", by = "site", overlay=T)
plot(area.density$residuals, area.density$fitted.values)
hist(area.density$residuals)

summary(area.density)



area.density.sites <- lapply(1:length(sites), function(x) {
  site_x <- sites[x]
  
  lm(sqrt(density) ~ poly(logRG, 2),
     data = uncaged %>% filter(site == site_x))
  
})

names(area.density.sites) <- sites


lapply(area.density.sites, summary)



ggplot(data = uncaged,
       aes(y = sqrt(oysters), 
           x = logRG,
           colour = site,
           fill = site,
           group = site)) +
  
  geom_jitter(width = 0.01, alpha = 0.8) +
  
  geom_smooth(method = 'lm', formula = y ~poly(x,2),
              size = 2, se = F, colour = 'white') +
  geom_smooth(method = 'lm', formula = y ~poly(x,2),
              size = 1, se = T) +
  
  theme_minimal()






# Plotting
par(mfrow=c(1,2))
visreg::visreg(area.count, "logRG", by = "site", overlay=T, ylim = c(0, 12)) # Effects of Rugosity on count data
visreg::visreg(area.density, "logRG", by = "site", overlay=T, ylim = c(0, 0.4)) # Area-independent effects on

# Models for single variables (logRG, fd and logHeight) per site, comparing treatment effects
# Hawkesbury
hr.rg.cg <- lm(sqrt(density) ~ poly(logRG, 2),
               data = data.hr %>% filter(treatment == 'caged'))
summary(hr.rg.cg)

hr.rg.uc <- lm(sqrt(density) ~ poly(logRG, 2),
               data = data.hr %>% filter(treatment == 'uncaged'))
summary(hr.rg.uc)


hr.rg.uc_vs_cg <- lm(sqrt(density) ~ poly(logRG, 2) * treatment,
            data = data.hr)
summary(hr.rg.uc_vs_cg)
anova(hr.rg.uc_vs_cg)
visreg::visreg(hr.rg.uc_vs_cg, "logRG", by = "treatment", overlay=T)
plot(hr.rg.uc_vs_cg$residuals, hr.rg.uc_vs_cg$fitted.values)
hist(hr.rg.uc_vs_cg$residuals)



hr.fd.cg <- lm(sqrt(density) ~ poly(fd, 2),
               data = data.hr %>% filter(treatment == 'caged'))
summary(hr.fd.cg)

hr.fd.uc <- lm(sqrt(density) ~ poly(fd, 2),
               data = data.hr %>% filter(treatment == 'uncaged'))
summary(hr.fd.uc)


hr.fd.uc_vs_cg <- lm(sqrt(density) ~ poly(fd, 2) * treatment,
            data = data.hr)
summary(hr.fd)
visreg::visreg(hr.fd, "fd", by = "treatment", overlay=T)
plot(hr.fd$residuals, hr.fd$fitted.values)
hist(hr.fd$residuals)



hr.height.cg <- lm(sqrt(density) ~ poly(logHeight, 2),
               data = data.hr %>% filter(treatment == 'caged'))
summary(hr.height.cg)

hr.height.uc <- lm(sqrt(density) ~ poly(logHeight, 2),
               data = data.hr %>% filter(treatment == 'uncaged'))
summary(hr.height.uc)


hr.height<- lm(sqrt(density) ~ poly(logHeight, 2) * treatment,
               data = data.hr)
summary(hr.height)
visreg::visreg(hr.height, "logHeight", by = "treatment", overlay=T)
plot(hr.height$residuals, hr.height$fitted.values)
hist(hr.height$residuals)

#Port Hacking

ph.rg.cg <- lm(sqrt(density) ~ poly(logRG, 2),
               data = data.ph %>% filter(treatment == 'caged'))
summary(ph.rg.cg)

ph.rg.uc <- lm(sqrt(density) ~ poly(logRG, 2),
               data = data.ph %>% filter(treatment == 'uncaged'))
summary(ph.rg.uc)

ph.rg <- lm(sqrt(density) ~ poly(logRG, 2) * treatment,
            data = data.ph)
summary(ph.rg)
anova(ph.rg)
visreg::visreg(ph.rg, "logRG", by = "treatment", overlay=T)
plot(ph.rg$residuals, ph.rg$fitted.values)
hist(ph.rg$residuals)


ph.fd.cg <- lm(sqrt(density) ~ poly(fd, 2),
               data = data.ph %>% filter(treatment == 'caged'))
summary(ph.fd.cg)

ph.fd.uc <- lm(sqrt(density) ~ poly(fd, 2),
               data = data.ph %>% filter(treatment == 'uncaged'))
summary(ph.fd.uc)

ph.fd <- lm(sqrt(density) ~ poly(fd, 2) * treatment,
            data = data.ph)
summary(ph.fd)
visreg::visreg(ph.fd, "fd", by = "treatment", overlay=T)
plot(ph.fd$residuals, ph.fd$fitted.values)
hist(ph.fd$residuals)



ph.height.cg <- lm(sqrt(density) ~ poly(logHeight, 2),
                   data = data.ph %>% filter(treatment == 'caged'))
summary(ph.height.cg)

ph.height.uc <- lm(sqrt(density) ~ poly(logHeight, 2),
                   data = data.ph %>% filter(treatment == 'uncaged'))
summary(ph.height.uc)


ph.height<- lm(sqrt(density) ~ poly(logHeight, 2) * treatment,
               data = data.ph)
summary(ph.height)
visreg::visreg(ph.height, "logHeight", by = "treatment", overlay=T)
plot(ph.height$residuals, ph.height$fitted.values)
hist(ph.height$residuals)

#Brisbane water


bw.rg.cg <- lm(sqrt(density) ~ poly(logRG, 2),
               data = data.bw %>% filter(treatment == 'caged'))
summary(bw.rg.cg)

bw.rg.uc <- lm(sqrt(density) ~ poly(logRG, 2),
               data = data.bw %>% filter(treatment == 'uncaged'))
summary(bw.rg.uc)


bw.rg <- lm(sqrt(density) ~ poly(logRG, 2) * treatment,
            data = data.bw)
summary(bw.rg)
visreg::visreg(bw.rg, "logRG", by = "treatment", overlay=T)
plot(bw.rg$residuals, bw.rg$fitted.values)
hist(bw.rg$residuals)



bw.fd.cg <- lm(sqrt(density) ~ poly(fd, 2),
               data = data.bw %>% filter(treatment == 'caged'))
summary(bw.fd.cg)

bw.fd.uc <- lm(sqrt(density) ~ poly(fd, 2),
               data = data.bw %>% filter(treatment == 'uncaged'))
summary(bw.fd.uc)


bw.fd <- lm(sqrt(density) ~ poly(fd, 2) * treatment,
            data = data.bw)
summary(bw.fd)
visreg::visreg(bw.fd, "fd", by = "treatment", overlay=T)
plot(bw.fd$residuals, bw.fd$fitted.values)
hist(bw.fd$residuals)



bw.height.cg <- lm(sqrt(density) ~ poly(logHeight, 2),
                   data = data.bw %>% filter(treatment == 'caged'))
summary(bw.height.cg)

bw.height.uc <- lm(sqrt(density) ~ poly(logHeight, 2),
                   data = data.bw %>% filter(treatment == 'uncaged'))
summary(bw.height.uc)


bw.height<- lm(sqrt(density) ~ poly(logHeight, 2) * treatment,
               data = data.bw)
summary(bw.height)
visreg::visreg(bw.height, "logHeight", by = "treatment", overlay=T)
plot(bw.height$residuals, bw.height$fitted.values)
hist(bw.height$residuals)


# Plotting into a panel with 9 plots

a <- visreg::visreg(hr.rg, "logRG", by = "treatment", overlay=T, gg=TRUE, 
                    points.par=list(size=1.5, alpha=0.6)) + #HR RG
  scale_colour_manual(values = c("#B50000", "#5F5F5F")) + ylim(-0.035, 0.55) +
  scale_fill_manual(values = alpha(c("#E87777", "#BDBCBC"), 0.5)) +
  theme_classic() + theme(legend.position = "none") 

b <- visreg::visreg(ph.rg, "logRG", by = "treatment", overlay=T, gg=TRUE,
                    points.par=list(size=1.5, alpha=0.6)) + # PH RG
  scale_colour_manual(values = c("#B50000", "#5F5F5F")) + ylim(-0.035, 0.55) +
  scale_fill_manual(values = alpha(c("#E87777", "#BDBCBC"), 0.5)) +
  theme_classic() + theme(legend.position = "none") 

c <- visreg::visreg(bw.rg, "logRG", by = "treatment", overlay=T, gg = TRUE,
                    points.par=list(size=1.5, alpha=0.6)) + # BW RG
  scale_colour_manual(values = c("#B50000", "#5F5F5F")) + ylim(-0.035, 0.55) +
  scale_fill_manual(values = alpha(c("#E87777", "#BDBCBC"), 0.5)) +
  theme_classic() + theme(legend.position = "none") 

d <- visreg::visreg(hr.fd, "fd", by = "treatment", overlay=T, gg=TRUE, #HR FD
                    points.par=list(size=1.5, alpha=0.6)) + 
  scale_colour_manual(values = c("#B50000", "#5F5F5F")) + ylim(-0.035, 0.55) +
  scale_fill_manual(values = alpha(c("#E87777", "#BDBCBC"), 0.5)) +
  theme_classic() + theme(legend.position = "none")

e <- visreg::visreg(ph.fd, "fd", by = "treatment", overlay=T, gg=TRUE, #PH FD
                    points.par=list(size=1.5, alpha=0.6)) + 
  scale_colour_manual(values = c("#B50000", "#5F5F5F")) + ylim(-0.035, 0.55) +
  scale_fill_manual(values = alpha(c("#E87777", "#BDBCBC"), 0.5)) +
  theme_classic() + theme(legend.position = "none")

f <-  visreg::visreg(bw.fd, "fd", by = "treatment", overlay=T, gg=TRUE, #BW FD
                     points.par=list(size=1.5, alpha=0.6)) + 
  scale_colour_manual(values = c("#B50000", "#5F5F5F")) + ylim(-0.035, 0.55) +
  scale_fill_manual(values = alpha(c("#E87777", "#BDBCBC"), 0.5)) +
  theme_classic() + theme(legend.position = "none") 

g <- visreg::visreg(hr.height, "logHeight", by = "treatment", overlay=T, gg=TRUE, #HR Height
                    points.par=list(size=1.5, alpha=0.6)) +
  scale_colour_manual(values = c("#B50000", "#5F5F5F")) + ylim(-0.035, 0.55) +
  scale_fill_manual(values = alpha(c("#E87777", "#BDBCBC"), 0.5)) +
  theme_classic() + theme(legend.position = "none")  

h <- visreg::visreg(ph.height, "logHeight", by = "treatment", overlay=T, gg=TRUE, # PH Height
                    points.par=list(size=1.5, alpha=0.6)) +
  scale_colour_manual(values = c("#B50000", "#5F5F5F")) + ylim(-0.035, 0.55) +
  scale_fill_manual(values = alpha(c("#E87777", "#BDBCBC"), 0.5)) +
  theme_classic() + theme(legend.position = "none")  

i <- visreg::visreg(bw.height, "logHeight", by = "treatment", overlay=T, gg=TRUE, # BW Height
                    points.par=list(size=1.5, alpha=0.6)) +
  scale_colour_manual(values = c("#B50000", "#5F5F5F")) + ylim(-0.035, 0.55) +
  scale_fill_manual(values = alpha(c("#E87777", "#BDBCBC"), 0.5)) +
  theme_classic() + theme(legend.position = "none")  

grid.arrange(a, b ,c, d, e, f, g, h, i, ncol=3, nrow=3)

## Model Selection
# Models for each site #
#Rugosity

rg.hr <- lm(sqrt(density) ~ poly(logRG, 2), data = data.hr %>% filter(treatment=="uncaged"))
anova(rg.hr) 
summary(rg.hr)
saveRDS(rg.hr, "models/rg_hr.rds")
AIC(rg.hr)

rg.ph <- lm(sqrt(density) ~ poly(logRG, 2), data = data.ph %>% filter(treatment=="uncaged"))
anova(rg.ph) 
summary(rg.ph)
saveRDS(rg.ph, "models/rg_ph.rds")
AIC(rg.ph)

rg.bw <- lm(sqrt(density) ~ poly(logHeight, 2) * poly(fd, 2), data = data.bw %>% filter(treatment=="uncaged"))
anova(rg.bw) 
summary(rg.bw)
saveRDS(rg.bw, "models/rg_bw.rds")
AIC(rg.bw)

# ∆Height

height.hr <- lm(sqrt(density) ~ poly(logHeight, 2), data = data.hr %>% filter(treatment=="uncaged"))
anova(height.hr) 
summary(height.hr)
saveRDS(height.hr, "models/height_hr.rds")
AIC(height.hr)

height.ph <- lm(sqrt(density) ~ poly(logHeight, 2), data = data.ph %>% filter(treatment=="uncaged"))
anova(height.ph) 
summary(height.ph)
saveRDS(height.ph, "models/height_ph.rds")
AIC(height.ph)

height.bw <- lm(sqrt(density) ~ poly(logHeight, 2), data = data.bw %>% filter(treatment=="uncaged"))
anova(height.bw) 
summary(height.bw)
saveRDS(height.bw, "models/height_bw.rds")
AIC(height.bw)

# FD

fd.hr <- lm(sqrt(density) ~ poly(fd, 2), data = data.hr %>% filter(treatment=="uncaged"))
anova(fd.hr) 
summary(fd.hr)
saveRDS(fd.hr, "models/fd_hr.rds")
AIC(fd.hr)

fd.ph <- lm(sqrt(density) ~ poly(fd, 2), data = data.ph %>% filter(treatment=="uncaged"))
anova(fd.ph) 
summary(fd.ph)
saveRDS(fd.ph, "models/fd_ph.rds")
AIC(fd.ph)

fd.bw <- lm(sqrt(density) ~ poly(fd, 2), data = data.bw %>% filter(treatment=="uncaged"))
anova(fd.bw) 
summary(fd.bw)
saveRDS(fd.bw, "models/fd_bw.rds")
AIC(fd.bw)

# Height x FD

fd_h.hr <- lm(sqrt(density) ~ poly(fd, 2) * poly(logHeight, 2), data = data.hr %>% filter(treatment=="uncaged"))
anova(fd_h.hr) 
summary(fd_h.hr)
saveRDS(fd_h.hr, "models/fd_height_hr.rds")
AIC(fd_h.hr)

fd_h.ph <- lm(sqrt(density) ~ poly(fd, 2) * poly(logHeight, 2), data = data.ph %>% filter(treatment=="uncaged"))
anova(fd_h.ph) 
summary(fd_h.ph)
saveRDS(fd_h.ph, "models/fd_height_ph.rds")
AIC(fd_h.ph)

fd_h.bw <- lm(sqrt(density) ~ poly(fd, 2) * poly(logHeight, 2), data = data.bw %>% filter(treatment=="uncaged"))
anova(fd_h.bw) 
summary(fd_h.bw)
saveRDS(fd_h.bw, "models/fd_height_bw.rds")
AIC(fd_h.bw)

# Height + FD

fdp_h.hr <- lm(sqrt(density) ~ poly(fd, 2) + poly(logHeight, 2), data = data.hr %>% filter(treatment=="uncaged"))
anova(fdp_h.hr) 
summary(fdp_h.hr)
saveRDS(fd_h.hr, "models/fd_plus_height_hr.rds")
AIC(fdp_h.hr)

fdp_h.ph <- lm(sqrt(density) ~ poly(fd, 2) + poly(logHeight, 2), data = data.ph %>% filter(treatment=="uncaged"))
anova(fdp_h.ph) 
summary(fdp_h.ph)
saveRDS(fdp_h.ph, "models/fd_plus_height_ph.rds")
AIC(fdp_h.ph)

fdp_h.bw <- lm(sqrt(density) ~ poly(fd, 2) + poly(logHeight, 2), data = data.bw %>% filter(treatment=="uncaged"))
anova(fdp_h.bw) 
summary(fdp_h.bw)
saveRDS(fdp_h.bw, "models/fd_plus_height_bw.rds")
AIC(fdp_h.bw)

## plotting Height and FD effects per site
# summarizing oyster reef data
oyster_reefs <- reefs %>% group_by(patch) %>% # summarising reef points by patch
  summarise(meanFD = mean(fd),
            seFD   = sd(fd)/sqrt(n()),
            meanlogHeight = mean(logHeight),
            seHeight = sd(logHeight)/sqrt(n()),
            sdFD = sd(fd), sdHeight = sd(logHeight))

# heatmaops per site (FD x HEIGHT)
heat_hr <- visreg2d(fd_h.hr, xvar = "fd", yvar = "logHeight", plot.type ="gg",
                 data = uncaged) +
  geom_contour(aes(z=z), color="grey30", linewidth = 0.2) +
  geom_point(data = oyster_reefs, aes (x = meanFD, y = meanlogHeight)) + 
  geom_errorbar(data = oyster_reefs, aes(x = meanFD, y = meanlogHeight, ymin = meanlogHeight - sdHeight, ymax = meanlogHeight + sdHeight, xmin = meanFD - seFD, xmax = meanFD + seFD), width=0.01) + 
  geom_errorbar(data = oyster_reefs, aes(x = meanFD, y = meanlogHeight, xmin = meanFD - sdFD, xmax = meanFD + sdFD), width=0.01) +
  scale_fill_continuous(type = "viridis", limits=c(0, 0.32), breaks=seq(0, 0.3,by=0.1)) ## change color and set similar scale

heat_hr

heat_ph <- visreg2d(fd_h.ph, xvar = "fd", yvar = "logHeight", plot.type ="gg",
                    data = uncaged) +
  geom_contour(aes(z=z), color="grey30", linewidth = 0.2) +
  geom_point(data = oyster_reefs, aes (x = meanFD, y = meanlogHeight)) + 
  geom_errorbar(data = oyster_reefs, aes(x = meanFD, y = meanlogHeight, ymin = meanlogHeight - sdHeight, ymax = meanlogHeight + sdHeight, xmin = meanFD - seFD, xmax = meanFD + seFD), width=0.01) + 
  geom_errorbar(data = oyster_reefs, aes(x = meanFD, y = meanlogHeight, xmin = meanFD - sdFD, xmax = meanFD + sdFD), width=0.01) + 
  scale_fill_continuous(type = "viridis", limits=c(0, 0.32), breaks=seq(0, 0.3,by=0.1)) ## change color and set similar scale

heat_ph

heat_bw <- visreg2d(fd_h.bw, xvar = "fd", yvar = "logHeight", plot.type ="gg",
                    data = uncaged) +
  geom_contour(aes(z=z), color="grey30", linewidth = 0.2) +
  geom_point(data = oyster_reefs, aes (x = meanFD, y = meanlogHeight)) + 
  geom_errorbar(data = oyster_reefs, aes(x = meanFD, y = meanlogHeight, ymin = meanlogHeight - sdHeight, ymax = meanlogHeight + sdHeight, xmin = meanFD - seFD, xmax = meanFD + seFD), width=0.01) + 
  geom_errorbar(data = oyster_reefs, aes(x = meanFD, y = meanlogHeight, xmin = meanFD - sdFD, xmax = meanFD + sdFD), width=0.01) +
  scale_fill_continuous(type = "viridis", limits=c(0, 0.32), breaks=seq(0, 0.3,by=0.1)) ## change color and set similar scale

grid.arrange(heat_hr, heat_ph ,heat_bw, ncol=3, nrow=1)


## Analysis of caging artifacts (cage controls x uncaged treatments) - cage control are cages with windowns that allow predator access

control_data <- read.csv("control.csv")
control.mod <- lm(sqrt(oysters) ~ treatment, data = control_data)
summary(control.mod) # no difference between cage control
hist(control.mod$residuals)
plot(control.mod$fitted.values, control.mod$residuals)



# 3D plots

library(plotly)

metrics_and_oyster_reefs <- bind_rows(metrics, oyster_reefs)

plot <- plot_ly(metrics, x = ~height, y = ~fd, z = ~rugosity,
                type = "scatter3d", mode = "markers",
                marker = list(size = 12, color = ~rugosity, colorscale = "Viridis"))

# THIS NEED RUGOSITY IN THE OYSTER_REEFS DATA TO WORK
# plot <- plot %>% add_trace(data = oyster_reefs %>%
#                              mutate(height = meanlogHeight,
#                                     fd = meanFD,
#                                     logRG = meanlogHeight)
#                            , colour = 'black')

# Show plot
plot
