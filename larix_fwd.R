# Code for deriving density (G) and size (MSD) parameters for fine woody debris
# (FWD, < 7 cm in diameter) of Larix cajanderi that can be used to calculate
# FWD biomass using the line-intersect method

## setting the working directory in which the excel files from
## Delcourt and Veraverbeke (2022) were saved
wdir <- ""

## loading required packages
library(tidyverse)
library(car)
library(emmeans)
library(multcompView)
library(openxlsx)

# Specific gravity (G) ----------------------------------------------------

## loading raw data from the water displacement method
sp.gr <- read.xlsx(paste0(wdir,"YA2019_fwd_specific_gravity.xlsx"),
                   sheet = "specific_gravity")
glimpse(sp.gr)
head(sp.gr)
sp.gr$size_class <- as.factor(sp.gr$size_class)
levels(sp.gr$size_class)

# density of the paraffin wax used to seal the surface of the
# woody pieces before immersion in water for volume determination
dens.paraffin <- 0.9

# constant equal to 1 when the mass is in grams and the volume is in
# cubic centimeters
k <- 1

## retrieving volume and specific gravity of each piece of wood using Eq. (6) and
# Eq. (7) from our paper
sp.gr <- sp.gr %>%
  mutate(diff_im = ((m_post_im - m_coated) / m_coated) * 100,
         m_paraffin = m_coated - m_0,
         v_0 = m_w_disp - (m_paraffin / dens.paraffin),
         G = (m_0 * k) / v_0) %>%
  filter(diff_im > -1 & diff_im < 1)
# the pieces whose difference in mass before and after immersion
# in water is +/- 1% were discarded 
glimpse(sp.gr)

## computing mean and standard deviation of specific gravity per
## diameter size class (last column of Table 3)
df.g <- sp.gr %>%
  group_by(size_class) %>%
  summarise(mean = mean(G, na.rm = TRUE),
                   sd = sd(G, na.rm = TRUE),
                   count = n()) %>%
  rename(g_mean = mean, g_sd = sd, g_count = count)

## examining significant differences in G between size classes
boxplot(G ~ size_class, data = sp.gr,
        xlab = "Diameter size class",
        ylab = "Specific gravity (g/cm3)")

# one-way ANOVA
res.aov.g <- aov(G ~ size_class, data = sp.gr)
summary(res.aov.g)

## checking ANOVA assumptions
# homogeneity of variances
plot(res.aov.g, 1)
leveneTest(G ~ size_class, data = sp.gr)
# normality
plot(res.aov.g, 2)
aov_residuals.g <- residuals(object = res.aov.g)
shapiro.test(x = aov_residuals.g)

# Tukey-Kramer test
TukeyHSD(res.aov.g)
multcompLetters(TukeyHSD(res.aov.g)$size_class[, 'p adj'])

# Mean squared diameter (MSD) ---------------------------------------------

## creating a function to compute mean squared diameter (cm2)
msd <- function(x){
  y <- sum(x^2)/length(x)
  return(y)
}

## loading field measurements of FWD diameters
m.sq.d <- read.xlsx(paste0(wdir,"YA2019_fwd_diameters.xlsx"),
                    sheet = "diameters")
glimpse(m.sq.d)
head(m.sq.d)

## retrieving diameter perpendicular to the transect line (d.lim) at the
## intersection point using diameter along the line (d_field) and the angle
## between the transect line and the piece of wood (angle) as measured in the field

## assigning of a size class to each piece using d.lim and the size class
## definition from McRae et al. (1979)

m.sq.d <- m.sq.d %>%
  mutate(d.lim = d_field * sin((angle * pi) / 180),
         size_class = ifelse(d.lim < 0.5, 1,
                             ifelse(d.lim >= 0.5 & d.lim < 1, 2,
                                    ifelse(d.lim >= 1 & d.lim < 3, 3,
                                           ifelse(d.lim >= 3 & d.lim < 5, 4,
                                                  ifelse(d.lim >=5 & d.lim < 7, 5, NA)))))) %>%
  filter(species=="LC" & d.lim > 0 & !is.na(size_class))
print(m.sq.d[196,])
m.sq.d[196,7] <- 3 # correcting one mistake from the size classification above
m.sq.d$size_class <- as.factor(m.sq.d$size_class)
glimpse(m.sq.d)
head(m.sq.d)

## computing size-class-specific values of MSD for each stand
## MSD values were divided by the square of the arithmetic
## class center (ACC) of the corresponding size class
df.msd <- m.sq.d %>%
  group_by(plotID, size_class) %>%
  summarise(msd.val = msd(d.lim)) %>%
  mutate(ACC = ifelse(size_class==1, 0.25,
                      ifelse(size_class==2, 0.75,
                             ifelse(size_class==3, 2,
                                    ifelse(size_class==4, 4,
                                           ifelse(size_class==5, 6, NA))))),
         msd.acc = msd.val / (ACC^2),
         val.log = log(msd.acc))
glimpse(df.msd)

## computing MSD per diameter size class (2nd column of Table 3)
df.msd.2 <- m.sq.d %>%
  group_by(size_class) %>%
  summarise(msd.val = msd(d.lim),
            msd.count = n())

## examining the effect of size classes on MSD relative to class midpoint
boxplot(msd.acc ~ size_class, data = df.msd,
        xlab = "Diameter size class",
        ylab = expression(paste("MSD/", "ACC"^"2")))

# one-way ANOVA
res.aov.msd <- aov(msd.acc ~ size_class, data = df.msd)
summary(res.aov.msd)

## checking ANOVA assumptions
# homogeneity of variances
plot(res.aov.msd, 1)
leveneTest(msd.acc ~ size_class, data = df.msd)
# normality
plot(res.aov.msd, 2)
aov_residuals.msd <- residuals(object = res.aov.msd)
shapiro.test(x = aov_residuals.msd)

# Kruskal-Wallis test and Wilcoxon rank sum test for pairwise
# comparisons with Benjamini-Hochberg corrections
kruskal.test(msd.acc ~ size_class, data = df.msd)
pairwise.wilcox.test(df.msd$msd.acc, df.msd$size_class,
                     p.adjust.method = "BH")

## saving Table 3 in the outputs folder
df.fwd <- merge(df.msd.2, df.g, by = "size_class", all = TRUE)
write.csv(mutate_if(df.fwd, is.numeric, ~ round(.,3)),
          paste0(wdir,"outputs/table3_fwd_summary.csv"),
          row.names = FALSE)

