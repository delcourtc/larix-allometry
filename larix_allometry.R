# Code for deriving the coefficients of site-specific allometric equations relating
# diameter at breast height (DBH) to stem, stem wood, stem bark, branches, foliage,
# and total aboveground biomass derived using linear and nonlinear regressions

## setting the working directory in which the excel file from Schepaschenko et al. (2017) is saved
wdir <- ''

## loading required packages
library(openxlsx)
library(tidyverse)
library(ggpubr)
library(broom)
library(car)

## inputs
schepa <- paste0(wdir,"Biomass_tree_DB.xlsx")

## loading data
tree.db <- read.xlsx(schepa, sheet = "Tree_db")

str(tree.db)
head(tree.db)

## renaming the tree density variable
colnames(tree.db)[which(startsWith(colnames(tree.db), "Tree.number"))] <- "Plot.density"
unique(tree.db$Species[startsWith(tree.db$Species,"Larix")])

## selecting only Cajander larch data and removing rows without biomass measurements
larix.db <- tree.db %>%
  filter(Species=="Larix cajanderi Mayr.") %>%
  filter_at(vars(Pst,Pbr,Pf,Pabo),all_vars(!is.na(.))) %>%
  select(ID, Tree.age, DBH, H.tree, Pst, Pbark, Pbr, Pf, Pabo,
         Latitude, Longitude, Plot.density, Reference, ID_Plot)

unique(larix.db$Reference)
unique(larix.db$ID_Plot) # ID_Plot is missing in some rows

unk.plt.dens <- unique(larix.db[is.na(larix.db$ID_Plot),]$Plot.density)
print(unk.plt.dens) # Plot.density is used to retrieve 3 plots without identifiers

## assigning a unique ID to each of these 3 plots
## computing stem wood biomass (Pbole)
larix.db <- larix.db %>%
  dplyr::mutate(Reference = ifelse(startsWith(Reference, "Moska"), "Moskalyuk", "Schepaschenko"),
                ID_Plot = ifelse(is.na(ID_Plot) &  Plot.density==unk.plt.dens[1], "U01",
                                 ifelse(is.na(ID_Plot) & Plot.density==unk.plt.dens[2], "U02",
                                        ifelse(is.na(ID_Plot) & Plot.density==unk.plt.dens[3], "U03",
                                               ID_Plot))),
                Pbole = Pst - Pbark)

str(larix.db)
rm(unk.plt.dens)


## graphical exploration

sites <- c("MG", "UY") #MG: Magadan Oblast; UY: Ust-Yansky district
lab.comp <- c("Pst", "Pbole", "Pbark", "Pbr", "Pf", "Pabo")
labels.comp <- c("Stem", "Stem wood", "Stem bark", "Branches", "Foliage", "Aboveground")

# Ust-Yansky site (Schepaschenko, 2015)

uy_ylims <- c(175,150,35,25,8.75,175)

g <- list()
for (i in seq(6)){
  
  g[[i]] <- ggplot(filter(larix.db, Reference=="Schepaschenko"),
                   aes_string(x = "DBH", y = lab.comp[i],
                              shape = "ID_Plot", color = "ID_Plot")) + 
    geom_point(size = 3.5)+
    coord_cartesian(clip="off") + 
    xlab("DBH (cm)") + 
    ylab(paste(labels.comp[i], "biomass (kg dry matter)")) +
    scale_x_continuous(limits = c(0,25), expand = c(0,0)) +
    scale_y_continuous(limits = c(0, uy_ylims[i]), expand = c(0,0)) + 
    scale_color_manual(name = "Ust-Yansky plots (Schepaschenko, 2015)",
                       values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a"),
                       breaks = c("2780", "2781", "2782", "2783"),
                       labels = c("Plot 2780", "Plot 2781", "Plot 2782", "Plot 2783")) +
    scale_shape_manual(name = "Ust-Yansky plots (Schepaschenko, 2015)",
                       values = c(16,17,15,3),
                       breaks = c("2780", "2781", "2782", "2783"),
                       labels = c("Plot 2780", "Plot 2781", "Plot 2782", "Plot 2783")) + 
    theme(panel.border = element_rect(colour = "black", fill=NA, size=0.6)) +
    theme(axis.text = element_text(colour = "black", size = 11),
          axis.title = element_text(colour = "black", size = 13)) +
    theme(legend.position = c(0.25,0.85),
          legend.text = element_text(colour = "black", size = 14),
          legend.title = element_text(colour = "black", size = 14)) +
    theme(aspect.ratio = 1)
}
ggarrange(plotlist = g, ncol = 3, nrow = 2,
          labels = c("a", "b", "c", "d", "e", "f"),
          common.legend = TRUE, legend = "bottom")

rm(g, i, uy_ylims)

# Magadan site (Moskalyuk, 2015)

mg_ylims <- c(1250,1125,175,150,25,1500)

p <- list()
for (i in seq(6)){
  
  p[[i]] <- ggplot(filter(larix.db, Reference=="Moskalyuk"),
                   aes_string(x = "DBH", y = lab.comp[i],
                              shape = "ID_Plot", color = "ID_Plot")) + 
    geom_point(size = 3.5)+
    coord_cartesian(clip="off") + 
    xlab("DBH (cm)") + 
    ylab(paste(labels.comp[i], "biomass (kg dry matter)")) +
    scale_x_continuous(limits = c(0,55), expand = c(0,0)) +
    scale_y_continuous(limits = c(0, mg_ylims[i]), expand = c(0,0)) + 
    scale_color_manual(name = "Magadan plots (Moskalyuk, 2015)",
                       values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a",
                                  "#66a61e", "#e6ab02", "#a6761d"),
                       breaks = c("1061", "1062", "1063", "7351",
                                  "U01", "U02", "U03"),
                       labels = c("Plot 1061", "Plot 1062", "Plot 1063", "Plot 7351",
                                  "Plot U01", "Plot U02", "Plot U03")) +
    scale_shape_manual(name = "Magadan plots (Moskalyuk, 2015)",
                       values = c(16,17,15,3,8,8,8),
                       breaks = c("1061", "1062", "1063", "7351",
                                  "U01", "U02", "U03"),
                       labels = c("Plot 1061", "Plot 1062", "Plot 1063", "Plot 7351",
                                  "Plot U01", "Plot U02", "Plot U03")) + 
    theme(panel.border = element_rect(colour = "black", fill=NA, size=0.6)) +
    theme(axis.text = element_text(colour = "black", size = 11),
          axis.title = element_text(colour = "black", size = 13)) +
    theme(legend.position = c(0.21,0.78),
          legend.text = element_text(colour = "black", size = 14),
          legend.title = element_text(colour = "black", size = 14)) +
    theme(aspect.ratio = 1)
}

ggarrange(plotlist = p, ncol = 3, nrow = 2,
          labels = c("a", "b", "c", "d", "e", "f"),
          common.legend = TRUE, legend = "bottom")

rm(p, i, mg_ylims)


## linear regression (ordinary least squares method after log-log transformation)

# plots 2780 and 2781 at the Ust-Yansky site (Schepaschenko, 2015) were removed from
# the analysis because the biomass-DBH relationships differed from those of plots 2782 and 2783.
# The latter two plots were sampled at the same location (similar lat and lon values).
larix.db <- larix.db %>%
  subset(Reference=="Moskalyuk" | ID_Plot %in% c("2782","2783")) %>%
  select(Pst, Pbole, Pbark, Pbr, Pf, Pabo, DBH,
         H.tree, Tree.age, Plot.density, Reference) %>%
  dplyr::mutate(Location = ifelse(Reference=="Moskalyuk", "MG", "UY")) %>%
  select(-Reference)

larix.db$Location = factor(larix.db$Location)
levels(larix.db$Location)

write.csv(larix.db, paste0(wdir,"outputs/Larch_tree_DB.csv"),
          row.names = FALSE)

## creating empty lists to be filled with linear regressions between biomass component
## and size parameters after log-log transformation by site and by biomass component

# size parameter: DBH
lr.D <- list()
lr.D[["MG"]] <- list()
lr.D[["UY"]] <- list()

# size parameter: DBH^2*H.tree
lr.D2H <- list()
lr.D2H[["MG"]] <- list()
lr.D2H[["UY"]] <- list()

for (site in sites){
  
  df <- subset(larix.db, Location == site)
  
    for (comp in lab.comp){
    
      dbh <- df[!is.na(df[comp]),][["DBH"]]
      h.tree <- df[!is.na(df[comp]),][["H.tree"]]
      comp.mass <- df[!is.na(df[comp]),][[comp]]
      
      lr.D[[site]][[comp]] <- lm(log(comp.mass)~log(dbh))
      lr.D2H[[site]][[comp]] <- lm(log(comp.mass)~I(log((dbh^2)*h.tree)))
  }
}
rm(site, comp, df, dbh, h.tree, comp.mass)

## Table C1 (Appendix C)
# Coefficients of site-specific equations relating diameter at breast height (DBH) to
# stem, stem wood, stem bark, branches, foliage, and total aboveground biomass derived
# using linear regressions.

df_lr <- as.data.frame(matrix(0,30,4))
colnames(df_lr) <- c(sites, "t.values", "p.values")
params <- c("_a", "_b", "_r2", "_rmse", "_CF")
rnames <- c()
for (i in lab.comp){
  for (j in params){
    rnames <- c(rnames, paste0(i, j))
  }
}
row.names(df_lr) <- rnames
rm(i, j, rnames)

for (i in seq(2)){

  df <- subset(larix.db, Location == sites[i])
  
  for (j in seq(6)){
    
    df_lr[5*j,i] <- exp((summary(lr.D[[i]][[j]])$sigma^2)/2) #CF
    df_lr[5*j-4,i] <- exp(coefficients(lr.D[[i]][[j]])[1])*df_lr[5*j,i] #a
    df_lr[5*j-3,i] <- coefficients(lr.D[[i]][[j]])[2] #b
    df_lr[5*j-2,i] <- summary(lr.D[[i]][[j]])$r.squared #R2
    df_lr[5*j-1,i] <- sqrt(mean((df[!is.na(df[lab.comp[j]]),][[lab.comp[j]]] - exp(predict(lr.D[[i]][[j]])))^2)) #RMSE
    
    # ANOVA with interaction term to show differences in a and b coefficients
    # among site-specific regressions
    dbh <- larix.db[!is.na(larix.db[lab.comp[j]]),][["DBH"]]
    comp.mass <- larix.db[!is.na(larix.db[lab.comp[j]]),][[lab.comp[j]]]
    loc <- larix.db[!is.na(larix.db[lab.comp[j]]),][["Location"]]
    fit <- lm(log(comp.mass)~log(dbh)*loc)
    
    df_lr[5*j-4,3] <- tidy(fit)[3,4] # t value, a
    df_lr[5*j-4,4] <- tidy(fit)[3,5] # p value
    df_lr[5*j-3,3] <- tidy(fit)[4,4] # t value, b
    df_lr[5*j-3,4] <- tidy(fit)[4,5] # p value
  }
}
rm(i, df, j, dbh, comp.mass, loc, fit)

df_lr <- df_lr %>%
  dplyr::mutate(t.values = ifelse(t.values == 0, NA, t.values),
                p.values = ifelse(t.values == 0, NA, p.values),
                p.sign = ifelse(p.values >= 0.05, "ns",
                                ifelse(p.values < 0.05 & p.values >= 0.01, "*",
                                       ifelse(p.values < 0.01 & p.values >= 0.001, "**",
                                              ifelse(p.values < 0.001, "***", NA)))))

write.csv(df_lr, paste0(wdir,"outputs/tableC1_ols_summary.csv"))

## Table C2 (Appendix C)
# Coefficient of determination (R2) from linear regressions between biomass components (Y)
# and size parameters (D) after log-log transformation.

df_r2 <- as.data.frame(matrix(0,12,4))
colnames(df_r2) <- c("MG.R2", "UY.R2", "MG.p.values", "UY.p.values")
row.names(df_r2) <- c("lr.Pst.D", "lr.Pst.D2H", "lr.Pbo.D", "lr.Pbo.D2H",
                      "lr.Pba.D", "lr.Pba.D2H", "lr.Pbr.D", "lr.Pbr.D2H",
                      "lr.Pf.D", "lr.Pf.D2H", "lr.Pabo.D", "lr.Pabo.D2H")

for (i in seq(2)){
  
  k <- seq(1,11,2)
  
  for (j in seq(6)){
    
    df_r2[k[j],i] <- summary(lr.D[[i]][[j]])$r.squared
    df_r2[k[j],i+2] <- tidy(lr.D[[i]][[j]])[2,5]

    df_r2[k[j]+1,i] <- summary(lr.D2H[[i]][[j]])$r.squared
    df_r2[k[j]+1,i+2] <- tidy(lr.D2H[[i]][[j]])[2,5]
  }
}
rm(i, j, k)

df_r2 <- df_r2 %>%
  dplyr::mutate(MG.p.sign = ifelse(MG.p.values >= 0.05, "ns",
                                     ifelse(MG.p.values < 0.05 & MG.p.values >= 0.01, "*",
                                            ifelse(MG.p.values < 0.01 & MG.p.values >= 0.001, "**",
                                                   ifelse(MG.p.values < 0.001, "***", NA)))),
                UY.p.sign = ifelse(UY.p.values >= 0.05, "ns",
                                   ifelse(UY.p.values < 0.05 & UY.p.values >= 0.01, "*",
                                          ifelse(UY.p.values < 0.01 & UY.p.values >= 0.001, "**",
                                                 ifelse(UY.p.values < 0.001, "***", NA)))))

write.csv(mutate_if(df_r2, is.numeric, ~ round(.,3)),
          paste0(wdir,"outputs/tableC2_ols_rsquared.csv"))

## nonlinear regression (weigthed least squares method)

# Weightings (c values) determination (2 methods)

for (site in sites){
  df <- subset(larix.db, Location == site)
  
  for (comp in lab.comp){
    print(paste(site, comp))
    dbh <- df[!is.na(df[comp]),][["DBH"]]
    comp.mass <- df[!is.na(df[comp]),][[comp]]
    
    # (1) approximation of the conditional variance of DBH (Picard et al., 2012)
    D <- quantile(dbh, (0:5)/5) ## dividing DBH values in 5 classes centered on DBHk, k={1,...,5}
    z <- findInterval(dbh, D, rightmost.closed = TRUE)
    sdB <- data.frame(D=(D[-1]+D[-6])/2, sdB=tapply(comp.mass,z,sd))
    ## fitting a linear regression between the standard deviation of biomass and the median DBH
    ## of each class k using a log-log transformation.
    # c is approximated as the slope of this regression
    print(paste("c =",coefficients(lm(log(sdB)~I(log(D)), data = sdB))[[2]]))
    
    # (2) graphical exploration (visual assessment of the plots of the weighted residuals
    # against the fitted values for c values between 0.5 and 4)
    start <- coef(lm(log(comp.mass)~I(log(dbh))))
    start[1] <- exp(start[1])
    names(start) <- c("a","b")
    
    g <- list()
    for (j in seq(8)){
      
      wls <- nls(comp.mass~a*(dbh^b), start = start, weights = 1/(dbh^j))
      x <- fitted(wls)
      y <- residuals(wls)/dbh^(j/2)
      df.wls <- data.frame(x = x, y = y)

      g[[j]] <- ggplot(df.wls, aes(x = x, y = y)) +
        geom_point(size = 2)+
        ggtitle(paste(site, paste0(comp, paste(", c =", as.character(j/2)))))+
        xlab("Fitted values")+
        ylab("Weighted residuals")
    }
    print(ggarrange(plotlist = g, ncol = 4, nrow = 2))
  }
}
rm(df, dbh, comp.mass, D, z, sdB, start, g, wls, x, y, df.wls)
dev.off()

# Based on these two methods, c values were selected for each site and biomass component as follows
n.trees <- list()
for (site in sites){
  df <- subset(larix.db, Location == site)
  for (j in seq(6)){
    n.trees[[site]][j] <- nrow(df[!is.na(df[lab.comp[j]]),])
  }
}

## Table B1 (Appendix B)
df.c <- data.frame(Bio.comp = labels.comp,
                   MG.n = n.trees[["MG"]],
                   MG.c = c(2,2,2.5,2,1.5,2),
                   UY.n = n.trees[["UY"]],
                   UY.c = c(2.5,2.5,2,1.5,1.5,2))

write.csv(df.c, paste0(wdir,"outputs/tableB1_wls_weightings.csv"),
          row.names = FALSE)
rm(n.trees, df, site, j)

## Table 5
# Allometric equations relating diameter at breast height (DBH) to stem, stem wood,
# stem bark, branches, foliage, and total aboveground biomass at Larix cajanderi
# sites in the Magadan Oblast and the Ust-Yansky district (Republic of Sakha).

df.wnlr <- list()

for (i in seq(2)){
  
  df <- subset(larix.db, Location == sites[i])
  df.wnlr[[i]] <- data.frame(matrix(ncol = 10, nrow = 0))
  colnames(df.wnlr[[i]]) <- c("loc", "bio.comp", "min.dbh", "max.dbh",
                              "n.trees", "a", "a.se", "b", "b.se", "rmse")
  
  for (j in seq(6)){
    
    dbh <- df[!is.na(df[lab.comp[j]]),][["DBH"]]
    comp.mass <- df[!is.na(df[lab.comp[j]]),][[lab.comp[j]]]
    start <- coef(lm(log(comp.mass)~I(log(dbh))))
    start[1] <- exp(start[1])
    names(start) <- c("a","b")
    
    wnlr <- nls(comp.mass~a*(dbh^b), start = start, weights = 1/(dbh^(2*df.c[j,2*i+1])))
    
    df.wnlr[[i]] <- rbind(df.wnlr[[i]],
                          list(loc = sites[i],
                               bio.comp = labels.comp[j],
                               min.dbh = min(dbh),
                               max.dbh = max(dbh),
                               n.trees = nrow(df[!is.na(df[lab.comp[j]]),]),
                               a = coefficients(wnlr)[[1]], # a
                               a.se = summary(wnlr)$coefficients[1,2], # a [SE]
                               b = coefficients(wnlr)[[2]], # b
                               b.se = summary(wnlr)$coefficients[2,2], # b [SE]
                               rmse = sqrt(mean((comp.mass - predict(wnlr))^2)))) # RMSE
    
    
    fit.st.Bt <- Boot(wnlr, method = 'case')
    # creating predictions of each bootstrapped model
    boot1_preds <- fit.st.Bt$t %>%
      as.data.frame() %>%
      drop_na() %>%
      mutate(iter = 1:n()) %>%
      group_by_all() %>%
      do(data.frame(DBH = seq(0, max(larix.db$DBH), length.out = 100))) %>%
      ungroup() %>%
      mutate(pred = a*(DBH^b))
    # calculating bootstrapped confidence intervals
    boot1_conf_preds <- group_by(boot1_preds, DBH) %>%
      summarise(conf_lower = quantile(pred, 0.025),
                conf_upper = quantile(pred, 0.975),
                .groups = 'drop')
    # this is used to draw confidence intervals on Figure 4
    write.csv(boot1_conf_preds,
              paste0(wdir,
                     paste0('outputs/Figures/Figure4_CIs/CIs_',
                            paste0(sites[i],
                                   paste0("_",
                                          paste0(lab.comp[j],".csv"))))),
              row.names = F)
  }
}

df.wnlr <- bind_rows(df.wnlr)
rm(i, j, df, dbh, comp.mass, start, wnlr)
rm(fit.st.Bt, boot1_conf_preds, boot1_preds)

write.csv(df.wnlr,
          paste0(wdir,"outputs/table5_wls_summary.csv"),
          row.names = FALSE)

