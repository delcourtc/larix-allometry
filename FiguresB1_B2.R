# Code for plotting Figure B1 and Figure B2

## setting the working directory in which the outputs from the larix_allometry.R script were saved
wdir <- ''

## loading required packages
library(tidyverse)

## inputs
larch <- paste0(wdir,"Larch_tree_DB.csv")
wght <- paste0(wdir, "tableB1_wls_weightings.csv")

## loading data
larix.db <- read.csv(larch)
df.c <- read.csv(wght)

sites <- c("MG", "UY") # MG: Magadan Oblast; UY: Ust-Yansky district
lab.comp <- c("Pst", "Pbole", "Pbark", "Pbr", "Pf", "Pabo")
labels.comp <- c("Stem", "Stem wood", "Stem bark",
                 "Branches", "Foliage", "Aboveground")
cols <- c("#020305", "#3B6797")
panels <- c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)")

for(i in seq(2)){
  
  png(paste0(wdir,paste0("Figures/FigureB",paste0(as.character(i), ".png"))),
      units = "cm", width = 15.7, height = 19.6, res = 300)
  par(mfrow = c(3,2), mgp = c(1.5,0.2,0), mar = c(3, 3.4, 1.5, 1),
      tck = 0.015, ps = 7)
  
  df <- subset(larix.db, Location == sites[i])
  
  for (j in seq(6)){
    
    dbh <- df[!is.na(df[lab.comp[j]]),][["DBH"]]
    comp.mass <- df[!is.na(df[lab.comp[j]]),][[lab.comp[j]]]
    start <- coef(lm(log(comp.mass)~I(log(dbh))))
    start[1] <- exp(start[1])
    names(start) <- c("a","b")
    
    wnlr <- nls(comp.mass~a*(dbh^b), start = start, weights = 1/(dbh^(2*df.c[j,2*i+1])))
    
    plot(fitted(wnlr), residuals(wnlr)/dbh^(df.c[j,2*i+1]),
         xlab = "Fitted values", ylab = "Weighted residuals",
         pch = 16, col = cols[i], cex = 1.1, cex.lab = 1.7, cex.axis = 1.6)
    
    text(x = 0.64*max(fitted(wnlr)), y = 0.9*max(residuals(wnlr)/dbh^(df.c[j,2*i+1])),
         labels = labels.comp[j], col = cols[i], cex = 1.9, pos = 4)
    text(x = 0.64*max(fitted(wnlr)), y = 0.78*max(residuals(wnlr)/dbh^(df.c[j,2*i+1])),
         labels = paste("c =", as.character(df.c[j,2*i+1])), col = cols[i], cex = 2, pos = 4)
    text(x = -0.15*(max(fitted(wnlr))-min(fitted(wnlr))),
         y = max(residuals(wnlr)/dbh^(df.c[j,2*i+1])),
         labels = panels[j], xpd = NA,
         cex = 2.5, font = 2)
    
  }
  dev.off()
}
