# Allometric equations and woody debris parameters for estimating aboveground biomass in Cajander larch forests of northeast Siberia
Python and R codes for the development of new allometric equations for *Larix cajanderi* trees in northeast Siberia

[![DOI](https://zenodo.org/badge/531168048.svg)](https://zenodo.org/badge/latestdoi/531168048)

This code uses tree biomass and size data retrieved from a comprehensive dataset of biomass measurements in Eurasia (Schepaschenko et al., 2017). The original excel file can be downloaded on https://doi.org/10.1594/PANGAEA.871491. We selected raw harvest data from *Larix cajanderi* stands sampled as part of two experiments conducted in the Magadan Oblast (MG) and the Ust-Yansky district (UY) (Republic of Sakha), Russia. These data were used to develop allometric equations relating tree diameter at breast height (DBH) to stem, stem wood, stem bark, branches, foliage, and total aboveground tree biomass.

To illustrate differences between our newly developed allometric relationships and existing biomass equations for *L. cajanderi*, we calculated aboveground tree biomass at 53 forest stands sampled near Yakutsk (Russia) using each allometric model. Field measurements collected within these stands are freely available from https://doi.org/10.5281/zenodo.7049450. This dataset was also used to determine density and size parameters for fine woody debris (FWD, < 7 cm in diameter) of *L. cajanderi* that can be used to calculate FWD biomass using the line-intersect method (Van Wagner, 1982; Nalder et al., 1999). The allometric equations and FWD parameters presented in our study can be used to refine estimates of aboveground biomass in Cajander larch forests of northeast Siberia.

The source codes for some display items in the article (Figures 4, 5, B1, B2, C1 and C2; Tables 3, 4, 5, A1, B1, C1, C2) are also included.

The larix_allometry.R script includes:
- preprocessing steps that generate the data for analysis
- graphical exploration of the relationships between tree components biomass and diameter at breast height (DBH) for each site
- development of site-specific allometric equations derived using linear regressions (ordinary least squares method) (including Table C1, Table C2)
- determination of the weightings used in nonlinear regressions following Picard et al. (2012)
- development of site-specific allometric equations derived using nonlinear regressions (weighted least squares method) (including Table 5, Table B1)

The larix_fwd.R script includes:
- determination of size-class-specific values of specific gravity (G)
- determination of size-class-specific values of mean squared diameter (MSD) (including Table 3)

## Workflow

1. Development of allometric equations (larix_allometry.R)
2. Figures related to allometric equations:
  - code for Figures B1 and B2 (FiguresB1_B2.R)
  - code for Figures 4 and C1 (Figures4_C1.py)
  - code for Figures 5 and C2 (Figures5_C2.py)
 3. Determination of FWD density and size parameters (larix_fwd.R)
 4. Figure related to FWD analysis (Figure3.py)
  
## Data requirements
  - Biomass measurements from Schepaschenko et al. (2017) (Biomass_tree_DB.xlsx)
  - Characteristics of the 53 forest stands sampled near Yakutsk in the summer of 2019 (YA2019_plots.xlsx)
  - Field measurements collected for estimates of aboveground tree biomass (YA2019_tree_transects.xlsx)
  - Raw data from the water displacement method required to derive size-class-specific values of specific gravity for FWD (YA2019_fwd_specific_gravity.xlsx)
  - Field measurements of FWD diameters for determination of size-class-specific values of mean squared diameter (MSD) (YA2019_fwd_diameters.xlsx)
  - Field inventory of fine woody debris (FWD) using the line-intersect method (YA2019_fwd_transects.xlsx)
  
## Software requirements

R code tested with R 4.0.3, required packages:
  - openxlsx 4.2.3
  - tidyverse 1.3.0
  - ggpubr 0.4.0
  - broom 0.7.5
  - car 3.0-10
  - emmeans 1.5.5-1
  - multcompView 0.1-8
  
Python code tested with Anaconda Python 3.8.5, required packages:
  - pandas 1.2.1
  - numpy 1.19.2
  - matplotlib 3.3.2

Please note that this code was not written by a professional software developer, so it may not be written in the most efficient way possible. Feel free to contact us if you have any comments, suggestions or questions regarding the code, data or the analysis in general.

## References

Delcourt, C. J. F. and Veraverbeke, S.: Field measurements for estimating aboveground and woody debris biomass in Cajander larch forests of northeast Siberia, Zenodo [data set], https://doi.org/10.5281/zenodo.7049450, 2022.

Nalder, I. A., Wein, R. W., Alexander, M. E., and de Groot, W. J.: Physical properties of dead and downed round-wood fuels in the Boreal forests of western and Northern Canada, Int. J. Wildland Fire, 9, 85–99, https://doi.org/10.1071/WF00008, 1999.

Picard, N., Saint-André, L., and Henry, M.: Manual for building tree volume and biomass allometric equations: from field measurements to prediction, Food and Agricultural Organization of the United Nations, Rome, Italy, and Centre de Coopération Internationale en Recherche Agronomique pour le Développement, Montpellier, France, 215 pp., E-ISBN: 9789251073476, 2012.

Schepaschenko, D., Shvidenko, A., Usoltsev, V., Lakyda, P., Luo, Y., Vasylyshyn, R., Lakyda, I., Myklush, Y., See, L., Mc-Callum, I., Fritz, S., Kraxner, F., and Obersteiner, M.: Biomass tree data base, PANGAEA [data set], https://doi.org/10.1594/PANGAEA.871491, 2017.

Van Wagner, C. E.: Practical aspects of the line intersect method, Canadian Forestry Service, Maritimes Forest Research Centre, Fredericton, New Brunswick, Information Report PI-X-12E, 11 pp., ISBN 0-662-11816-2, 1982.
