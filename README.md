# Allometric equations for estimating aboveground tree biomass in Cajander larch forests of northeast Siberia
Python and R codes for the development of new allometric equations for Larix cajanderi trees in northeast Siberia

This code uses tree biomass and size data retrieved from a comprehensive dataset of biomass measurements in Eurasia (Schepaschenko et al., 2017). The original excel file can be downloaded on https://doi.org/10.1594/PANGAEA.871491. We selected raw harvest data from L. cajanderi stands sampled as part of two experiments conducted in the Magadan Oblast (MG) and the Ust-Yansky district (UY) (Republic of Sakha), Russia. The source codes for some display items in the article (Figures 4, C1, B1 and B2; Tables 5, B1, C1, C2) are also included.

The larix_allometry.R script includes:
- preprocessing steps that generate the data for analysis
- graphical exploration of the relationships between tree components biomass and diameter at breast height (DBH) for each site
- development of site-specific allometric equations derived using linear regressions (ordinary least squares method) (including Table C1, Table C2)
- determination of the weightings used in nonlinear regressions following Picard et al. (2012)
- development of site-specific allometric equations derived using nonlinear regressions (weighted least squares method) (including Table 5, Table B1)

## Workflow

1. Development of allometric equations (larix_allometry.R)
2. Figures
  - code for Figures B1 and B2 (FiguresB1_B2.R)
  - code for Figures 4 and C1 (Figures4_C1.py)

Please note that this code was not written by a professional software developer, so it may not be written in the most efficient way possible. Feel free to contact us if you have any comments, suggestions or questions regarding the code, data or the analysis in general.

## References

Picard, N., Saint-André, L., and Henry, M.: Manual for building tree volume and biomass allometric equations: from field measurements to prediction, Food and Agricultural Organization of the United Nations, Rome, Italy, and Centre de Coopération Internationale en Recherche Agronomique pour le Développement, Montpellier, France, 215 pp., E-ISBN: 9789251073476, 2012.

Schepaschenko, D., Shvidenko, A., Usoltsev, V., Lakyda, P., Luo, Y., Vasylyshyn, R., Lakyda, I., Myklush, Y., See, L., Mc-Callum, I., Fritz, S., Kraxner, F., and Obersteiner, M.: Biomass tree data base, PANGAEA [data set], https://doi.org/10.1594/PANGAEA.871491, 2017.
