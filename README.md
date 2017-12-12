# geospat
Material for the Course "Spatial Data Analysis" (CDSS, Mannheim, Spring 2018)


# Course Content

The first law of geography states that “everything is related to everything else, but near things are more related than distant things” (Tobler 1970: 236). In the Social Sciences, geographic data and spatial analyses offer rich insights into a variety of relevant research questions (Franzese and Hays 2008). The course covers crucial concepts involved in spatial analysis, introduces a toolbox of statistical models and pays particular attention to the accessible implementation of spatial analysis in free software (working with R and packages for spatial analysis). This implies that participants should bring their Laptop with R (required) and RStudio (recommended) installed. 
 

# Parts

**I. Concepts

The first part of the course deals with the concepts involved in spatial analysis. We will spend some time getting to know “W”, the connectivity matrix which defines spatial dependencies (Neumayer and Plümper 2016). Alternative conceptions of neighborhood and the weighting of connections are discussed. For example, spatial proximity does not necessarily imply geographic proximity, as for instance trade or information exchange can bring distance things closely together. Further concepts handled include the geo-referencing of data, regarding both its use and generation.

** II. Models

To test hypotheses based on spatial data, tailored statistical tools are needed. The second part of the course is dedicated to spatial correlation coefficients (such as Moran’s I), varieties of spatial regression models (variants of spatial lags; categorical, count and duration specifications), spatio-temporal models as well as extensions to the multilevel case. The options offered by a Bayesian approach to spatial data analysis are also discussed.

** III. Implementation

As mentioned, one focus of the course is the accessible implementation of spatial analysis. To this end, the free statistical software R is used (https://www.r-project.org/). Its advantages in addition to the open source character are the provision of user-written packages, including several on geodata and spatial analysis (such as sp, maptools or spdep) as well as powerful graphical capabilities (see Bivand et al. 2013). Recommended is the combination of R with the RStudio editor/environment (https://www.rstudio.com/). For the Bayesian models, we rely on JAGS (http://mcmc-jags.sourceforge.net/) or Stan (http://mc-stan.org/). During this part of the course, geographic data are handled, autocorrelations and spatial regression models are estimated, and everything is thoroughly visualized.

** IV. Applications

The last part of the course is devoted to the participants’ own applications in the field of geographic data and spatial analysis. Anticipated throughout the course, they will look for and handle their own geographic data and we will jointly identify adequate spatial models to test hypotheses. The final meeting is dedicated to the presentation and discussion of the paper outlines, with a focus on a (preliminary) spatial analysis of geographic data.


# References

Bivand, R. S., E. Pebesma and V. Gómez-Rubio (2013): Applied Spatial Data Analysis with R (2nd Edition). New York: Springer.

Franzese, R. J. and J. C. Hays (2008): Interdependence in Comparative Politics. Substance, Theory, Empirics, Substance. Comparative Political Studies 41(4/5): 742-780.

Neumayer, E. and T. Plümper (2016): W. Political Science Research and Methods 4(1): 175-193.

Tobler, W. R. (1970): A Computer Movie Simulating Urban Growth in the Detroit Region. Economic Geography 46(Suppl.): 234-340.

