# ml4bd: Microbial communities and environmental conditions in Bangladesh
### Machine Learning for Life Sciences - Autumn 2026 (Lovisa, Justine, Abu)

## In this repo...
* In [scripts](https://github.com/abu85/ml4bd/tree/main/scripts), you can find all scripts used for analysis and visualization of the data
* In [figures](https://github.com/abu85/ml4bd/tree/main/figures), all figures produced during analyses are stored

## Source dataset
- 17 samples (16S rDNA ampliSeq, freshwater) 
- 3 locations (districts) 
- 4041 taxa (ASV)
- Environmental variables measured at sampling locations: pH, DO (Dissolved Oxygen) concentrations, Temperature, Ammonia concentration
- Sequencing artifacts effects or noise was removed (we removed ASVs that have less than 5 counts over all samples)

## Shiny App
Please access our results through the Shiny app: [https://lovisalindquist.shinyapps.io/ml4bd/](https://lovisalindquist.shinyapps.io/ml4bd/)

## Used packages
- phyloseq 1.44.0
- vegan 2.6-10
- mixOmics version 6.24.0
- ggplot2 version 4.0.3
- 

# Members contributions
- All members contributed to the choice of methods
- All members worked on the preprocessing of the data
- Justine ran the NMDS/PERMANOVA and PLS-DA
- Abu ran the sPLS-DA
- Lovisa ran the hierarchical clustering, built the Shiny App and deployed it
- All members participated in the communication of the results (powerpoint, narrative, github...)
