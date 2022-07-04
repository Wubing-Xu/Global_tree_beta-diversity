# Global_tree_beta-diversity
This repository contains data and code necessary to reproduce the results shown in the manuscript" Current global tree beta-diversity is shaped by Quaternary climate change".

## Data
Tree species occurrences came from Serra-Diaz et al. (2018). These occurrences were compiled from five major comprehensive biodiversity infrastructures, including the [Global Biodiversity Information Facility](https://www.gbif.org), the [Botanical Information and Ecological Network](http://bien.nceas.ucsb.edu/bien), the [Latin American Seasonally Dry Tropical Forest Floristic Network]( http://www.dryflor.info), the [RAINBIO](https://gdauby.github.io/rainbio/index.html) database and the [Atlas of Living Australia](https://www.ala.org.au). All species occurrences have been deposited in the [BIEN] (https://bien.nceas.ucsb.edu/bien/) database. 

Species range maps were then constructed from occurrences based on the alpha-hull method using different values of alpha (2, 4, 6, and 10). Species ranges were then rasterized to grid-cells in a resolution of 200 km with an equal-area Behrmann projection. Files in folder `data/tree_pam` provide matrixes with grid-cells as rows and species as columns containing all tree species or the species used in this study. The alpha-hull range maps can be found in [this repository]( https://github.com/wyeco/TC_conservation).

The phylogeny and imputed trait data were provided in the folder ` data/ traits_phylogeny `.  The raw trait data are not available due to data privacy and sharing restrictions of TRY, but they can be obtained by submitting a data request to [TRY]( https://www.try-db.org/) database. 

## R analysis files
These scripts were used to prepare the data, calculate beta-diversity, run null models, fit statistical models, and generate figures and tables.

Please note that some of the code in this repository was written to run on a HPC cluster.

**Contact:** wbingxu@gmail.com (Wubing Xu)
