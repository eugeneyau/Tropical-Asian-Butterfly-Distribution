# Spatial occurrence records and distributions of tropical Asian butterflies

### Eugene Yu Hin Yau, Emily E. Jones, Toby Pak Nok Tsang, Shuang Xing, Richard Corlett, Patrick Roehrdanz, David J Lohman, Adam Kai Chi Lee, Catherine Wai Ching Hai, Shawan Chowdhury, Jane K. Hill, Jade A. T. Badon, Gerard Chartier, Cheong Weei Gan, Yves Basset, I-Ching Chen, Suzan Benedick, Anuj Jain, Tiffany Ki, Krushnamegh Kunte, Lien Van Vu, Sarah Scriven, Alice C. Hughes, Timothy C. Bonebrake*
#### *E-mail: tbone@hku.hk

While regional and global distributions are available for a range of species, insects remain poorly documented biogeographically – this knowledge gap is particularly acute for tropical insects. Intensive research on the butterflies of tropical Asia in recent decades, combined with increasingly available records from citizen science in the region, have provided an opportunity to map the distributions of the butterflies of Southeast Asia. We collated records from online databases (e.g. GBIF) and extracted records from the literature to compile a dataset on occurrences. We found 303,176 occurrence records for 3,270 butterfly species across tropical Asia. Using species distribution models based on these records, we modeled the distribution of 1,000 species and identified areas of limited sampling (e.g. the Philippines and New Guinea) and also areas of high predicted diversity (Peninsular Malaysia and Sabah). The dataset provided here can be leveraged for a range of future studies on Asian butterflies including 1) research on biogeography, 2) conservation planning, 3) sampling prioritization, and 4) biodiversity threat mapping. We encourage the use and development of this database as a means for establishing a data foundation for tropical insect conservation planning in tropical Asia. For more details on our dataset please visit our [`published data paper`](https://doi.org/link)

<br>

[![](https://img.shields.io/badge/Citation-Scientific%20Data-blue)](https://doi.org/link)
[![](https://img.shields.io/badge/Archive-figshare/192596-blue)](https://doi.org/link)
[![](https://img.shields.io/badge/License-CC%20BY-blue)](https://creativecommons.org/licenses/by/4.0/)

# Table of Contents

- Cleaned butterfly occurrence records from [GBIF](https://www.gbif.org/), [Borneo Butterfly Distribution Database (B2D2)](https://www-users.york.ac.uk/~jkh6/), [Spatial occurrence data for the animals of Bangladesh derived from Facebook](https://doi.pangaea.de/10.1594/PANGAEA.948104), and other credible published sources. 
  - [`Tropical Asian Butterfly Occurrence Database`](https://drive.google.com/file/d/17MxkXPFb8Z_BJiF1wKfRxd8aGZRYcZ51/view?usp=sharing)
- Distribution maps of tropical Asian butterflies as predicted by species distribution models(SDMs).
  - [`Predicted Individual Species Distribution Maps`](https://figshare.com/LINK)
- R script used to construct SDMs.
  - [`Code/SDM`](https://github.com/eugeneyau/Tropical-Asian-Butterfly-Distribution/tree/main/Code/SDM)
     - [`R script`](https://github.com/eugeneyau/Tropical-Asian-Butterfly-Distribution/blob/main/Code/SDM/SEA_Bfy_SDM.R)
     - [`R Markdown file`](https://github.com/eugeneyau/Tropical-Asian-Butterfly-Distribution/blob/main/Code/SDM/RMD_SEA_Bfy_SDM.Rmd)
- R scripts used to match species names / clean dataset, map diversity patterns, and extract SDM variable importance.
  - [`Code/Supplementary`](https://github.com/RhettRautsaw/VenomMaps/tree/master/code)
     - [`Match species names / clean dataset`](https://github.com/eugeneyau/Tropical-Asian-Butterfly-Distribution/blob/main/Code/Supplementary/SEA_Bfy_Data-Update_sp_name.R)
     - [`Map diversity patterns by stacking single species distributions`](https://github.com/eugeneyau/Tropical-Asian-Butterfly-Distribution/blob/main/Code/Supplementary/SEA_Bfy_Data-Plot_zoned_alpha_diversity.R)
     - [`SDM variable importance`](https://github.com/eugeneyau/Tropical-Asian-Butterfly-Distribution/blob/main/Code/Supplementary/SEA_Bfy_Data-Var_imp.R)
     - [`SDM performance evaluation`](https://github.com/eugeneyau/Tropical-Asian-Butterfly-Distribution/blob/main/Code/Supplementary/SEA_Bfy_Data-Model_eval_summary.R)
- Files essential for running our R scripts, please download all of them before running our codes.
  - [`Code/Files`](https://github.com/eugeneyau/Tropical-Asian-Butterfly-Distribution/tree/main/Code/Files)
- JavaScript code used in Google Earth Engine to extract and filter Landsat data for use as SDM variable (NDVImean), view in [`Google Earth Engine`](https://code.earthengine.google.com/7e1c649f06f22536419886e34a14d830) or download code from here:
  - [`Code/Variables`](https://github.com/RhettRautsaw/VenomMaps/tree/master/code)


