# Spatial occurrence records and distributions of tropical Asian butterflies

### Eugene Yu Hin Yau, Emily E. Jones, Toby Pak Nok Tsang, Shuang Xing, Richard T. Corlett, Patrick Roehrdanz, David J. Lohman, Adam Kai Chi Lee, Catherine Wai Ching Hai, Shawan Chowdhury, Jane K. Hill, Jade A. T. Badon, Cheong Weei Gan, Yves Basset, I-Ching Chen, Suzan Benedick, Anuj Jain, Tiffany L.T. Ki, Krushnamegh Kunte, Akihiro Nakamura, Lien Van Vu, Sarah A. Scriven, Alice C. Hughes, Timothy C. Bonebrake*
#### *E-mail: tbone@hku.hk

<img align="right" src="https://github.com/eugeneyau/Tropical-Asian-Butterfly-Distribution/blob/main/md_images/Sampling%20density%20with%20zones1.png" width=380> 

Insect biogeography is poorly documented globally, particularly in the tropics. Recent intensive research in tropical Asia, combined with increasingly available records from citizen science, provides an opportunity to map the distributions of tropical Asian butterflies. We compiled a dataset of 724,247 occurrences of 3,581 tropical Asian butterfly species by aggregating records from GBIF (651,285 records), published literature (21,271), published databases (37,695), and unpublished data (13,993). Here, we present this dataset and single-species distribution maps of 1,520 species. Using these maps, along with records of the 2,071 remaining species, we identified areas of limited sampling (e.g., the Philippines, Myanmar, and New Guinea) and predicted areas of high diversity (Peninsular Malaysia and Borneo). This dataset can be leveraged for a range of studies on Asian and tropical butterflies, including 1) species biogeography, 2) sampling prioritization to fill gaps, 3) biodiversity hotspot mapping, and 4) conservation evaluation and planning. We encourage the continued development of this dataset and the associated code as a tool for the conservation of tropical Asian insects. For more details on our dataset and SDMs, please refer to our [`data paper`](https://doi.org/link). 


<br>

[![](https://img.shields.io/badge/Citation-Scientific%20Data-blue)](https://doi.org/link)
[![](https://img.shields.io/badge/Archive-figshare/25037645-blue)](https://doi.org/10.6084/m9.figshare.25037645)
[![](https://img.shields.io/badge/License-CC%20BY%204.0-blue)](https://creativecommons.org/licenses/by/4.0/)

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
  - [`Code/Supplementary`](https://github.com/eugeneyau/Tropical-Asian-Butterfly-Distribution/tree/main/Code/Supplementary)
     - [`Match species names / clean dataset`](https://github.com/eugeneyau/Tropical-Asian-Butterfly-Distribution/blob/main/Code/Supplementary/SEA_Bfy_Data-Update_sp_name.R)
     - [`Map diversity patterns by stacking single species distributions`](https://github.com/eugeneyau/Tropical-Asian-Butterfly-Distribution/blob/main/Code/Supplementary/SEA_Bfy_Data-Plot_zoned_alpha_diversity.R)
     - [`SDM variable importance`](https://github.com/eugeneyau/Tropical-Asian-Butterfly-Distribution/blob/main/Code/Supplementary/SEA_Bfy_Data-Var_imp.R)
     - [`SDM performance evaluation`](https://github.com/eugeneyau/Tropical-Asian-Butterfly-Distribution/blob/main/Code/Supplementary/SEA_Bfy_Data-Model_eval_summary.R)
- Files essential for running our R scripts, please download all of them before running our codes.
  - [`Code/Files`](https://github.com/eugeneyau/Tropical-Asian-Butterfly-Distribution/tree/main/Code/Files)
- JavaScript code used in Google Earth Engine to extract and filter Landsat data for use as SDM variable (NDVImean), view in [`Google Earth Engine`](https://code.earthengine.google.com/7e1c649f06f22536419886e34a14d830) or download code from here:
  - [`Code/Variables`](https://github.com/eugeneyau/Tropical-Asian-Butterfly-Distribution/blob/main/Code/Variables/GEE_NDVImean.txt)


