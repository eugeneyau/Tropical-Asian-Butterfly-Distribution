# Spatial occurrence records and distributions of tropical Asian butterflies

### Eugene Yu Hin Yau, Emily E. Jones, Toby Pak Nok Tsang, Shuang Xing, Richard T. Corlett, Patrick Roehrdanz, David J. Lohman, Adam Kai Chi Lee, Catherine Wai Ching Hai, Shawan Chowdhury, Jane K. Hill, Jade A. T. Badon, Cheong Weei Gan, Yves Basset, I-Ching Chen, Suzan Benedick, Anuj Jain, Tiffany L.T. Ki, Krushnamegh Kunte, Akihiro Nakamura, Lien Van Vu, Sarah A. Scriven, Alice C. Hughes, Timothy C. Bonebrake*
#### *E-mail: tbone@hku.hk

<img align="right" src="https://github.com/eugeneyau/Tropical-Asian-Butterfly-Distribution/blob/main/md_images/Figure%201.jpg" width=420> 

Insect biogeography is poorly documented globally, particularly in the tropics. Recent intensive research in tropical Asia, combined with increasingly available records from citizen science, provides an opportunity to map the distributions of tropical Asian butterflies. We compiled a dataset of 730,190 occurrences of 3,752 tropical Asian butterfly species by aggregating records from GBIF (651,285 records), published literature (27,217), published databases (37,695), and unpublished data (13,993). Here, we present this dataset and single-species distribution maps of 1,576 species. Using these maps, along with records of the 2,176 remaining species, we identified areas of limited sampling (e.g., the Philippines, Myanmar, and New Guinea) and predicted areas of high diversity (Peninsular Malaysia and Borneo). This dataset can be leveraged for a range of studies on Asian and tropical butterflies, including 1) species biogeography, 2) sampling prioritization to fill gaps, 3) biodiversity hotspot mapping, and 4) conservation evaluation and planning. We encourage the continued development of this dataset and the associated code as a tool for the conservation of tropical Asian insects. For more details on our dataset and SDMs, please refer to our [`data paper`](https://doi.org/10.32942/X2C904). 


<br>

[![](https://img.shields.io/badge/Citation-Scientific%20Data-blue)](https://doi.org/link)
[![](https://img.shields.io/badge/Archive-figshare/25037645-blue)](https://doi.org/10.6084/m9.figshare.25037645)
[![](https://img.shields.io/badge/License-CC%20BY%204.0-blue)](https://creativecommons.org/licenses/by/4.0/)

# Table of Contents

- Cleaned butterfly occurrence records from [GBIF](https://doi.org/10.15468/dd.nvw5wr), [Borneo Butterfly Distribution Database (B2D2)](https://www-users.york.ac.uk/~jkh6/), [Spatial occurrence data for the animals of Bangladesh derived from Facebook](https://doi.pangaea.de/10.1594/PANGAEA.948104), and other credible published sources. 
  - [`Tropical Asian Butterfly Occurrence Database`](https://doi.org/10.6084/m9.figshare.25037645)
- Distribution maps of tropical Asian butterflies as predicted by species distribution models(SDMs) can be downloaded as separate raster files or one single PDF file from our [`Figshare repository`](https://doi.org/10.6084/m9.figshare.25037645).
- R script used to construct SDMs:
  - [`Code/SDM`](https://github.com/eugeneyau/Tropical-Asian-Butterfly-Distribution/tree/main/Code/SDM)
     - [`R Markdown file`](https://github.com/eugeneyau/Tropical-Asian-Butterfly-Distribution/blob/main/Code/SDM/RMD_TropicalAsia_Bfy_SDM.Rmd)
     - [`R script`](https://github.com/eugeneyau/Tropical-Asian-Butterfly-Distribution/blob/main/Code/SDM/TropicalAsia_Bfy_SDM.R)
- Please download files essential for running our R scripts from our [`Figshare repository (SDMsupp_files.zip)`](https://doi.org/10.6084/m9.figshare.25037645) before running our codes.
- Additional R scripts used to clean our dataset, prepare for distribution modeling, analyze SDM outputs, and buffer occurrence points of unmodelled species:
  - [`Code/Supplementary`](https://github.com/eugeneyau/Tropical-Asian-Butterfly-Distribution/tree/main/Code/Supplementary)
     - [`Harmonize species names`](https://github.com/eugeneyau/Tropical-Asian-Butterfly-Distribution/blob/main/Code/Supplementary/Bfy_Data_supp_update_sp_name.R)
     - [`Clean dataset family names`](https://github.com/eugeneyau/Tropical-Asian-Butterfly-Distribution/blob/main/Code/Supplementary/Bfy_Data_supp_update_family_name.R)
     - [`Identify possible biogeographic range of dispersal for each species`](https://github.com/eugeneyau/Tropical-Asian-Butterfly-Distribution/blob/main/Code/Supplementary/Bfy_Data_supp_id_landmass_mask.R)
     - [`Calculate spatial sampling effort`](https://github.com/eugeneyau/Tropical-Asian-Butterfly-Distribution/blob/main/Code/Supplementary/Bfy_Data_supp_get_dens_ras.R)
     - [`Map diversity patterns by stacking single species distributions`](https://github.com/eugeneyau/Tropical-Asian-Butterfly-Distribution/blob/main/Code/Supplementary/Bfy_Data_supp_plot_alpha_diversity.R)
     - [`SDM performance evaluation`](https://github.com/eugeneyau/Tropical-Asian-Butterfly-Distribution/blob/main/Code/Supplementary/Bfy_Data_supp_eval_summary%20(PO).R)
     - [`SDM variable importance`](https://github.com/eugeneyau/Tropical-Asian-Butterfly-Distribution/blob/main/Code/Supplementary/Bfy_Data_supp_var_imp_analysis.R)
     - [`Buffer occurrence points of unmodelled species`](https://github.com/eugeneyau/Tropical-Asian-Butterfly-Distribution/blob/main/Code/Supplementary/Bfy_Data_supp_unmodelled%20species_point_richness.R)
     - [`Compare SDM results with that of Daru (2024)`](https://github.com/eugeneyau/Tropical-Asian-Butterfly-Distribution/blob/main/Code/Supplementary/Bfy_Data_supp_validate_daru2024_map.R)
- JavaScript code used in Google Earth Engine to extract and filter Landsat data for use as SDM variable (NDVImean), view in [`Google Earth Engine`](https://code.earthengine.google.com/7e1c649f06f22536419886e34a14d830) or download code from here:[`Code/Variables`](https://github.com/eugeneyau/Tropical-Asian-Butterfly-Distribution/blob/main/Code/Variables/GEE_NDVImean.txt)


