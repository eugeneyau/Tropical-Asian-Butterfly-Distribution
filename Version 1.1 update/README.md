# Occurrence Record Version 1.1

Occurrence Records of Tropical Asian Butterflies: 1970-2024 v1.1 addresses multiple taxonomic errors identified in the GBIF records in our original dataset which stemmed from issues with taxonomic harmonization (see #1 below) and incorporates 17,635 GBIF records omitted from v1 (see #2 below).

V1.1 now includes 747,825 occurrence records for 3,738 species. Records of Nymphalidae (1,358 spp.; 327,643 records) comprise 44% of the dataset, followed by Lycaenidae (1,076 spp; 156,162 records), Papilionidae (265 spp.; 100,904), Pieridae (391 spp.; 93,376 records), Hesperiidae (638 spp.; 65,613 records), and Riodinidae (33 spp.; 4,127 records). Updated dataset metadata is available at v1.1 Metadata for Occurrence Records of Tropical Asian Butterflies 1970-2024_23July2026EJ. Updated SDMs incorporating these changes are now available in SDMsupp_files.zip. 

_Please note: Data in v1.1 is reordered from v1. Index column data are not unique identifiers for data in v1 and v1.1; they only denote the row number. The datasetKey data for all non-GBIF data may still be used as record identifiers._

<br>

# 📰 News
* **June 2026**: Released occurrence dataset v1.1 (Occurrence Records of Tropical Asian Butterflies: 1970-2024 v1.1) that addresses multiple taxonomic errors identified in the GBIF records in our original dataset, which stemmed from issues with taxonomic harmonization and incorporates 17,635 GBIF records omitted from v1.
* **June 2025**: Occurrence dataset v1 released with our [`data paper`](https://doi.org/10.1038/s41597-025-05333-w).

# Table of Contents

- Cleaned butterfly occurrence records from [GBIF](https://doi.org/10.15468/dd.nvw5wr), [Borneo Butterfly Distribution Database (B2D2)](https://www-users.york.ac.uk/~jkh6/), [Spatial occurrence data for the animals of Bangladesh derived from Facebook](https://doi.pangaea.de/10.1594/PANGAEA.948104), and other credible published sources. 
  - [`Tropical Asian Butterfly Occurrence Database`](https://doi.org/10.6084/m9.figshare.25037645)
  - [`Version 1.1 update backbone taxonomy`](https://github.com/eugeneyau/Tropical-Asian-Butterfly-Distribution/tree/main/Version%201.1%20update)
     - [`R script to correct GBIF binomial synonym harmonization`](https://github.com/eugeneyau/Tropical-Asian-Butterfly-Distribution/blob/main/Version%201.1%20update/Supp_update_sp_name--correct_GBIF.R)
     - [`Updated harmonization list for GBIF records`](https://github.com/eugeneyau/Tropical-Asian-Butterfly-Distribution/blob/main/Version%201.1%20update/GBIF_harmonization_YFL.csv)
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


