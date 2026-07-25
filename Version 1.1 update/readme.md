# Occurrence Record Version 1.1

Occurrence Records of Tropical Asian Butterflies: 1970-2024 v1.1 addresses multiple taxonomic errors identified in the GBIF records in our original dataset, which stemmed from issues with taxonomic harmonization (see #1 below), and incorporates 17,635 GBIF records omitted from v1 (see #2 below).

V1.1 now includes 747,825 occurrence records for 3,738 species. Records of Nymphalidae (1,358 spp.; 327,643 records) comprise 44% of the dataset, followed by Lycaenidae (1,076 spp; 156,162 records), Papilionidae (265 spp.; 100,904), Pieridae (391 spp.; 93,376 records), Hesperiidae (638 spp.; 65,613 records), and Riodinidae (33 spp.; 4,127 records). Updated dataset metadata is available at v1.1 Metadata for Occurrence Records of Tropical Asian Butterflies 1970-2024_23July2026EJ. Updated SDMs incorporating these changes are now available in SDMsupp_files.zip. 

_Please note: Data in v1.1 is reordered from v1. Index column data are not unique identifiers for data in v1 and v1.1; they only denote the row number. The datasetKey data for all non-GBIF data may still be used as record identifiers._
<br>

### Description of major changes
- We bypassed the GBIF harmonization pipeline and harmonized GBIF’s “verbatimScientificName” data (available in file 0169126-240321170329656.csv) directly using Lamas’ Catalogue of Butterflies (2015). Data from “verbatimScientificName” that could not be harmonized using this catalogue were harmonized to other reputable sources (i.e., Catalogue of Life, etc.). GBIF data that produced inconsistent binomials across databases were omitted. A full list of these updated harmonizations (GBIF_harmonization_YFL.csv produced by Ling Yuet Fung) as well as the R script for implementing them is available [here](https://github.com/eugeneyau/Tropical-Asian-Butterfly-Distribution/tree/main/Version%201.1%20update). Frequent updates to GBIF backbone taxonomy are sourced from the Catalogue of Life. For information on GBIF’s data processing and backbone, see [`GBIF documentation`](https://techdocs.gbif.org/en/data-processing/).
- Due to changes described in #1, we were able to include 17,635 GBIF entries omitted from the former dataset due to missing data in the “species” column.
