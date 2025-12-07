# Marine-fungi-response-to-Orange-County-oil-spill-2021

This repository contains scripts, supplementary spreadsheets, and databases for "Subtle changes in community composition of marine fungi following the Orange County, California, oil spill"

# Folder: Spreadsheets
Spreadsheet_A.1(S1). Metadata and number of reads per subsample at different steps of the sequence data processing. Sequence processing was carried out using the AMPtk pipeline (Choreno_Parra_etal_script_oilspill_AMPtk.sh), and filtering for fungi on the resulting OTU table (Choreno_Parra_etal_script_oilspill_Diversity.R).

Spreadsheet_A.2(S2). Taxonomy and number of counts per fungal OTU in each sample. This table results from filtering for "Fungi" in the OTU table originally obtained using the AMPtk pipeline (Choreno_Parra_etal_script_oilspill_AMPtk.sh). The table contains the data prior to filtering for marine fungal OTUs (Choreno_Parra_etal_script_oilspill_Diversity.R). 

Spreadsheet_A.3(S3). Direction and significance of shifts in relative abundance of fungal genera across oiling degrees and time after the oil spill. Obtained using Choreno_Parra_etal_script_oilspill_Diversity.R

Spreadsheet_A.4(S4). Taxonomy and significance of indicator OTUs across oiling degrees and time after the oil spill. Obtained using Choreno_Parra_etal_script_oilspill_Diversity.R

# Folder: Environmental_databases
noaa_46253_mwd.txt. Database containing Mean Wave Direction data. Retrieved from the National Data Buoy Center, station 46253: https://www.ndbc.noaa.gov/.

jws_hb_buoy_water_temp.csv. Database containing Surface Water Temperature data. Available on JWS HB Buoy operated by California State University, Long Beach: https://sensors.ioos.us/#metadata/135130/station/data.

sccoos_newport_salinity.csv. Database containing Salinity data. Retrieved from the Southern California Coastal Ocean Observing System (SCCOOS) collected at the Newport Beach Pier Automated Shore Station: https://sccoos.org/autoss/.

sccoos_newport_chlorophyll.csv. Database containing Chlorophyll a data. Retrieved from the Southern California Coastal Ocean Observing System (SCCOOS) collected at the Newport Beach Pier Automated Shore Station: https://sccoos.org/autoss/.

CUTI_daily.csv. Database containing CUTI. Available on Jacox et al. (2018): https://mjacox.com/upwelling-indices/.

BEUTI_daily.csv. Database containing BEUTI. Available on Jacox et al. (2018): https://mjacox.com/upwelling-indices/.

bco_dmo_pom.csv. Database containing POC and PON data. Available on Gerace et al. (2025): https://www.bco-dmo.org/dataset/564351.

# Folder: Scripts
Choreno_Parra_etal_script_oilspill_AMPtk.sh. This script contains code to process raw ITS2 sequences obtained from surface seawater samples collected after the Orange County oil spill. Sequences are available from NCBI’s SRA under BioProject ID PRJNA1344880: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1344880/.

Choreno_Parra_etal_script_oilspill_Diversity.R. This script contains R code for the analysis of alpha (richness) and beta diversity, trend analysis, comparison of relative abundances of genera, and analysis of indicator OTUs. Use with Spreadsheet_A.1(S1) and Spreadsheet_A.2(S2).

Choreno_Parra_etal_script_oilspill_Environment.R. This script contains R code to process environmental data. Use with files in the folder “Environmental_databases”.


