setwd("C:/Users/krzysztof.jurdzins/OneDrive - KTH/Skrivbordet/ENVGEN_ndb/EnvPredict_github/EnvPredict/code")


library(tidyverse)

## Define infiles

metadata_file = "../env_data/combined/physical_chemical_processed.tsv"
long_phys_chem_data_file = "../env_data/combined/physical_chemical_data_2019-2020_2022-03-16.txt"
some_shark_data_file = "../env_data/combined/zooplankton_2015_2020_fromSHARK.txt"

## Read the data

metadata = read.table(metadata_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
long_phys_chem = read.table(long_phys_chem_data_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
some_shark = read.table(some_shark_data_file, sep = "\t",
                        header = TRUE, stringsAsFactors = FALSE)

## Get the station ids

stations = unique(long_phys_chem$station_name)

setdiff(unique(metadata$station_name), stations)

station_ids = c()
for (station in stations) {
  ix = which(long_phys_chem$station_name == station)[1]
  id = long_phys_chem$station_id[ix]
  station_ids = c(station_ids, id)
}
names(station_ids) = stations

## Get BY2 ARKONA station id - missing in physical_chemical data

ix = which(metadata$station_name == "BY2 ARKONA")
id = str_split(metadata$sample_id[ix], '_')[[length(ix)]][2]
station_ids[names(station_ids) == "BY2 ARKONA"] = id

ix = match(metadata$station_name, names(station_ids))
metadata$station_id = station_ids[ix]

## Check that the same BY2 Arkona code is used in an SHARK dataset
some_shark$station_id[some_shark$station_name == 'BY2 ARKONA'][1] == id

## Reorder metadata

metadata = metadata[,c(1,2, ncol(metadata), 3:(ncol(metadata)-1))]

## Get the station_id_date

metadata$station_id_date = paste(metadata$station_id, gsub('-', '', metadata$date), sep = "_")

## Reorder
iy_date = which(colnames(metadata) == "date")
ncol(metadata)
metadata = metadata[,c(1:iy_date, ncol(metadata), (iy_date+1):(ncol(metadata)-1))]

## Save

write.table(metadata,
            "../env_data/combined/physical_chemical_processed_translation.tsv",
            sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)

