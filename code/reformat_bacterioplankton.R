library(dplyr)
library(tidyr)
raw_dataframe = read.delim('bacterioplankton_2015_2020_fromSHARK.txt', dec= ',')


dataframe = raw_dataframe %>%
  select(c('Stationsnamn', 'Nationellt.övervakningsstations.ID','Provtagningsdatum..start.',
           'Provets.övre.djup','Parameter', 'Mätvärde', 'Mätenhet')) %>%
  mutate(Sample_name = paste0(Nationellt.övervakningsstations.ID,'_', Provtagningsdatum..start.)) %>%
  select(-c('Nationellt.övervakningsstations.ID', 'Provtagningsdatum..start.')) %>%
  mutate(Parameter = paste0(Parameter, ' (',Mätenhet,')')) %>%
  mutate(Stationsnamn = if_else(Stationsnamn == 'RA1', 'RÅNEÅ-1', Stationsnamn)) %>%
  filter(Provets.övre.djup <= 10) %>%
  filter(!(Stationsnamn == 'RÅNEÅ-1' & Provets.övre.djup > 5)) %>%
  select(-c('Provets.övre.djup', 'Stationsnamn')) %>%
  group_by(Parameter, Sample_name) %>%
  summarize(Average_Mätvärde = mean(Mätvärde)) %>%
  ungroup() %>%
  pivot_wider(names_from = Sample_name, values_from = Average_Mätvärde)

write.table(dataframe, "bacterioplankton_processed.tsv", sep = "\t", quote = FALSE, row.names = TRUE)
