library(dplyr)
library(tidyr)
raw_dataframe = read.delim('picoplankton_2015_2020_fromSHARK.txt', dec= ',', encoding = 'UTF-8')


dataframe = raw_dataframe %>%
  select(c('Stationsnamn', 'Nationellt.övervakningsstations.ID','Provtagningsdatum..start.','Vetenskapligt.namn',
           'Provets.nedre.djup','Parameter', 'Mätvärde', 'Mätenhet')) %>%
  mutate(Sample_name = paste0(Nationellt.övervakningsstations.ID,'_', Provtagningsdatum..start.)) %>%
  select(-c('Nationellt.övervakningsstations.ID', 'Provtagningsdatum..start.')) %>%
  mutate(Parameter = paste0(Parameter, ' (',Mätenhet,')')) %>%
  mutate(Stationsnamn = if_else(Stationsnamn == 'RA1', 'RÅNEÅ-1', Stationsnamn)) %>%
  filter(Provets.nedre.djup <= 10) %>%
  filter(!(Stationsnamn == 'RÅNEÅ-1' & Provets.nedre.djup > 5)) %>%
  select(-c('Provets.nedre.djup', 'Stationsnamn')) %>%
  group_by(Parameter, Sample_name) %>%
  summarize(Average_Mätvärde = mean(Mätvärde)) %>%
  ungroup() %>%
  pivot_wider(names_from = Sample_name, values_from = Average_Mätvärde)

write.table(dataframe, "picoplankton_processed.tsv", sep = "\t", quote = FALSE, row.names = TRUE)


# add method in the name of the parameter