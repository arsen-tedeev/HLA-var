#
#
# processing of files containing
# structural data
#

mir_atomdist <- read.table("./data/mir_files/mir_atomdist.txt", sep = "", header = TRUE)
head(mir_atomdist)

require(dplyr)
# removing columns with atom names
# and distance value
mir_atomdist <- select(mir_atomdist, -c(dist, atom.to, atom.from))

mir_general <- read.table("./data/mir_files/mir_general.txt", sep = "", header = TRUE, 
                          na.strings = NA, fill = TRUE)
head(mir_general)

mir_general %>% 
  select(pdb.id, chain.supertype) %>% 
  filter(chain.supertype == "MHCI" | chain.supertype == "MHCII") %>% 
  distinct() %>% 
  select(chain.supertype) %>% 
  table() %>% 
  prop.table() %>% 
  `*`(100) %>% 
  round(2) # the ration between MHC class I and MHC class II

unique(mir_general$complex.species) # data has both human and mouse MHC

length(unique(mir_general$pdb.id)) # 234 name in total

# removing PDB-structures with 
# peptide length equal to 0
mir_general <- mir_general %>% 
  filter(!(chain.supertype == "PEPTIDE" & seq.length == 0))

# leaving only alpha-chain from MHC class I
pdb_names <- mir_general %>% 
  filter(chain.supertype == "MHCI" & chain.type == "MHCa") %>% 
  select(1) %>% 
  pull()

length(pdb_names) # 161 name in total

mir_general %>% 
  filter(chain.supertype == "MHCI") %>% 
  select(1) %>% 
  pull() %>% 
  unique() %>% 
  length() # 161 > 173 => there were some PDB-structures of MHC class I 
# that did not have alpha-chain (already filtrated)

mir_general %>% 
  group_by(pdb.id) %>% 
  select(pdb.id, chain.id) %>% 
  distinct() %>%
  summarise(n = n()) %>% 
  select(n) %>% 
  unique() # there are some PDB-structure that have 4 chains

mir_resmarkup <- read.table("./data/mir_files/mir_resmarkup.txt", sep = "", header = TRUE, fill = TRUE)

head(mir_resmarkup)

mir_resmarkup <- mir_resmarkup %>% 
  select(-residue.aa) # removing of redundant column

# changing column order
mir_atomdist_inv <- mir_atomdist %>% 
  relocate(residue.index.to, residue.index.from, .after = pdb.id) %>% 
  relocate(chain.id.to, chain.id.from, .after = pdb.id) 

# binding original mir_atomdist and the one with 
# changed column order together, leaving only MHC class I 
atomdist_bind <- rbind(mir_atomdist, mir_atomdist_inv) %>% 
  filter(pdb.id %in% pdb_names)

# filtration of the rest of data
mir_general <- mir_general %>% 
  filter(pdb.id %in% pdb_names)
mir_resmarkup <- mir_resmarkup %>% 
  filter(pdb.id %in% pdb_names)

# after filtration each PDB-structure have 5 chains
mir_general %>% 
  group_by(pdb.id) %>% 
  select(pdb.id, chain.id) %>% 
  distinct() %>%
  summarise(n = n()) %>% 
  select(n) %>% 
  unique()

# sample is greatly skewed towards HLA-A2 allele
mir_general %>% 
  filter(chain.supertype == "MHCI" & chain.type == "MHCa") %>% 
  select(allele.info) %>% 
  table() %>% 
  prop.table() %>% 
  `*`(100) %>% 
  round(2)

general.from <- mir_general %>% 
  select(pdb.id, chain.type, chain.id) %>% 
  rename(chain.type.from = chain.type, chain.id.from = chain.id)

general.to <- mir_general %>% 
  select(pdb.id, chain.type, chain.id) %>% 
  rename(chain.type.to = chain.type, chain.id.to = chain.id)

# assigning biological names according to chain ID 
# and leaving only contacts that belong to alpha-chain of MHC
assign.data <- atomdist_bind %>% 
  merge(general.to, by = c("pdb.id", "chain.id.to")) %>% 
  merge(general.from, by = c("pdb.id", "chain.id.from")) %>% 
  relocate(chain.type.from, chain.type.to, .after = pdb.id) %>% 
  filter(!(chain.type.from == "MHCb" | chain.type.to == "MHCb")) %>% 
  filter(chain.type.from == "MHCa") %>% 
  merge(distinct(select(mir_general, pdb.id, complex.species)), by.x = "pdb.id") # adding organism information

head(assign.data)

# eliminating unnecessary dataframes  
rm(list = c("general.to", "general.from", "mir_atomdist_inv", "atomdist_bind"))

mir_general <- mir_general %>% 
  select(pdb.id, complex.species, chain.type, chain.id) %>% 
  filter(chain.type == "MHCa") %>% 
  relocate(chain.id, chain.type, .after = pdb.id)

# dataframe containing alpha-chain residue index information 
# for each PDB-structure  
index.data <- merge(mir_general, mir_resmarkup, by.x = c("pdb.id", "chain.id")) %>% 
  arrange(pdb.id, residue.index, region.type)

head(index.data)

# alpha-chain sequence "assembly"
hla.pdb.seq <- index.data %>% 
  group_by(pdb.id) %>% 
  summarise(sequence = paste(residue.ins.code, collapse = "")) %>% 
  merge(distinct(select(mir_general, pdb.id, complex.species)), by.x = "pdb.id") # adding organism information

head(hla.pdb.seq)

# vector of chain lengths
pdb_len <- lapply(hla.pdb.seq$sequence, nchar) %>% 
  unlist()

require(ggplot2)
# plot depicting distribution of alpha-chain lengths
ggplot(mapping = aes(pdb_len))+
  geom_histogram(binwidth = 1, fill = "grey", color = "black")+
  theme_bw()+
  labs(x = "Length", y = "Number of chains", title = "\u03b1-chain length distribution (PDB-structures)")

ggsave("pdb_len_distr.png", path = "./output/plots", width = 4.2, height = 3, units = "in")

# elimination of mouse alleles
dfList <- list(index.data_1 = index.data_1, hla.pdb.seq_1 = hla.pdb.seq_1) %>% 
  lapply(function(x) {
    x <- x %>% 
      filter(!(complex.species == "Mouse")) %>% 
      select(-complex.species)
    return(x)
  }) %>% 
  list2env(envir = .GlobalEnv) # writing back to the global environment
