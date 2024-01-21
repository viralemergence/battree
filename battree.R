
library(tidyverse)
library(vroom)
library(ggtreeExtra)
library(MoMAColors)
library(MetBrewer)

setwd("~/Documents/Github/")

##### 1. Get component datasets

vir <- vroom("virion/virion/virion.csv.gz")
vir %>% filter(HostOrder == 'chiroptera') %>% pull(Host) %>% unique() -> data

steph <- c("Artibeus jamaicensis", "Carollia perspicillata", "Eidolon helvum", "Hypsignathus monstrosus", "Rhinolophus alcyone", "Rousettus aegyptiacus")

dan <- read_csv("battree/Dan.csv") %>% pull(Species) %>% unique()
dan2 <- c("Rhinolophus mehelyi",
          "Rhinolophus euryale",
          "Rhinolophus ferrumequinum",
          "Rhinolophus hipposiderus",
          "Miniopterus schreibersii",
          "Plecotus austriacus",
          "Myotis myotis",
          "Myotis blythii",
          "Myotis emarginatus",
          "Myotis bechsteinii",
          "Myotis capaccinii",
          "Pteropus alecto",
          "Pipistrellus bodenheimeri")
dan <- c(dan, dan2) %>% unique()

adam <- read_csv("battree/Adam.csv") %>% pull(Taxon) %>% unique()

##### 2. Get phylogeny and build the data frame itself

tree <- ape::read.nexus("batgap/01_data processing/MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre")
taxa <- read.csv("batgap/01_data processing/taxonomy_mamPhy_5911species.csv",header=T) %>% filter(ord == "CHIROPTERA", `extinct.` == 0) 
tree <- keep.tip(tree, taxa$tiplabel)

df <- taxa %>% select(Species_Name, tiplabel) %>% 
  mutate(species = gsub("_", " ", Species_Name)) %>% select(species, tiplabel)

df %>% 
  mutate(`Databases` = as.numeric(str_to_lower(species) %in% data),
         `Field sampling` = as.numeric(species %in% dan),
         `Museum sampling` = as.numeric(species %in% adam),
         `Cell lines` = as.numeric(species %in% steph)) %>%
  pivot_longer(cols = c('Databases', 'Field sampling', 'Museum sampling', 'Cell lines'), names_to = 'Data sources') %>%
  select(-species) -> tipdf

tipdf %>% mutate(`Data sources` = factor(`Data sources`, 
                                         levels = c('Databases', 'Field sampling', 'Museum sampling', 'Cell lines'))) -> tipdf

vir %>% filter(Host == 'homo sapiens') %>% pull(Virus) %>% unique() -> zoo
 vir %>%
   filter(HostOrder == 'chiroptera') %>%
   select(Host, Virus) %>% distinct() %>%
   mutate(iszoo = as.numeric(Virus %in% zoo)) %>% 
   group_by(Host) %>% summarize(nzoo = sum(iszoo)) -> nzoo

df %>% 
  mutate(species = str_to_lower(species)) %>%
  dplyr::left_join(nzoo, by = join_by('species' == 'Host')) %>%
  mutate(nzoo = replace_na(nzoo, 0)) -> dfvir

dfvir %>% mutate(ID = tiplabel) %>% select(ID, nzoo) %>%
  mutate(any = as.factor(nzoo>0)) %>% mutate(any = replace_na(any, 'FALSE')) %>%
  rename(Zoonoses = nzoo) -> dfvir

##### 3. Make a figgy

p <- ggtree(tree, layout="fan", size=0.15, open.angle=5); p

p %<+% dfvir + geom_point(
  mapping=aes(size = Zoonoses), alpha = 0.3, color = MetBrewer::met.brewer("Hokusai2")[4]) + 
  #scale_alpha_discrete(range = c(0, 0.5), guide = 'none') + 
  scale_size_continuous(range = c(-1, 8), breaks = c(0,1,10,20)) -> p; p

p + 
  geom_fruit(data = tipdf, geom = geom_tile,
             mapping=aes(y = tiplabel, x = `Data sources`, alpha = value, fill = `Data sources`),
             color = "white", offset = 0.04, size = 0.005, pwidth = 0.3) +
  scale_alpha_continuous(range=c(0, 1),guide = 'none') +
  scale_fill_manual(values=rev(MetBrewer::met.brewer("Egypt")[1:4])) -> p; p 

ggsave(filename = "BatTree.jpg", path = "~/Documents/Github/BatTree", width = 9, height = 6, device='jpg', dpi=700)
