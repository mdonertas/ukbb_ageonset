library(tidyverse)
library(ggpubr)
library(ggridges)
library(ggrepel)
library(RColorBrewer)
library(scales)
library(igraph)
theme_set(theme_pubr(base_size = 10, legend = 'bottom'))
pntnorm <- (1/0.352777778)
sexcolors <- setNames(c('rosybrown4', 'slategray'), c('Female', 'Male'))
disCoding <- read_tsv('./data/raw/ukbb/datacoding/coding6.tsv')
nodelabels=setNames(disCoding$meaning,disCoding$node_id)
nodelabels=c(nodelabels,`0`='Top')
disTreedf=data.frame(node=unname(nodelabels[as.character(disCoding$node_id)]),parent=unname(nodelabels[as.character(disCoding$parent_id)]))
disTree <- graph_from_data_frame(select(disTreedf,parent,node),directed = T)
disTreecl <- reshape::melt(sapply(neighbors(disTree,'Top')$name,function(nm){
  subcomponent(disTree,nm,'out')$name
}))%>%
  rename(node=value,
         cluster=L1)
disTreecl=setNames(disTreecl$cluster,disTreecl$node)
disTree_level1=neighbors(disTree,'Top','out')$name
discatcolors <- setNames(c(brewer.pal(8,'Set2'),brewer.pal(4,'Pastel1')),disTree_level1)
discatcolors['unclassifiable'] <- 'gray80' 
discatcolors['Top'] <- 'gray40' 
