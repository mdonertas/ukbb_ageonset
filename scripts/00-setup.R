library(tidyverse)
library(ggforce)
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
mhcchr=6
mhcstart=28477797
mhcend=33448354

# ageonsetcolors = setNames(rev(brewer.pal(4,'Oranges')),1:4)
ageonsetcolors=setNames(c('#D55E00','#0072B2','#F0E442','#CC79A7','#E69F00','#56B4E9','gray80',"#C6E2FF"),c('1','2','3','1-2','1-3','2-3','1-2-3','4'))
