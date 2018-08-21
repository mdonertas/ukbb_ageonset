library(tidyverse)
library(igraph)
library(GGally)
library(intergraph)

disCoding <- read_tsv('/Volumes/research1/thornton/ukbb_ageonset/data/raw/ukbb/datacoding/coding6.tsv')
nms=setNames(disCoding$meaning,disCoding$node_id)
nms=c(nms,`0`='Top')
nmsdf=data.frame(node=unname(nms[as.character(disCoding$node_id)]),parent=unname(nms[as.character(disCoding$parent_id)]))
disTree <- graph_from_data_frame(select(nmsdf,parent,node),directed = T)

data.frame(id=c('Top',sapply(all_simple_paths(disTree,'Top',V(disTree)$name),function(nmx)paste(nmx$name,collapse='.'))),
           value=1)%>%
        arrange(id)%>%
        write_csv('~/Desktop/diseases.csv')

cardnodes=subcomponent(disTree,'cardiovascular',mode = 'out')$name
cardTree=induced_subgraph(disTree,cardnodes)

data.frame(id=c('cardiovascular',sapply(all_simple_paths(cardTree,'cardiovascular',V(cardTree)$name),function(nmx)paste(nmx$name,collapse='.'))),
           value=1)%>%
        arrange(id)%>%
        write_csv('~/Desktop/cardDiseases.csv')
