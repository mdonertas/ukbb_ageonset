library(tidyverse)
library(igraph)

disCoding <- read_tsv('./data/raw/ukbb/datacoding/coding6.tsv')
disTree <- graph_from_data_frame(select(disCoding,parent_id,node_id),directed = T)

cardnodes=subcomponent(disTree,'1071',mode = 'out')$name
cardTree=induced_subgraph(disTree,cardnodes)
nms=unname(setNames(disCoding$meaning,disCoding$node_id)[V(cardTree)$name])
V(cardTree)$desc=nms
pdf('./results/UKBB_disease_EDA/cardioTree.pdf',width=30,height = 15)

plot(cardTree, 
     layout=layout_as_tree,
     vertex.size=1.5,
     vertex.frame.color=NA,
     vertex.label=nms,
     vertex.label.cex = .75,             # node label size
     vertex.label.family = "Helvetica", # node label family
     vertex.label.color = 'gray25',    
     asp=0.5, 
     edge.arrow.size=0.2,
     edge.color='gray80',
     vertex.color='gray70')
dev.off()
