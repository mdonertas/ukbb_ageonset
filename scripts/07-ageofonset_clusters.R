library(tidyverse)
library(RFRlib)
library(dendextend)
library(RColorBrewer)
library(ggdendro)
library(factoextra)
library(jpeg)
theme_set(theme_rfr())

SRdisease <- readRDS('./data/processed/traits_clean/SRdisease_baseline_propagated.rds')
traits <- readRDS('./data/processed/traits_clean/traitData_baseline.rds')
traits <- traits %>%
  mutate(Sex = c('Female','Male')[Sex+1])%>%
  mutate(BMI = Weight/((`Standing height`/100)^2))
# prevDF <- readRDS('./data/processed/traits_clean/SRdisease_prevdf.rds')
disSet <- readRDS('./data/processed/traits_clean/SRdiseaseSet.rds')

SRdisease <- SRdisease %>% 
  filter(Disease %in% disSet$Disease)

SRdisease <- traits %>%
  select(eid, Sex) %>%
  right_join(SRdisease)

totagex=sapply(seq(0,70,by=1),function(x){
  sum(traits$`Age when attended assessment centre`>=x)
})

agedist <- tapply(SRdisease$Age,SRdisease$Disease,FUN = function(x){
  # predict(smooth.spline(table(cut(x[complete.cases(x)],seq(0,71,by=1)))/totagex),seq(0,71,by=0.1))$y
  predict(smooth.spline(density(x[complete.cases(x)])),seq(0,71,by=0.1))$y
},simplify = T)

agedist=t(sapply(agedist,c))

hc <- hclust(dist(agedist))

dend <- as.dendrogram(hc)
owd=getwd()
system('mkdir -p ./results/UKBB_disease_EDA/ageonset')
for(hcut in 2:25){
  system(paste('mkdir -p /nfs/research1/thornton/ukbb_ageonset/results/UKBB_disease_EDA/ageonset/',hcut,'cluster',sep=''))
  setwd(paste('/nfs/research1/thornton/ukbb_ageonset/results/UKBB_disease_EDA/ageonset/',hcut,'cluster',sep=''))
  print(hcut)
  if(hcut<=1){
    discl=cutree(hc,h=hcut)
  } else if(hcut>1){
    discl=cutree(hc,k=hcut)
  }
  
  par(mar=c(20, 4, 4, 2) + 0.1)
  colsx=c(brewer.pal(8,'Dark2'),brewer.pal(8,'Set2'),brewer.pal(8,'Set1'))[1:(luniq(discl))]
  
  dend <- dend %>% 
    set("branches_k_color", k = luniq(discl),
        value = colsx) %>% 
    set("labels_colors", k = luniq(discl),
        value = colsx) %>% 
    set("branches_lwd", 2)
  agedf <- lapply(unique(discl[hc$order]),function(i){
    dismat <- agedist[names(which(discl==i)),]
    data.frame(mnsx=if(class(dismat)=='matrix'){colMeans(dismat)}else{dismat},
               sdx=if(class(dismat)=='matrix'){apply(dismat,2,sd)}else{rep(0,length(seq(0,71,by=0.1)))},
               cl=i,
               numDiseases=luniq(names(which(discl==i))),
               age=seq(0,71,0.1))
  })
  
  agedf=reshape2::melt(agedf,id.vars=colnames(agedf[[1]]))%>%
    rename(clorder=L1)
  names(colsx)=unique(discl[hc$order])
  agedistplots <- lapply(unique(agedf$clorder),function(i){
    filter(agedf,clorder==i) %>%
      ggplot(aes(x=age))+
      geom_segment(aes(xend=age,y=mnsx-sdx,yend=mnsx+sdx),color='gray70')+
      geom_line(aes(y=mnsx),color='black',size=2)+
      theme_minimal()+
      xlab('')+ylab('')+
      theme(axis.text.y = element_blank(),
            panel.grid.major.x = element_line(color='gray70',size = 1),
            panel.grid.major.y = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color=colsx[i],linetype='dashed',size=2,fill=NA))+
      scale_x_continuous(breaks = seq(0,70,by=10),limits = c(0,71))
  })
  
  midpoints=sapply(unique(discl[hc$order]),function(i){
    mean(unname(which(discl[hc$order]==i)))
  })
  
  sapply(1:length(agedistplots),function(i){
    ggsave(paste(i,'.jpeg',sep=''),
           agedistplots[[i]],device = 'jpeg',width = 3,height = 2)
  })
  
  myimgdat=data.frame(midpoints=midpoints,
                      xstart=midpoints-3,
                      xend=midpoints+3,
                      ystart=rep(-0.65,luniq(discl)),
                      yend=rep(-0.85,luniq(discl)),
                      filename=paste(1:luniq(discl),'.jpeg',sep=''))
  
  mynewimgdat=myimgdat[1,]
  for(i in 2:nrow(myimgdat)){
    if((myimgdat$xstart[i]<mynewimgdat$xend[i-1]) & (myimgdat$ystart[i]==mynewimgdat$ystart[i-1])){
      newadd=myimgdat[i,]
      newadd$ystart=newadd$ystart-0.2
      newadd$yend=newadd$yend-0.2
      mynewimgdat=rbind(mynewimgdat,newadd)
    } else{
      mynewimgdat=rbind(mynewimgdat,myimgdat[i,])
    }
  }
  myimgdat=mynewimgdat
  
  p=fviz_dend(hc, k = luniq(discl), # Cut in four groups
              cex = 1.1, # label size
              k_colors = colsx,
              color_labels_by_k = T,
              rect = TRUE, # Add rectangle around groups
              rect_border = colsx,
              rect_fill = F)+
    ylim(-1,1)+
    theme(axis.text = element_blank(),
          axis.ticks = element_blank())+
    ylab('')+
    ggtitle('')
  
  for(i in 1:nrow(myimgdat)){
    p = p +
      annotation_raster(readJPEG(as.character(myimgdat$filename[i])),
                        xmin = myimgdat$xstart[i],
                        xmax = myimgdat$xend[i],
                        ymin = myimgdat$ystart[i],
                        ymax = myimgdat$yend[i],interpolate = T)
  }
  ggsave(paste(hcut,'.pdf',sep=''),p,width = 35,height = 18)
}
setwd(owd)
