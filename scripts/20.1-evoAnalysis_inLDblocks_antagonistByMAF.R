antaplot = list()
k=0
signifSNPs = signifSNPs %>%
  mutate(snppos = paste(CHR,BP,sep='_'))
mafxx = signifSNPs %>% 
  mutate(betagr = c(1:10/10)[cut(abs(signifSNPs$BETA),breaks = quantile(abs(signifSNPs$BETA),probs = seq(0,1,length.out = 11)),include.lowest = T)])
# mafxx = signifSNPs %>%
#   mutate(MAF = ifelse(UKBB_RAF<=0.5,UKBB_RAF,1-UKBB_RAF)) %>% group_by(snppos) %>% summarise(MAF = min(MAF))
for(mafx in c(10:1/10)){
  k=k+1
  print(mafx)
  ld_cl1cl2_anta = apply(LDblocks,1,function(x){
    st=x['start']
    en=x['stop']
    ch=x['CHR']
    snpsexc = unique(filter(mafxx,betagr>=mafx)$snppos)
    xx = signifSNPs %>%
      select(-disID,-BETA,-P_BOLT_LMM_INF) %>%
      filter(CHR==ch & BP<=en & BP>=st) %>%
      filter(snppos %in% snpsexc) %>%
      filter(cl1_cl2 == ' antagonist' & !grepl('agonist',cl1_cl3) & !grepl('agonist', cl2_cl3)) %>%
      unique() %>%
      mutate(cluster = factor(cluster)) %>% 
      select(SNP,cluster,UKBB_RAF,ALL_RAF,AFR_RAF,EAS_RAF,AMR_RAF,SAS_RAF,EUR_RAF) %>% 
      gather(key = 'pop', value ='raf',-SNP,-cluster) %>%
      mutate(pop = factor(gsub('_RAF','',pop),levels = c('UKBB','ALL','AFR','AMR','EAS','EUR','SAS'))) %>%
      na.omit() 
    if(nrow(xx)>0){
      xx = xx %>% group_by(cluster,pop) %>% summarise(n=length(unique(SNP))) %>% ungroup() %>% right_join(xx)
      xx %>%
        group_by(cluster,pop,n) %>% 
        summarise(RAF_mean = mean(raf,na.rm=T),
                  RAF_min = min(raf,na.rm=T),
                  RAF_max = max(raf,na.rm=T),
                  RAF_sd = sd(raf,na.rm=T),
                  RAF_q1 = quantile(raf,probs = 0.25, na.rm=T),
                  RAF_q3 = quantile(raf,probs = 0.75, na.rm=T),
                  RAF_median = quantile(raf,probs = 0.75, na.rm=T)) %>%
        mutate(LDblock=x['block'])
    } else {NULL}
  }) 
  ld_cl1cl2_anta=ld_cl1cl2_anta[!sapply(ld_cl1cl2_anta,is.null)]
  ld_cl1cl2_anta = reshape2::melt(ld_cl1cl2_anta,id.vars = colnames(ld_cl1cl2_anta[[1]])) 
  ld_cl1cl2_anta = ld_cl1cl2_anta %>% group_by(cluster,pop) %>% summarise(nsum = sum(n),nmean=mean(n)) %>% ungroup() %>% right_join(ld_cl1cl2_anta)
  antaplot[[k]] = ld_cl1cl2_anta %>% 
    mutate(cluster=as.factor(cluster)) %>%
    filter(pop == 'UKBB') %>%
    ggplot(aes(y = RAF_median, x= cluster, fill = cluster)) +
    geom_violin()+
    geom_sina(size=1,alpha=1)+
    geom_boxplot(width=0.1,fill='white') +
    scale_fill_manual(values = ageonsetcolors) +
    stat_summary(fun.y = 'median', aes(label = round(..y..,2)), geom = "label",fill='white')  +
    guides(fill = F) +
    xlab('Age of onset cluster') + ylab('Risk Allele Frequency')+
    scale_y_continuous(limits=c(0,1.1),breaks = seq(0,1,by=0.2))+
    stat_compare_means(label = 'p.format',label.y=0.01,label.x = 0.5,hjust=0,size=6/pntnorm) +
    geom_text(data=unique(filter(ld_cl1cl2_anta,pop=='UKBB' & cluster==1)%>%select(cluster,nsum,nmean)),aes(x=0.5,label=paste('n_total=',nsum,', n_mean=',nmean,sep='')),y=1.1,size=6/pntnorm,hjust=0) +
    ggtitle(paste(mafx*100 - 10,' - ',mafx*100,' % BETA',sep=''))
  print(antaplot[[k]])
}
p1=ggarrange(antaplot[[1]],antaplot[[2]],antaplot[[3]],antaplot[[4]],antaplot[[5]],antaplot[[6]],antaplot[[7]],antaplot[[8]],antaplot[[9]],antaplot[[10]],ncol=5,nrow=2,labels='auto')
ggsave('./results/evoAnalysis/cl12_antagonist_maffiletrs.pdf',p1, width = 28,height = 14, units = 'cm',useDingbats=F)
ggsave('./results/evoAnalysis/cl12_antagonist_maffiletrs.png',p1, width = 28,height = 14, units = 'cm')
