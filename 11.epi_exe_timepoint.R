# compare EXE vs EPI and TE vs ICM
#
#


library(openxlsx)
library(tidyverse)
library(Mfuzz)
library(RColorBrewer)
library(ComplexHeatmap)
library(patchwork)
library(ggpubr)
library(ggrepel)


load("./mouse.uniprot.function.rda")

df.id <- fun$protein %>% mutate(ID = protein)

sample <- read.xlsx("Supplementary Table 1.xlsx",sheet = "FileInformation")
sample$category <- factor(sample$category,levels = unique(sample$category),ordered = T)






sample <- sample %>% arrange(category,time)
sample$sample <- factor(sample$sample,levels = sample$sample,ordered = T) 
sample$class <- factor(sample$class,levels = unique(sample$class),ordered = T) 

df.raw <- read.xlsx("Supplementary Table 1.xlsx",sheet = "RawIntensity")
df.nor <- read.xlsx("Supplementary Table 1.xlsx",sheet = "NormalizedIntensity") 
df.scale <-  read.xlsx("Supplementary Table 1.xlsx",sheet = "ZscoreInteisty") 

df.q <- df.nor %>%  gather(sample,int,-ID) 
df.q <-  sample %>% left_join(df.q,by = "sample")
df.q$sample <- factor(df.q$sample,levels = unique(df.q$sample),ordered = T)

df.sc <- df.scale %>%  gather(sample,int,-ID) 
df.sc <- sample %>% left_join(df.sc,by = "sample")
df.sc$sample <- factor(df.sc$sample,levels = unique(df.sc$sample),ordered = T)




df.sc.mean <- df.sc %>% group_by(ID,class) %>% summarise(scale = mean(int))
df.scale.mean <- data.frame(class =unique(sample$class)) %>% left_join(df.sc.mean)
df.scale.mean$class <- factor(df.scale.mean$class,levels = unique(df.scale.mean$class),ordered = T)


sample <- read.xlsx("Supplementary Table 1.xlsx",sheet = "FileInformation")
sample$category <- factor(sample$category,levels = unique(sample$category),ordered = T)


annot_sample <- sample %>% select(category,time,class)
row.names(annot_sample) <- sample$sample

# global color setting for sample class
sample2col <- sample %>% arrange(time)
col.samplecategory <- brewer.pal(3, "Set2")
names(col.samplecategory) <- unique(sample$category)

col.time <- brewer.pal(length(unique(sample$time)), "YlOrRd")
names(col.time) <- unique(sample$time)

col.sampleclass <- c("grey",brewer.pal(10, "Paired"))
names(col.sampleclass) <- c(unique(sample2col$class))
annotcol <- list(class = col.sampleclass,
                 category = col.samplecategory,
                 time = col.time)




#-------------------------------------
#
# all protein 2 kegg
#
#-------------------------------------

coma <- "ICM-EPI"
comb <- "TE-EXE"
comlist <- c("3.5d-TE:3.5d-ICM",
             "4.5d-TE:4.5d-ICM",
             "5.0d-EXE:5.0d-EPI",
             "5.5d-EXE:5.5d-EPI",
             "6.5d-EXE:6.5d-EPI"
             )



result <- data.frame()
result.scale <- data.frame()

fc <- 1.5
for(com in comlist) {
  
  coma <- strsplit(com,":")[[1]][1]
  comb <- strsplit(com,":")[[1]][2]
  to.stat <- df.scale %>% select(ID,sample$sample) %>% 
    gather(sample,intensity,-ID) %>% 
    left_join(sample)  %>%
    filter(class == coma | class == comb)  
  

  stat.p <- to.stat  %>%
    mutate(intensity = if_else(sample == unique(to.stat$sample)[1],intensity +0.1,intensity) )%>% 
    group_by(ID) %>% summarise(pvalue.t.test = t.test(intensity ~ class)$p.value,
                               pvalue.w.test = wilcox.test(intensity ~ class)$p.value) %>% 
    mutate(qvalue.t.test = p.adjust(pvalue.t.test,n = length(unique(to.stat$ID))),
           qvalue.w.test = p.adjust(pvalue.w.test,n = length(unique(to.stat$ID))),
           )
  
  stat.0 <- to.stat %>% group_by(ID,class) %>% summarise(mean = mean(intensity)) %>% 
    spread(class,mean) %>% select(ID,coma,comb)
  colnames(stat.0) <- c("ID","case","ctrl")
  stat <- stat.0 %>% mutate(ratio= case - ctrl) %>% left_join(stat.p)  %>% mutate(comparison = paste0(coma,"/",comb))
  
  
  stat.reg <- stat %>% mutate(regulated = if_else(pvalue.t.test < 0.05 & ratio >0,"TE-EXE",
                                                  if_else(pvalue.t.test < 0.05 & ratio < 0,"ICM-EPI","none")))
  
  
  
  stat.sig <- stat.reg%>% filter(regulated != "none")  %>% arrange(ratio)
  sig.pro <- stat.sig %>% select(ID) %>% 
    left_join(to.stat) %>% select(ID,sample,intensity) 
  
  sig.pro$ID <- factor(sig.pro$ID,levels = unique(sig.pro$ID),ordered = T)
  to.scale <- sig.pro %>% select(ID,sample,intensity) %>% spread(ID,intensity) 
  sig.scale <- to.scale %>% select(-sample) %>% scale() %>% bind_cols(sample = to.scale$sample) %>% 
    gather(ID,scale,-sample) %>% left_join(sample) %>% 
    select(ID,class,sample,scale) %>%  mutate(comparison = paste0(coma,"/",comb))
  
  result <- result %>% bind_rows(stat.reg)
  result.scale <- result.scale %>% bind_rows(sig.scale)
  
  
  
}


# volcano

result <- result %>%   mutate(comparison = str_replace_all(comparison,"d.+","d"))  %>% 
  left_join(fun$protein,by = c("ID"="protein")) 

wb <- createWorkbook()
df.sig <- result %>% 
  filter(regulated != "none")
addWorksheet(wb,"sig")
writeData(wb,"sig",df.sig)
addWorksheet(wb,"all")
writeData(wb,"all",result)
saveWorkbook(wb,file = "../result/3.two.stage/timepoint.comparison.xlsx",overwrite = T)


dbar <- result %>% filter(regulated != "none") %>% 
  mutate(log = ratio) %>% 
  group_by(comparison) %>% summarise(max = max(log),min= min(log)) %>% 
  select(comparison,max,min)


dsig <- result %>% filter(regulated != "none") %>% 
  mutate(log = ratio )%>% 
  mutate(lab = if_else(!is.na(tf_category),"TF",
                       if_else(!is.na(kinase_family),"Kinase","Others")))

dkey <- dsig %>% filter(lab != "Others")
dlab <- dkey %>% 
  group_by(comparison,regulated) %>% top_n(n = 20,wt = - pvalue.t.test)
  
ggplot() +
  geom_col(data = dbar,aes(comparison,y = max),fill = "#dcdcdc",alpha = 0.5,width = 0.7) +
  geom_col(data = dbar,aes(comparison,y = min),fill = "#dcdcdc",alpha = 0.5,width = 0.7) +
  geom_jitter(data = dsig,aes(comparison,log,size = -log10(pvalue.t.test)),width = 0.35,col = "grey",alpha = 0.75) +
  geom_jitter(data = dkey,aes(comparison,log,size = -log10(pvalue.t.test),col = lab),width = 0.35) +
  

  geom_tile(data = dbar,aes(comparison,y=0,fill = comparison),width = 0.7,height = 0.5,show.legend = F) +
  geom_text(data = dbar,aes(comparison,y=0,label = comparison),size =4,col = "white") +
  geom_text_repel(data = dlab,min.segment.length = Inf,
                  aes(comparison,y=log,label = genename,col = lab),size = 3) +
#  annotate(geom = "text",x = "3.5d",y = 12,label = "ICM-EPI") +
#  annotate(geom = "text",x = "3.5d",y = -12,label = "TE-EXE") +
  ggsci::scale_fill_npg() +
  scale_size(range = c(0.5,3)) +
  scale_color_brewer(palette = "Set1",direction = -1) +
  labs(y = "Average log2FC",size = "-log10(p)",col = "Type") +
  theme_minimal() +
  theme(legend.position = "top",
        panel.grid = element_blank(),
        axis.line.y = element_line(size = 0.5),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())
  

ggsave("../result/3.two.stage/two.leneage.vol.png",width = 8,height = 5)
ggsave("../result/3.two.stage/two.leneage.vol.pdf",width = 8,height = 5)


# venny

sig.list <- result %>% filter(regulated != "none") %>%  
  group_by(comparison) %>% summarise(sig.pro = list(ID))

sig.pro <- sig.list$sig.pro
names(sig.pro)<- sig.list$comparison

ggVennDiagram::ggVennDiagram(x = sig.pro,label = "count",label_alpha = 0) +
  scale_fill_distiller(palette = "Paired",direction = 1) +
  scale_x_continuous(expand = c(0.2,0.2)) +
  theme(legend.position = "top")

ggsave("../result/3.two.stage/two.leneage.overlap.png",width = 5,height = 5)



df.sig <- result %>% filter(regulated != "none")

ggplot(df.sig,aes(ratio,-log10(pvalue.t.test),col = regulated)) +
  geom_point() +
  facet_grid(comparison ~.,scales = "free_y")








# heatmap
for( i in 1:5 ){
  c <- unique(result.scale$comparison)[i]
  to.plot <- result.scale %>% filter(comparison == c)
  to.plot$ID <- factor(to.plot$ID,levels = unique(to.plot$ID),ordered = T)
  p <- ggplot(to.plot,aes(sample,ID,fill = scale))+
    geom_tile() +
    facet_wrap( comparison ~.,scales = "free") +
    scale_fill_gradient2(high = "orange",mid = "black",low = "blue4",midpoint = 0) +
    theme_void() +
    theme(legend.position = "left")
          #axis.text.x =element_text(angle = 45))
  
  if (i == 1) {
    p0 <- p
  } else {
    p0 <- p0/p
  }
}

ggsave(p0,file = "../result/3.two.stage/two.leneage.png",width = 3,height = 10)



##-------------
## lineage
##-------------

wb<- createWorkbook()

to.pre <- df.sig %>% 
  filter(comparison == "3.5d" | 
         comparison == "4.5d" #| 
           #comparison == "5.0d"  
           ) %>% 
  filter(regulated == "ICM-EPI")


df.pre.icp <- to.pre %>% 
  group_by(ID,genename) %>% summarise(n=n()) %>% filter(n>1) %>% mutate(type = "ICM-EPI") 



to.pre <- df.sig %>% 
  filter(comparison == "3.5d" | 
         comparison == "4.5d" #| 
         #comparison == "5.0d" 
         ) %>% 
  filter(regulated == "TE-EXE")

df.pre.te <- to.pre %>% 
  group_by(ID,genename) %>% summarise(n=n()) %>% filter(n>1)%>% mutate(type = "TE-EXE")



to.heat <- df.pre.icp %>% 
  bind_rows(df.pre.te) %>% 
  left_join(df.scale)  %>% 
  ungroup() %>% 
  filter(!is.na(genename)) %>% 
  left_join(df.id)

addWorksheet(wb,"ICM_TE")
writeData(wb,"ICM_TE",to.heat)


to.sample <- sample %>% filter(time == "3.5" |
                               time == "4.5" #| 
                               #time == "5.0"
                               ) %>% 
  arrange(category)
heat <- to.heat %>% select(to.sample$sample) %>% as.data.frame()
row.names(heat) <- to.heat$genename
#heat[heat>2] <-2

p <- pheatmap::pheatmap(heat,scale = "row",
         annotation_col = annot_sample,
         annotation_colors = annotcol,
         show_colnames = F,
         cluster_rows = T,
         cluster_cols = F)



ggsave(p,file = "../result/3.1.lineage_cell_type/ICM_TE_E3.4_E4.pdf",height = 15,width = 6)


p <- pheatmap::pheatmap(heat,scale = "row",
                        clustering_method = "ward.D",
                        annotation_col = annot_sample,
                        annotation_colors = annotcol,
                        treeheight_row = 10,
                        show_colnames = F,
                        show_rownames = F,
                        cluster_rows = T,
                        cluster_cols = F)



ggsave(p,file = "../result/3.1.lineage_cell_type/ICM_TE_E3.4_E4.png",height = 6,width = 6)




to.post.epi <- df.sig %>% 
  filter(comparison == "5.5d" | 
           comparison == "6.5d" #| 
         #comparison == "5.0d"  
  ) %>% 
  filter(regulated == "ICM-EPI")


df.post.epi <- to.post.epi %>% 
  group_by(ID,genename) %>% summarise(n=n()) %>% filter(n>1) %>% mutate(type = "ICM-EPI")



to.post.exe <- df.sig %>% 
  filter(comparison == "5.5d" | 
           comparison == "6.5d" #| 
         #comparison == "5.0d" 
  ) %>% 
  filter(regulated == "TE-EXE")

df.post.exe <- to.post.exe %>% 
  group_by(ID,genename) %>% summarise(n=n()) %>% filter(n>1)%>% mutate(type = "TE-EXE")



to.heat <- df.post.epi %>% 
  bind_rows(df.post.exe) %>% 
  left_join(df.scale)  %>% 
  ungroup() %>% 
  filter(!is.na(genename)) %>% 
  left_join(df.id)

addWorksheet(wb,"ICM_TE")
writeData(wb,"ICM_TE",to.heat)


to.sample <- sample %>% filter(time == "5.5" |
                                 time == "6.5" #| 
                               #time == "5.0"
) %>% 
  arrange(category)
heat <- to.heat %>% select(to.sample$sample) %>% as.data.frame()
row.names(heat) <- to.heat$genename
#heat[heat>2] <-2

p <- pheatmap::pheatmap(heat,scale = "row",
                        annotation_col = annot_sample,
                        annotation_colors = annotcol,
                        treeheight_row = 10,
                        show_colnames = F,
                        cluster_rows = T,
                        cluster_cols = F)


ggsave(p,file = "../result/3.1.lineage_cell_type/EPI_EXE_E5.5_E6.5.pdf",height = 15,width = 6)


p <- pheatmap::pheatmap(heat,scale = "row",
                        annotation_col = annot_sample,
                        annotation_colors = annotcol,
                        treeheight_row = 10,
                        show_colnames = F,
                        show_rownames = F,
                        cluster_rows = T,
                        cluster_cols = F)


ggsave(p,file = "../result/3.1.lineage_cell_type/EPI_EXE_E5.5_E6.5.png",height = 6,width = 6)


saveWorkbook(wb,"../result/3.1.lineage_cell_type/significant.xlsx")








#-------------------------------
#
# ICM -> EPI
#
#-------------------------------


 







#---------- plot pathway ------
library(tidytext)
library(ggsci)

count.com <- result %>% ungroup() %>% filter(regulated != "none") %>% 
  group_by(comparison,regulated) %>% summarise(count.com = n())

num.bg <- fun$bgnum %>% filter(database == "kegg" ) %>% select(bg) %>% as.numeric()
to.fun <- result %>% filter(regulated != "none") %>% 
  mutate(term = kegg_pathway) %>% 
  filter(!is.na(genename)) %>% 
  filter(!is.na(term)) %>% 
  separate_rows(term,sep = ";") %>% 
  separate(term ,into =c("l1","l2","l3"),sep = "\\|\\|") %>% 
  separate(l3,into = c("termid","term"),sep = "::") %>% 
  group_by(comparison,regulated,term) %>% summarise(count.pro = n(),sig.pro = paste(genename,collapse = ";")) %>% 
  left_join(count.com) %>% left_join(fun$term$kegg_pathway) %>% 
  mutate(count.term = bgcount,count.bg = num.bg) %>% 
  group_by(comparison,regulated,level1,level2,termid,term,sig.pro,count.pro,count.term,count.com,count.bg) %>% 
  summarise(pvalue = phyper(count.pro,count.term,count.bg - count.term,count.com,lower.tail = F,log.p = F)) 

to.fun$qvalue <-  p.adjust(to.fun$pvalue,n = nrow(to.fun),method = "BH")
out.fun <- to.fun %>% arrange(comparison,regulated,pvalue)
write.xlsx(out.fun,"../result/3.two.stage/to.pathway.xlsx",overwrite = T)













#--------------------------
# combine up and down
#--------------------------

to.plot<- openxlsx::read.xlsx("../result/3.two.stage/to.pathway.combine.xlsx",sheet = "select") %>% 
  mutate(time = str_replace_all(comparison,"d.+","d"))





to.label.1 <- to.plot %>% mutate(label=term,type = "2") %>% 
  mutate(label = str_replace(term,"endoplasmic reticulum","E.R."))
to.label.2 <- to.plot %>% mutate(pvalue = 1,label =sig.pro,type = "1")  %>% 
  mutate(label = str_replace(label,
                             "Myh11;Itgb1;Fn1;Cfl1;Tmsb4x;Ezr;Nckap1;Map2k1;Itgav;Arpc4;Actb;Mapk1;Rock2;Mapk3;Actn1;Actr3",
                             "Myh11;Itgb1;Fn1;Cfl1;Tmsb4x;Ezr;Nckap1;Map2k1;Itgav;Arpc4...")) %>% 
  mutate(label = str_replace(label,
                             "Myh11;Itgb1;Ezr;Nedd4;Cgn;Arpc4;Actb;Rock2;Myl6;Actn1;Dlg1;Map2k7;Actr3;Cldn3",
                             "Myh11;Itgb1;Ezr;Nedd4;Cgn;Arpc4;Actb;Rock2;Myl6;Actn1;Dlg1...")) %>% 
  mutate(label = str_replace(label,
                             "Rpsa;Rpl5;Rpl36;Rpl27;Rpl37a;Rps23;Rpl23a;Rps24;Rps28;Fau;Rpl30;Rps9;Rpl35;Rpl22l1",
                             "Rpsa;Rpl5;Rpl36;Rpl27;Rpl37a;Rps23;Rpl23a;Rps24;Rps28;Fau...")) %>% 
  
  mutate(label = str_replace(label,
                             "Rbm25;Dhx15;Srsf3;Isy1;Snrnp200;Hnrnpa3;Tcerg1;Thoc1;Prpf3;Wbp11;Srsf9;Bcas2;Acin1",
                             "Rbm25;Dhx15;Srsf3;Isy1;Snrnp200;Hnrnpa3;Tcerg1;Thoc1;Prpf3...")) %>%
  mutate(label = str_replace(label,
                             "Arid1a;Cox5a;Rps6ka1;Actb;Gnas;Smarcc1;Smarcd1;Kdm3a;Sdha;Atp5pb;Akt1s1;Coa3",
                             "Arid1a;Cox5a;Rps6ka1;Actb;Gnas;Smarcc1;Smarcd1;Kdm3a;Sdha...")) 

  


to.label <- to.label.1 %>% bind_rows(to.label.2)
#to.label$term <- factor(to.label$term,levels = unique(to.label$term),ordered = T)





col <- c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF")

col.i <- col[1]
i <- to.label %>% filter(time == "3.5d")
p1 <- 
  ggplot(i,aes(x = -log10(pvalue),
                   y = reorder_within(term,-pvalue,comparison)))+
  geom_bar(alpha = 0.6,stat = "identity",aes(fill = type),
           position = "dodge") +
  geom_text(aes(x= 0,label = label,col = type,size = type),hjust = 0,position = position_dodge(width = 0.9))+

  facet_wrap(time~., nrow = 1,scale = "free") +
    scale_size_manual(values = c(1.5,3)) +
  scale_fill_manual(values = c("grey",col.i)) +
  scale_color_manual(values = c(col.i,"black")) +
  labs(x = expression(-log[10]*'('*italic(P)*value*')')) +
  theme_bw() +
  theme(legend.position = "none",panel.grid = element_blank(),
        strip.background.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

col.i <- col[2]
i <- to.label %>% filter(time == "4.5d")
p2 <-
  ggplot(i,aes(x = -log10(pvalue),
               y = reorder_within(term,-pvalue,comparison)))+
  geom_bar(alpha = 0.6,stat = "identity",aes(fill = type),
           position = "dodge") +
  geom_text(aes(x= 0,label = label,col = type,size = type),hjust = 0,position = position_dodge(width = 0.9))+
  
  facet_wrap(time~., nrow = 1,scale = "free") +
  scale_size_manual(values = c(1.5,3)) +
  scale_fill_manual(values = c("grey",col.i)) +
  scale_color_manual(values = c(col.i,"black")) +
  labs(x = expression(-log[10]*'('*italic(P)*value*')')) +
  theme_bw() +
  theme(legend.position = "none",panel.grid = element_blank(),
        strip.background.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())



col.i <- col[3]
i <- to.label %>% filter(time == "5.0d")
p3 <- 
  ggplot(i,aes(x = -log10(pvalue),
               y = reorder_within(term,-pvalue,comparison)))+
  geom_bar(alpha = 0.6,stat = "identity",aes(fill = type),
           position = "dodge") +
  geom_text(aes(x= 0,label = label,col = type,size = type),hjust = 0,position = position_dodge(width = 0.9))+
  
  facet_wrap(time~., nrow = 1,scale = "free") +
  scale_size_manual(values = c(1.5,3)) +
  scale_fill_manual(values = c("grey",col.i)) +
  scale_color_manual(values = c(col.i,"black")) +
  labs(x = expression(-log[10]*'('*italic(P)*value*')')) +
  theme_bw() +
  theme(legend.position = "none",panel.grid = element_blank(),
        strip.background.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())


col.i <- col[4]
i <- to.label %>% filter(time == "5.5d") 
p4 <- 
  ggplot(i,aes(x = -log10(pvalue),
               y = reorder_within(term,-pvalue,comparison)))+
  geom_bar(alpha = 0.6,stat = "identity",aes(fill = type),
           position = "dodge") +
  geom_text(aes(x= 0,label = label,col = type,size = type),hjust = 0,position = position_dodge(width = 0.9))+
  
  facet_wrap(time~., nrow = 1,scale = "free") +
  scale_size_manual(values = c(1.5,3)) +
  scale_fill_manual(values = c("grey",col.i)) +
  scale_color_manual(values = c(col.i,"black")) +
  labs(x = expression(-log[10]*'('*italic(P)*value*')')) +
  theme_bw() +
  theme(legend.position = "none",panel.grid = element_blank(),
        strip.background.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())




col.i <- col[5]
i <- to.label %>% filter(time == "6.5d")
p5 <- 
  ggplot(i,aes(x = -log10(pvalue),
               y = reorder_within(term,-pvalue,comparison)))+
  geom_bar(alpha = 0.6,stat = "identity",aes(fill = type),
           position = "dodge") +
  geom_text(aes(x= 0,label = label,col = type,size = type),hjust = 0,position = position_dodge(width = 0.9))+
  
  facet_wrap(time~., nrow = 1,scale = "free") +
  scale_size_manual(values = c(1.5,3)) +
  scale_fill_manual(values = c("grey",col.i)) +
  scale_color_manual(values = c(col.i,"black")) +
  labs(x = expression(-log[10]*'('*italic(P)*value*')')) +
  theme_bw() +
  theme(legend.position = "none",panel.grid = element_blank(),
        strip.background.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

p1 +p2+p3+p4+p5 +plot_layout(nrow = 1)
ggsave("../result/3.two.stage/stage.function.kegg.png",width = 12,height = 3)




#---------------------------------
# separate up and down
#--------------------------------


to.plot<- openxlsx::read.xlsx("../result/3.two.stage/pathway.xlsx",sheet = "up_down") %>% 
  mutate(time = str_replace_all(comparison,"d.+","d"))





to.label.1 <- to.plot %>% mutate(label=term,type = "2") %>% 
  mutate(label = str_replace(label,"endoplasmic reticulum","E.R.")) %>% 
  mutate(label = str_replace(label,"metabolism","metab."))
to.label.2 <- to.plot %>% mutate(pvalue = 1,label =sig.pro,type = "1")  %>% 
  mutate(label = str_replace(label,
                             "Rpsa;Rpl36;Rpl27;Rpl37a;Rps23;Rpl23a;Rps24;Rps28;Fau;Rpl30;Rpl35;Rpl4",
                             "Rpsa;Rpl36;Rpl27;Rpl37a;Rps23;Rpl23a;Rps24;Rps28...")) %>% 
  mutate(label = str_replace(label,
                             "Rbm25;Dhx15;Srsf3;Hnrnpa3;Tcerg1;Thoc1;Prpf3;Wbp11;Bcas2;Acin1",
                             "Rbm25;Dhx15;Srsf3;Hnrnpa3;Tcerg1;Thoc1;Prpf3...")) %>% 
  mutate(label = str_replace(label,
                             "Capn2;Ufd1;Eif2s1;Ganab;Map2k7;Ubqln1;Dnajc3;Derl1;Dnajb1;Cul1",
                             "Capn2;Ufd1;Eif2s1;Ganab;Map2k7;Ubqln1;Dnajc3...")) %>% 
  mutate(label = str_replace(label,
                             "Arid1a;Cox5a;Actb;Gnas;Smarcd1;Kdm3a;Sdha;Atp5pb;Akt1s1",
                             "Arid1a;Cox5a;Actb;Gnas;Smarcd1;Kdm3a;Sdha;Atp5pb...")) %>% 
  mutate(label = str_replace(label,
                             "Ezr;Nedd4;Cgn;Arpc4;Actb;Actn1;Dlg1;Map2k7;Actr3;Cldn3", 
                             "Ezr;Nedd4;Cgn;Arpc4;Actb;Actn1;Dlg1;Map2k7;Actr3..." ))




to.label <- to.label.1 %>% bind_rows(to.label.2)
#to.label$term <- factor(to.label$term,levels = unique(to.label$term),ordered = T)





col <- c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF")

col.i <- col[1]
i <- to.label %>% filter(time == "3.5d")
p1 <- 
  ggplot(i,aes(x = -log10(pvalue),
               y = reorder_within(term,-pvalue,comparison)))+
  geom_bar(alpha = 0.6,stat = "identity",aes(fill = type),
           position = "dodge") +
  geom_text(aes(x= 0,label = label,col = type,size = type),hjust = 0,position = position_dodge(width = 0.9))+

  facet_grid(regulated ~ time,scale = "free_y",space = "free_y") +
  scale_size_manual(values = c(2,3)) +
  scale_fill_manual(values = c("grey",col.i)) +
  scale_color_manual(values = c(col.i,"black")) +
  labs(x = expression(-log[10]*'('*italic(P)*value*')')) +
  theme_bw() +
  theme(legend.position = "none",panel.grid = element_blank(),
        strip.background.x = element_blank(),
        strip.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

col.i <- col[2]
i <- to.label %>% filter(time == "4.5d")
p2 <- 
  ggplot(i,aes(x = -log10(pvalue),
               y = reorder_within(term,-pvalue,comparison)))+
  geom_bar(alpha = 0.6,stat = "identity",aes(fill = type),
           position = "dodge") +
  geom_text(aes(x= 0,label = label,col = type,size = type),hjust = 0,position = position_dodge(width = 0.9))+
  
  facet_grid(regulated ~ time,scale = "free_y",space = "free_y") +
  scale_size_manual(values = c(2,3)) +
  scale_fill_manual(values = c("grey",col.i)) +
  scale_color_manual(values = c(col.i,"black")) +
  labs(x = expression(-log[10]*'('*italic(P)*value*')')) +
  theme_bw() +
  theme(legend.position = "none",panel.grid = element_blank(),
        strip.background.x = element_blank(),
        strip.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())


col.i <- col[3]
i <- to.label %>% filter(time == "5.0d")
p3 <- 
  ggplot(i,aes(x = -log10(pvalue),
               y = reorder_within(term,-pvalue,comparison)))+
  geom_bar(alpha = 0.6,stat = "identity",aes(fill = type),
           position = "dodge") +
  geom_text(aes(x= 0,label = label,col = type,size = type),hjust = 0,position = position_dodge(width = 0.9))+
  
  facet_grid(regulated ~ time,scale = "free_y",space = "free_y") +
  scale_size_manual(values = c(2,3)) +
  scale_fill_manual(values = c("grey",col.i)) +
  scale_color_manual(values = c(col.i,"black")) +
  labs(x = expression(-log[10]*'('*italic(P)*value*')')) +
  theme_bw() +
  theme(legend.position = "none",panel.grid = element_blank(),
        strip.background.x = element_blank(),
        strip.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())




col.i <- col[4]
i <- to.label %>% filter(time == "5.5d")
p4 <- 
  ggplot(i,aes(x = -log10(pvalue),
               y = reorder_within(term,-pvalue,comparison)))+
  geom_bar(alpha = 0.6,stat = "identity",aes(fill = type),
           position = "dodge") +
  geom_text(aes(x= 0,label = label,col = type,size = type),hjust = 0,position = position_dodge(width = 0.9))+
  
  facet_grid(regulated ~ time,scale = "free_y",space = "free_y") +
  scale_size_manual(values = c(2,3)) +
  scale_fill_manual(values = c("grey",col.i)) +
  scale_color_manual(values = c(col.i,"black")) +
  labs(x = expression(-log[10]*'('*italic(P)*value*')')) +
  theme_bw() +
  theme(legend.position = "none",panel.grid = element_blank(),
        strip.background.x = element_blank(),
        strip.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())




col.i <- col[5]
i <- to.label %>% filter(time == "6.5d")
p5 <- 
  ggplot(i,aes(x = -log10(pvalue),
               y = reorder_within(term,-pvalue,comparison)))+
  geom_bar(alpha = 0.6,stat = "identity",aes(fill = type),
           position = "dodge") +
  geom_text(aes(x= 0,label = label,col = type,size = type),hjust = 0,position = position_dodge(width = 0.9))+
  
  facet_grid(regulated ~ time,scale = "free_y",space = "free_y") +
  scale_size_manual(values = c(2,3)) +
  scale_fill_manual(values = c("grey",col.i)) +
  scale_color_manual(values = c(col.i,"black")) +
  labs(x = expression(-log[10]*'('*italic(P)*value*')')) +
  theme_bw() +
  theme(legend.position = "none",panel.grid = element_blank(),
        strip.background = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

p1 + p2 + p3 + p4 + p5 +plot_layout(nrow = 1)
ggsave("../result/3.two.stage/stage.function.kegg.up.down.png",width = 12,height = 4)
ggsave("../result/3.two.stage/stage.function.kegg.up.down.pdf",width = 12,height = 4)
