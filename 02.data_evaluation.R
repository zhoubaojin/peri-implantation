library(tidyverse)
library(openxlsx)
library(RColorBrewer)
library(ggpubr)
library(scatterpie)


load("./mouse.uniprot.function.rda")

sample <- read.xlsx("Supplementary Table 1.xlsx",sheet = "FileInformation") %>% 
  arrange(time)
sample$class <- factor(sample$class,levels = unique(sample$class),ordered = T)
col.samplecategory <- brewer.pal(3, "Set2")
names(col.samplecategory) <- unique(sample$category)

col.time <- brewer.pal(length(unique(sample$time)), "YlOrRd")
names(col.time) <- unique(sample$time)

col.sampleclass <- c("grey",brewer.pal(10, "Paired"))
names(col.sampleclass) <- c(unique(sample$class))
annotcol <- list(class = col.sampleclass,
                 category = col.samplecategory,
                 time = col.time)





#------------------------------------------------------------
#
# id
#
#-------------------------------------------------------------


df.q.t <- read.xlsx("Supplementary Table 1.xlsx",sheet = "RawIntensity") %>% 
  gather(sample,int,-ID) %>% 
  filter(int > 0)%>% 
  left_join(sample)

df.count <- df.q.t %>%
  group_by(category,time,class,sample) %>% summarise(n =n()) %>% arrange(category,time)
df.count$sample <- factor(df.count$sample,levels = unique(df.count$sample),ordered = T)
df.count$category <- factor(df.count$category,levels = unique(df.count$category),ordered = T)


df.count.mean <- df.count %>% group_by(class) %>% summarise(mean = mean(n),sd = sd(n))

#pid <-
ggplot(df.count.mean,aes(class,mean,fill= class)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=.2) +
  annotate(geom = "text",y = 1200,x = "5.5d-EPI",label = "average: 895",hjust = 0)+
  geom_point(data = df.count,aes(class,n,fill = class),shape = 21,size =3) +
  geom_hline(yintercept = mean(df.count$n),lty  =3) +
  scale_y_continuous(labels = scales::comma,limits = c(0,1800)) +
  scale_fill_manual(values = col.sampleclass) +
  scale_color_manual(values = col.sampleclass) +
  theme_bw() +
  labs(y = "No. protein groups")  +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 60,hjust =1),
        axis.title.x = element_blank())

ggsave("../result/1.global/1.id.pdf",width =4.5,height = 3.5)  
ggsave("../result/1.global/1.id.png",width =4.5,height = 3.5)  

#--------------------------------------
#
# overlap
#
#--------------------------------------

for(c in unique(sample$class)) {
  df.1 <- df.q.t %>% filter(class == c) %>% 
    group_by(sample)%>% summarise(protein = list(ID))
  x <- df.1$protein
  names(x) <- df.1$sample
  p <- ggVennDiagram::ggVennDiagram(x = x) +
    scale_x_continuous(expand = c(0.1,0.1)) +
    scale_fill_gradient(low = "white",high = "white") +
    theme(legend.position = "none")
  ggsave(p,file = paste0("../result/1.global/1.id.overlap.",c,".png"),width = 3,height = 3)
}




for(i in 1:11) {
  message(i) 
  c.i <- unique(sample$class)[i]
  df.i <- df.q.t %>% filter(class == as.character(c.i))  %>% 
    group_by(ID,class) %>% summarise(mean = log10(mean(int))) %>% arrange(desc(mean)) 
  df.i$rank <- 1:nrow(df.i)
  if (i == 1) {
    p <- ggplot(df.i,aes(rank,mean)) +
      annotate("segment",x = -50 + 50 * (i -1) , y = 4 + 1 *(i-1) , xend = 2500 +50 * (i -1), yend = 4 +1 * (i -1)) +
      annotate("segment",x = -50 + 50 * (i -1) , y = 4 + 1 *(i-1) , xend =  -50 + 50 * (i -1), yend = 10 +1 * (i -1)) +
      theme_void()
  } else {
    df.i.c <- df.i %>% mutate(newrank = rank + 50 * (i -1), newmean = mean + 1 * (i -1))
    p <- p + 
      annotate("segment",x = -50 + 50 * (i -1) , y = 4 + 1 *(i-1) , xend = 2500 +50 * (i -1), yend = 4 +1 * (i -1)) +
      annotate("segment",x = -50 + 50 * (i -1) , y = 4 + 1 *(i-1) , xend =  -50 + 50 * (i -1), yend = 10 +1 * (i -1)) 
  }
}


for(i in 11:1) {

  c.i <- unique(sample$class)[i]
  message(c.i) 
  df.i <- df.q.t %>% filter(class == as.character(c.i))  %>% 
    group_by(ID,class) %>% summarise(mean = log10(mean(int))) %>% arrange(desc(mean)) 
  df.i$rank <- 1:nrow(df.i)
 
  df.i.c <- df.i %>% mutate(newrank = rank + 50 * (i -1), newmean = mean + 1 * (i -1))
  p <- p + geom_point(data = df.i.c, aes(newrank,newmean),col = col.sampleclass[i],size =2) 
  
}

p 
ggsave("../result/1.global/1.id.rank.3d.png",width = 8,height = 6)



#
# overlap 2 cell type
#

df.q.t <- read.xlsx("Supplementary Table 1.xlsx",sheet = "RawIntensity") %>% 
  gather(sample,int,-ID) %>% 
  filter(int > 0)%>% 
  left_join(sample)



list.cell <-  df.q.t %>%
  group_by(ID,class) %>% summarise(mean = mean(int) )%>% ungroup() %>% 
  mutate(is = 1) %>% select(ID,class,is) %>% spread(class,is) %>% 
  select(ID,rev(as.character(unique(sample$class)))) %>% 
  as.data.frame()



list.cell[,-1][is.na(list.cell[,-1])] <- 0



library(UpSetR)


png("../result/1.global/1.id.overlap.upset.png",res = 300,units = "cm",height =12,width = 15)
upset(  list.cell,  point.size = 5,
        sets = rev(as.character(unique(sample$class))),keep.order = T,
        nintersects = 20,
        order.by = "freq",
   
      
        line.size = 1.5, # 线宽 
        mainbar.y.label = "Overlap proteins",  
        sets.x.label = "Protein count",
        mb.ratio =   c(0.5,0.5),
        text.scale = c(    1.5, # y axis title size
                           1.1, # y axis text size     
                           1.5, # x axis title size    
                           1.1, # x axis text size
                           1.5, # dataset name    
                           1.0  # bar size  
        ),  
        sets.bar.color = rev(col.sampleclass), 
        matrix.color = "#FB9A99", 
        main.bar.color = "#1f77b4", 
        shade.color = "#d0e5f5",  
        shade.alpha = 0.5,  
)
dev.off()



pdf("../result/1.global/1.id.overlap.upset.pdf",height = 5,width =6)
upset(  list.cell,  point.size = 5,
        sets = rev(as.character(unique(sample$class))),keep.order = T,
        nintersects = 20,
        order.by = "freq",
        
        
        line.size = 1.5, # 线宽 
        mainbar.y.label = "Overlap proteins",  
        sets.x.label = "Protein count",
        mb.ratio =   c(0.5,0.5),
        text.scale = c(    1.5, # y axis title size
                           1.1, # y axis text size     
                           1.5, # x axis title size    
                           1.1, # x axis text size
                           1.5, # dataset name    
                           1.0  # bar size  
        ),  
        sets.bar.color = rev(col.sampleclass), 
        matrix.color = "#FB9A99", 
        main.bar.color = "#1f77b4", 
        shade.color = "#d0e5f5",  
        shade.alpha = 0.5,  
)
dev.off()

#--------------------------------
#
# rank
#
#--------------------------------


df.raw <- read.xlsx("Supplementary Table 1.xlsx",sheet = "RawIntensity") %>% 
  gather(sample,int,-ID) %>% 
  filter(int > 0) %>%
  left_join(sample)  #%>% mutate(class = as.character(class))

df.1 <- df.raw %>% 
  filter(class == "Morula") %>% 
  group_by(ID,class) %>% summarise(mean = mean(int)) %>% 
  arrange(desc(mean))
df.1$rank <- 1:nrow(df.1)

df.2 <- df.raw %>% 
  filter(class == "3.5d-TE") %>%
  group_by(ID,class) %>% summarise(mean = mean(int)) %>% 
  arrange(desc(mean))
df.2$rank <- 1:nrow(df.2)

df.3 <- df.raw %>% 
  filter(class == "4.5d-TE") %>%
  group_by(ID,class) %>% summarise(mean = mean(int)) %>% 
  arrange(desc(mean))
df.3$rank <- 1:nrow(df.3)

df.4 <- df.raw %>% 
  filter(class == "5.0d-EXE") %>%
  group_by(ID,class) %>% summarise(mean = mean(int)) %>% 
  arrange(desc(mean))
df.4$rank <- 1:nrow(df.4)

df.5 <- df.raw %>% 
  filter(class == "5.5d-EXE") %>%
  group_by(ID,class) %>% summarise(mean = mean(int)) %>% 
  arrange(desc(mean))
df.5$rank <- 1:nrow(df.5)

df.6 <- df.raw %>% 
  filter(class == "6.5d-EXE") %>%
  group_by(ID,class) %>% summarise(mean = mean(int)) %>% 
  arrange(desc(mean))
df.6$rank <- 1:nrow(df.6)








df.22 <- df.raw %>% 
  filter(class == "3.5d-ICM") %>%
  group_by(ID,class) %>% summarise(mean = mean(int)) %>% 
  arrange(desc(mean))
df.22$rank <- 1:nrow(df.22)

df.33 <- df.raw %>% 
  filter(class == "4.5d-ICM") %>%
  group_by(ID,class) %>% summarise(mean = mean(int)) %>% 
  arrange(desc(mean))
df.33$rank <- 1:nrow(df.33)

df.44 <- df.raw %>% 
  filter(class == "5.0d-EPI") %>%
  group_by(ID,class) %>% summarise(mean = mean(int)) %>% 
  arrange(desc(mean))
df.44$rank <- 1:nrow(df.44)

df.55 <- df.raw %>% 
  filter(class == "5.5d-EPI") %>%
  group_by(ID,class) %>% summarise(mean = mean(int)) %>% 
  arrange(desc(mean))
df.55$rank <- 1:nrow(df.55)

df.66 <- df.raw %>% 
  filter(class == "6.5d-EPI") %>%
  group_by(ID,class) %>% summarise(mean = mean(int)) %>% 
  arrange(desc(mean))
df.66$rank <- 1:nrow(df.66)


df.rank <- df.1 %>% 
  bind_rows(df.2) %>%  bind_rows(df.22) %>%
  bind_rows(df.3) %>%  bind_rows(df.33) %>% 
  bind_rows(df.4) %>%  bind_rows(df.44) %>% 
  bind_rows(df.5) %>%  bind_rows(df.55) %>% 
  bind_rows(df.6) %>%  bind_rows(df.66) 






df.rank$class <- factor(df.rank$class,levels = unique(df.rank$class),ordered = T)
ggplot(df.rank) +
  geom_point(aes(rank,log10(mean),col = class)) +
  scale_color_manual(values = col.sampleclass) +
  labs(col = "",x = "Rank proteins",y = expression(log[10]*'(mean intenisty)')) +
  
  theme_bw() +
  scale_x_continuous(labels = scales::comma) +
  guides(col = guide_legend(ncol = 3, byrow = TRUE))+
  theme(legend.position = c(0.6,0.75) ,
        panel.grid  = element_blank(),
        legend.background = element_blank(),
        legend.direction = "horizontal")

ggsave("../result/1.global/1.id.rank.6.2.png",width = 4.5,height = 3.5)
ggsave("../result/1.global/1.id.rank.6.2.pdf",width = 4.5,height = 3.5)



top <- df.rank %>% group_by(class) %>% summarise(count= n()) %>% mutate(cut = count * 0.05)  %>% 
  left_join(df.rank,by = "class" ) %>% mutate(delta = cut - rank) %>% 
  filter(rank <= 20)







df.morula <- df.rank %>% filter(class == "Morula")  %>% mutate(log = log10(mean))

top5 <- df.rank %>% 
  filter(class != "Morula") %>% 
  group_by(class) %>% top_n(n = 5, wt = mean) %>% 
  group_by(ID) %>% summarise(cell = paste(class,collapse = ";"),ncell = n()) %>% 
  left_join(fun$protein, by = c("ID" = "protein")) %>% 
  select(ID,genename,cell,ncell)  %>% 
  inner_join(df.morula)  %>% 
  separate_rows(cell,sep = ";") 


to.plot <- top5 %>% 
  mutate(count   = 1 ) %>% spread(cell,count)



t10 <- df.morula %>% ungroup() %>%  filter(rank == 10) %>% select(log) %>% as.numeric()
t50 <- df.morula %>% ungroup() %>% filter(rank == 50) %>% select(log) %>% as.numeric()


library(PieGlyph)
library(ggrepel)
ggplot()+
  geom_point(data = df.morula,aes(x = rank,y = log),col = "grey",alpha=0.1) +
  geom_hline(yintercept = t10,lty = 3) +
  annotate(geom = "text",x = 1500,y = t10,label = "top 10 proteins in morula" ,vjust = -1.2 )+
  
  labs(x = "Rank",y = expression(log[10]*'(Intensity)'),size = "Frequency of none-morula top5" )+
  
  geom_pie_glyph(data = to.plot, aes(x = rank, y = log),slices = unique(sample$class)[-1])+
  geom_text_repel(data = to.plot, aes(x = rank, y = log,label = genename)) +
  scale_fill_manual(values = col.sampleclass) +
  scale_x_continuous(labels = scales::comma)+
  labs(fill = "Cell type Top 5",x = "Rank morula proteins",y = expression(log[10]*'(Mean intensity)')) +
  theme_bw() +
  theme(legend.position = c(0.5,0.2),
        legend.key.size = unit(0.2, 'cm'),
    
        panel.grid = element_blank(),
        legend.title.position = "top",
        legend.text = element_text(size = 6),
    legend.direction = "horizontal")


ggsave("../result/1.global/1.id.rank.morula.pdf",width = 4,height = 3.5)
ggsave("../result/1.global/1.id.rank.morula.png",width = 4,height = 3.5)

















#----------------------------------------------------
#
# pca
# 
#-----------------------------------------------------
df.nor.final <- read.xlsx("Supplementary Table 1.xlsx",sheet = "NormalizedIntensity") %>% 
  gather(sample,int,-ID) %>% 
  mutate(log = log10(int)) %>% 
  left_join(sample) 



to.pca <- df.nor.final %>%
  select(sample,log,ID )%>% spread(ID,log) %>%
  left_join(sample)

df <- to.pca %>% select(-category,-time,-class,-file,-sample) %>% as.matrix()

# tsne <- Rtsne(df, initial_config = NULL, k = 2, initial_dims = 30, perplexity = 10,
#               max_iter = 1000, min_cost = 0, epoch_callback = NULL, whiten = TRUE,
#               epoch=100)
# 
# score <- tsne$Y %>% as.data.frame() %>% bind_cols(data.frame(class = to.pca$class))
# ggscatter(score, x = "V1", y = "V2",
#           fill = "class",shape = 21,
#           ellipse = TRUE, ellipse.type = "convex",size =1) +
#   
#   geom_point(aes(fill = class),size =2,shape = 21) +
#   geom_vline(xintercept = 0,lty = 3,col = "grey")+
#   geom_hline(yintercept = 0,lty = 3,col = "grey")+
#   scale_color_manual(values = col.sampleclass) +
#   scale_fill_manual(values = col.sampleclass) 
# 


rownames(df) <- to.pca$sample
pca_res <- prcomp(df, scale. = TRUE,center = T)
toscore <- score <- data.frame(sample = rownames(pca_res$x),pca_res$x )
score <- sample %>% left_join(toscore,by ="sample")
pc <- summary(pca_res)$importance


#---- pca score -----
#p2 <-

score.c <- score %>% 
  group_by(class) %>% summarise(PC1 = mean(PC1),PC2 = mean(PC2)) %>% ungroup() %>% mutate(sample = class)
                                                                                          

library(ggrepel)
score.c %>% bind_rows(score) %>% 
  ggplot(aes(PC1,PC2)) +

  geom_vline(xintercept = 0,lty = 3,col = "grey")+
  geom_hline(yintercept = 0,lty = 3,col = "grey")+
  stat_ellipse(geom = "polygon",aes(fill = class),alpha = 0.2) +
  geom_point(data = score,shape =21,size =3,col = "black",aes(fill = class))+
  geom_text_repel(data = score.c,aes(label = sample)) +
  
  scale_color_manual(values = col.sampleclass) +
  scale_fill_manual(values = col.sampleclass) +
  labs(x = paste0("PC1 (",round(pc[2,"PC1"] *100,2),"%)"),
       y = paste0("PC2 (",round(pc[2,"PC2"] *100,2),"%)"),
       col = "",
       fill = "")+
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank()) 



ggsave("../result/1.global/1.pca.pdf",width = 4,height = 4)
ggsave("../result/1.global/1.pca.png",width = 4,height = 4)








#--------------------------------------------
#
# heatmap 
#
#--------------------------------------------
library(pheatmap)
out.nor <- read.xlsx("Supplementary Table 1.xlsx",sheet = "NormalizedIntensity")
toheat <- out.nor %>% gather(sample,int,-ID) 
toheat <- sample %>% arrange(category) %>%  left_join(toheat)  

toheat$sample <- factor(toheat$sample,levels = unique(toheat$sample),ordered = T)
heat <- toheat %>% select(ID,sample,int) %>% spread(sample,int) %>% 
  select(contains("Morula"),contains("ICM"),contains("EPI")) %>% 
  as.data.frame() 

sampleclass <- sample %>% select(class,time,category) %>% as.data.frame()
rownames(sampleclass) <- sample$sample



df.gap <- sample %>% group_by(class) %>% summarise(n = n()) %>% mutate(sum = cumsum(n))

heatcol <- colorRampPalette(c("blue","green","yellow","orange","red"))(100)
p <- pheatmap(log10(heat),
              clustering_method = "complete",
              annotation_col = sampleclass,
              annotation_names_col = F,
              annotation_colors  = annotcol,
              #     gaps_col = df.gap$sum,
              
              show_rownames = F,
              show_colnames = F,
              treeheight_row = 20,
              treeheight_col = 20,
              cluster_cols = F,
              color = heatcol)

p1 <- ggplotify::as.ggplot(p)


heat <- toheat %>% select(ID,sample,int) %>% spread(sample,int) %>% 
  select(contains("Morula"),contains("TE"),contains("EXE")) %>% 
  as.data.frame() 

sampleclass <- sample %>% select(class,time,category) %>% as.data.frame()
rownames(sampleclass) <- sample$sample



df.gap <- sample %>% group_by(class) %>% summarise(n = n()) %>% mutate(sum = cumsum(n))

heatcol <- colorRampPalette(c("blue","green","yellow","orange","red"))(100)
p <- pheatmap(log10(heat),
              clustering_method = "complete",
              annotation_col = sampleclass,
              annotation_names_col = F,
              annotation_colors  = annotcol,
              #     gaps_col = df.gap$sum,
              
              show_rownames = F,
              show_colnames = F,
              treeheight_row = 20,
              treeheight_col = 20,
              cluster_cols = F,
              color = heatcol)

p2 <- ggplotify::as.ggplot(p)


p1 +p2

ggsave(file = "../result/1.global/1.id.heatmap.png",width = 12,height = 6)





heat <- toheat %>% select(ID,sample,int) %>% spread(sample,int) %>% 
  select(contains("Morula"),contains("ICM"),contains("EPI"),contains("TE"),contains("EXE")) %>% 
  as.data.frame() 

sampleclass <- sample %>% select(class,time,category) %>% as.data.frame()
rownames(sampleclass) <- sample$sample



df.gap <- sample %>% group_by(class) %>% summarise(n = n()) %>% mutate(sum = cumsum(n))

heatcol <- colorRampPalette(c("blue","white","red"))(100)
p <- pheatmap(log10(heat),scale = "row",
              clustering_method = "complete",
              annotation_col = sampleclass,
              annotation_names_col = F,
              annotation_colors  = annotcol,
              #     gaps_col = df.gap$sum,
              
              show_rownames = F,
              show_colnames = F,
              treeheight_row = 20,
              treeheight_col = 20,
              cluster_cols = F,
              color = heatcol)

p1 <- ggplotify::as.ggplot(p)

ggsave(file = "../result/1.global/1.id.heatmap.scale.png",width = 6,height = 5)
ggsave(file = "../result/1.global/1.id.heatmap.scale.pdf",width = 6,height = 5)


#------------------------------
#
# global tf
#
#-----------------------------

load("./mouse.uniprot.function.rda")

df.id <- fun$protein


df.raw <- read.xlsx("Supplementary Table 1.xlsx",sheet = "RawIntensity")  %>% 
  gather(sample,int,-ID) %>% 
  filter(int >0) %>% 
  left_join(sample) %>% 
  group_by(class,ID) %>% summarise(mean = mean(int))


count.class <- df.raw %>% 
  group_by(class) %>% summarise(n=n())


type.tf <- df.raw %>% left_join(df.id,by= c("ID"="protein")) %>%  filter(!is.na(tf_family)) %>% 
  group_by(class) %>% summarise(count =n()) %>% mutate(type = "TF")

type.pk <- df.raw %>% left_join(df.id,by= c("ID"="protein")) %>% filter(!is.na(kinase_category)) %>%
  group_by(class) %>% summarise(count =n()) %>% mutate(type = "Kinase")


type.re <- df.raw %>% left_join(df.id,by= c("ID"="protein")) %>%
  mutate(is.rec = str_match(description,"receptor")) %>%
  filter(!is.na(is.rec) ) %>% 
  group_by(class) %>% summarise(count =n()) %>% mutate(type = "Receptor")


count.type <- type.tf %>% bind_rows(type.pk) %>% bind_rows(type.re) %>% 
  left_join(count.class) %>% mutate(per = count/n)


library(ggsci)


ggplot(count.type,aes(per,class,fill = type)) +
  geom_bar(stat = "identity",alpha = 0.8) +
  labs(x = "proportion",y = "",fill = "") +
  scale_fill_npg()  +
  scale_x_continuous(labels = scales::percent)+
  theme_bw() +
  theme(legend.position = "top",
        panel.grid = element_blank()
        )
ggsave("../result/1.global/3.pktf.propotion.png",width = 4,height = 4)
ggsave("../result/1.global/3.pktf.propotion.pdf",width = 4,height = 4)



#------------
#
# TF identification
#
#------------
df.raw <- read.xlsx("Supplementary Table 1.xlsx",sheet = "RawIntensity")  %>% 
  gather(sample,int,-ID) %>% 
  filter(int >0) %>% 
  left_join(sample) 



count.sample<- df.raw %>% 
  group_by(class,sample) %>% summarise(n=n())

type.tf <- df.raw %>% left_join(df.id,by= c("ID"="protein")) %>%  filter(!is.na(tf_family)) %>% 
  group_by(sample) %>% summarise(count =n()) %>% mutate(type = "TF")

type.pk <- df.raw %>% left_join(df.id,by= c("ID"="protein")) %>% filter(!is.na(kinase_category)) %>%
  group_by(sample) %>% summarise(count =n()) %>% mutate(type = "Kinase")

type.re <- df.raw %>% left_join(df.id,by= c("ID"="protein")) %>%
  mutate(is.rec = str_match(description,"receptor")) %>%
  filter(!is.na(is.rec) ) %>% 
  group_by(sample) %>% summarise(count =n()) %>% mutate(type = "Receptor")




count.type <- type.tf %>% bind_rows(type.pk) %>% bind_rows(type.re) %>% 
  left_join(count.sample) %>% mutate(per = count/n)

count.class <- count.type %>% 
  group_by(type,class) %>% summarise(mean = mean(per),sd = sd(per))

ggplot(count.type,aes(class,per,fill = class)) +
  geom_bar(data = count.class,aes(class,mean),stat = "identity") +
  geom_point(shape = 21,size = 2) +
  facet_wrap( type~. ,ncol = 1) +
  scale_fill_manual(values = col.sampleclass) +
  scale_y_continuous(labels = scales::percent,expand = expansion(mult = c(0.1,0.2))) +
  stat_compare_means(label.y = 0.05,label.x.npc = 0.1,size = 3) +
  stat_compare_means(comparisons = list(c(unique(sample$class)[2],unique(sample$class)[3]),
                                        c(unique(sample$class)[4],unique(sample$class)[5]),
                                        c(unique(sample$class)[6],unique(sample$class)[7]),
                                        c(unique(sample$class)[8],unique(sample$class)[9]),
                                        c(unique(sample$class)[10],unique(sample$class)[11])
  ),label.y = 0.04,size =3,method = "t.test",label = "p.format",hide.ns = T) +
  labs(y = "Propotion",x = "") +
  theme_bw() +
  theme(legend.position = "none",
        strip.background.x = element_blank(),
    axis.text.x = element_text(angle = 45,hjust = 1))

ggsave("../result/1.global/3.pktf.id.per.png",width = 4,height = 8)


# count number
count.class <- count.type %>% 
  group_by(type,class) %>% summarise(mean = mean(count),sd = sd(count))

ggplot(count.type,aes(class,count,fill = class)) +
  geom_bar(data = count.class,aes(class,mean),stat = "identity") +
  geom_point(shape = 21,size = 2) +
  facet_wrap( type~. ,ncol = 1) +
  scale_fill_manual(values = col.sampleclass) +
  scale_y_continuous(expand = expansion(mult = c(0.1,0.2))) +
  stat_compare_means(label.y = 50,label.x.npc = 0.1,size = 3) +
  stat_compare_means(comparisons = list(c(unique(sample$class)[2],unique(sample$class)[3]),
                                        c(unique(sample$class)[4],unique(sample$class)[5]),
                                        c(unique(sample$class)[6],unique(sample$class)[7]),
                                        c(unique(sample$class)[8],unique(sample$class)[9]),
                                        c(unique(sample$class)[10],unique(sample$class)[11])
  ),label.y = 40,size =3,method = "t.test",label = "p.format",hide.ns = T) +
  labs(y = "Protein count",x = "") +
  theme_bw() +
  theme(legend.position = "none",
        strip.background.x = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1))










#------------------------------
#
# global tf
#
#-----------------------------
df.nor <- read.xlsx("Supplementary Table 1.xlsx",sheet = "NormalizedIntensity") 
df.scale <-  read.xlsx("Supplementary Table 1.xlsx",sheet = "ZscoreInteisty") 


df.q <- df.nor %>%  gather(sample,int,-ID) 
df.q <-  sample %>% left_join(df.q,by = "sample")
df.q$sample <- factor(df.q$sample,levels = unique(df.q$sample),ordered = T)

df.sc <- df.scale %>%  gather(sample,int,-ID) 
df.sc <- sample %>% left_join(df.sc,by = "sample")
df.sc$sample <- factor(df.sc$sample,levels = unique(df.sc$sample),ordered = T)



df.tf <- df.id %>% filter(!is.na(tf_family)) %>% 
  filter(protein %in% df.nor$ID) %>% 
  mutate(ID = protein,type = "TF")  %>% 
  select(type,ID,genename)
df.pk <- df.id %>% filter(!is.na(kinase_category)) %>%
  filter(protein %in% df.nor$ID) %>% 
  mutate(ID = protein,type = "Kinase") %>% 
  select(type,ID,genename)

df.re <- df.id %>% mutate(is.rec = str_match(description,"receptor")) %>%
  filter(!is.na(is.rec) ) %>% 
  filter(protein %in% df.nor$ID)  %>% 
  mutate(ID = protein,type = "Receptor") %>% 
  anti_join(df.tf,by = "ID") %>% 
  anti_join(df.pk,by ="ID") %>% 
  
  select(type,ID,genename)




df.tf.p <- df.tf %>% 
  bind_rows(df.pk,df.re) %>% 
  inner_join(df.sc)

df.tf.class <- df.tf.p %>% group_by(type,class,ID) %>% summarise(mean = mean(int))
to.plot <- sample %>% arrange(time)%>% 
  
  left_join(df.tf.class)
to.plot$class <- factor(to.plot$class,levels = unique(to.plot$class),ordered = T)





p <- ggplot(to.plot,aes(class,mean)) +
  geom_violin(aes(fill = class),col = NA) +
  geom_boxplot(outliers = F,width = 0.2) +
 # geom_jitter(size = .2) +
  #facet_grid(type ~. ) +
  facet_wrap(.~ type,ncol = 1) +
  stat_compare_means(label.y = 5,hjust = -0.1,size =3) +  
  stat_compare_means(comparisons = list(c(unique(sample$class)[2],unique(sample$class)[3]),
                                        c(unique(sample$class)[4],unique(sample$class)[5]),
                                        c(unique(sample$class)[6],unique(sample$class)[7]),
                                        c(unique(sample$class)[8],unique(sample$class)[9]),
                                        c(unique(sample$class)[10],unique(sample$class)[11])
  ),label.y = 3,size =3,label = "p.format",hide.ns = T) +
  scale_y_continuous(expand = c(0.2,0.2)) +
  scale_fill_manual(values = col.sampleclass) +
  
  labs(x = "", y = "Z-score") +
  theme_bw() +
  theme(legend.position = "none",
        strip.text.x = element_text(colour = "white"),
        axis.text.x = element_text(angle = 30,hjust = 1))




pp <- ggplot_gtable(ggplot_build(p))
strips <- which(grepl('strip-', pp$layout$name))

col.i <- ggsci::pal_npg()(3)

for (i in seq_along(strips)) {
  k <- which(grepl('rect', pp$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  pp$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- col.i[i]
  
}

ppp <- grid::grid.draw(pp)

ggsave(pp,file= "../result/1.global/3.pktf.png",width = 4,height = 4)
ggsave(pp,file= "../result/1.global/3.pktf.pdf",width = 4,height = 4)




ggplot(to.plot,aes(mean,class,fill = class)) +
  geom_boxplot(outliers = F) +
  #geom_jitter(size = .2,shape =21) +
  facet_wrap(type ~.,) +
  stat_compare_means(label.x = 5,hjust = -0.1,size =3) +
  stat_compare_means(comparisons = list(c(unique(sample$class)[2],unique(sample$class)[3]),
                                        c(unique(sample$class)[4],unique(sample$class)[5]),
                                        c(unique(sample$class)[6],unique(sample$class)[7]),
                                        c(unique(sample$class)[8],unique(sample$class)[9]),
                                        c(unique(sample$class)[10],unique(sample$class)[11])
  ),label.x = 3,size =3,label = "p.format",hide.ns = T) +
  scale_x_continuous(expand = c(0.2,0.2)) +
  scale_color_manual(values = col.sampleclass) +
  scale_fill_manual(values = col.sampleclass) +
    
  labs(x = "Z-score", y = "") +
  theme_bw() +
  theme(legend.position = "none")




