library(tidyverse)
library(openxlsx)
library(RColorBrewer)
library(ComplexHeatmap)

load("./mouse.uniprot.function.rda")
df.id <- fun$protein %>% mutate(ID = protein)

df.name <- df.id %>% select(ID,genename)


#
# sample
#
sample <- read.xlsx("D:/project/proteome_embr/manuscript/Table1.Proteome.xlsx",sheet = "FileInformation")
sample$category <- factor(sample$category,levels = unique(sample$category),ordered = T)




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


sample <- sample %>% arrange(category,time)
sample$sample <- factor(sample$sample,levels = sample$sample,ordered = T) 
sample$class <- factor(sample$class,levels = unique(sample$class),ordered = T) 



#
# zscore matrix
#

df.raw <- read.xlsx("D:/project/proteome_embr/manuscript/Table1.Proteome.xlsx",sheet = "RawIntensity")
df.nor <- read.xlsx("D:/project/proteome_embr/manuscript/Table1.Proteome.xlsx",sheet = "NormalizedIntensity") 
df.scale <-  read.xlsx("D:/project/proteome_embr/manuscript/Table1.Proteome.xlsx",sheet = "ZscoreInteisty") 

df.q <- df.nor %>%  gather(sample,int,-ID) 
df.q <-  sample %>% left_join(df.q,by = "sample")
df.q$sample <- factor(df.q$sample,levels = unique(df.q$sample),ordered = T)

df.sc <- df.scale %>%  gather(sample,int,-ID) 
df.sc <- sample %>% left_join(df.sc,by = "sample")
df.sc$sample <- factor(df.sc$sample,levels = unique(df.sc$sample),ordered = T)


df.sc.mean <- df.sc %>% group_by(ID,class) %>% summarise(scale = mean(int))
df.scale.mean <- data.frame(class =unique(sample$class)) %>% left_join(df.sc.mean) %>% 
  mutate(class = str_replace(class,"ICM$","ICM-EPI")) %>% 
  mutate(class = str_replace(class,"d-EPI","d-ICM-EPI")) %>% 
  mutate(class = str_replace(class,"TE$","TE-EXE")) %>% 
  mutate(class = str_replace(class,"d-EXE","d-TE-EXE")) 
df.scale.mean$class <- factor(df.scale.mean$class,levels = unique(df.scale.mean$class),ordered = T)



df.zscore <- df.sc %>% spread(sample,int)




#
# cluster information
#

all.cluster <- read.xlsx("D:/project/proteome_embr/manuscript/Table2.Cluster.version202503.xlsx") %>% 
  arrange(cluster)
df.cluster <- all.cluster %>% 
  group_by(genename) %>%top_n(n = 1,wt = membership) %>% ungroup() %>% 
  filter(!is.na(genename)) %>%
  filter(membership > 0.1)





to.heat <- df.cluster %>% 
  select(ID) %>% left_join(df.scale)

heat <-to.heat%>% select(sample$sample) %>% as.data.frame()
row.names(heat ) <- to.heat$ID

limit <- min(max(heat),abs(min(heat)))
heat[heat > limit] <- limit
heat[heat < -limit] <- -limit

annot.sample <- sample %>% select(category,time,class) %>% as.data.frame()
row.names(annot.sample) <- sample$sample


df.gene <- df.cluster %>% group_by(cluster) %>% summarise(n=n()) %>% 
  mutate(cl = paste0(str_replace(paste0("C",cluster+100),"C1","C")," (n=",n,")"))  %>% 
  left_join(df.cluster)


annot.gene <- df.gene %>% 
  select(cl,membership) %>% as.data.frame()

row.names(annot.gene) <- df.gene$genename




col.membership <- circlize::colorRamp2(c( 0,0.5,1), c("green", "yellow","red")) 

heatcol <- colorRampPalette(c("blue","white","red"))(100)


#col1 <- colorRampPalette(c( "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5"))
col1 <- colorRampPalette(RColorBrewer::brewer.pal(9,"Spectral"))
col.cluster <- col1(10)
names(col.cluster) <- unique(annot.gene$cl)



#
# begin complexheatmap
#



df.cl.fun <- read.xlsx("Supplementary Table 1.xlsx",sheet = "ClusterGO")
df.cell <- read.xlsx("./cell_type.xlsx")

df.stemcell <- df.cl.fun %>% mutate(is = str_match(term,"stem cell")) %>% 
  filter(!is.na(is)) %>% 
  separate_rows(gene,sep = ";") %>% 
  group_by(gene) %>% summarise(n=n(),term = paste(term,collapse = ";"))


df.col.cl <- data.frame(cluster = 1:10,
                        col = col.cluster)

df.mark <- df.cluster %>% mutate(pos = 1:nrow(df.cluster) ) %>% 
  select(genename,cluster,pos) %>% 
  left_join(df.col.cl) %>% 
  filter(cluster != 1) %>% 
  inner_join(df.stemcell,by = c("genename" = "gene"))


mark.gene  = rowAnnotation(
  foo = anno_mark(at = df.mark$pos, 
                  labels = df.mark$genename,
                  labels_gp = gpar(col= df.mark$col,fontsize = 8)),
  annotation_name_side = "bottom"
)


row_annot_obj <- rowAnnotation(Cluster = anno_simple(annot.gene$cl, col = col.cluster,border = T),
                               Membership = anno_simple(annot.gene$membership, col = col.membership,border = T),
                               annotation_name_side = "top"
                               
                               )

heat.a.cl = Legend(labels = unique(annot.gene$cl), title = "cluster",
                   legend_gp = gpar(fill = col.cluster,ncol = 3))
heat.a.mem = Legend(title = "membership",
                    labels = c("0", "0.5", "1"),at = c(0,0.5,1),col_fun = col.membership)




col_annot_obj <- columnAnnotation(df = annot.sample,
                                  col = annotcol )





p <-
  ComplexHeatmap::Heatmap(heat, 
                          heatmap_legend_param=list(title="z-score"),
                          
                          # column setting
                          top_annotation = col_annot_obj,
                          column_split = sample$category,
                          show_column_names = F,
                          cluster_columns = F,
                          column_title_gp=grid::gpar(fontsize=10),
                          ## row setting
                          left_annotation= row_annot_obj,
                 #         right_annotation = mark.gene,
                          clustering_method_rows = "complete",
                          row_split = annot.gene$cl,
                          row_title = NULL,
                          show_row_names = F,
                          
                        #  row_title_rot = F,
                        #  row_title_side = "left",
                        #  row_title_gp=grid::gpar(fontsize=7.5),
                          
                          cluster_rows = F,
                          border = TRUE
  )



heat.a.cl = Legend(labels = unique(annot.gene$cl), title = "cluster",
                    legend_gp = gpar(fill = col.cluster),ncol = 2)

png("../result/2.cluster/cluster.heatmap.png",units = "cm",res = 300,height = 25,width = 20)

draw(p, 
     
     annotation_legend_list = list(heat.a.cl,
                                   heat.a.mem),
     merge_legend = TRUE,
     heatmap_legend_side = "left", annotation_legend_side = "left",
     padding = unit(c(2, 3, 2, 5), "mm"))

dev.off()



#-------------------------------------------
#
# boxplotplot of interesting proteins
#
#-------------------------------------------
library(ggsci)

test <- df.cluster %>% 
  inner_join(df.stemcell,by = c("genename" = "gene")) %>% 
  
  left_join(df.sc) %>% filter(!is.na(time)) 
  
test.0 <- test %>% filter(category != "Morula")
test.1 <- test %>% filter(category == "Morula") %>% 
  mutate(category = "ICM-EPI")
test.2 <- test %>% filter(category == "Morula") %>% 
  mutate(category = "TE-EXE")

to.plot <- test.0 %>% 
  bind_rows(test.1) %>% 
  bind_rows(test.2  ) %>% 
  mutate(time = as.numeric(time))




t1 <- data.frame(genename = c("Dppa2","Pcm1",
                              "Cnot1","Kdm2b",
                              "Ociad1","Zc3h13",
                              "Klf10","Mapk3",
                              "Anxa6","Dnmt3l",
                              "Lin28a","Nsun2"
                              )) %>% 
  left_join(to.plot)
t1$genename <- factor(t1$genename,levels = unique(t1$genename),ordered = T)




ggplot(t1,aes(time,int,col  = category)) +
  geom_point() +
  geom_smooth(span = 1) +
  facet_wrap(genename ~.,ncol = 2) +
  labs(y = "Zscore",col = "",x = "Time (day)")+
  scale_color_jco() +
  theme_bw() +
  theme(legend.position = "top",
       # legend.key.size = unit(0.1,"cm"),
        strip.background.x = element_blank()) 






ggsave("../result/2.cluster/cluster.heatmap.bp.pdf",width = 4,height =9 )


#---------------------
#  target gene 
#---------------------

test <- df.cluster %>% 
  
  left_join(df.sc) %>% filter(!is.na(time)) 

test.0 <- test %>% filter(category != "Morula")
test.1 <- test %>% filter(category == "Morula") %>% 
  mutate(category = "ICM-EPI")
test.2 <- test %>% filter(category == "Morula") %>% 
  mutate(category = "TE-EXE")

to.plot <- test.0 %>% 
  bind_rows(test.1) %>% 
  bind_rows(test.2  ) %>% 
  mutate(time = as.numeric(time))





t1 <- data.frame(genename = c("Smarcad1","Smarcc2","Smarcal1",
                              "Smarce1","Smarca5",
                              "Smarca4","Smarcd1",
                              "Smarcc1"
)) %>% 
  left_join(to.plot)
t1$genename <- factor(t1$genename,levels = unique(t1$genename),ordered = T)




ggplot(t1,aes(time,int,col  = category)) +
  geom_point() +
  geom_smooth(span = 1) +
  facet_wrap(cluster ~ genename ,ncol = 2) +
  labs(y = "Zscore",col = "",x = "Time (day)")+
  scale_color_jco() +
  theme_bw() +
  theme(legend.position = "top",
        # legend.key.size = unit(0.1,"cm"),
        strip.background.x = element_blank()) 

dir.create("../result/5.interesting_protein")

ggsave("../result/5.interesting_protein/smarc.png",width = 4,height =9 )





# yap1

t1 <- data.frame(genename = c("Yap1"
)) %>% 
  left_join(to.plot)
t1$genename <- factor(t1$genename,levels = unique(t1$genename),ordered = T)




ggplot(t1,aes(time,int,col  = category)) +
  geom_point() +
  geom_smooth(span = 1) +
  facet_wrap(cluster ~ genename ,ncol = 2) +
  labs(y = "Zscore",col = "",x = "Time (day)")+
  scale_color_jco() +
  theme_bw() +
  theme(legend.position = "top",
        # legend.key.size = unit(0.1,"cm"),
        strip.background.x = element_blank()) 




# 3.5d
t1 <- data.frame(genename = c("Eed","Naa30",
                              "Uhrf1","Kdm2b","Arid1a","Smarce1","Smarcd1",
                              "Ep400","Mta3",
                              "Hdac6","Jarid2","Rere",
                              "Ncor2","Chd7",#"p300",
                              "Kmt2d","Iws1","Smarca4","Smarcd1",
                              "Kdm3a","Ctcf","Rad50",
                              "Sirt1","Dnmt3l",
                              "Klf10",
                              "Sub1","Rest","Wdhd1",
                              "Mta2","Baz2a","Dpf1","Setd1a"
                        
                              

)) %>% 
  left_join(to.plot)
t1$genename <- factor(t1$genename,levels = unique(t1$genename),ordered = T)



for (i in 1:length(unique(t1$genename))){
  g <- unique(t1$genename)[i] %>% as.character()
  p <- t1 %>% filter(genename == g)
pp <- ggplot(p,aes(time,int,col  = category)) +
  geom_point() +
  geom_smooth(span = 1) +
  facet_wrap(  genename ~. ) +
  labs(y = "Zscore",col = "",x = "Time (day)")+
  scale_color_jco() +
  theme_bw() +
  theme(legend.position = "none",
        # legend.key.size = unit(0.1,"cm"),
        strip.background.x = element_blank()) 



  ggsave(pp,file = paste0("../result/5.interesting_protein/timepoint.",i,"-",g,".png"),width = 2,height =2 )
}
