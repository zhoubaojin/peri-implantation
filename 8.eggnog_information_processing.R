library(tidyverse)
library(openxlsx)
library(RColorBrewer)
library(ComplexHeatmap)
library(ggsci)

#--------------------------
# target function 
# 1. Dr Wu.
# 2. Epigenome
# 3. Mitochrondria
# 4. other transcription factor
#--------------------------

load("./mouse.uniprot.function.rda")
df.id <- fun$protein %>% mutate(ID = protein)

df.name <- df.id %>% select(ID,genename)

#
# sample
#
sample <- read.xlsx("Supplementary Table 1.xlsx",sheet = "FileInformation")
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
df.scale.mean <- data.frame(class =unique(sample$class)) %>% left_join(df.sc.mean) %>% 
  mutate(class = str_replace(class,"ICM$","ICM-EPI")) %>% 
  mutate(class = str_replace(class,"d-EPI","d-ICM-EPI")) %>% 
  mutate(class = str_replace(class,"TE$","TE-EXE")) %>% 
  mutate(class = str_replace(class,"d-EXE","d-TE-EXE")) 
df.scale.mean$class <- factor(df.scale.mean$class,levels = unique(df.scale.mean$class),ordered = T)



df.zscore <- df.sc %>% spread(sample,int)

all.cluster <- read.xlsx("Supplementary Table 2.xlsx") %>% 
  arrange(cluster)
df.cluster <- all.cluster %>% 
  group_by(genename) %>%top_n(n = 1,wt = membership) %>% ungroup() %>% 
  filter(!is.na(genename)) %>%
  filter(membership > 0.1)




#--------------------------------------------------------
#
# dataset global genetic
#
#--------------------------------------------------------


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




df.m <- df.cluster %>% 
  #  filter(cluster != "1") %>% 
  filter(eggnog_category == "Information storage and processing")  %>% 
  mutate(EggNOG = str_replace_all(eggnog_term,",.+","")) %>% 
  mutate(EggNOG = str_replace_all(EggNOG," and dynamics",""))  %>% 
  mutate(EggNOG = str_replace_all(EggNOG," and modification",""))  


df.m.c <- df.m %>% 
  mutate(lineage = if_else(cluster == 1,"Morula",
                           if_else(cluster >1 & cluster <= 5,"Multiple TPs","Single TP"))) %>% 
  group_by(lineage,EggNOG) %>% summarise(n=n()) %>% ungroup()



# p <- data.frame(EggNOG = c("Replication",
#                            "Transcription",
#                            "RNA processing",
#                            "Chromatin structure",
#                            "Translation"
# )) %>% 
#   left_join(df.m.c)

p <- data.frame(EggNOG = c("Translation",
                           "Chromatin structure",
                           "RNA processing",
                           "Transcription",
                           "Replication"
),i = c(1:5)) %>% 
  left_join(df.m.c)

p <- data.frame(lineage = c("Multiple TPs","Single TP","Morula")) %>% left_join(p)
p$EggNOG <- factor(p$EggNOG,levels = unique(p$EggNOG),ordered = T)
p$lineage <- factor(p$lineage,levels = unique(p$lineage),ordered = T)


pos <- p %>% 
  group_by(EggNOG) %>% summarise(sum  = sum(n)) %>% 
  left_join(p) %>% mutate(per= round(n/sum,2)) %>% 
  mutate(cum = cumsum(per)) %>% mutate(ymin = i-cum)  

pos$tmp <- c(1,pos$ymin[-nrow(pos)])
pos <- pos %>% mutate(ymax = if_else(lineage == "Multiple TPs",1,tmp)) %>% 
  mutate(pos= ymin + (ymax -ymin)/2)

ggplot(p,aes(n,EggNOG)) +
  geom_bar(aes(fill = lineage),stat = "identity",position="fill")  +
  geom_text(data= pos,aes(pos,EggNOG,label = per),size = 2.5) +
  scale_fill_npg() +
  scale_x_continuous(labels = scales::percent) +
  
  labs(y = "",fill = "",x = "") +
  theme_minimal() +
  theme(legend.position = "top")

ggsave("../result/3.dna_rna_epi/information.term.png",width = 4.5,height = 2.2)
ggsave("../result/3.dna_rna_epi/information.term.pdf",width = 4.5,height = 2.2)



t1 <- data.frame(eggnog_term = c("Replication, recombination and repair",
                                 "Transcription",
                                 "Chromatin structure and dynamics",
                                 "RNA processing and modification",
                                 "Translation, ribosomal structure and biogenesis"))  %>%
  
  left_join(to.plot) %>%  
  mutate(EggNOG = str_replace_all(eggnog_term,",.+","")) %>% 
  mutate(EggNOG = str_replace_all(EggNOG," and dynamics",""))  %>% 
  mutate(EggNOG = str_replace_all(EggNOG," and modification","")) %>%   
  mutate(lineage = if_else(cluster == 1,"Morula",
                           if_else(cluster >1 & cluster <= 5,"Multiple TPs","Single TP"))) 

t1$EggNOG <- factor(t1$EggNOG,levels = unique(t1$EggNOG),ordered = T)
t1$genename <- factor(t1$genename,levels = unique(t1$genename),ordered = T)






to.heat <- t1 %>% 
  mutate(de = if_else(class == "Morula" & category == "ICM-EPI","yes","no" ))  %>% 
  filter(de  == "no") %>% 
  select(genename,sample,int,EggNOG,lineage) %>% spread(sample,int) #%>% 



annot.sample <- sample %>% select(time,class,category) %>% as.data.frame()
row.names(annot.sample) <- sample$sample


heat <- to.heat %>% select(-genename,-EggNOG,-lineage) %>% as.matrix()
row.names(heat) <- to.heat$genename



to.annot.g <- to.heat %>% 
  group_by(EggNOG ) %>% summarise(n = n()) %>% mutate(type = paste0(EggNOG," (n=",n,")")) %>% 
  left_join(to.heat)
annot.g <- to.annot.g %>%
  mutate(exp.type = lineage) %>% 
  select(type,exp.type) %>% as.data.frame()

row.names(annot.g) <- to.annot.g$genename


col.g <- ggsci::pal_jco(alpha = 0.7)(length(unique(annot.g$type)))
names(col.g) <- unique(annot.g$type)
col.tp <- ggsci::pal_npg()(3)
names(col.tp) <- c("Morula","Single TP","Multiple TPs")


col_annot_obj <- columnAnnotation(df = annot.sample,col = annotcol)
row_annot_obj <- rowAnnotation(type = anno_simple(as.character(annot.g$type),col = col.g),
                               exp.type = anno_simple(as.character(annot.g$exp.type),col = col.tp),
                               annotation_name_side = "top")




heat.a.type = Legend(labels = names(col.g), title = "EggNOG",
                     legend_gp = gpar(fill = col.g),gap = unit(1, "cm"))
heat.a.tp = Legend(labels = names(col.tp), title = "Exp.type",
                   legend_gp = gpar(fill = col.tp))




p <- 
  Heatmap(heat,name = "Z-score",
          column_split = sample$category,column_title = NULL,
          row_split =  annot.g$type,row_title = NULL,
          column_title_gp=grid::gpar(fontsize=10),
          top_annotation = col_annot_obj,
          left_annotation = row_annot_obj,
          cluster_rows = F,cluster_columns = F,
          show_column_names = F,show_row_names = F)

png("../result/3.dna_rna_epi/all_info.png",units = "cm",width = 17,height = 22.5,res = 300)
draw(p, 
     
     annotation_legend_list = list(heat.a.type,
                                   heat.a.tp),
     merge_legend = TRUE,
     heatmap_legend_side = "bottom", annotation_legend_side = "bottom",
     
     padding = unit(c(2, 2, 2, 2), "mm"))
dev.off()


pdf("../result/3.dna_rna_epi/all_info.pdf",width = 6,height = 9)
draw(p, 
     
     annotation_legend_list = list(heat.a.type,
                                   heat.a.tp),
     merge_legend = TRUE,
     heatmap_legend_side = "bottom", annotation_legend_side = "bottom",
     
#     padding = unit(c(2, 2, 2, 2), "mm")
)
dev.off()









data.frame(a=c("a","b"),
           b=c(1,2),
           c=c("Transcription co-regulators","Transcription factor")) %>%
  
  ggplot(aes(a,c,fill = c)) +
  geom_bar(stat ="identity") +
  scale_fill_jco() +
  labs(fill= "") +
  theme_void()
ggsave("../result/3.dna_rna_epi/annotate.pdf",width = 3,height = 3)



#--------------------------------------------------------
#
# dataset from PhD Wu
#
#--------------------------------------------------------


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





t1 <- data.frame(genename = c("Eed","Naa30",
                              "Uhrf1","Kdm2b","Arid1a","Smarce1",
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

target$wu <- data.frame(genename = unique(t1$genename))


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





# heeat 


t1 <- data.frame(genename = c("Eed","Naa30",
                              "Kdm2b","Arid1a",
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


to.heat <- t1 %>% 
  mutate(de = if_else(class == "Morula" & category == "ICM-EPI","yes","no" ))  %>% 
  filter(de  == "no") %>% 
  select(genename,sample,int) %>% spread(sample,int) %>% 
  select(genename,
         contains("35"),
         contains("45"),
         contains("50"),
         contains("55"),
         contains("65")
         )

heat <- to.heat %>% select(-genename) %>% as.matrix()
row.names(heat) <- to.heat$genename


file.g <- read.xlsx("../result/5.interesting_protein/epi.xlsx")
annot.g <- file.g %>% select(type)



row.names(annot.g) <- file.g$genename
heat[heat >3] <- 3

# heatcol <- colorRampPalette(c("steelblue4","steelblue2","white","pink","red2","red3"))(100)
# p <- pheatmap::pheatmap(heat,#color = heatcol,border_color = F,
#                    cluster_cols = F,cluster_rows = F,
#                    annotation_row = annot.g,
#                    annotation_names_row = F,
#                    annotation_col = annot.sample,
#                    annotation_colors = annotcol,
#                    show_colnames = F,
#                    gaps_col = c(6,12,18,24),
#                    gaps_row = c(2,4,6,9,18,21)
#                    )
# ggplotify::as.ggplot(p)
# 
# ggsave()


to.sample <- data.frame(sample = colnames(heat)) %>% left_join(sample)
annot.sample <- to.sample %>% select(time,class,category) %>% as.data.frame()
row.names(annot.sample) <- to.sample$sample

annotcol.dem <- list()
annotcol.dem$class <- annotcol$class[-1]
annotcol.dem$category <- annotcol$category[-1]
annotcol.dem$time <- annotcol$time[-1]



col_annot_obj <- columnAnnotation(df = annot.sample,col = annotcol.dem)



#col.g <- brewer.pal(length(unique(file.g$type)),"Set1")



col.g <- ggsci::pal_npg(alpha = 0.7)(length(unique(file.g$type)))
names(col.g) <- unique(file.g$type)

row_annot_obj <- rowAnnotation(type = anno_simple(file.g$type,col = col.g))
                               

heat.a.type = Legend(labels = unique(file.g$type), title = "type",
                   legend_gp = gpar(fill = col.g))




p <- Heatmap(heat,name = "Z-score",
             column_split = to.sample$time,#column_title = NULL,
             column_title_gp=grid::gpar(fontsize=10),
             row_split = c(rep("3.5",4),
                           rep("4.5",5),
                           rep("5.0",9),
                           rep("5.5",3),
                           rep("6.5",7)),row_title = NULL,
        top_annotation = col_annot_obj,right_annotation = row_annot_obj,
        cluster_rows = F,cluster_columns = F,show_column_names = F)

png("../result/5.interesting_protein/epi.png",units = "cm",width = 15,height = 20,res = 300)
draw(p, 
     
     annotation_legend_list = list(heat.a.type),
     merge_legend = TRUE,
     heatmap_legend_side = "bottom", annotation_legend_side = "bottom",
     padding = unit(c(2, 2, 2, 2), "mm"))
dev.off()












#--------------------------------------------------------
#
# dataset global epi
#
#--------------------------------------------------------


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




df.m <- df.cluster %>% 
  mutate(is = str_match(description,"methyltrans")) %>% 
  filter(!is.na(is)) %>%
  mutate(is1 = str_match(description,"DNA")) %>% 
  mutate(is2 = str_match(description,"RNA")) %>% 
  filter(genename != "Bhmt" & genename != "Shmt2") %>% 
  #mutate(type = "Writer") %>% 
  mutate(type = if_else(!is.na(is1),"Writer-DNA-m",
                           if_else(!is.na(is2),"Writer-RNA-m","Writer-Pro-m"))) %>% 
  arrange(type,genename)


df.ac <- df.cluster %>% 
  mutate(is = str_match(description,"cetyltrans")) %>% 
  filter(!is.na(is)) %>% 
  filter(genename != "Acat2" & genename != "Acat1" & genename != "Dlat") %>% 
  mutate(type = "Writer-Pro-a")

df.dem <- df.cluster %>% 
  mutate(is = str_match(description,"demeth")) %>% 
  filter(!is.na(is)) %>% 
  mutate(type = "Eraser-m")

df.dea <- df.cluster %>% 
  mutate(is = str_match(description,"deace")) %>% 
  filter(!is.na(is)) %>% 
  filter(genename != "Ydjc") %>% 
  mutate(type = "Eraser-a")


t1 <- df.m %>% 
  bind_rows(df.ac) %>% 
  bind_rows(df.dem) %>% 
  bind_rows(df.dea) %>% 
  select(genename,type) %>% 
  left_join(to.plot)
t1$genename <- factor(t1$genename,levels = unique(t1$genename),ordered = T)




#target$epi <- data.frame(genename = unique(t1$genename))



no.heat <- t1 %>% 
  mutate(de = if_else(class == "Morula" & category == "ICM-EPI","yes","no" ))  %>% 
  filter(de  == "no") %>% 
  select(genename,sample,int) %>% spread(sample,int) 

to.heat <- data.frame(genename = c("Prmt5",
                                   "Mettl25",
                                   "Mettl18",
                                   "Kmt2d",
                                   "Pcmt1",
                                   "Setd1a",
                                   "Prmt1",
                                   "Nsd1",
                                   
                                   "Esco1",
                                   "Crebbp",
                                   "Naa30",
                                   "Naa50",
                                   "Naa15",
                                   "Hat1",
                                   "Naa25",
                                   
                                   "Dmap1",
                                   "Dnmt3a",
                                   "Dnmt3b",
                                   "Dnmt1",
                                   "Dnmt3l",
                               
                                   
                                   "Nsun2",
                                   "Rnmt",
                                   "Mettl5",
                                   "Trmt1",
                                   "Trmt2a",
                                   "Ftsj3",
                                   "Fbl",
                                   "Emg1",
                                   "Nop2",
                                   "Cmtr1",
                                   
                                   
                                   
                                   "Hdac1",
                                   "Sirt1",
                                   "Hdac6",
                                   "Sap18",
                                   
                                   "Kdm2b",
                                   "Kdm3a",
                                   "Kdm3b"
                                   
                                   
                                  
                           
                              
                                  
                                
                             
                                   
                           
                                   
                             )) %>% left_join(no.heat)

annot.sample <- sample %>% select(time,class,category) %>% as.data.frame()
row.names(annot.sample) <- sample$sample


heat <- to.heat %>% select(-genename) %>% as.matrix()
row.names(heat) <- to.heat$genename


file.g <- t1 %>% group_by(genename,type) %>% summarise() %>% ungroup() 
to.annot.g <- to.heat %>% select(genename) %>% left_join(file.g)
annot.g <-  to.annot.g %>% select(type) %>% as.data.frame()
row.names(annot.g) <- to.annot.g$genename

gap <- to.heat %>% select(genename) %>% left_join(file.g)
#heat[heat >3] <- 3


col_annot_obj <- columnAnnotation(df = annot.sample,col = annotcol)


col.g <- ggsci::pal_npg()(length(unique(file.g$type)))
names(col.g) <- unique(file.g$type)

row_annot_obj <- rowAnnotation(type = anno_simple(to.annot.g$type,col = col.g))


heat.a.type = Legend(labels = unique(file.g$type), title = "type",
                     legend_gp = gpar(fill = col.g))




p <- 
  Heatmap(heat,name = "Z-score",
          row_names_gp = gpar(col = c(rep("black",7),"red",
                                      rep("black",6),"red",
                                      rep("black",4),"red",
                                      rep("black",9),"red",
                                      rep("black",6),"red"
                                      )
                                      ),
             column_split = sample$category,column_title = NULL,
          row_split = gap$type,row_title = NULL,
             column_title_gp=grid::gpar(fontsize=10),
             top_annotation = col_annot_obj,right_annotation = row_annot_obj,
             cluster_rows = T,cluster_columns = F,show_column_names = F)

png("../result/3.dna_rna_epi/epi.all.png",units = "cm",width = 15,height = 20,res = 300)
draw(p, 
     
     annotation_legend_list = list(heat.a.type),
     merge_legend = TRUE,
     heatmap_legend_side = "right", annotation_legend_side = "right",
     padding = unit(c(2, 2, 2, 2), "mm"))
dev.off()


pdf("../result/3.dna_rna_epi/epi.all.pdf",width = 6,height = 8)
draw(p, 
     
     annotation_legend_list = list(heat.a.type),
     merge_legend = TRUE,
     heatmap_legend_side = "right", annotation_legend_side = "right",
     padding = unit(c(2, 2, 2, 2), "mm"))
dev.off()











#--------------------------------------------------------
#
# dataset global mitochondia
#
#--------------------------------------------------------


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




df.m <- df.cluster %>% 
  mutate(is = str_match(description,"itochon")) %>% 
  filter(!is.na(is)) %>%
  filter(eggnog_category == "Information storage and processing")  %>% 
  arrange(genename)



t1 <- df.m %>% 
  mutate(type  = str_replace(eggnog_term,",.+","")) %>% 
  mutate(type  = str_replace(type,"processing and ","")) %>% 
  select(genename,type)  %>% left_join(to.plot)
t1$genename <- factor(t1$genename,levels = unique(t1$genename),ordered = T)


target$mito <- data.frame(genename = unique(t1$genename))




no.heat <- t1 %>% 
  mutate(de = if_else(class == "Morula" & category == "ICM-EPI","yes","no" ))  %>% 
  filter(de  == "no") %>% 
  select(genename,sample,int) %>% spread(sample,int) 

to.heat <- data.frame(genename = c("Mrpl53",
                                   "Mrpl45",
                                   "Gfm1",
                                   "Mtrfr",
                                   "Mrps17",
                                   "Mrpl39",
                                   "Mrpl40",
                                   
                                   "Mrpl50",
                                   "Mrpl12",
                                   "Mrps30",
                                   "Timm50",
                                   "Mrpl18",
                                   "Mrps18c",
                                   "Mrps16",
                                   "C1qbp",
                                   "Tfam",
                                   "Ssbp1",
                                   "Tufm",
                                   "Mrps5",
                                   "Atp5if1",
                                   "Yars2",
                                   "Mrps34",
                                   "Tars2",
                                   "Lrpprc",
                                   "Mrpl11",
                                   "Mrpl43",
                                   "Slirp",
                                   "Mrps28",
                                   "Iars2",
                                   "Mrpl41",
                                   "Mrps23"
                                   
                                   )) %>% 
  left_join(no.heat)
  

annot.sample <- sample %>% select(time,class,category) %>% as.data.frame()
row.names(annot.sample) <- sample$sample


heat <- to.heat %>% select(-genename) %>% as.matrix()
row.names(heat) <- to.heat$genename


file.g <- t1 %>% group_by(genename,type) %>% summarise() %>% ungroup() 
to.annot.g <- to.heat %>% select(genename) %>% left_join(file.g)
annot.g <-  to.annot.g %>% select(type) %>% as.data.frame()
row.names(annot.g) <- to.annot.g$genename






col_annot_obj <- columnAnnotation(df = annot.sample,col = annotcol)



#col.g <- brewer.pal(length(unique(file.g$type)),"Set1")



col.g <- ggsci::pal_jco(alpha = 0.7)(length(unique(file.g$type)))
names(col.g) <- unique(file.g$type)

row_annot_obj <- rowAnnotation(type = anno_simple(to.annot.g$type,col = col.g))


heat.a.type = Legend(labels = unique(file.g$type), title = "type",
                     legend_gp = gpar(fill = col.g))




p <- 
  Heatmap(heat,name = "Z-score",
          column_split = sample$category,column_title = NULL,
          column_title_gp=grid::gpar(fontsize=10),
          # row_split = c(rep("3.5",4),
          #               rep("4.5",5),
          #               rep("5.0",9),
          #               rep("5.5",3),
          #               rep("6.5",7)),row_title = NULL,
          top_annotation = col_annot_obj,right_annotation = row_annot_obj,
          cluster_rows = F,cluster_columns = F,show_column_names = F)

png("../result/3.dna_rna_epi//mitochondira.png",units = "cm",width = 15,height = 20,res = 300)
draw(p, 
     
     annotation_legend_list = list(heat.a.type),
     merge_legend = TRUE,
     heatmap_legend_side = "right", annotation_legend_side = "right",
     padding = unit(c(2, 2, 2, 2), "mm"))
dev.off()



pdf("../result/3.dna_rna_epi//mitochondira.pdf",width = 6,height = 8)
draw(p, 
     
     annotation_legend_list = list(heat.a.type),
     merge_legend = TRUE,
     heatmap_legend_side = "right", annotation_legend_side = "right",
     padding = unit(c(2, 2, 2, 2), "mm"))
dev.off()








