#--------------------------------------------------
# function network
# 1. 10 clusters
# 2. combine key bp and kegg, choose key pathway and genes
# 3. add tf and regulator information
#---------------------------------------------------

library(tidyverse)
library(RColorBrewer)
load("./mouse.uniprot.function.rda")

df.id <- fun$protein %>% mutate(ID = protein)

sample <- openxlsx::read.xlsx("Supplementary Table 1.xlsx",sheet = "FileInformation")
sample$category <- factor(sample$category,levels = unique(sample$category),ordered = T)


col1 <- colorRampPalette(brewer.pal(8, "Spectral"))
col.cluster <- col1(10)
#"#D53E4F" "#ED6245" "#F99153" "#FDBE6E" "#FBE28C" 、
#"#E8F296" "#BEE5A0" "#8CD1A4" "#5AB5AA" "#3288BD"

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



#------------------
#
# 2. cluster sub-cellular
#
#------------------
out.cluster <- read.xlsx("Supplementary Table 2.xlsx")


to.fun <-  out.cluster %>% 
  mutate(term = subcellularlocation) %>%
  filter(!is.na(term)) %>% 
  mutate(cluster = str_replace(paste0("C",cluster+100),"C1","C")) %>% 
  group_by(cluster,term) %>% summarise(num.target=n(),gene = paste0(genename,collapse = ";")) 


to.plot <- out.cluster %>%
  mutate(cluster = str_replace(paste0("C",cluster+100),"C1","C")) %>% 
  group_by(cluster ) %>% summarise(num.cl=n()) %>% 
  left_join(to.fun,  ) %>% mutate(per = num.target /num.cl) %>% 
  mutate(label = paste0(round(per*100,1),"%"))



ggplot(to.plot,aes(cluster,per,fill = term)) +
  geom_bar(width = 0.9, stat="identity",col = "white") + 
  coord_polar(theta = "y") +
  xlab("") + ylab("") +
  scale_fill_brewer(palette = "Paired")+
  geom_text(data = to.plot, hjust = 1, size = 3,
            aes(x = cluster, y = 0, label = cluster)) +
  geom_text(aes(label = label),position=position_stack(vjust = 0.5),size = 1) +
  # guides(fill=guide_legend(nrow=3)) +
  theme_void()+
  labs(fill = "") +
  theme(legend.position = "right"#,
        #  legend.key.size = unit(0.1,"cm")
  )



ggsave("../result/2.cluster/10.cluster/cluster.sub.png",width = 5,height = 4)
ggsave("../result/2.cluster/10.cluster/cluster.sub.pdf",width = 5,height = 4)
#-----
#
# overview of kegg
#
#-----

df.kegg <- openxlsx::read.xlsx("Supplementary Table 2.xlsx",sheet = "ClusterKEGG")

# dotplot
kegg.num <- df.kegg %>%  
  mutate(cluster = str_replace(paste0("C",cluster +100),"C1","C"))  %>% 
  filter(pvalue < 0.01) %>% 
  select(term,cluster,num.target) 
 
kegg.heat <- df.kegg %>% mutate(log = -log10(pvalue)) %>% 
  filter(level1 != "Human Diseases") %>% 
  filter(pvalue <0.01) %>%  
  mutate(cluster = str_replace(paste0("C",cluster +100),"C1","C"))  %>% 
  select(level1,level2,termid,term,log,cluster) %>% 
  spread(cluster,log)  %>% 
  gather(cluster,log,C01:C10) %>% left_join(kegg.num,by = c("term","cluster")) %>% 
  arrange(desc(termid))
kegg.heat$term<- factor(kegg.heat$term,levels = unique(kegg.heat$term),ordered = T)

t <- "Metabolism"
df.p1 <- kegg.heat %>% filter(level1 == t) 
p1 <-
ggplot(df.p1) +
  geom_point(aes(cluster,term,col = log,size = num.target),alpha = 0.5) +
  facet_grid(level1 ~.,scale= "free_y",space = "free_y") +
  theme_minimal() +
  labs(size = "Protein count",col = expression(-log[10]*italic(P)*val),title = t)+
  scale_color_distiller(palette = "Spectral",direction = -1) +
  theme(axis.title = element_blank(),strip.text = element_blank())


t <- "Genetic Information Processing"
df.p2 <- kegg.heat %>% filter(level1 == t) 
p2 <-
  ggplot(df.p2) +
  geom_point(aes(cluster,term,col = log,size = num.target),alpha = 0.5) +
  facet_grid(level1 ~.,scale= "free_y",space = "free_y") +
  theme_minimal() +
  labs(size = "Protein count",col = expression(-log[10]*italic(P)*val),title = t)+
  scale_color_distiller(palette = "Spectral",direction = -1) +
  theme(axis.title = element_blank(),strip.text = element_blank())



t <- "Environmental Information Processing"
df.p3 <- kegg.heat %>% filter(level1 == t) 
p3 <-
  ggplot(df.p3) +
  geom_point(aes(cluster,term,col = log,size = num.target),alpha = 0.5) +
  facet_grid(level1 ~.,scale= "free_y",space = "free_y") +
  theme_minimal() +
  labs(size = "Protein count",col = expression(-log[10]*italic(P)*val),title = t)+
  scale_color_distiller(palette = "Spectral",direction = -1) +
  theme(axis.title = element_blank(),strip.text = element_blank())

t <- "Cellular Processes"
df.p4 <- kegg.heat %>% filter(level1 == t) 
p4 <-
  ggplot(df.p4) +
  geom_point(aes(cluster,term,col = log,size = num.target),alpha = 0.5) +
  facet_grid(level1 ~.,scale= "free_y",space = "free_y") +
  theme_minimal() +
  labs(size = "Protein count",col = expression(-log[10]*italic(P)*val),title = t)+
  scale_color_distiller(palette = "Spectral",direction = -1) +
  theme(axis.title = element_blank(),strip.text = element_blank())


t <- "Organismal Systems"
df.p5 <- kegg.heat %>% filter(level1 == t) 
p5 <-
  ggplot(df.p5) +
  geom_point(aes(cluster,term,col = log,size = num.target),alpha = 0.5) +
  facet_grid(level1 ~.,scale= "free_y",space = "free_y") +
  theme_minimal() +
  labs(size = "Protein count",col = expression(-log[10]*italic(P)*val),title = t)+
  scale_color_distiller(palette = "Spectral",direction = -1) +
  theme(axis.title = element_blank(),strip.text = element_blank())




library(patchwork)
p1 /p2/p3/p4 /p5 +plot_layout(heights = c(length(unique(df.p1$term)),
                                      length(unique(df.p2$term)),
                                      length(unique(df.p3$term)),
                                      length(unique(df.p4$term)),
                                      length(unique(df.p4$term))))

ggsave("./../result/2.cluster/cluster.kegg.dot.png",width = 8,height = 24)
ggsave("./../result/2.cluster/cluster.kegg.dot.pdf",width = 8,height = 24)


# heatmap
kegg.heat <- df.kegg %>% mutate(log = -log10(pvalue)) %>% 
  filter(pvalue <0.01) %>% 
  mutate(cluster = str_replace(paste0("C",cluster +100),"C1","C")) %>% 
  select(level1,level2,termid,term,log,cluster) %>% 
  spread(cluster,log) 


to.heat <- fun$term$kegg_pathway %>% arrange(level1id,level2id) %>% ungroup() %>% 
  select(termid) %>% left_join(kegg.heat) %>%filter(!is.na(level1)) %>% 
  filter(level1 != "Human Diseases")

to.heat$level1 <- factor(to.heat$level1,levels = unique(to.heat$level1),ordered = T)
to.heat$level2 <- factor(to.heat$level2,levels = unique(to.heat$level2),ordered = T)

gap <- to.heat %>% group_by(level1) %>% summarise(n = n()) %>% mutate(sum = cumsum(n))
annot.path <- to.heat %>% select(level1,level2) %>% as.data.frame()
row.names(annot.path) <- to.heat$term

heat <- to.heat %>% select(-level1,-level2,-term,-termid) %>% as.matrix()
row.names(heat) <- to.heat$term
heat[is.na(heat)] <-0
heat[heat >10] <- 10

heatcol <- colorRampPalette(brewer.pal(6,"Paired"))(100)
col.level1 <- brewer.pal(length(unique(to.heat$level1)),"Set1")
names(col.level1) <- unique(to.heat$level1)
col.level2 <- colorRampPalette(brewer.pal(8,"Set2"))(length(unique(to.heat$level2)))
names(col.level2) <- unique(to.heat$level2)

col.annol <- list(level1 = col.level1,
                  level2 = col.level2)
p <- pheatmap::pheatmap(heat,cluster_rows = F,cluster_cols = F,color = heatcol,
                   annotation_colors = col.annol,
         annotation_row = annot.path,
         gaps_row = gap$sum)
pp <-ggplotify::as.ggplot(p)
ggsave(pp,file = "../result/2.cluster/cluster.pathway.overview.png",width = 10,height = 10)





for(path.c in unique(to.heat$level1) ){
  message(path.c)
  to.heat.i <- to.heat %>% filter(level1 == path.c)

  
  annot.path <- to.heat.i %>% select(level2) %>% as.data.frame()
  row.names(annot.path) <- to.heat.i$term
  
  heat <- to.heat.i %>% select(-level1,-level2,-term,-termid) %>% as.matrix()
  row.names(heat) <- to.heat.i$term
  heat[is.na(heat)] <-0
  heat[heat >10] <- 10
  
  heatcol <- colorRampPalette(brewer.pal(6,"Paired"))(100)
  names(col.level1) <- unique(to.heat.i$level1)
  col.level2 <- colorRampPalette(brewer.pal(8,"Set2"))(length(unique(to.heat.i$level2)))
  names(col.level2) <- unique(to.heat.i$level2)
  
  col.annol <- list(level2 = col.level2)
  p <- pheatmap::pheatmap(heat,cluster_rows = F,cluster_cols = F,
                          annotation_names_row = F,
                          annotation_colors = col.annol,
                          annotation_row = annot.path)
  pp <-ggplotify::as.ggplot(p)
  
  ggsave(pp,file = paste0("../result/2.cluster/cluster.pathway.",path.c,".png"),
                          width = 10,height = 2+ nrow(heat) *0.1)
  
}

#-----
#
# network. to cytoscape
#
#------

df.cl <-  openxlsx::read.xlsx("Supplementary Table 2.xlsx")
df.go <- openxlsx::read.xlsx("Supplementary Table 2.xlsx",sheet = "ClusterGO")
df.kegg <- openxlsx::read.xlsx("Supplementary Table 2.xlsx",sheet = "ClusterKEGG")



# out
library(openxlsx)
wb<- openxlsx::createWorkbook()

out.kegg.1 <- df.kegg %>% 
  filter(pvalue < 0.05) %>% 
  dplyr::select(cluster,level1,level2,term,pvalue) %>% 
  spread(cluster,pvalue)

out.kegg.2 <- df.kegg %>% 
  filter(pvalue < 0.05) %>% 
  dplyr::select(cluster,level1,level2,term,num.target) %>% 
  spread(cluster,num.target)

addWorksheet(wb,"Pvalue")
writeData(wb,"Pvalue",out.kegg.1)

addWorksheet(wb,"ProteinNumber")
writeData(wb,"ProteinNumber",out.kegg.2)

#saveWorkbook(wb,file = "../result/2.cluster/10.cluster/kegg.information.xlsx")


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

df.type <- df.tf %>% bind_rows(df.pk) %>% bind_rows(df.re) 
df.type.cluster <- df.type %>% select(type,ID) %>% inner_join(df.cl,by = "ID" ) %>% 
  arrange(cluster)




df.c1.k <- df.kegg %>% 
  filter(term == "Peroxisome" |
           term == "Pyruvate metabolism" |
           term == "Amino sugar and nucleotide sugar metabolism" 
           ) %>% 
  filter(cluster == 1 | cluster ==2)

df.c2.k <-df.kegg %>% 
  filter(term == "Glycolysis / Gluconeogenesis"  |
           term == "Citrate cycle (TCA cycle)" 
  ) %>% filter(cluster == 2 )

df.c2.g <- df.go %>% 
  filter(term == "positive regulation of embryonic development") %>% 
  filter(cluster == 2)


df.c3.k <- df.kegg %>% 
  filter(term == "Ribosome" |
           term == "Spliceosome" |
           term == "Proteasome") %>% 
  filter(cluster == 3  |cluster == 4 |cluster == 5 |cluster == 10)
  





df.c4.k <- df.kegg %>% 
  filter(term == "Thermogenesis" |
           term == "Spliceosome" | 
           term == "Oxidative phosphorylation" |
           term == "Focal adhesion") %>% 
  filter(cluster == 4)
 
df.c4.g <- df.go %>% 
  filter(term == "aerobic electron transport chain") %>% 
  filter(cluster == 4)


df.c6.k <-df.kegg %>% 
  filter(
         term == "Citrate cycle (TCA cycle)" 
  ) %>% 
  filter(cluster == 6)
df.c6.g <-df.go %>% 
  filter(
    term == "carbohydrate catabolic process" 
  ) %>% 
  filter(cluster == 6)


# Nucleocytoplasmic transport
# mRNA surveillance pathway
#

# df.c6.k <-df.kegg %>% 
#   filter(
#     term == "Citrate cycle (TCA cycle)" 
#   )

df.c7.k <-df.kegg %>% 
  filter(
    term == "Cell cycle"  |
      term == "Cysteine and methionine metabolism"
  ) %>% 
  filter(cluster == 7)

df.c7.g <-df.go %>% 
  filter(
      term == "response to light stimulus" |
      term == "steroid hormone mediated signaling pathway"
  ) %>% 
  filter(cluster == 7)

df.c8.k <-df.kegg %>% 
  filter(
    term == "HIF-1 signaling pathway"  |
      term == "Regulation of actin cytoskeleton" |
      term == "Oxytocin signaling pathway"
  )%>% 
  filter(cluster == 8)


df.c9.g <- df.go %>% 
  filter(term == "maternal placenta development" |
           term== "Decidualization") %>% 
  filter(cluster == 9)


df.c9.k <-df.kegg %>% 
  filter(
    term == "Tight junction"
  ) %>% 
  filter(cluster == 9)



df.c10.k <-df.kegg %>% 
  filter(
    term == "HIF-1 signaling pathway"
  )%>% 
  filter(cluster == 10)

df.c10.g <- df.go %>% 
  filter(
    term == "spindle assembly" |
    term == "mitotic sister chromatid segregation"
  ) %>% 
  filter(cluster == 10)


test <- df.go %>% filter(cluster == 10) %>%
arrange(pvalue)



#------------------------------------------
# out
#------------------------------------------


target.term <- df.c1.k %>% 
  bind_rows(df.c2.k) %>% 
  bind_rows(df.c2.g) %>% 
  bind_rows(df.c3.k) %>% 
  bind_rows(df.c4.k) %>% 
  bind_rows(df.c4.g) %>% 
  bind_rows(df.c6.k) %>% 
  bind_rows(df.c7.g) %>% 
  bind_rows(df.c7.k) %>% 
  bind_rows(df.c7.g) %>% 
  bind_rows(df.c8.k) %>% 
  
  bind_rows(df.c9.k) %>% 
  bind_rows(df.c9.g) %>% 
  bind_rows(df.c10.k)   %>% 
  bind_rows(df.c10.g)

  
# net 1, cluster to term
net.1 <- target.term %>%
  mutate(source = str_replace(paste0("C",cluster +100),"C1","C"),target = term) %>% 
  group_by(source,target) %>% summarise()


# net 2, cluster to gene
net.2 <- target.term %>% separate_rows(gene,sep = ";") %>% 
  mutate(source = str_replace(paste0("C",cluster +100),"C1","C"),target = gene) %>% 
  group_by(source,target) %>% summarise()
  


# net 3, term to gene
net.3 <- target.term %>% separate_rows(gene,sep = ";") %>% 
  mutate(source = term,target = gene) %>% 
  group_by(source,target) %>% summarise()

net <- net.1 %>% bind_rows(net.2,net.3)  



# node
# color: cluster information
# shape: cluster, go, kegg, protein (kinase, receptor, tf)
# size: 
#   cluster: 100, 
#   term: -log10(pvalue)
#   protein: 50
# 

node.term <- target.term %>% group_by(term) %>% top_n(n = 1,wt = -pvalue)  %>% ungroup() %>%
  mutate(key = term,
         size = -log10(pvalue),
         col = str_replace(paste0("C",cluster +100),"C1","C"),
         type = if_else(is.na(level1),"GO","KEGG")) %>% 
  select(key,col,size,type)



#df.type <- read_delim("../result/")
#df.pro <- net.2 %>% group_by(source,target) %>% summarise() %>% ungroup() %>% 
#  left_join(fun$protein,by = c("target" = "genename")) %>% 
#  mutate(type = if_else())

node.pro <- net.2 %>% 
  mutate(key = target,col = source) %>% 
  group_by(key,col) %>% summarise() %>% ungroup() %>% 
  mutate(type = "protein",size = 2)



node.cluster <- target.term %>%
  mutate(key = str_replace(paste0("C",cluster +100),"C1","C")) %>% 
  group_by(key) %>% summarise() %>% ungroup()%>% 
  mutate(col = key,
         size = 0.5,
         type = "cluster")

node <- node.pro %>% bind_rows(node.term,node.cluster)

write.table(x = node,file = "../result/2.cluster/10.cluster/function/network/to.cytocape.node.txt",quote = F,sep = "\t",row.names = F)
write.table(x = net,file = "../result/2.cluster/10.cluster/function/network/to.cytocape.net.txt",quote = F,sep = "\t",row.names = F)


node.fig <- data.frame(key = str_replace(paste0("C",1:10 +100),"C1","C")) %>% 
  mutate(fig = key)


node.label.col <- node.pro %>% 
  select(key) %>% 
  inner_join(df.type.cluster ,by = c("key"= "genename")) %>% 
  mutate(procol= type) %>%
  select(key,procol)

write.table(x = node.fig,file = "../result/2.cluster/10.cluster/function/network/to.cytocape.node.fig.txt",quote = F,sep = "\t",row.names = F)
write.table(x = node.label.col,file = "../result/2.cluster/10.cluster/function/network/to.cytoscale.node.procol.txt",quote = F,sep = "\t",row.names = F)

