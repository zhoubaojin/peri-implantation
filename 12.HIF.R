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
df.nor <- read.xlsx("Supplementary Table 1.xlsx",sheet = "NormalizedIntensity") 

df.q <- df.nor %>%  gather(sample,int,-ID) 
df.q <-  sample %>% left_join(df.q,by = "sample")
df.q$sample <- factor(df.q$sample,levels = unique(df.q$sample),ordered = T)

#--------------------------------------------------------
#
# dataset HIF
#
#--------------------------------------------------------


test <- df.q %>% filter(!is.na(time)) 

test.0 <- test %>% filter(category != "Morula")
test.1 <- test %>% filter(category == "Morula") %>% 
  mutate(category = "ICM-EPI")
test.2 <- test %>% filter(category == "Morula") %>% 
  mutate(category = "TE-EXE")

to.plot <- test.0 %>% 
  bind_rows(test.1) %>% 
  bind_rows(test.2  ) %>% 
  mutate(time = as.numeric(time))




df.m <- data.frame(genename = c("Aldoa",
                                "Eno3",
                                "Cul2",
                                "Plcg1",
                                "Pdha1",
                                "Pdhb",
                                "Mapk1",
                                "Map2k1"))


t1 <- df.m %>% left_join(df.name) %>% 
  left_join(to.plot)
t1$genename <- factor(t1$genename,levels = unique(t1$genename),ordered = T)



to.sc <- t1 %>% 
  mutate(de = if_else(class == "Morula" & category == "ICM-EPI","yes","no" ))  %>% 
  filter(de  == "no") %>% 
  select(genename,sample,int) %>% spread(genename,int) 

to.heat <- to.sc %>%select(-sample) %>% scale()%>% bind_cols(sample = to.sc$sample) %>% 
  gather(genename,int,-sample) %>% spread(sample,int)


annot.sample <- sample %>% select(time,class,category) %>% as.data.frame()
row.names(annot.sample) <- sample$sample


heat <- to.heat %>% select(-genename) %>% as.matrix()
heat[heat > 2] <- 2
row.names(heat) <- to.heat$genename




col_annot_obj <- columnAnnotation(df = annot.sample,col = annotcol)



p <- 
  Heatmap(heat,name = "Z-score",
          column_split = sample$category,column_title = NULL,
          column_title_gp=grid::gpar(fontsize=10),
          top_annotation = col_annot_obj,
          #right_annotation = row_annot_obj,
          cluster_rows = F,cluster_columns = F,show_column_names = F)

png("../result/3.two.stage//HIF.png",units = "cm",width = 20,height = 6,res = 300)
draw(p, 
     merge_legend = TRUE,
     heatmap_legend_side = "right", annotation_legend_side = "right",
     padding = unit(c(2, 2, 2, 2), "mm"))
dev.off()

pdf("../result/3.two.stage/HIF.pdf",width = 8,height = 2)
draw(p, 
     merge_legend = TRUE,
     heatmap_legend_side = "right", annotation_legend_side = "right",
     padding = unit(c(2, 2, 2, 2), "mm"))
dev.off()


#---- mrna ----

df <- openxlsx::read.xlsx("../raw/WT_EPI EXE samples_FPKM_matrix.xlsx") 

df.q.rna <- df %>% 
  gather(sample,fpkm,-X1) %>% 
  mutate(class = str_replace_all(sample,"-\\d$","")) %>% 
  mutate(time = as.numeric(str_match(class,"\\d\\.\\d"))) %>% 
  mutate(cate = str_match(class,"^\\w\\w"))  %>% 
  mutate(cate = if_else(cate == "IC","ICM-EPI",
                        if_else(cate == "Ep","ICM-EPI","TE-EXE")))

unique(df.q.rna$class)


df.rna.invivo <- df.q.rna %>% 
  filter(class != "Epi-5.0") %>% 
  filter(class != "Epi-5.25") %>% 
  filter(class != "Epi-5.5_invitro") %>% 
  filter(class != "ExE-5.5_invitro")
  

unique(df.rna.invivo$class)
unique(df.rna.invivo$cate)

target <- df.rna.invivo %>% filter(X1 == "Aldoa" | X1 == "Eno3")


p1 <- ggplot(target,aes(time,log2(fpkm),col  = cate)) +
  geom_point() +
  geom_smooth(span = 0.5) +
  facet_wrap(  X1 ~. ,scales = "free") +
  geom_vline(xintercept = 4.8,lty=3,col = "red") +
  labs(y = expression(log[2]*FPKM),col = "",x = "Time (day)")+
  scale_color_jco() +
  theme_bw() +
  theme(legend.position = "top",
        # legend.key.size = unit(0.1,"cm"),
        strip.background.x = element_blank()) 

#ggsave("../result/6.hif_val/mrna.png",width = 4,height = 3)




target <- df.q.rna %>% 
  mutate(int = log2(fpkm)) %>% 
  filter(X1 == "Aldoa" |
           X1 == "Eno3") %>% 
  filter(class == "Epi-5.0_invivo" |
           class == "Epi-5.0" |
           class == "Epi-5.25_invivo" | 
           class == "Epi-5.25" |
           class == "Epi-5.5_invitro" |
           class == "Epi-5.5_invivo" ) %>% 
  mutate(class = if_else(class == "Epi-5.0","Epi-5.0_invitro",
                         if_else(class == "Epi-5.25","Epi-5.25_invitro",class))) %>% 
  mutate(exp = str_match(class,"in\\w+"))



library(ggpubr)
library(rstatix)

bxp <- ggboxplot(
  target, x = "time", y = "int", facet.by = "X1",
  color = "exp", palette = "jco"
) 

# Add p-values onto the box plots
stat.test <- target %>%
  group_by(X1,time) %>%
  t_test(int ~ exp) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")

stat.test <- stat.test %>%
  add_xy_position(x = "time", dodge = 0.8)

p2 <- bxp + stat_pvalue_manual(
  stat.test,  label = "p", tip.length = 0
) +
  facet_wrap(  X1 ~. ,scales = "free") +
  scale_y_continuous(expand = expansion(mult = c(0.1,0.2))) +
  labs(color = "",y = expression(log[2]*FPKM),x = "Time (day)") +
  theme_bw() +
  theme(legend.position = "top",
        # legend.key.size = unit(0.1,"cm"),
        strip.background.x = element_blank()) 
library(patchwork)
p1 +p2

ggsave("../result/6.hif_val/mrna.png",width = 8,height = 3.5)
ggsave("../result/6.hif_val/mrna.pdf",width = 8,height = 3.5)





