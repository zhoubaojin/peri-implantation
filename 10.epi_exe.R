# 1. recalculated protein in two linage
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
library(ggsci)


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


#--------------------------------
#
# all protein 2 kegg
#
#--------------------------------

coma <- "ICM-EPI"
comb <- "TE-EXE"

to.stat <- df.nor %>% select(ID,sample$sample) %>% 
  gather(sample,intensity,-ID) %>% 
  left_join(sample) %>% 
  filter(category == coma | category == comb) %>% 
  group_by(ID) %>% summarise(pvalue.t.test = t.test(intensity ~ category)$p.value,
                             pvalue.w.test = wilcox.test(intensity ~ category)$p.value)


to.vol <- df.nor %>% select(ID,sample$sample) %>% 
  gather(sample,intensity,-ID) %>% 
  left_join(sample) %>% 
  filter(category == coma | category == comb) %>% 
  group_by(ID,category) %>% summarise(mean = mean(intensity)) %>% spread(category,mean) %>% ungroup() %>% 
  mutate(ratio = `TE-EXE` / `ICM-EPI`) %>% 
  left_join(to.stat) %>% 
  mutate(regulated = if_else(pvalue.t.test < 0.01 & ratio > 1 ,comb,
                             if_else(pvalue.t.test < 0.01 & ratio < 1,coma,"None"))) %>% 
  left_join(df.id)


to.plot <- to.vol %>% group_by(regulated) %>% summarise(n = n()) %>% mutate(reg = paste0(regulated,"(n=",n,")") ) %>% 
  left_join(to.vol)

df.label <- to.vol  %>% filter(regulated != "None")
df.label.tf <- to.vol %>% filter(regulated != "None") %>%  filter(!is.na(tf_category) |
                                                                    !is.na(kinase_category))

ggplot(to.plot,aes(log2(ratio),-log10(pvalue.t.test))) +
  geom_point(aes(fill = reg,size = reg),shape = 21) +
  geom_text_repel(aes(label = genename),data = df.label,size = 3) +
  geom_point(data = df.label.tf,size = 3,shape = 21,col = "red") +
  scale_size_manual(values = c(2,0.1,2)) +
  scale_fill_manual(values = c("#FC8D62","grey","#8DA0CB")) +
  geom_vline(xintercept = c(-1,1),lty = 3,col = "red") +
  geom_hline(yintercept = -log10(0.01),lty = 3,col  = "red") +
  labs(x= expression(log[2]*'('*Fold*' '*change*')'),y = expression(-log[10]*'('*italic(Pvalue)*')'),
       fill = "",size = "",shape = "")+
  theme_bw() +
  theme(legend.position = "top")

#ggsave("../result/3.EPI_EXE.vol.png",width = 4,height = 4)







#tarppt#target 
df.target <- df.label %>% filter(regulated == coma) %>% 
  top_n(n =2,wt = -pvalue.t.test) %>% 
  select(ID,pvalue.t.test) %>% left_join(df.nor) %>% left_join(df.id)  %>% arrange(pvalue.t.test) 

df.target$genename <- factor(df.target$genename,levels = unique(df.target$genename),ordered = T)

p1 <- df.target %>% ggplot(aes(category,int)) +
  geom_boxplot(outliers = F) +
  geom_jitter(aes(col = time)) +
  facet_wrap(genename~.,scales = "free_y") +
  scale_y_continuous(expand = c(0.2,0.2) ) +
  scale_color_manual(values = col.time) +
  stat_compare_means(comparisons = list(c(coma,comb)),method = "t.test") +
  labs(x="",y = "Intensity") +
  theme_bw() +
  theme(legend.position = "none"
  )


df.target <- df.label %>% filter(regulated == comb) %>% 
  top_n(n =2,wt = -pvalue.t.test) %>% 
  select(ID,pvalue.t.test) %>% left_join(df.nor) %>% left_join(df.id)  %>% arrange(pvalue.t.test) 

df.target$genename <- factor(df.target$genename,levels = unique(df.target$genename),ordered = T)

p2 <- df.target %>% ggplot(aes(category,int)) +
  geom_boxplot(outliers = F) +
  geom_jitter(aes(col = time)) +
  facet_wrap(genename~.,scales = "free_y") +
  scale_y_continuous(expand = c(0.2,0.2) ) +
  stat_compare_means(comparisons = list(c(coma,comb)),method = "t.test") +
  scale_color_manual(values = col.time) +
  labs(x="",y = "Intensity") +
  theme_bw() +
  theme(legend.position = "none"
  )

p1 /p2
ggsave("../result/3.EPI_EXE.boxplot.png",width = 5,height = 5)



wb <- createWorkbook()
addWorksheet(wb,"DEP")
writeData(wb,"DEP",to.vol)




#------------------------
# kegg function
#------------------------

#---- fun -------

num.kegg <- fun$bgnum %>% filter(database == "kegg" ) %>% select(bg) %>% as.numeric()
num.reg <- df.label %>% group_by(regulated) %>% summarise(num.sig = n()) %>% ungroup() 

df.fun <- df.label %>% mutate(term = kegg_pathway) %>% 
  separate_rows(term,sep=";") %>% 
  filter(!is.na(term)) %>% 
  group_by(term,regulated) %>% summarise(num.target = n(),sig.protein = paste(genename,collapse = ";"))  %>% ungroup() %>% 
  separate(term,into =c("l1","l2","l3"),sep = "\\|\\|") %>% separate(l3,into= c("termid","term"),sep = "\\:\\:") %>% 
  left_join(num.reg) %>% 
  select(regulated,termid,sig.protein,num.target,num.sig) %>% 
  left_join(fun$term$kegg_pathway,by = "termid") %>% 
  group_by(regulated,level1id,level1,level2id,level2,termid,term,sig.protein,num.target) %>% 
  summarise(pvalue = phyper(lower.tail = F,log.p = F,num.target -1,bgcount,(num.kegg - bgcount),num.sig))

to.plot <- df.fun %>% filter(pvalue < 0.05) %>% 
  filter(level1 != "Human Diseases") %>% 
  mutate(log = if_else(regulated == coma,log10(pvalue),-log10(pvalue)),
         just =if_else(regulated == coma,1,0 )) %>% 
  arrange(regulated,log) 

to.plot$term<- factor(to.plot$term,levels = unique(to.plot$term),ordered = T)
ggplot(to.plot,aes(log,term)) +
  geom_vline(xintercept = c(-log10(0.05),log10(0.05)),lty = 3,col = "red") +
  geom_bar(aes(fill = regulated),stat = "identity") +
  geom_text(aes(x = 0,y = term,label = sig.protein,hjust = just),size = 3)+
  scale_x_continuous(labels = abs) +
  geom_vline(xintercept = c(-log10(0.05),log10(0.05)),lty = 3,col = "red") +
  scale_fill_manual(values = col.samplecategory) +
  labs(fill = "",y = "",x = expression(-log[10]*'('*italic(Pvalue)*')')) +
  theme_bw()+
  theme(legend.position = "top")

ggsave("../result/3.two.stage/3.EPI_EXE.fun.kegg.png",width = 6,height = 5)



addWorksheet(wb,"KEGG")
writeData(wb,"KEGG",df.fun)





#-------- go ---------
num.go <- fun$bgnum %>% filter(database == "go" ) %>% select(bg) %>% as.numeric()
num.reg <- df.label %>% group_by(regulated) %>% summarise(num.sig = n()) %>% ungroup() 

df.fun <- to.vol %>% 
  filter(regulated != "None") %>% mutate(term = go_biologicalprocess) %>% 
  separate_rows(term,sep=";") %>% 
  filter(!is.na(term)) %>% 
  group_by(term,regulated) %>% summarise(num.target = n(),sig.protein = paste(genename,collapse = ";"))  %>% ungroup() %>% 
  #separate(term,into =c("l1","l2","l3"),sep = "\\|\\|") %>% 
  separate(term,into= c("termid","term"),sep = "\\:\\:") %>% 
  left_join(num.reg) %>% 
  select(regulated,termid,sig.protein,num.target,num.sig) %>% 
  left_join(fun$term$go,by = "termid") %>% 
  group_by(regulated,termid,term,sig.protein,num.target,bgcount) %>% 
  summarise(pvalue = phyper(lower.tail = F,log.p = F,num.target ,bgcount,(num.kegg - bgcount),num.sig))



addWorksheet(wb,"GO")
writeData(wb,"GO",df.fun)
saveWorkbook(wb,file = "../result/3.two.stage/1.EPI_EXE.xlsx",overwrite = T)

df.fun.sig <- df.fun %>% filter(pvalue < 0.05)

to.plot <- df.fun.sig %>% 
  filter(term == "regulation of transepithelial transport" |
           term == "regulation of cellular component size" |
           term == "centrosome duplication" |
           term == "cell morphogenesis" |
           term == "actin filament organization" |
           term == "stem cell development" |
           term == "cerebral cortex development" |
           term == "Okazaki fragment processing involved in mitotic DNA replication" |
           term == "cell morphogenesis involved in differentiation" |
           term == "positive regulation of nodal signaling pathway involved in determination of lateral mesoderm left/right asymmetry" |
           term == "nucleobase-containing compound metabolic process" |
           term == "regulation of vitamin D receptor signaling pathway" 
         
  ) %>% 
  mutate(term = str_replace_all(term,
                                "positive regulation of nodal signaling pathway involved in determination of lateral mesoderm left/right asymmetry",
                                "pos. reg. of nodal signaling path. invol.\ndetermination of lateral mesoderm left/right asymmetry")) %>% 
  mutate(term = str_replace_all(term,
                                "Okazaki fragment processing involved in mitotic DNA replication",
                                "Okazaki frag. processing invol. mitotic DNA replication")) %>% 
  
  mutate(log = if_else(regulated == coma,-log10(pvalue),-log10(pvalue))) %>% 
#         just =if_else(regulated == coma,0,1 )) %>% 
  mutate(fil = if_else(regulated == "ICM-EPI" & term == "centrosome duplication",1,0)) %>%  filter(fil == 0) %>% 
  mutate(fil = if_else(regulated == "TE-EXE" & term == "cellular nitrogen compound biosynthetic process",1,0)) %>%  filter(fil == 0) %>% 
  arrange(regulated,log)  
  # mutate(sig.protein = str_replace(sig.protein,
  #                                  "Umps;Celf1;Lig1;Gmps;Pcbp2;Smad2;Trim24;Kdm3b;Rrp1b;Igf2bp3;Xrcc3",
  #                                  "Umps;Celf1;Lig1;Gmps;Pcbp2;Smad2;\nTrim24;Kdm3b;Rrp1b;Igf2bp3;Xrcc3")) %>%
  # mutate(sig.protein = str_replace(sig.protein,
  #                                  "Plxnb2;Krt8;Actb;Slc9a3r1;Spr;Strip1;Shtn1;Uchl1;Coro1b",
  #                                  "Plxnb2;Krt8;Actb;Slc9a3r1;Spr;\nStrip1;Shtn1;Uchl1;Coro1b"))

to.plot$term<- factor(to.plot$term,levels = unique(to.plot$term),ordered = T)
ggplot(to.plot,aes(log,term)) +
  geom_bar(aes(fill = regulated),stat = "identity",alpha = 0.5) +
  geom_text(aes(x = 0,y = term,label = term),hjust = 0,size = 3)+
  scale_fill_manual(values = col.samplecategory) +
  scale_color_manual(values = col.samplecategory) +
  facet_grid(regulated ~. ,scale = "free_y",space = "free_y") +
  labs(fill = "",col = "",y = "",x = expression(-log[10]*'('*italic(Pvalue)*')')) +
  theme_bw()+
  theme(legend.position = "none",panel.grid = element_blank(),
        strip.background = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y  = element_blank())

ggsave("../result/3.two.stage/1.EPI_EXE.bp.png",width = 4,height = 4)
ggsave("../result/3.two.stage/1.EPI_EXE.bp.pdf",width = 4,height = 4)




#-------------------------------------------
#
# heatmap
#
#-------------------------------------------
library(ComplexHeatmap)

df.sig <- to.vol %>% filter(regulated != "None") %>%
  arrange(regulated,pvalue.t.test)

  gather(sample,int,-ID) %>% inner_join(df.sig,by = "ID") %>% 
  select(genename,sample,int) %>% spread(genename,int) %>% 
  left_join(sample,by = "sample") %>% 
  filter(class != "Morula" ) %>% arrange(category)

to.heat <- sig.to.scale %>% select(df.sig$genename) %>% scale() %>% as.data.frame()
rownames(to.heat) <- sig.to.scale$sample




# annot sample
sample2col <- sample %>% arrange(category,time) %>% filter(category != "Morula" )
col.samplecategory <- brewer.pal(3,"Set2")[-1]
names(col.samplecategory) <- unique(sample2col$category)

col.time <- brewer.pal(length(unique(sample2col$time)), "YlOrRd")
names(col.time) <- unique(sample2col$time)

col.sampleclass <- c(brewer.pal(10, "Paired"))
names(col.sampleclass) <- c(sort(as.character(unique(sample2col$class))))
annotcol <- list(class = col.sampleclass,
                 category = col.samplecategory,
                 time = col.time)

annot.sample <- data.frame(sample = colnames(to.heat)) %>% left_join(sample) %>% select(category,time,class) %>% as.data.frame()
row.names(annot.sample) <- colnames(to.heat)

row_annot_obj <- rowAnnotation(category = anno_simple(as.character(sample2col$category),col = col.samplecategory),
                               time = anno_simple(as.character(sample2col$time),col = col.time),
                               class = anno_simple(as.character(sample2col$class),col = col.sampleclass),
                               annotation_name_side = "top" )



col.type <- ggsci::pal_jco()(2)
# annot gene
col_annot_obj <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = col.type),
                                                    labels = c("ICM-EPI enriched (n=21)", "TE-EXE enriched (n=35)"), 
                                                    labels_gp = gpar(col = "white", fontsize = 10)))



col.pvalue <- circlize::colorRamp2(c( 1, 2, 3, 4), c("green", "yellow", "orange","red")) 

heat.a.p = Legend(title = "P-value", col_fun = col.pvalue, at = c(1, 2, 3, 4), 
                  labels = c("0.1", "0.01", "0.001", "0.0001"))
heat.a.cate = Legend(labels = unique(sample2col$category), title = "Lineage",
                    legend_gp = gpar(fill = col.samplecategory),ncol = 1)
heat.a.time = Legend(labels = unique(sample2col$time), title = "Time",
                      legend_gp = gpar(fill = col.time),ncol = 1)
heat.a.class = Legend(labels = sort(unique(as.character(sample2col$class))), title = "Class",
                      legend_gp = gpar(fill = col.sampleclass),nrow = 2)

heat.annot.bottom  = HeatmapAnnotation(
  `pvalue` = anno_simple(-log10(df.sig$pvalue.t.test),col = col.pvalue),
  annotation_name_side = "left"
)



heat.1 <- ComplexHeatmap::Heatmap(to.heat,
                                  cluster_rows = F,cluster_columns = F,
                                  column_split = df.sig$regulated,
                                  column_title = NULL,
                                  show_row_names = F,
                                  left_annotation= row_annot_obj,
                                  row_split = sample2col$category,row_title = NULL,
                                  top_annotation = col_annot_obj,
                                  bottom_annotation = heat.annot.bottom,
                                  heatmap_legend_param=list(title="z-score"))


png("../result/3.two.stage/1.two.leneage.heatmap.png",units = "cm",res = 300,height = 15,width = 30)
draw(heat.1, 
     annotation_legend_list = list(heat.a.p,
                                   heat.a.cate,
                                   heat.a.time,
                                   heat.a.class),
     merge_legend = TRUE,heatmap_legend_side = "bottom", annotation_legend_side = "bottom",
     padding = unit(c(3, 3, 3, 3), "mm"))


dev.off()
pdf("../result/3.two.stage/1.two.leneage.heatmap.pdf",height = 5,width = 10)
draw(heat.1, 
     annotation_legend_list = list(heat.a.p,
                                   heat.a.cate,
                                   heat.a.time,
                                   heat.a.class),
     merge_legend = TRUE,heatmap_legend_side = "bottom", annotation_legend_side = "bottom",
     padding = unit(c(3, 3, 3, 3), "mm"))


dev.off()



#------------------------
# validation with mRNA
#------------------------
df.sig.pro <- to.vol %>% filter(regulated != "None")


df <- openxlsx::read.xlsx("../raw/WT_EPI EXE samples_FPKM_matrix.xlsx") 
df.q.rna <- df %>% 
  gather(sample,fpkm,-X1) %>% 
  mutate(class = str_replace_all(sample,"-\\d$","")) %>% 
  mutate(cate = str_match(class,"^\\w\\w"))  %>% 
  mutate(cate = if_else(cate == "IC","ICM-EPI",
                        if_else(cate == "Ep","ICM-EPI","TE-EXE")))
  

df.ratio.rna <- 
  df.q.rna  %>% 
  #filter(class != "Epi-4.75" & class != "Epi-5.0_invivo" & class != "Epi-5.25" & class != "Epi-5.5_invitro") %>% 
  #filter(class != "ExE-5.5_invitro" & class != "ExE-4.75" & class != "Epi-5.25_invivo" )
  group_by(X1,cate) %>% summarise(mean = mean(fpkm)) %>% 
  spread(cate,mean) %>% mutate(ratio.rna = `TE-EXE` / `ICM-EPI`)


ratio.two <- df.ratio.rna %>% 
  inner_join(df.sig.pro,by= c("X1" = "genename")) %>% 
  filter(X1 != "Uchl1") %>% 
  mutate(rna = log2(ratio.rna),protein = log2(ratio)) %>% 
  mutate(type = if_else(protein >0 & rna > 0, "both.up",
                            if_else(protein <0 & rna <0,"both.down","reverse")))

#p1 <-

p <- ratio.two %>% group_by(type ) %>% summarise(n = n()) %>% 
  mutate(lab = paste0(type," (n=",n,")")) %>% 
  left_join(ratio.two)

ggscatter(p,x = "protein",y = "rna",add = "reg.line")+

  geom_point(aes(col = lab)) +
  geom_text_repel(aes(label = X1)) +
  stat_regline_equation(label.y.npc = 0.8) +
  stat_cor(label.y.npc = 0.9)+
  geom_vline(xintercept = 0,lty = 3) +
  geom_hline(yintercept = 0,lty = 3) +
  labs(x = "Protein level",y = "mRNA level",col = "") +
  theme_bw() +
  scale_color_npg() +
  theme(panel.grid = element_blank(),
        legend.background = element_blank(),
        legend.position = c(0.7,0.2))

ggsave("../result/3.two.stage/1.EPI_EXE.protein.vs.rna.png",width = 4,height = 4)
ggsave("../result/3.two.stage/1.EPI_EXE.protein.vs.rna.pdf",width = 4,height = 4)
# 
# 
# df.ratio.rna <- 
#   df.q.rna  %>% 
#   filter(# class != "Epi-4.75" & 
#          #  class != "Epi-5.0_invivo" & 
#          #    class != "Epi-5.25" & 
#              class != "Epi-5.5_invitro") %>%
#   filter(class != "ExE-5.5_invitro" #& 
#            #class != "ExE-4.75" & 
#            #class != "Epi-5.25_invivo"
#            ) %>% 
#   group_by(X1,cate) %>% summarise(mean = mean(fpkm)) %>% 
#   spread(cate,mean) %>% mutate(ratio.rna = `TE-EXE` / `ICM-EPI`)
# 
# 
# ratio.two <- df.ratio.rna %>% 
#   inner_join(df.sig.pro,by= c("X1" = "genename")) %>% 
#   filter(X1 != "Uchl1") %>% 
#   mutate(rna = log2(ratio.rna),protein = log2(ratio))
# 
# #p2 <- 
# ggscatter(ratio.two,x = "protein",y = "rna",add = "reg.line")+
#   geom_text_repel(aes(label = X1)) +
#   stat_regline_equation(label.y.npc = 0.8) +
#   stat_cor(label.y.npc = 0.9)+
#   geom_vline(xintercept = 0,lty = 3) +
#   geom_hline(yintercept = 0,lty = 3) +
#   theme_bw() +
#   
#   theme(panel.grid = element_blank())
# 
# 
# p1 +p2



# up
sig.p.up <- to.vol %>% filter(regulated == "TE-EXE") %>% arrange(genename) %>% 
  dplyr::select(ID,genename)

to.heat <- sig.p.up %>% left_join(df.scale)
heat <- to.heat %>% dplyr::select(-ID,-genename,-contains("Morula")) %>% as.data.frame()
name <- to.heat %>% left_join(sig.p.up)

rownames(heat) <- name$genename

p <- pheatmap::pheatmap(heat,cluster_cols = F,cluster_rows = F,show_colnames = F,
                        annotation_col = annot.sample,annotation_colors = annotcol)

pp <- ggplotify::as.ggplot(p) + labs(title = "Proteome")


sig.rna <- df.q.rna %>% 
  filter(X1 %in% sig.p.up$genename ) #%>% 
  #filter(class != "Epi-4.75" & class != "Epi-5.0_invivo" & class != "Epi-5.25" & class != "Epi-5.5_invitro") %>% 
  #filter(class != "ExE-5.5_invitro" & class != "ExE-4.75" & class != "Epi-5.25_invivo" )



to.heat <- sig.rna %>% dplyr::select(X1,sample,fpkm) %>% 
  spread(sample,fpkm) %>% 
  dplyr::select(X1,contains("ICM"),contains("Epi"),contains("TE"),contains("ExE"))

heat <- to.heat %>% dplyr::select(-X1) %>% as.data.frame()
row.names(heat) <- to.heat$X1


annot.
p <- pheatmap::pheatmap(heat,scale = "row",cluster_cols = F,cluster_rows = F)
ppt <- ggplotify::as.ggplot(p) + labs(title = "Transcriptome")

pp +ppt

ggsave("../result/3.two.stage/exe.pro_rna.png",width = 12,height = 6)












# down
sig.p.up <- to.vol %>% filter(regulated == "ICM-EPI")
to.heat <- df.scale %>% filter(ID %in% sig.p.up$ID) 

heat <- to.heat %>% dplyr::select(-ID,-contains("Morula")) %>% as.data.frame()
name <- to.heat %>% left_join(sig.p.up)
rownames(heat) <- name$genename

p <- pheatmap::pheatmap(heat,cluster_cols = F)
pp <- ggplotify::as.ggplot(p) + labs(title = "Proteome")


sig.rna <- df.q.rna %>% 
  filter(X1 %in% sig.p.up$genename ) %>% 
  filter(class != "Epi-5.0_invivo" & class != "Epi-5.25") %>% 
  filter(class != "ExE-5.5_invitro")



to.heat <- sig.rna %>% dplyr::select(X1,sample,fpkm) %>% 
  spread(sample,fpkm) %>% 
  dplyr::select(X1,contains("ICM"),contains("Epi"),contains("TE"),contains("ExE"))

heat <- to.heat %>% dplyr::select(-X1) %>% as.data.frame()
row.names(heat) <- to.heat$X1


p <- pheatmap::pheatmap(heat,scale = "row",cluster_cols = F)
ppt <- ggplotify::as.ggplot(p) + labs(title = "Transcriptome")

pp +ppt
ggsave("../result/3.two.stage/epi.pro_rna.png",width = 12,height = 6)






