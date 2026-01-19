library(tidyverse)
library(openxlsx)
library(pheatmap)
library(patchwork)
library(ggpubr)
library(ggrepel)


#----------------------------------------------------
#
# figure1. boxplot of raw and normalized intensity discritubtion
# figure2. pearson correlation of inner cell type 
# figure3. pearson correlation matrix of all cell type
# figure4. cross0analysis fo wang's study. chose one morula cell
# figure5. cross-analysis of cell's study. chose one morula cell
# figure6. cross-analysis of cell's study. chose tf of morula/blastocyst
# figure7. chose selected tf
#
#---------------------------------------------------
load("./mouse.uniprot.function.rda")
upper.panel <- function(x, y){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y,use = "na.or.complete",method = "pearson" ), 3)
  txt <- paste0("", format(r,nsmall = 3))
  cex.cor <- 0.5/strwidth(txt)
  hcolor = heat.colors(1000)
  rect(xleft = 0,ybottom = 0,xright = 1,ytop = 1,col=hcolor[(1-r)*1000])
  text(0.5, 0.5, txt, cex = cex.cor * r +0.5)
}

lower.panel<-function(x, y){
  df <- data.frame(x,y)
  c <- densCols(x,y, colramp=colorRampPalette(c("black", "white")))
  df$dens <- col2rgb(c)[1,] + 1L
  cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F",
                              "#FCFF00", "#FF9400", "#FF3100"))(256)
  df$col <- cols[df$dens]
  points(y~x,data=df[order(df$dens),],pch=20,col=col,cex=0.1)
  abline(lm (y~x),lty=3,col="red")
}


sample <- read.xlsx("Supplementary Table 1.xlsx",sheet = "FileInformation") %>% arrange(time)
sample$class <- factor(sample$class,levels = unique(sample$class),ordered = T)

library(RColorBrewer)
col.samplecategory <- brewer.pal(3, "Set2")
names(col.samplecategory) <- unique(sample$category)

col.time <- brewer.pal(length(unique(sample$time)), "YlOrRd")
names(col.time) <- unique(sample$time)

col.sampleclass <- c("grey",brewer.pal(10, "Paired"))
names(col.sampleclass) <- c(unique(sample$class))
annotcol <- list(class = col.sampleclass,
                 category = col.samplecategory,
                 time = col.time)


#----------------------------------------------------
#
#
# chapter1  normalization
#
#
#---------------------------------------------------


file <- read.xlsx("Supplementary Table 1.xlsx",sheet = "RawIntensity") 
df <- file %>% gather(sample,intensity,-ID) %>% left_join(sample) %>% arrange(class) %>% 
  filter(intensity >0) %>% mutate(log = log10(intensity)) %>% select(-intensity) 

df$sample <- factor(df$sample,levels = unique(df$sample),ordered = T)


# raw intensity distribution
mean <-  df  %>% group_by(class,ID) %>% summarise(int = mean(na.omit(log)))

library(ggridges)
ggplot(mean,aes(int,class,fill= class)) +
  geom_density_ridges() +
  scale_fill_manual(
    values = col.sampleclass
  ) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.y = element_text(hjust=0))
#ggsave("../result/1.global/0.raw.png",width = 4,height = 4)

  
to.nor <- df %>% mutate(type = "Raw intensity")


to.cor <-  df %>%select(ID,sample,log) %>%  spread(sample,log) %>% select(-ID)
cor.raw <- cor(to.cor,use = "na.or.complete") %>% as.data.frame()
df.cor <- cor.raw %>% as.data.frame() %>%  bind_cols(source = row.names(cor.raw)) %>% 
  gather(target,cor,-source) %>% mutate(is = if_else(source == target,1,0)) %>% filter(is == 0) %>% 
  left_join(sample,by = c("target" = "sample")) %>%  mutate(target.class = class) %>% 
  select(source,target,target.class,cor) %>% 
  left_join(sample,by = c("source" = "sample")) %>%  mutate(is = if_else(target.class == class,1,0)) %>% filter(is == 1) 

df.cor$i <- 1:nrow(df.cor)
df.cor.top <- df.cor %>% group_by(i) %>% summarise(relation = paste(sort(c(target,source)),collapse = "-")) %>% ungroup() %>%
  left_join(df.cor) %>% group_by(relation) %>% top_n(n= 1,wt  = i) %>% mutate(type = "Raw intensity")
to.cor.type <- df.cor.top



cor.raw[cor.raw < 0.3] <- 0.3


annot.class <- sample %>% select(class) %>% as.data.frame()
row.names(annot.class) <- sample$sample

p1 <- pheatmap(cor.raw,cluster_rows = F,cluster_cols = F,annotation_names_row = F,annotation_names_col = F,
               annotation_row = annot.class,annotation_col = annot.class,
               show_rownames = F,show_colnames = F,
               annotation_colors = annotcol)




file <- read.xlsx("Supplementary Table 1.xlsx",sheet = "NormalizedIntensity")
df <- file %>% gather(sample,intensity,-ID) %>% left_join(sample) %>% arrange(class) %>% 
  filter(intensity >0) %>% mutate(log = log10(intensity)) %>% select(-intensity)  %>% 
  mutate(type = "Normalized intensity")
df$sample <- factor(df$sample,levels = unique(df$sample),ordered = T)



# plot1 boxplot distribution
to.nor <- to.nor %>% bind_rows(df) %>% 
  arrange(class)
to.nor$sample <- factor(to.nor$sample,levels = unique(to.nor$sample),ordered = T)
to.nor$type <- factor(to.nor$type,levels = unique(to.nor$type),ordered = T)
ggplot(to.nor) +
  geom_boxplot(aes(sample,log,col = class),outliers = F) +
  facet_grid(type ~.) +
  scale_color_manual(values =col.sampleclass ) +
  labs(y = expression(log[10]*'(Intensity)')) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())

ggsave("../result/1.global/QC/normalized.intensity.png",width = 4,height = 3)





to.cor <- df %>% select(sample,log,ID )%>%   spread(sample,log) %>% select(-ID)
cor.raw <- cor(to.cor,use = "na.or.complete") %>% as.data.frame()

p2 <- pheatmap(cor.raw,cluster_rows = F,cluster_cols = F,annotation_names_row = F,annotation_names_col = F,
               annotation_row = annot.class,annotation_col = annot.class,
               show_rownames = F,show_colnames = F,
               annotation_colors = annotcol)

# plot2 cormatrix
pp1 <- ggplotify::as.ggplot(p1)
pp2 <- ggplotify::as.ggplot(p2)
pp1 +pp2
ggsave("../result/1.global/QC/correlation.matrix.png",width = 14,height = 5.5)



df.cor <- cor.raw %>% as.data.frame() %>%  bind_cols(source = row.names(cor.raw)) %>% 
  gather(target,cor,-source) %>% mutate(is = if_else(source == target,1,0)) %>% filter(is == 0) %>% 
  left_join(sample,by = c("target" = "sample")) %>%  mutate(target.class = class) %>% 
  select(source,target,target.class,cor) %>% 
  left_join(sample,by = c("source" = "sample")) %>%  mutate(is = if_else(target.class == class,1,0)) %>% filter(is == 1) 

df.cor$i <- 1:nrow(df.cor)
df.cor.top <- df.cor %>% group_by(i) %>% summarise(relation = paste(sort(c(target,source)),collapse = "-")) %>% ungroup() %>%
  left_join(df.cor) %>% group_by(relation) %>% top_n(n= 1,wt  = i) %>% mutate(type = "Normalized intensity")






# plot3 cor inner cell type
to.cor.type <- to.cor.type %>% bind_rows(df.cor.top)
to.cor.type$type <- factor(to.cor.type$type,levels = unique(to.cor.type$type),ordered = T)

mean.r <- to.cor.type %>% group_by(type) %>% summarise(mean = mean(cor))

ggbarplot(to.cor.type,x = "class",y = "cor",add = "mean_se",fill = "class",facet.by = "type" ) +
  geom_jitter(aes(class,cor,fill = class),shape=21,size = 1) +
  geom_hline(aes(yintercept = mean ),data = mean.r,lty = 3,col = "tomato")+
  geom_text(data= mean.r,aes(label = paste("mean r:",round(mean,2)),x = "Morula", y = 1),hjust =0,vjust = -0.1,col = "black") +
  facet_grid(type ~.) +
  scale_y_continuous(limits = c(0,1.1),n.breaks = 5) +
  scale_fill_manual(values = col.sampleclass) +
  labs(y = "Pearson correlation ") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())

ggsave("../result/1.global/QC/correlation.r.png",width = 4,height = 3)











##--------------------------------------------
##
## target
##
##-------------------------------------------
df.target <- data.frame(genename = c("Hmga1","Hmga2","Hmgb1","Hmgb2","Lin28a","Sox15","Stat3","Carhsp1"))



df.cell  <- read.xlsx("../raw/publish/mmc1.xlsx",sheet = "Mouse_log2_normalized_withQC") %>% 
  gather(sample,log,-Gene_Name) %>% 
  mutate(int = 2 ^ log) %>% 
  mutate(genename = Gene_Name) %>% 
  filter(genename %in% df.target$genename) %>% 
  mutate(class =str_replace_all(sample,"_Rep\\d+$","")) %>% 
  mutate(data = "Zhu et al.")


to.scale <- df.cell %>% 
  filter(class == "Mor" |class == "Blast") %>% 
  select(genename,sample,int) %>%spread(genename,int) 
df.cell.scale <- to.scale %>% select(-sample) %>% scale() %>% bind_cols(sample = to.scale$sample) %>% 
  gather(genename,scale,-sample) %>% mutate(class =str_replace(sample,"_\\w+",""))  %>% 
  mutate(data = "Zhu et al.")


df.cell.per <- df.cell %>% group_by(genename) %>% summarise(max = max(int)) %>% 
  left_join(df.cell) %>% mutate(per = int/max)


df <- read.xlsx("Supplementary Table 1.xlsx",sheet = "NormalizedIntensity") %>% 
  gather(sample,int ,-ID) %>%
  left_join(sample) %>% 
  filter(class == "Morula" | 
           class == "3.5d-ICM" | 
           class == "3.5d-TE" | 
           class == "4.5d-ICM" | 
           class == "4.5d-TE" 
           ) %>% 
  left_join(fun$protein,by = c("ID" = "protein")) %>% 
 # mutate(class = str_replace(class,"3.5d-ICM","ICM")) %>% 
#  mutate(class = str_replace(class,"3.5d-TE","TE")) %>% 
  filter(genename %in% df.target$genename) %>% 
  mutate(data = "This study")


to.scale <- df %>% 
  filter(class == "Morula" | class == "3.5d-ICM" |class == "3.5d-TE") %>% 
  select(genename,sample,int) %>%spread(genename,int) 

df.scale <- to.scale %>% select(-sample) %>% scale() %>% bind_cols(sample = to.scale$sample) %>% 
  gather(genename,scale,-sample) %>% left_join(sample) %>% mutate(data = "This study")

df.per <- df %>% group_by(genename) %>% summarise(max = max(int)) %>% 
  left_join(df) %>% mutate(per = int/max)
  
  
df.bind <- df.cell.per %>% bind_rows(df.per)  
#to.plot <- data.frame(class = c("Mor","Blast","Morula","ICM","TE") ) %>% left_join(df.bind)

select <-  c("Mor","Blast","Morula","3.5d-ICM","3.5d-TE")




tp <- df.scale %>% bind_rows(df.cell.scale) 


p <- data.frame(class = select) %>% left_join(tp)
p <- df.target %>% left_join(p)
p$genename <- factor(p$genename,levels = unique(p$genename),ordered = T)
p$class <- factor(p$class,levels = unique(p$class),ordered = T)


ggplot(p,aes(class,scale,fill = data)) +
  geom_boxplot() +
  geom_jitter(shape =21) +
  facet_wrap( . ~ genename,scales = "free",nrow = 2) +
  labs(y = "Z-score",x  = "",fill = "") +
  scale_y_continuous(expand = expansion(mult = c(0.1,0.1))) +
  # ggpubr::stat_compare_means(
  #   method = "t.test",label = "p.signif",hide.ns = T,
  #   comparisons = list(c("Mor","Blast"),
  #                      c("Morula","3.5d-ICM"),
  #                      c("Morula","3.5d-TE")))+  
  ylim(-1.5,2) +
  ggsci::scale_fill_jco() +
  theme_bw() +
  theme(legend.position = "right",
        panel.grid = element_blank(),
        strip.background.x = element_blank(),
        axis.text.x = element_text(angle = 30,hjust = 1))

ggsave("../result/1.global/QC/comparision.cell.morula.zscore.pdf",width = 8,height = 4)
ggsave("../result/1.global/QC/comparision.cell.morula.zscore.png",width = 8,height = 4)
