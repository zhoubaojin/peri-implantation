library(tidyverse)
library(openxlsx)
library(RColorBrewer)
library(ComplexHeatmap)
library(ggsci)

#--------------------------
# target function 
# 1. TF
# 2. 
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







#----------------------------
# check tf version:
# 1. animal TF
# 2. tflink
#----------------------------

# tf format from tflink
# df.tf.link <- read_delim("D:/software/genome/function/TFLink_Mus_musculus/TFLink_Mus_musculus_interactions_All_simpleFormat_v1.0.tsv")
# df.tf.l <- df.tf.link %>% 
#   group_by(Name.TF) %>% summarise() %>% 
#   filter(Name.TF != "-")
# 
# 
# df.tf.mouse <- df.id %>% filter(!is.na(tf_category))
# 
# 
# 
# # tf format from animal tf
df.tf <- df.cluster %>% select(cluster,ID) %>%
  left_join(df.id,by = c("ID" = "protein")) %>%
#  filter(!is.na(tf_category)) %>%
  filter(tf_category != "Chromatin remodeling factors") %>% 
  select(cluster,tf_category,tf_family,eggnog_term,ID,genename,description)

# 
# df.tf.2 <- df.cluster %>% 
#   filter(genename %in% unique(df.tf.l$Name.TF)) %>% 
#   select(cluster,eggnog_term,ID,genename,description)
# 
# 
# ggVennDiagram::ggVennDiagram(list(animalTF = df.tf.mouse$genename,
#                                   tflink = df.tf.l$Name.TF))
# 
# ggVennDiagram::ggVennDiagram(list(animalTF = df.tf$genename,
#                                   tflink = df.tf.2$genename))



#---- to heatmap ----
df.type.cluster <- df.tf %>% 
  mutate(type = tf_category) %>% 
#  filter(!is.na(type)) %>% 
  left_join(df.sc) %>% 
  select(type,ID,genename,cluster,sample,int,sample) %>% 
  spread(sample,int) %>% 
  arrange(cluster,type)





heat <- df.type.cluster %>% select(sample$sample) %>% as.data.frame()
row.names(heat ) <- df.type.cluster$genename

limit <- min(max(heat),abs(min(heat)))
heat[heat > limit] <- limit
heat[heat < -limit] <- -limit

annot.sample <- sample %>% select(category,time,class) %>% as.data.frame()
row.names(annot.sample) <- sample$sample

annot.gene <- df.type.cluster %>% 
  mutate(cluster = str_replace(paste0("C",cluster+100),"C1","C")) %>% 
  select(cluster,type) %>% as.data.frame()
row.names(annot.gene) <- df.type.cluster$genename





col.type <- ggsci::pal_jco()(length(unique(annot.gene$type)))
names(col.type) <- unique(annot.gene$type)


col1 <- colorRampPalette(brewer.pal(8, "Spectral"))
col.cluster <- col1(num.c)
names(col.cluster) <- unique(sort(annot.gene$cluster))



# corr

to.cor <-  df.type.cluster  %>%
  select(genename,`Morula-1`:`E65EXE-3`) %>% 
  gather(sample,int,-genename) %>% spread(genename,int)

df.cor <- data.frame()
for(i in 2:(ncol(to.cor)-1)){
  for(j in (i+1):ncol(to.cor)){
    df.cor.i <- data.frame(from = colnames(to.cor)[i],
                           to = colnames(to.cor)[j],
                           cor = as.numeric(cor(to.cor[,i],to.cor[,j],method = "pearson")))
    df.cor <- df.cor %>% bind_rows(df.cor.i)
  }
}

add.index <- df.type.cluster %>% select(genename,cluster)
add.index$index <- 1:nrow(add.index)



col.reg <- ggsci::pal_npg(alpha = 0.2)(2)

cor.cut <- 0.7
cor.filter <- df.cor %>% filter(cor > cor.cut | cor < -cor.cut) %>% 
  mutate(color = if_else(cor > 0.5,col.reg[1],col.reg[2])) %>% 
  mutate(genename = from) %>% left_join(add.index) %>% mutate(from.i = index,from.c = cluster) %>% select(-genename,-index,-cluster) %>% 
  mutate(genename = to) %>% left_join(add.index) %>% mutate(to.i = index,to.c = cluster) %>% select(-genename,-index,-cluster) %>% 
  mutate(i =  from.c - to.c) %>% filter(i != 0)




library(circlize) 
split <- str_replace_all(paste0("C",df.type.cluster$cluster+100),"C1","C")
split <- factor(split,levels = unique(split))





do.circle <- T
a1 <- 1.5
a3 <- 8

b1 <- 1.2
b3 <- 6

c1 <- 1.0
c3 <- 3

do.circle <- T
if(do.circle) {
  pdf("../result/3.dna_rna_epi/cluster.tf.pdf",width = 7,height = 7)
  circos.clear()
  
  num.c <- 10
  circos.par(start.degree = 87, gap.after = c(rep(2,num.c-1),9))
  
  
  
  #----- track 1 -----
  heat.type <- df.type.cluster %>% select(type) %>% as.data.frame()
  row.names(heat.type) <- df.type.cluster$genename
  
  circos.heatmap(heat.type,  cluster = FALSE,
                 split = split,
                 col = col.type,
                 track.height = 0.02,
                 rownames.side = "outside", # gene name
                 rownames.cex = 0.7,
                 show.sector.labels = T#, # splite name
           #      bg.border = NA, bg.lwd = 1, bg.lty = 1
  )
  
  
  
  
  #----- track 2 : molular -------
  col_fun1 = colorRamp2(c(-limit, 0, limit), c("blue", "white", "red"))
  heat.1 <- heat %>% select(1:3)
  circos.heatmap(heat.1,split = split,col = col_fun1,
                 show.sector.labels = F, # splite name
                 track.height = 0.05#, 
             #    bg.border = NA ,bg.lwd = 1, bg.lty = 1
  )
  
  
  circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
    
    
    if(CELL_META$sector.numeric.index == num.c) { # the last sector
      
      # type
      col_type = col.samplecategory[1]
      circos.rect(CELL_META$cell.xlim[2] + convert_x(a1, "mm"), 0,
                  CELL_META$cell.xlim[2] + convert_x(a1 +(a3-a1)/3, "mm"), 3,
                  col = col_type,border = NA)
      
      
      col_class = col.sampleclass[1]  
      
      
      circos.rect(CELL_META$cell.xlim[2] + convert_x(a1 +(a3-a1)/3, "mm"), 0,
                  CELL_META$cell.xlim[2] + convert_x(a3 -(a3-a1)/3, "mm"), 3,
                  col = col_class,border = NA)
      
      col_time = col.time[1]  
      
      circos.rect(CELL_META$cell.xlim[2] + convert_x(a3 -(a3-a1)/3, "mm"), 0,
                  CELL_META$cell.xlim[2] + convert_x(a3, "mm"), 3,
                  col = col_time,border = NA)
      
      
    }
  },bg.border = NA)
  
  
  #----- track 3 : exe -------
  heat.2 <- heat %>% select(4:18)
  circos.heatmap(heat.2,split = split,col = col_fun1,
                 show.sector.labels = F, # splite name
                 track.height = 0.25
  )
  
  circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
    
    
    if(CELL_META$sector.numeric.index == num.c) { # the last sector
      
      # type
      col_type = col.samplecategory[2]
      circos.rect(CELL_META$cell.xlim[2] + convert_x(b1, "mm"), 0,
                  CELL_META$cell.xlim[2] + convert_x(b1+(b3-b1)/3, "mm"), 15,
                  col = col_type,border = NA)
      
      
      for(i in 1:5) {
        if(i == 1 ) { col_class = col.sampleclass[10]  }
        if(i == 2 ) { col_class = col.sampleclass[8]  }
        if(i == 3 ) { col_class = col.sampleclass[6]  }
        if(i == 4 ) { col_class = col.sampleclass[4]  }
        if(i == 5 ) { col_class = col.sampleclass[2]  }
        
        circos.rect(CELL_META$cell.xlim[2] + convert_x(b1+(b3-b1)/3, "mm"), (i-1)*3,
                    CELL_META$cell.xlim[2] + convert_x(b3-(b3-b1)/3, "mm"), i*3,
                    col = col_class,border = NA)
        
      }
      
      for(i in 1:5) {
        if(i == 1 ) { col_time = col.time[6]  }
        if(i == 2 ) { col_time = col.time[5]  }
        if(i == 3 ) { col_time = col.time[4]  }
        if(i == 4 ) { col_time = col.time[3]  }
        if(i == 5 ) { col_time = col.time[2]  }
        
        circos.rect(CELL_META$cell.xlim[2] + convert_x(b3-(b3-b1)/3, "mm"), (i-1)*3,
                    CELL_META$cell.xlim[2] + convert_x(b3, "mm"), i*3,
                    col = col_time,border = NA)
        
      }
    }
  },bg.border = NA)
  
  
  #------ track 4 : icm -------
  heat.3 <- heat %>% select(19:33)
  circos.heatmap(heat.3,split = split,col = col_fun1,
                 show.sector.labels = F, # splite name
                 track.height = 0.25, 
                 bg.border = NA#, 
        #         bg.lwd = 1, 
        #         bg.lty = 1
  )
  
  circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
    
    
    if(CELL_META$sector.numeric.index == num.c) { # the last sector
      
      # type
      col_type = col.samplecategory[3]
      circos.rect(CELL_META$cell.xlim[2] + convert_x(c1, "mm"), 0,
                  CELL_META$cell.xlim[2] + convert_x(c1+(c3-c1)/3, "mm"), 15,
                  col = col_type,border = NA)
      
      
      for(i in 1:5) {
        if(i == 1 ) { col_class = col.sampleclass[11]  }
        if(i == 2 ) { col_class = col.sampleclass[9]  }
        if(i == 3 ) { col_class = col.sampleclass[7]  }
        if(i == 4 ) { col_class = col.sampleclass[5]  }
        if(i == 5 ) { col_class = col.sampleclass[3]  }
        
        circos.rect(CELL_META$cell.xlim[2] + convert_x(c1+(c3-c1)/3, "mm"), (i-1)*3,
                    CELL_META$cell.xlim[2] + convert_x(c3-(c3-c1)/3, "mm"), i*3,
                    col = col_class,border = NA)
        
      }
      
      for(i in 1:5) {
        if(i == 1 ) { col_time = col.time[6]  }
        if(i == 2 ) { col_time = col.time[5]  }
        if(i == 3 ) { col_time = col.time[4]  }
        if(i == 4 ) { col_time = col.time[3]  }
        if(i == 5 ) { col_time = col.time[2]  }
        
        circos.rect(CELL_META$cell.xlim[2] + convert_x(c3-(c3-c1)/3, "mm"), (i-1)*3,
                    CELL_META$cell.xlim[2] + convert_x(c3, "mm"), i*3,
                    col = col_time,border = NA)
        
      }
    }
  },bg.border = NA)
  
  
  
  

  
  
  # correlation
  
  for (i in 1:nrow(cor.filter) ) {
    source <- cor.filter[i,5]
    target <- cor.filter[i,7]
    color <- cor.filter[i,4]
    
    
    circos.heatmap.link(source,target,col = color)
    
  }
  
  
  dev.off()
  
  
  
}

