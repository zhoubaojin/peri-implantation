library(openxlsx)
library(tidyverse)
library(Mfuzz)
library(RColorBrewer)
library(ComplexHeatmap)
library(patchwork)
library(ggpubr)
library(ggrepel)


#------------------
#
# 1. loading data and prepare for mfuzzy cluster
#
#------------------


load("./mouse.uniprot.function.rda")

df.id <- fun$protein %>% mutate(ID = protein)

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







#------------------
#
# 2. cluster test
#
#------------------



library(Mfuzz)
to.mfuzz <- df.scale.mean %>% spread(class,scale) %>%  select(ID,sample$class)
write.table(x = to.mfuzz,file= "../result/processing//to.fuzzy.xls",row.names = F,sep = "\t",quote = F)

cluster.input <- table2eset("../result/processing/to.fuzzy.xls")


for (n in 5:25){
  para.cluster = n
  para.m <- mestimate(cluster.input)
  para.membership <- 0.5
  
  cluster.output <- mfuzz(eset = cluster.input,
                          c=para.cluster,
                          m=para.m)  #m = fuzzification parameter
  
  
  save(cluster.output,file = paste0("../result/processing/cluster.",n,".morula"))
  
  load( paste0("../result/processing/cluster.",n,".morula"))
  
  
  
  # if you want to filter some protein with low membership
  
  df.cluster.membership <- cluster.output$membership %>% as.data.frame() %>% 
    bind_cols(data.frame(ID = row.names(cluster.output$membership))) %>% 
    gather(class,m,-ID)  %>% 
    group_by(ID) %>% summarise(membership = max(m)) %>% 
    filter(membership > 0.5)
  
  
  df.cluster.zcore <- data.frame(ID =  names(cluster.output$cluster),
                                 cluster =  as.character(cluster.output$cluster )) %>% 
    inner_join(df.cluster.membership,by = "ID")  %>%   # ilter low membership feature
    left_join(df.scale.mean,by = "ID") %>% 
    arrange(membership)
  df.cluster.zcore$ID <- factor(df.cluster.zcore$ID,levels = unique(df.cluster.zcore$ID),ordered = T)
  
  
  
  num.cluster <- df.cluster.zcore %>% group_by(cluster,ID) %>% summarise() %>% ungroup()%>% 
    group_by(cluster) %>% summarise(n=n()) %>% mutate(cluster.num = paste0("C",cluster," (n=",n,")")) %>% 
    select(cluster.num,cluster)
  
  
  
  
  # prepare two data 
  
  df.cluster.center <- cluster.output$centers %>% as.data.frame() %>% 
    bind_cols(data.frame(cluster = row.names(cluster.output$centers))) %>% 
    gather(class,center,-cluster)  %>% 
    left_join(num.cluster)
  
  
  
  to.plot <- df.cluster.zcore %>% left_join(num.cluster)
  p1 <- ggplot(to.plot) +
    geom_line(aes(class,scale,group = ID,col = membership)) +
    facet_grid(cluster.num ~.,) +
    geom_line(data= df.cluster.center,aes(class,center,group = "1"),size =1,lty = 2) +
    scale_color_gradientn(colours = c("red","yellow","green","blue","purple"),
                          values = c(1.0,0.8,0.6,0.4,0.2,0)) +
    labs(y = "Z-score",color = "Membership") +
    theme_bw() +
    theme(
      legend.position = "top",
      strip.text.y = element_text(angle = 0,hjust = 0),
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 45,hjust = 1),
      panel.grid = element_blank())
  
  
  
  # heat
  to.heat <- data.frame(ID =  names(cluster.output$cluster),
                        cluster =  as.character(cluster.output$cluster )) %>% 
    inner_join(df.cluster.membership,by = "ID")  %>% 
    left_join(df.sc,by = "ID") 
  
  
  
  p2 <- 
    ggplot(to.heat,aes(sample,ID)) +
    geom_tile(aes(fill = int)) +
    facet_grid(cluster ~ .,scales = "free") +
    scale_fill_gradientn(colours = c("red4","red","white","blue","blue4"),
                         values = c(1.0,0.7,0.6,0.3,0.1,0))  +
    labs(fill = "Z-score")+
    theme(legend.position = "top",
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          strip.text.y = element_blank())
  
  p1 + p2
  ggsave(paste0("../result/processing//cluster.num.",para.cluster,".png"),width = 8,height = 8)
  
}


#------------------
#
# 3. chose cluster
#
#------------------




num.c <- 10
load( paste0("../result/2.cluster/backup/25.to.18/processing/cluster.",num.c,".morula"))



# if you want to filter some protein with low membership

df.cluster.membership <- cluster.output$membership %>% as.data.frame() %>% 
  bind_cols(data.frame(ID = row.names(cluster.output$membership))) %>% 
  gather(class,m,-ID)  %>% 
  group_by(ID) %>% summarise(membership = max(m)) #%>% 
  #filter(membership > 0.5)


df.cluster.zcore <- data.frame(ID =  names(cluster.output$cluster),
                               cluster =  as.character(cluster.output$cluster )) %>% 
  inner_join(df.cluster.membership,by = "ID")  %>%   # ilter low membership feature
  left_join(df.scale.mean,by = "ID") %>% 
  arrange(membership)
df.cluster.zcore$ID <- factor(df.cluster.zcore$ID,levels = unique(df.cluster.zcore$ID),ordered = T)


# for changing



new.c <- data.frame(cluster = as.character(1:10),
                    newcluster = c("08","03","02","07","04","01","05","06","09","10"))
num.cluster <- df.cluster.zcore %>% group_by(cluster,ID) %>% summarise() %>% ungroup()%>% 
 left_join(new.c,by = "cluster") %>% 
  group_by(newcluster,cluster) %>% summarise(n=n()) %>% mutate(cluster.num = paste0("C",newcluster," (n=",n,")")) %>% 
  select(cluster.num,cluster,newcluster)




# prepare two data 



to.plot <- df.cluster.zcore %>% left_join(num.cluster) %>% mutate(newcluster = as.numeric(newcluster)) %>% arrange(newcluster)

to.plot$cluster.num <- factor(to.plot$cluster.num,levels = unique(to.plot$cluster.num),ordered = T)

df.cluster.center <- cluster.output$centers %>% as.data.frame() %>% 
  bind_cols(data.frame(cluster = row.names(cluster.output$centers))) %>% 
  gather(class,center,-cluster)  %>% 
  left_join(num.cluster) %>% 
  mutate(newcluster = as.numeric(newcluster)) %>% arrange(newcluster)

df.cluster.center$cluster.num <- factor(df.cluster.center$cluster.num,levels = unique(df.cluster.center$cluster.num),ordered = T)



#p1 <-
ggplot(to.plot,aes(class,scale)) +
  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf,
           alpha = .5,fill = col.samplecategory[1])+
  annotate("rect", xmin = 1.5, xmax = 6.5, ymin = -Inf, ymax = Inf,
           alpha = .5,fill = col.samplecategory[2])+
  annotate("rect", xmin = 6.5, xmax = 11.5, ymin = -Inf, ymax = Inf,
           alpha = .5,fill = col.samplecategory[3])+

  geom_line(aes(group = ID,col = membership),alpha =.2) +
  facet_wrap(cluster.num ~.,ncol = 5) +
  geom_line(data= df.cluster.center,aes(class,center,group = "1"),linewidth = 1,lty = 3) +
  scale_color_gradientn(colours = c("red","yellow","green"),
                        values = c(1.0,0.75,0.5,0.25,0)) +
  labs(y = "Z-score",color = "Membership") +
  theme_bw() +
  theme(
    legend.position = "top",
    strip.text.y = element_text(angle = 0,hjust = 0),
    strip.background = element_blank(),
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank())

ggsave("../result/2.cluster/10.cluster/cluster.line.png",width = 9,height = 4)




# cytoscape
for(i in 1:length(unique(to.plot$cluster.num))){
  to.plot.i <- to.plot %>% filter(cluster.num %in% unique(to.plot$cluster.num)[i] )
  df.center.i <- df.cluster.center %>% filter(cluster.num %in% unique(to.plot$cluster.num)[i] )
  ggplot(to.plot.i,aes(class,scale)) +
    annotate("rect", xmin = 0, xmax = 1.5, ymin = -Inf, ymax = Inf,
             alpha = .5,fill = col.samplecategory[1])+
    annotate("rect", xmin = 1.5, xmax = 6.5, ymin = -Inf, ymax = Inf,
             alpha = .5,fill = col.samplecategory[2])+
    annotate("rect", xmin = 6.5, xmax = 12, ymin = -Inf, ymax = Inf,
             alpha = .5,fill = col.samplecategory[3])+
    
    geom_line(aes(group = ID,col = membership),alpha =.2) +
    facet_wrap(cluster.num ~.,ncol = 5) +
    geom_line(data= df.center.i,aes(class,center,group = "1"),linewidth = 1,lty = 3) +
    scale_color_gradientn(colours = c("red","yellow","green"),
                          values = c(1.0,0.75,0.5,0.25,0)) +
    scale_y_continuous(breaks =  c(3,0,-3)) +
    labs(y = "Z-score",color = "Membership") +
    theme_void() +
    theme(
      legend.position = "none",
      strip.text.y = element_text(angle = 0,hjust = 0),
      strip.background = element_rect(fill = col.cluster[i],color = NA),
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank())
  ggsave(paste0("../result/2.cluster/10.cluster/function/network/c",i,".png"),width = 1.5,height = 1)  
}






to.heat <- to.plot %>% 
  group_by(newcluster,ID) %>% summarise() %>% ungroup() %>% 
  inner_join(df.cluster.membership,by = "ID")  %>% 
  left_join(df.sc,by = "ID") %>% 
  select(newcluster,ID,membership,sample,int) %>% spread(sample,int) 

out.cluster <- to.heat %>% left_join(df.id,by = "ID") %>% mutate(cluster = newcluster) %>% 
  select(cluster,ID,genename,description,everything()) %>% select(-newcluster)












#--------------------------------
#
# 4. cluster function 
#
#--------------------------------
out.cluster <- read.xlsx("D:/project/proteome_embr/manuscript/Table2.Cluster.version202503.xlsx")


to.fun <-  out.cluster %>% 
  mutate(term = go_biologicalprocess) %>%
  separate_rows(term,sep =";") %>%
  filter(!is.na(term)) %>%
  group_by(cluster,term) %>% summarise(num.target=n(),gene = paste0(genename,collapse = ";"))  %>%
  separate(term,into = c("termid","term"),sep = "::") %>% select(-termid) %>%
  left_join(fun$term$go,by = "term" ) %>% mutate(num.term = bgcount)


num.all <- fun$bgnum %>% filter(database == "go") %>% select(bg) %>% as.numeric()
cluster.num <- out.cluster %>% 
  mutate(term = go_biologicalprocess) %>%
  separate_rows(term,sep =";") %>%
  filter(!is.na(term)) %>%
  group_by(cluster,ID) %>% summarise() %>% 
  group_by(cluster) %>% summarise(num.cluster = n())

fun.pvalue.go <- to.fun %>% 
  left_join(cluster.num) %>% 
  group_by(cluster,termid,term,gene,num.target,num.term,num.cluster) %>% 
  summarise(pvalue = phyper(num.target-1,num.term,num.all - num.term,num.cluster,lower.tail = F,log.p = F )) %>% 
  mutate(qvalue = p.adjust(pvalue,method = "BH",n = num.cluster)) %>% ungroup()


wb <- createWorkbook()
addWorksheet(wb,"ClusterProtein")
writeData(wb,"ClusterProtein",out.cluster)
addWorksheet(wb,"ClusterGO")
writeData(wb,"ClusterGO",fun.pvalue.go)



to.fun <-  out.cluster %>% 
  mutate(term = kegg_pathway) %>%
  separate_rows(term,sep =";") %>%
  filter(!is.na(term)) %>%
  group_by(cluster,term) %>% summarise(num.target=n(),gene = paste0(genename,collapse = ";"))  %>%
  separate(term,into = c("l1","l2","l3"),sep = "\\|\\|") %>% 
  separate(l3,into = c("termid","term"),sep = "::") %>% select(-termid) %>%
  left_join(fun$term$kegg_pathway,by = "term" ) %>% mutate(num.term = bgcount)


num.all <- fun$bgnum %>% filter(database == "kegg") %>% select(bg) %>% as.numeric()
cluster.num <- out.cluster %>% 
  mutate(term = kegg_pathway) %>%
  separate_rows(term,sep =";") %>%
  filter(!is.na(term)) %>%
  group_by(cluster,ID) %>% summarise() %>% 
  group_by(cluster) %>% summarise(num.cluster = n())

fun.pvalue.kegg <- to.fun %>% 
  left_join(cluster.num) %>% 
  group_by(cluster,level1,level2,termid,term,gene,num.target,num.term,num.cluster) %>% 
  summarise(pvalue = phyper(num.target-1,num.term,num.all - num.term,num.cluster,lower.tail = F,log.p = F )) %>% 
  mutate(qvalue = p.adjust(pvalue,method = "BH",n = num.cluster)) %>% 
  ungroup()

addWorksheet(wb,"ClusterKEGG")
writeData(wb,"ClusterKEGG",fun.pvalue.kegg)

saveWorkbook(wb,file = "../result/2.cluster/10.cluster/Table2.Cluster.xlsx",overwrite = T)













