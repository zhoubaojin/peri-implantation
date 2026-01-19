library(tidyverse)
library(ggrepel)
library(RColorBrewer)
library(patchwork)
library(pheatmap)
library(ropls)
library(pROC)
library(ggpubr)
library(ggcorrplot)
library(impute)
library(readxl)
library(VennDiagram)
library(ggVennDiagram)
library(networkD3)
library(htmlwidgets)
library(visNetwork)
library(circlize)
library(igraph)
library(ggupset)
library(GGally)
library(Mfuzz)
library(plotly)
library(openxlsx)
library(ggpubr)
cols <- colorRampPalette(brewer.pal(8, "Spectral"))


load("../bin/mouse.uniprot.function.rda")
df.id <- fun$protein


##---------------------------------------------------
# preparation
## input:
## 1 protein profiling
## 2 peptide profiling
#----------------------------------------------------


## global setting 1, protein
pro <-read_excel("./proteinGroups.txt") %>% 
  filter(is.na(`Potential contaminant`)) %>% 
  filter(is.na(Reverse)) %>% 
  mutate(filter.reverse = str_match(`Protein IDs`,"^REV__")) %>% filter(is.na(filter.reverse)) %>%
  mutate(filter.cont = str_match(`Protein IDs`,"^CON__")) %>% filter(is.na(filter.cont)) %>%
  filter(`Q-value` < 0.01)

sample <- read_excel("sample.xlsx")









##################################################################################################
##
## part 1 identification stat
##
##################################################################################################


df.q <- pro %>% mutate(ID = str_replace_all(`Protein IDs`,";.+","")) %>% 
  mutate(ID = str_replace_all(str_replace_all(ID,"^\\w+\\|",""),"\\|\\w+$","")) %>%
  select(ID,contains("Intensity ")) 


df.id.proteingroup <- pro %>%  mutate(ID = str_replace_all(`Protein IDs`,";.+","")) %>% 
  mutate(ID = str_replace_all(str_replace_all(ID,"^\\w+\\|",""),"\\|\\w+$","")) %>%
  mutate(ProteinGroups = `Protein IDs`) %>%
  select(ID,ProteinGroups,`Q-value`,Score,`Sequence coverage [%]`)
  
  
df.q.t <- df.q %>% gather(sample,int,-ID) %>% 
  mutate(file = str_replace(file, "Intensity ",""),int =as.numeric(int)) %>%
  filter(int >0) %>% mutate(log = log10(int)) %>% 
  inner_join(sample)









# out 
out.raw <- sample %>% select(file) %>% left_join(df.q.t,by = "file") 
out.raw$sample <-factor(out.raw$sample,levels = unique(out.raw$sample),ordered = T)
out.raw <- out.raw %>% select(sample,int,ID) %>% spread(sample,int)



# count protein number
df.count <- df.q.t %>%
  group_by(category,time,class,sample) %>% summarise(n =n()) %>% arrange(category,time)
df.count$sample <- factor(df.count$sample,levels = unique(df.count$sample),ordered = T)
df.count$category <- factor(df.count$category,levels = unique(df.count$category),ordered = T)


pid <- ggplot(df.count,aes(class,n,fill= class)) +

  geom_boxplot() +
  geom_point(size = 3,aes(fill= class),shape = 21) +
  scale_y_continuous(labels = scales::comma,limits = c(0,1800)) +
  scale_fill_manual(values = col.sampleclass) +
  scale_color_manual(values = col.sampleclass) +
  theme_classic() +
  labs(y = "No. protein groups") +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank() )




# count missing value

count.pro <- df.q.t %>% group_by(ID) %>% summarise(sn = n()) %>% mutate(per = sn / length(unique(df.q.t$file)))
count.pro.q.1 <- count.pro %>% filter(sn >2)
count.pro.q.25 <- count.pro %>% filter(per >=0.25)
count.pro.q.50 <- count.pro %>% filter(per >=0.50)
count.pro.q.100 <- count.pro %>% filter(per ==1)



# count missing in group
count.pro.group <- df.q.t %>% left_join(sample) %>% 
  group_by(class,ID) %>% summarise(n=n()) %>% ungroup() %>% 
  #filter(n > 1) %>% group_by(ID) %>% summarise()
  group_by(ID) %>% summarise(sn = max(n)) %>% 
  filter(sn >=2)

pmiss <- 
  ggplot(count.pro) +
  geom_histogram(aes(per),bins = 100,col = "white",fill = "steelblue4") +
  theme_classic() +
  scale_x_continuous(labels = scales::percent) +
  scale_y_continuous(labels = scales::comma,limits = c(0,nrow(count.pro)*1.1)) +
  geom_hline(yintercept = nrow(count.pro),lty = 3) +
  annotate(geom = "text",size = 3,x= 0  , y=nrow(count.pro),hjust = 0,vjust = -0.2,label = paste0(nrow(count.pro), " identified proteins")) +
    
  geom_hline(yintercept = nrow(count.pro.q.1),lty = 3) +
  annotate(geom = "text",size = 3,x= 0  , y=nrow(count.pro.q.1),hjust = 0,vjust = -0.2,label = paste0(nrow(count.pro.q.1), ">=3 sample")) +
    
    geom_hline(yintercept = nrow(count.pro.group),lty = 3) +
    annotate(geom = "text",size = 3,x= 0.2  , y=nrow(count.pro.group),hjust = 0,vjust = -0.2,label = paste0(nrow(count.pro.group), ">=2 sample in each phenotype")) +
    
  annotate(geom = "rect",xmin = 0.25,xmax = 1,ymin = 0,ymax = nrow(count.pro.q.25),fill = "green4",alpha =0.2 )+
  annotate(geom = "text",size = 3,x= 0.25  , y=nrow(count.pro.q.25),hjust = 0,vjust = 0,label = paste0(nrow(count.pro.q.25), " proteins in >= 25% samples"),col = "green4") +
  annotate(geom = "text",size = 3,x= 0.5, y=nrow(count.pro.q.50),hjust = 0,vjust = 1,label = paste0(nrow(count.pro.q.50), " proteins in >= 50% samples"),col = "red4") +
  annotate(geom = "rect",xmin = 0.5,xmax = 1,ymin = 0,ymax = nrow(count.pro.q.50),fill = "red",alpha =0.2 )+
  geom_hline(yintercept = nrow(count.pro.q.100),lty = 3,col="black") +
  annotate(geom = "text",size = 3,x= 1, y=nrow(count.pro.q.100),hjust = 1,vjust = 1,label = paste0(nrow(count.pro.q.100), " proteins in all samples"),col = "black") +
  
  labs(x = "Proportion of quantified proteins in all samples",y = "Protein count")



##################################################################################################
##
## part 2 quantification  preprocessing
##
##################################################################################################





# impute 3


count.pro.q <- count.pro.group
df.demissing  <- count.pro.q %>% select(ID) %>% left_join(df.q.t,by = "ID")  %>% 
  inner_join(sample) %>%  arrange(sample)

df.min.group <- df.demissing %>% group_by(ID,class) %>% summarise(min = min(int) *0.5)
df.min.sample <- df.demissing %>% group_by(ID) %>% summarise(min = min(int) *0.1)

df.demissing.pre <- df.demissing %>% select(ID,sample,int) %>% spread(ID,int)
toimput <- sample %>% left_join(df.demissing.pre)

for( k in 6:ncol(toimput)){
  for (l in 1:nrow(toimput)){
    i.sample <- as.character(toimput[l,5])
    i.class <- sample %>% filter(sample %in% i.sample)  %>% select(class) 
    i.pro    <- colnames(toimput[,k])
    
    if (is.na(toimput[l,k])){
      message(paste0("protein number: ",k-5, "; protein: ",i.pro))
      
      i.min.class <- df.min.group %>% filter(ID %in% i.pro) %>% filter(class %in% i.class$class)
      
      if(nrow(i.min.class) == 0) {                            
      #  i.impute <- df.min.sample %>% filter(ID %in% i.pro) 
      #  toimput[l,k] <- i.impute$min
      } else {
    
        toimput[l,k] <- i.min.class$min
      }
      
      
    }
  }
}



to.nor <- toimput %>% gather(ID,int,6:ncol(toimput)) %>% select(ID,sample,int)


p1 <- to.nor %>% group_by(sample) %>% summarise(sum = sum(na.omit(int))) %>% 
  ggplot()+geom_bar(aes(sample,sum),stat= "identity") +
  theme(axis.title = element_blank(),
        axis.text.x = element_blank())





# choose norm method
df.nor.r <- to.nor %>% group_by(sample) %>% summarise(sn = sum(na.omit(int))) 
df.nor.pf  <- to.nor %>% left_join(df.nor.r,by = "sample") %>%
  mutate(nor = int/sn *mean(df.nor.r$sn),log= log10(nor)) 

p2 <- df.nor.pf %>% group_by(sample) %>% summarise(sum = sum(na.omit(nor))) %>% 
  ggplot()+geom_bar(aes(sample,sum),stat= "identity")+
  theme(axis.title = element_blank(),
        axis.text.x = element_blank())





# inpute2
toimput.2 <- sample %>% left_join(df.nor.pf)%>% select(sample,class,ID,nor) %>% spread(ID,nor)
for( i in 3:ncol(toimput.2)){
  toimput.2[,i][is.na(toimput.2[,i])] <- min(na.omit(toimput.2[,i])) * 0.2
}


to.nor.2 <- toimput.2 %>% gather(ID,int,3:ncol(toimput.2)) %>% select(ID,sample,int)
p3 <- to.nor.2 %>% group_by(sample) %>% summarise(sum = sum(na.omit(int))) %>% 
  ggplot()+geom_bar(aes(sample,sum),stat= "identity")+
  theme(axis.title = element_blank(),
        axis.text.x = element_blank())



df.nor.r.2 <- to.nor.2 %>% group_by(sample) %>% summarise(sn = sum(na.omit(int))) 
df.nor.final  <- to.nor.2 %>% left_join(df.nor.r.2,by = "sample") %>%
  mutate(nor = int/sn *mean(df.nor.r.2$sn),log= log10(nor)) 

p4 <- df.nor.final %>% group_by(sample) %>% summarise(sum = sum(na.omit(nor))) %>%
  ggplot()+geom_bar(aes(sample,sum),stat= "identity")+
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1))


p1 /p2 /p3 /p4
ggsave("../result/0.norm.png",width = 8,height = 8)




out.nor <- sample %>% left_join(df.nor.final)
out.nor$sample <- factor(out.nor$sample,levels = unique(out.nor$sample),ordered = T)
out.nor <- out.nor %>%  select(sample,ID,nor) %>% spread(sample,nor) 


out.scale <- toheat %>% select(ID,sample,int) %>% spread(sample,int) 
out.sample <- data.frame(category = c("Morula","ICM-EPI","TE-EXE")) %>% 
  left_join(sample)
wb<- createWorkbook()

addWorksheet(wb,"FileInformation")
writeData(wb,"FileInformation",out.sample)
addWorksheet(wb,"RawIntensity")
writeData(wb,"RawIntensity",out.raw)
addWorksheet(wb,"IdentificationCount")
writeData(wb,"IdentificationCount",count.pro)
addWorksheet(wb,"NormalizedIntensity")
writeData(wb,"NormalizedIntensity",out.nor)
addWorksheet(wb,"ZscoreInteisty")
writeData(wb,"ZscoreInteisty",out.scale)

saveWorkbook(wb,file = "../result/Table1.Proteome.xlsx",overwrite = T)


