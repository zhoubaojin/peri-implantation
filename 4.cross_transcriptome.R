library(tidyverse)
library(openxlsx)
file <- read.xlsx("../raw/WT_EPI EXE samples_FPKM_matrix.xlsx")


sample <- read.xlsx("Supplementary Table 1.xlsx",sheet = "FileInformation")

sample.rna <- data.frame(c = c("ICM-3.5","Epi-4.5","Epi-5.0_invivo","ExE-5.5_invivo","TE-3.5","TE-4.5","ExE-4.75","ExE-5.5_invivo"),
                         class = c("3.5d-ICM","4.5d-ICM","5.0d-EPI","5.5d-EPI","3.5d-TE","4.5d-TE","5.0d-EXE","5.5d-EXE"))




df.rna <- file %>% 
  gather(sample,fpkm,-X1) %>% 
  #filter(fpkm > 0 ) %>% 
  mutate(c = str_replace_all(sample,"-\\d$","")) %>% 
  inner_join(sample.rna,by = "c") %>% 
  group_by(class,X1) %>% summarise(fpkm= log2(mean(fpkm +1)))



load("../bin/mouse.uniprot.function.rda")
df.pro <- read.xlsx("Supplementary Table 1.xlsx",sheet = "NormalizedIntensity") %>% 
  gather(sample,int,-ID) %>% 
  left_join(sample)  %>% 
  filter(class %in% sample.rna$class) %>%
  left_join(fun$protein,by = c("ID" = "protein")) %>% 
  filter(!is.na(genename)) %>% 
  group_by(genename,class) %>% summarise(intensity = log10(mean(int)) )

ggVennDiagram::ggVennDiagram(list(rna = unique(df.rna$X1),
                                  pro = unique(df.pro$genename)))


df.inner <- data.frame(genename = unique(df.pro$genename))%>% filter(genename %in% df.rna$X1)

df.cor <- data.frame()
for(i in df.inner$genename){

  i.pro <- df.pro %>% filter(genename == i)
  df.cross <- df.rna %>% filter(X1 == i) %>% 
    left_join(i.pro,by = "class")
  df.i <- data.frame(genename = i,cor = cor(df.cross$fpkm,df.cross$intensity,method = "pearson"))
  df.cor <- df.cor %>% bind_rows(df.i)
}


to.cor.rna <- df.rna %>%
  filter(X1 %in% df.inner$genename) %>% 
  spread(class,fpkm)  %>% select(-X1)

colnames(to.cor.rna) <- paste0("mRNA-",colnames(to.cor.rna))

to.cor.pro<- df.pro %>% ungroup( )%>% 
  filter(genename %in% df.inner$genename) %>% 
  spread(class,intensity) %>% select(-genename)


colnames(to.cor.pro) <- paste0("protein-",colnames(to.cor.pro))

to.cor <- to.cor.rna %>% bind_cols(to.cor.pro)


df <- cor(to.cor) %>% as.data.frame() %>% 
  select(1:8) %>% 
  slice(9:16)


df %>% bind_cols(proteome = row.names(df)) %>% 
  
  gather(transcriptome,cor,-proteome) %>% 
  mutate(transcriptome = str_replace_all(transcriptome,"mRNA-","")) %>% 
  mutate(proteome = str_replace_all(proteome,"protein-","")) %>% 
  ggplot(aes(proteome,transcriptome)) +
  geom_tile(aes(fill = cor)) +
  geom_text(aes(label = round(cor,2)),size = 3) +
  labs(fill = "Correlation") +
  scale_fill_gradient(low = "white",high = "steelblue") +
  theme(axis.text.x = element_text(angle = 90,hjust = 1),
        legend.position = "top") 

ggsave("../result/1.global/2.mrna.png",width = 3.5,height = 4)
