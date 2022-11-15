library(DESeq2)
library(apeglm)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library("edgeR")
library('pheatmap')
library('RColorBrewer')
library(RUVSeq)
library(EDASeq)

######Start with the featurecount table###########

files <- list.files(pattern = '.txt')    
for(i in files) {
  x <- read_delim(i,  delim = "\t", escape_double = FALSE, 
                  col_types = cols(Chr = col_skip(), Start = col_skip(), 
                                   End = col_skip(), Strand = col_skip(), 
                                   Length = col_skip()), trim_ws = TRUE, 
                  skip = 1) %>% 
    rename_at(vars(contains('GRCm39_107_')), list( ~ gsub('GRCm39_107_', '', .))) %>%  
    rename_at(vars(contains('_trimmed_finalAligned.sortedByCoord.out.bam')), list( ~ gsub('_trimmed_finalAligned.sortedByCoord.out.bam', '', .))) %>%  
    rename_at(vars(contains('-')), list( ~ gsub('-', '_', .)))
  assign(i,x)  
}


Metodikai_antisense_NM <- featureCounts_GRCm39.107_s2_no_multi.txt %>% 
  setnames(., old = c("smpl01", "smpl02", "smpl03", "smpl04",  "smpl05",  
                      "smpl06",  "smpl07", "smpl08",  "smpl09",  "smpl10",  "smpl11",  "smpl12"),
           new = c( "CPU300_1", "CPU30_2", "CPU300_3", "CPU30_4",
                    "CPU300_5", "CPU30_6", "CPU300_7", "CPU30_8", "MS_330_9", "MS_330_10", "MS_330_11", "MS_330_12"))

Annotiationtable107 = getBM(attributes = c("ensembl_gene_id","external_gene_name", "description"), 
                            filters = c('ensembl_gene_id') ,
                            values = Metodikai_NS_NM$Geneid,
                            mart = ensembl107)

Metodikai_antisense_NM_biomart <- inner_join(Metodikai_antisense_NM, Annotiationtable107, by = c("Geneid" = "ensembl_gene_id")) %>% 
  dplyr::select(Geneid,external_gene_name, description, everything())


Metodikai_antisense_NM_cpm <- as.data.frame(cpm(Metodikai_antisense_NM)) %>% 
  rownames_to_column(., var='Geneid')


Metodikai_antisense_NM_biomart_cpm <- Metodikai_antisense_NM_biomart %>% 
  inner_join(Metodikai_antisense_NM_cpm_jav, by = c("Geneid" = "cpm_Geneid"))


data_met_as_nm <- Metodikai_antisense_NM %>% 
  data.frame(row.names = 1)

#####stack all aligned bar plot#####

dfstack_all <- data.frame(all_stack)
dlstack_all <-data.frame(stack_points)


metastack_all<-full_join(dfstack_all, dlstack_all, by = c('Name'='Name')) 
dfstack_all$Name = factor(dfstack_all$Name,
                       levels = c("CPU30", "CPU300", "MS330"))

dfstack_all$Type = factor(dfstack_all$Type,
                          levels = c("Total", "Aligned","Assigned"))
ggplot()+
  geom_bar(data = dfstack_all, 
           aes(x=mean, y = Name , fill=Type),
           position="stack",stat="identity", width = 0.8, colour="black", size=0.3)+
  labs(x=NULL, y=NULL)+
  theme_classic()+
  scale_x_continuous(expand=c(0, 0), limits = c(0, 50), n.breaks =10)+
  theme(legend.position="bottom", legend.direction="horizontal")+
  geom_point(data = Astack_points, aes(x= Aligned , y=Name), position = position_dodge2(width = 0.3),
            show.legend = FALSE, shape = 1, color= "red", size = 0.8, stroke = 0.5)+
  geom_point(data = Tstack_points, aes(x= Aligned , y= Name), position = position_dodge2(width = 0.3),
             show.legend = FALSE, shape = 1, color= "purple", size = 0.8, stroke = 0.5)+
  geom_point(data = assigned_points, aes(x= Aligned , y= Name), position = position_dodge2(width = 0.3),
             show.legend = FALSE, shape = 1, color= "black", size = 0.8, stroke = 0.5)+
  geom_errorbar(data = dfstack_all, 
                aes(x= Error, y = Name, xmin=Error-SEM, xmax=Error+SEM), width=0.2, size = 0.3, position = "identity")+
  scale_fill_manual(values = c("#FFCC00", "#0000FE","#98CB00"))+
  theme(axis.ticks = element_line(colour = "black"),
                                        axis.text = element_text(colour = "black")) +labs(x = "Reads (M)")+
  theme(text = element_text(size = 5)) + theme(axis.text = element_text(size = 6),
                                               axis.text.x = element_text(colour = "black"),
                                               axis.text.y = element_text(colour = "black"))+
  theme(legend.key.size = unit(2, 'mm'))+
  theme(legend.text = element_text(size=5))+
  theme(legend.title = element_text(size=5)) + theme(axis.line = element_line(size = 0.2)) + theme(axis.ticks = element_line(size = 0.2))
ggsave("total_stack.pdf", width = 80 , height = 40 , units = "mm")

######Transcriptome coverage#####

Count_as_nm <- as.data.frame(apply(t(data_met_as_nm), 1, function(x) length(x[x>=5])))

colnames(Count_as_nm)
coverage_as_nm <- Count_as_nm %>% 
  rename(Gene_number = 'apply(t(data_met_as_nm), 1, function(x) length(x[x >= 5]))') %>% 
  rownames_to_column(var='Name') %>% 
  arrange(Name) %>% 
  mutate(Type = c("CPU 30", "CPU 30", "CPU 30", "CPU 30", 
                  "CPU 300", "CPU 300", "CPU 300", "CPU 300",
                  "MS 330", "MS 330", "MS 330", "MS 330"))

pointNDG <- coverage_as_nm %>% 
  select(Type, Gene_number) %>% 
  rename(Name = Type)

meansem_NDG <- coverage_as_nm %>% 
  pivot_wider(., names_from = Type, values_from = Gene_number) %>%
  mutate(Type = c("CPU 30", "CPU 30", "CPU 30", "CPU 30", 
                  "CPU 300", "CPU 300", "CPU 300", "CPU 300",
                  "MS 330", "MS 330", "MS 330", "MS 330")) %>% 
  select("CPU 30", "CPU 300", "MS 330") %>%
  t() %>%
  data.frame() %>% 
  mutate(Mean = rowMeans(., na.rm = TRUE)) %>%
  mutate(SEM = rowSds(as.matrix(.), na.rm = TRUE)/sqrt(4)) %>%
  data.frame() %>% 
  rownames_to_column(., var = "Name") %>% 
  select(Name, Mean, SEM)




dfNDG <- data.frame(meansem_NDG)	
dlNDG <-data.frame(pointNDG)	

metaNDG<-inner_join(dfNDG, dlNDG, by = c('Name'='Name')) 

dfNDG$Name = factor(dfNDG$Name,
                    levels = c("MS 330","CPU 300", "CPU 30"))


ggplot(data = dfNDG, 
       aes(x=Name, y = Mean, 		 
           ymin=Mean-SEM, ymax=Mean+SEM, 
           fill=Name))+ 
  geom_col(show.legend = FALSE, colour = "black", size = 0.2)+	
  labs(x=NULL, y="Number of detected genes \n (Read>=5)")+		
  theme_classic()+				
  scale_y_continuous(expand=c(0, 0), limits = c(0, 16500), n.breaks = 10)+ 
  scale_fill_manual(values = c("#0000FE", "#FFCC00", "#98CB00"))+
  geom_errorbar(width=0.1, size = 0.2)+	
  geom_jitter(data = metaNDG, aes(x= Name, y=Gene_number), position = position_dodge2(width = 0.3),
              show.legend = FALSE, shape = 1, size = 0.7, stroke = 0.25)+	
  theme(axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black"))+
  theme(axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black")) +labs(x = "Reads (M)")+
  theme(text = element_text(size = 5)) + theme(axis.text = element_text(size = 6),
                                               axis.text.x = element_text(colour = "black"),
                                               axis.text.y = element_text(colour = "black"))+
  theme(legend.key.size = unit(2, 'mm'))+
  theme(legend.text = element_text(size=5))+
  theme(legend.title = element_text(size=5)) + theme(axis.line = element_line(size = 0.2)) + theme(axis.ticks = element_line(size = 0.2))
ggsave("NGD.pdf", width = 50 , height = 35.375 , units = "mm")   




####coeficient of variance####

CV_filter_table<- Metodikai_antisense_NM_cpm_jav%>% 
  mutate(Mean = rowMeans(.[2:13])) %>% 
  arrange(-Mean) %>% 
  slice_head(., n=1000)

###CPU300CV#####
CV_data_300 <- CV_filter_table %>% 
  select(cpm_Geneid, cpm_CPU300_1, cpm_CPU300_3, cpm_CPU300_5, cpm_CPU300_7) %>%
  data.frame(row.names = 1) %>% 
  t() %>% 
  data.frame()

CV_data_300_all <- as.data.frame(sapply(CV_data_300, function(x) sd(x) / mean(x) * 100))
colnames(CV_data_300_all) <- c("Name","CV")
CV_data_300_all <- CV_data_300_all %>% 
  rownames_to_column(., var = "Geneid")

CV_data_300_all <- CV_data_300_all %>%
  mutate(Name = "CPU300_CV") %>% 
  select(Name, `CV of CPU300`)

###CPU30CV#####
CV_data_30 <- CV_filter_table %>% 
  select(cpm_Geneid, cpm_CPU30_2, cpm_CPU30_4, cpm_CPU30_6, cpm_CPU30_8) %>%
  data.frame(row.names = 1) %>% 
  t() %>% 
  data.frame()

CV_data_30_all <- as.data.frame(sapply(CV_data_30, function(x) sd(x) / mean(x) * 100))
colnames(CV_data_30_all) <- c("Name","CV")
CV_data_30_all <- CV_data_30_all %>% 
  rownames_to_column(., var = "Geneid")

CV_data_30_all <- CV_data_30_all %>%
  mutate(Name = "CPU30_CV") %>% 
  select(Name, `CV of CPU30`)

###MS330CV#####
CV_data_330 <- CV_filter_table %>% 
  select(cpm_Geneid, cpm_MS_330_9, cpm_MS_330_10, cpm_MS_330_11, cpm_MS_330_12) %>%
  data.frame(row.names = 1) %>% 
  t() %>% 
  data.frame()

CV_data_330_all <- as.data.frame(sapply(CV_data_330, function(x) sd(x) / mean(x) * 100))
colnames(CV_data_330_all) <- c("Name","CV")

CV_data_330_all <- CV_data_330_all %>% 
  rownames_to_column(., var = "Geneid")

CV_data_330_all <- CV_data_330_all %>%
  mutate(Name = "MS330_CV") %>% 
  select(Name, `CV of MS330`) %>% 
  
CV_data_all <- CV_data_300_all %>% 
  bind_rows(., CV_data_30_all) %>% 
  bind_rows(., CV_data_330_all)


####Violin-plot for CVs#####
df <- data.frame(meansem)
dl <-data.frame(CV_data_all)

dl$Name = factor(dl$Name,
                 levels = c("MS330_CV","CPU300_CV", "CPU30_CV"))
ggplot(data = dl, 
       aes(x=Name, y = CV, fill=Name))+
  geom_violin(adjust = 2, show.legend = FALSE, size=0.1)+
  geom_dotplot(binaxis = "y", binwidth = 0.6, stackdir = "centerwhole", show.legend = FALSE, stroke = 0.1)+
  scale_fill_manual(values = c("#0000FE", "#FFCC00", "#98CB00"))+
  labs(x=NULL, y=NULL)+
  theme_classic()+
  stat_summary(fun = mean, geom = "crossbar", alpha = 0.6, show.legend = FALSE, size=0.05)+
  scale_y_continuous(expand=c(0, 0), limits = c(0, 190), n.breaks = 10)+
  theme(text = element_text(size = 6))+  theme(axis.ticks = element_line(colour = "black"),
                                               axis.text = element_text(colour = "black"))+
  theme(axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black"))+
  theme(text = element_text(size = 5)) + theme(axis.text = element_text(size = 5),
                                               axis.text.x = element_text(colour = "black"),
                                               axis.text.y = element_text(colour = "black"))+
  theme(legend.key.size = unit(2, 'mm'))+
  theme(legend.text = element_text(size=5))+
  theme(legend.title = element_text(size=5)) + theme(axis.line = element_line(size = 0.2)) + theme(axis.ticks = element_line(size = 0.2))+
  scale_x_discrete(labels=c("MS 330 CV", "CPU 300 CV", "CPU 30 CV")) + theme(axis.ticks = element_line(colour = "black"),
                                                                             axis.text = element_text(colour = "black")) + theme(axis.line = element_line(size = 0.12),
                                                                                                                                 axis.ticks = element_line(size = 0.12))
ggsave("violin_50x35.375.pdf", width = 50 , height = 35.375 , units = "mm")

#####Create the PCA plot with DESeq2######

infomet_as_nm = data.frame(Names=colnames(data_met_as_nm))
rownames(infomet_as_nm) = infomet_as_nm$Names
infomet_as_nm$type = 1:12
infomet_as_nm$type = 'CPU30'
infomet_as_nm$type[c(1,3,5,7)] = 'CPU300'
infomet_as_nm$type[c(9,10,11,12)] = 'MS330'

infomet_as_nm$type = factor(infomet_as_nm$type)

dds_as_nm <- DESeqDataSetFromMatrix(countData =data_met_as_nm, colData = infomet_as_nm, design = ~type)

rlog_as_nm = rlog(dds_as_nm)
pcadata_as_nm <-plotPCA(rlog_as_nm, intgroup = "type", returnData=TRUE)
percentVar_as_nm <- round(100 * attr(pcadata_as_nm, "percentVar"))
ggplot(pcadata_as_nm, aes(PC1, PC2, color=type, label=name)) +
  xlab(paste0("PC1: ",percentVar_as_nm[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_as_nm[2],"% variance")) + 
  coord_fixed()+
  geom_jitter(size= 0.1)+
  geom_text_repel(size= 1)+
  scale_color_manual(values = c( "#98CB00","#FFCC00","#0000FE"))+
  theme_bw()+
  stat_ellipse(level = 0.8, size= 0.1)+
  theme(legend.position="bottom", legend.direction="horizontal")+
  theme(text = element_text(size = 5)) + theme(axis.text = element_text(size = 6),
                                               axis.text.x = element_text(colour = "black"),
                                               axis.text.y = element_text(colour = "black"))+
  theme(legend.key.size = unit(2, 'mm'))+
  theme(legend.text = element_text(size=5))+
  theme(legend.title = element_text(size=5))+
  scale_y_continuous(limits = c(-60,38), n.breaks = 8)+
  scale_x_continuous(limits = c(-50,70), n.breaks = 10)+
  ggtitle("My PCA Graph") + theme(axis.title = element_text(size = 6)) + theme(panel.grid.major = element_line(size = 0.1),
                                                                               panel.grid.minor = element_line(size = 1)) + theme(panel.grid.minor = element_line(size = 0.1),
                                                                                                                                  panel.background = element_rect(size = 0.1),
                                                                                                                                  plot.background = element_rect(size = 0.1)) + theme(axis.line = element_line(size = 0.1),
                                                                                                                                                                                      axis.ticks = element_line(colour = "black",
                                                                                                                                                                                                                size = 0.1))+ theme(panel.border = element_rect(size = 0.1))
ggsave("PCA.pdf", width = 60 , height = 60 , units = "mm")

####heatmap CPU30 vs CPU300####	
test <- data.frame(met1013_top1000, row.names = 1)	#input: top 1000 most abundant gene in cpu30 and cpu300
set.seed(100)
pheatmap(test,					
         scale = "row", 					
         color = colorRampPalette(rev(brewer.pal(n =11, name = "RdBu")))(200),	
         kmeans_k = 500,					
         cellwidth = 20,					
         cellheight = 1,					
         fontsize = 12,					
         border_color = NA,				        
         cluster_cols = FALSE,			
         clustering_method = "average",						
         show_rownames = FALSE,			
         angle_col = 45,					
         gaps_col = c(4),				
         filename = "1013_top1000.pdf")
#####Human vs mouse top100 receptor####

Mouse_receptors <- receptors

ensembl_human_all = getBM(attributes = c("ensembl_gene_id","external_gene_name", "description"), 
                          filters = c('ensembl_gene_id') ,
                          values = human_cpm_ensembl$...1,
                          mart = ensemblhuman99)

Human_str_teljes_cpm <- Human_str_teljes %>%
  select(-symbol) %>%
  unique() %>% 
  data.frame(row.names = 1) %>% 
  cpm() %>%
  data.frame() %>% 
  rownames_to_column(., var = "Ensembl")

Human_str_teljes_cpm <- Human_str_teljes_cpm %>% 
  inner_join(Human_str_teljes, by = c("Ensembl" = "ens")) %>% 
  select(Ensembl, symbol, "ChIN..21","ChIN..22")

Human_receptors_str <- Human_all_receptors %>% 
  inner_join(Human_str_teljes_cpm, by = c("Gene_symbol" = "symbol")) %>%
  select(Gene_symbol,"ChIN..21","ChIN..22") %>%
  filter(!grepl("OR", Gene_symbol )) %>% 
  mutate(Mean_hC = rowMeans(.[,2:3])) %>% 
  arrange(-Mean_hC) %>%
  unique() %>% 
  mutate(gene_name = tolower(Gene_symbol)) %>% 
  inner_join(human_mouse_converter, by = c("Gene_symbol" = "Human")) %>% 
  slice_head(., n = 100)

Mouse_receptors_human <- Mouse_receptor_ion %>% 
  inner_join(Metodikai_antisense_NM_biomart_cpm, by =c("symbol" = "external_gene_name")) %>% 
  select(symbol, cpm_CPU300_1, cpm_CPU300_3, cpm_CPU300_5, cpm_CPU300_7) %>% 
  mutate(Mean_mm = rowMeans(.[2:5])) %>% 
  arrange(-Mean_mm) %>%
  mutate(gene_name = tolower(symbol)) %>% 
  slice_head(., n= 100)

Közös_lista_MH <- Mouse_receptors_human %>% 
  inner_join(Human_receptors_str, by = c("symbol" = "Mouse")) %>% 
  select(symbol, Gene_symbol, Mean_mm, Mean_hC) %>% 
  mutate(Melyik = "Közös")

Csak_human_receptor <- Human_receptors_str %>% 
  anti_join(Mouse_receptors_human, c("Mouse" = "symbol" )) %>% 
  select(Gene_symbol, gene_name, Mean_hC, Mouse)

Csak_human_receptor_mouse <- Csak_human_receptor %>% 
  inner_join(Mouse_receptors_human_all, by= c("Mouse" = "symbol" )) %>% 
  select(Gene_symbol, Mouse, Mean_hC, Mean_mm) %>% 
  mutate(Melyik = "Csak humán")

Csak_human_receptor_mouse_neg <- Csak_human_receptor %>% 
  anti_join(Csak_human_receptor_mouse, by= "Gene_symbol")


colnames(Csak_human_receptor_mouse)

Csak_mouse_receptor <- Mouse_receptors_human %>% 
  anti_join(Human_receptors_str, by = c("symbol" = "Mouse" )) %>% 
  select(symbol, gene_name, Mean_mm)

Csak_mouse_receptor_human <- Csak_mouse_receptor %>% 
  inner_join(Human_receptors_str_all, by= c("symbol" = "Mouse" )) %>% 
  select(symbol, symbol, Mean_mm, Mean_hC ) %>% 
  mutate(Melyik = "Csak egér")

Csak_mouse_receptor_human_neg <- Csak_mouse_receptor %>% 
  anti_join(Csak_mouse_receptor_human, by= "symbol")

TopHuman_mouse_receptor_list <- Közös_lista_MH %>% 
  bind_rows(., Csak_human_receptor_mouse) %>% 
  bind_rows(., Csak_mouse_receptor_human) %>% 
  unique() %>% 
  select(Gene_symbol, symbol, Mean_hC, Mean_mm, Melyik)


Human_receptors_str_source <- Human_all_receptors %>% 
  inner_join(Human_str_teljes_cpm, by = c("Gene_symbol" = "symbol")) %>%
  filter(!grepl("OR", Gene_symbol )) %>% 
  mutate(Mean = rowMeans(.[,3:4])) %>% 
  arrange(-Mean) %>%
  unique() %>% 
  mutate(gene_name = tolower(Gene_symbol)) %>% 
  inner_join(human_mouse_converter, by = c("Gene_symbol" = "Human")) %>% 
  slice_head(., n = 100) %>% 
  select(Ensembl, Gene_symbol, ChIN..21, ChIN..22, Mean) %>% 
  rename("Ensembl ID" = Ensembl) %>% 
  rename("ChIN #21" = ChIN..21) %>% 
  rename("ChIN #22" = ChIN..22)

Mouse_receptors_human_source <- Mouse_receptor_ion %>% 
  inner_join(Metodikai_antisense_NM_biomart_cpm, by =c("symbol" = "external_gene_name")) %>% 
  select(Geneid, symbol, cpm_CPU300_1, cpm_CPU300_3, cpm_CPU300_5, cpm_CPU300_7) %>% 
  mutate(Mean = rowMeans(.[3:6])) %>% 
  arrange(-Mean) %>%
  mutate(gene_name = tolower(symbol)) %>% 
  slice_head(., n= 100) %>%
  select(Geneid, symbol, cpm_CPU300_1, cpm_CPU300_3, cpm_CPU300_5, cpm_CPU300_7, Mean) %>% 
  rename("Ensembl ID" = Geneid) %>%
  rename_at(vars(contains('cpm_')), list( ~ gsub('cpm_', '', .))) %>%
  rename("Gene_symbol" = symbol)

####Analysis with RUV_seq and DESeq2#####

data_met_as_nm_300vs330 <- data_met_as_nm %>% 
  dplyr::select(!c(CPU30_2, CPU30_4, CPU30_6, CPU30_8))
data_met_as_nm_300vs330 <- data300vs330
data_met_as_nm_300vs330 <- data.frame(row300vs330, row.names = 1)

filter <- apply(data_met_as_nm_300vs330, 1, function(x) length(x[x>5])>=2)
filtered <- data_met_as_nm_300vs330[filter,]
x <- as.factor((rep(c("CHAT300", "MS330"), each=4)))
set <- newSeqExpressionSet(as.matrix(filtered),
                           phenoData = data.frame(x, row.names=colnames(filtered)))


library(RColorBrewer)
colors <- brewer.pal(3, "Set2")
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set, col=colors[x], cex=1.2)

set <- betweenLaneNormalization(set, which="upper")
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set, col=colors[x], cex=1.2)

design <- model.matrix(~x, data=pData(set))
y <- DGEList(counts=counts(set), group=x)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)
top <- topTags(lrt, n=nrow(set))$table
empirical <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1:5000]))]


set1 <- RUVg(set, empirical, k=1)
pData(set1)

plotRLE(set1, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set1, col=colors[x], cex=1.2)


dds <- DESeqDataSetFromMatrix(countData = counts(set1),
                              colData = pData(set1),
                              design = ~ W_1 + x)
dds <- DESeq(dds)
res <- results(dds, alpha =0.05)
summary(res)

resultsNames(dds)
res_slfc <- lfcShrink(dds, coef = "x_MS330_vs_CHAT300", type="apeglm")

rest_slfc = as.data.frame(res_slfc)
rest_slfc = rownames_to_column(rest_slfc, var='ensgene')
rest = as.data.frame(res)
rest = rownames_to_column(rest, var='ensgene')
rest = as.data.frame(res)
rest = rownames_to_column(rest, var='ensgene')

res__met_as_nm_300vs330 = getBM(attributes = c('ensembl_gene_id',"external_gene_name", "description"), 
                                filters = 'ensembl_gene_id',
                                values = rest_slfc$ensgene,
                                mart = ensembl107)

rest_slfc_with_all = left_join(rest_slfc, res__met_as_nm_300vs330,
                               by= c("ensgene" = "ensembl_gene_id"))

resrest_slfc_with_all__met_as_nm_300vs330 = rest_slfc_with_all %>% 
  dplyr::select(ensgene, external_gene_name,description, everything()) %>%
  inner_join(Metodikai_antisense_M_cpm_jav, by = c("ensgene" = "cpm_Geneid")) %>% 
  dplyr::select(!c(baseMean, lfcSE, cpm_CPU30_2, cpm_CPU30_4, cpm_CPU30_6, cpm_CPU30_8)) 



####Volcano-plot######

Envolcano_as_mn <- resrest_slfc_with_all__met_as_nm_300vs330 %>% 
  na.omit() %>%
  dplyr::select(ensgene, external_gene_name, log2FoldChange, padj)

asd <- resrest_slfc_with_all__met_as_nm_300vs330 %>%
  dplyr::select(ensgene, external_gene_name, log2FoldChange, padj)

asd_withlo10padj <- asd %>% 
  mutate(log10padj = log10(padj))

asd1 <- asd_withlo10padj %>%
  filter(log2FoldChange < -12 |
           log2FoldChange > 12 |
           log10padj < -150)

selected_genes <- asd1$external_gene_name

asd2 <- inner_join(asd1, Envolcano_as_mn, by = "ensgene")

p1 <- EnhancedVolcano(Envolcano_as_mn,
                      lab = Envolcano_as_mn$external_gene_name,
                      x = 'log2FoldChange',
                      y = 'padj',
                      selectLab = selected_genes,
                      labSize = 1,
                      shape= c(16,16,16,16),
                      axisLabSize = 1,
                      labCol = 'black',
                      labFace = 'bold',
                      colAlpha = 1,
                      title='', subtitle='',
                      pCutoff = 0.05,
                      FCcutoff = 1, col = c('grey30', 'grey30', 'royalblue', 'red2'),
                      pointSize = 0.2,
                      drawConnectors = TRUE,
                      widthConnectors = 0.1,
                      lengthConnectors = unit(0.01, "npc"),
                      arrowheads = FALSE,
                      cutoffLineWidth = 0.1,
                      hlineWidth = 0.1,
                      vlineWidth = 0.1,
                      borderWidth = 0.3,
                      boxedLabels = FALSE,
                      max.overlaps = 60,
                      legendLabels = c('NS', expression(Log[2]~FC), 'p-value', expression(p-adj.~and~log[2]~FC)),
                      legendPosition = 'top',
                      colConnectors = "black",
                      directionConnectors = "both",
                      legendLabSize = 0.5,
                      legendIconSize = 0.5,
                      caption = bquote(~Log[2]~"fold change cutoff: 1; adjusted p-value cutoff: 0.05"))

p1+ggplot2::scale_x_continuous(
  breaks=seq(-15,15, 5))+ 
  scale_y_continuous(
    breaks=seq(0,300, 50))+
  theme_bw()+
  theme(axis.text = element_text(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))+
  theme(legend.position="top")+theme(text = element_text(size = 5)) + theme(axis.text = element_text(size = 6),
                                                                            axis.text.x = element_text(colour = "black"),
                                                                            axis.text.y = element_text(colour = "black"))+
  theme(legend.key.size = unit(2, 'mm'))+
  theme(legend.text = element_text(size=5))+
  theme(legend.title = element_text(size=5)) + theme(axis.line = element_line(size = 0.2)) + theme(axis.ticks = element_line(size = 0.2))
ggsave("volcano_plot.pdf", width = 82 , height = 78 , units = "mm")





####heatmap MS330 vs CPU300####	

test <- data.frame(met1013_DE_sign, row.names = 1)	#input: counts with significant changes in ms330 vs. cpu300
test <- rlog(as.matrix(test), blind = FALSE)
set.seed(99)
pheatmap(test,					
         scale = "row", 					
         color = colorRampPalette(rev(brewer.pal(n =11, name = "PuOr")))(200),	
         kmeans_k = 500,					
         cellwidth = 20,					
         cellheight = 1,					
         fontsize = 12,					
         border_color = NA,				
         cluster_cols = FALSE,			
         clustering_method = "average",			
         show_rownames = FALSE,			
         angle_col = 45,					
         gaps_col = c(4),				
         filename = "1013_DE_rlog.pdf")
#####Significant transcription factors#####

TF_list <- TF_list_KEGG_brite %>%
  inner_join(Annotiationtable107, by = c("Gene symbol" = "external_gene_name")) %>%
  select(ensembl_gene_id)



TF_list_desc = getBM(attributes = c("ensembl_gene_id","description"), 
                     filters = c('ensembl_gene_id') ,
                     values = TF_list$Ensembl,
                     mart = ensembl104)

TF_list <- TF_list %>% 
  inner_join(TF_list_desc, by =c("Ensembl" = "ensembl_gene_id"))

write_xlsx(TF_list, "C:\\Users\\gocz.balazs\\Documents\\R\\codes_plots\\TF_list_KEGG_brite.xlsx")


TF_met <- semi_join(resrest_slfc_with_all__met_as_nm_300vs330, TF_list, by = c("ensgene" = "ensembl_gene_id"))

######Significant transporters#####

transporter_list_szig <- semi_join(resrest_slfc_with_all__met_as_nm_300vs330, transporter, by = c("ensgene" = "ensembl_gene_id"))

######Significant ion channels######

ion_channel_list_szig <- semi_join(resrest_slfc_with_all__met_as_nm_300vs330, Ion_channel, by = c("ensgene" = "ensembl_gene_id"))

######Significant receptors######

receptor_list_szig <- semi_join(resrest_slfc_with_all__met_as_nm_300vs330, receptors, by = c("ensgene" = "ensembl_gene_id"))


#####Boxplots for significant categories####

#########ion channels########

IC_selected_genes_up<-ion_channel_list_szig %>%
  filter(padj<0.05) %>% 
  mutate(MeanCPU = rowMeans(.[4:7])) %>% 
  mutate(MeanMS = rowMeans(.[8:11])) %>% 
  filter(MeanCPU >20 | MeanMS >20) %>% 
  .[1:5,]

IC_selected_genes_down <-ion_channel_list_szig %>% 
  filter(padj<0.05) %>% 
  mutate(MeanCPU = rowMeans(.[4:7])) %>% 
  mutate(MeanMS = rowMeans(.[8:11])) %>% 
  filter(MeanCPU >20 | MeanMS >20) %>% 
  arrange(log2FoldChange) %>% 
  .[1:5,]

IC_selected_genes <- IC_selected_genes_up %>% 
  bind_rows(., IC_selected_genes_down) %>% 
  select(external_gene_name, cpm_CPU300_1:cpm_MS_330_12) #%>% 
#data.frame(row.names = 1)

IC_selected_genes <-data.frame(IC_selected_genes, row.names = 1)



sg <- split(IC_selected_genes, rownames(IC_selected_genes))

for (i in rownames(IC_selected_genes)) {
  assign(paste0("sg_", i), as.data.frame(t(sg[[i]])) %>%
           mutate(names = c((rep(c("CPU300", "MS330"), each=4)))) %>% 
           rename(points = i))
}


for (i in rownames(IC_selected_genes)) {
  assign(paste0("sg_meanse_", i), as.data.frame(t(sg[[i]])) %>%
           mutate(names = c((rep(c("CPU300", "MS330"), each=4)))) %>% 
           rename(points = i) %>% 
           group_by(names) %>% 
           summarise( 
             n=n(),
             mean=mean(points),
             sd=sd(points))%>%
           mutate(se=sd/sqrt(n)))
}

boxlista <- purrr::map2(lista, meanlista, dplyr::inner_join, by = "names")

for (i in rownames(IC_selected_genes)) {
  assign(paste0("sg_box_", i), as.data.frame(boxlista[[i]]))
}



meanlista <- list("Ano3" = sg_meanse_Ano3,
                  "Asic4" = sg_meanse_Asic4, 
                  "Cacna1e" = sg_meanse_Cacna1e, 
                  "Cacna2d3" = sg_meanse_Cacna2d3,
                  "Kcnk2" = sg_meanse_Kcnk2,
                  "Kcnq5" = sg_meanse_Kcnq5,
                  "Kcnt2" = sg_meanse_Kcnt2,
                  "Scn4b" = sg_meanse_Scn4b,
                  "Trpc3" = sg_meanse_Trpc3,
                  "Trpc6" = sg_meanse_Trpc6)


mutate(SEM = std_mean(lista[[i]]$points[1:4] & std_mean(lista[[i]]$points[5:8]))) %>%
  mutate(Mean1 = mean(lista[[i]]$points[1:4])) %>% 
  mutate(Mean2 = mean(lista[[i]]$points[5:8])) %>%
  
  lista <- list("Ano3" = sg_Ano3,
                "Asic4" = sg_Asic4, 
                "Cacna1e" = sg_Cacna1e, 
                "Cacna2d3" = sg_Cacna2d3,
                "Kcnk2" = sg_Kcnk2,
                "Kcnq5" = sg_Kcnq5,
                "Kcnt2" = sg_Kcnt2,
                "Scn4b" = sg_Scn4b,
                "Trpc3" = sg_Trpc3,
                "Trpc6" = sg_Trpc6)



lista_box <- list("Ano3" = sg_box_Ano3,
                  "Asic4" = sg_box_Asic4, 
                  "Cacna1e" = sg_box_Cacna1e, 
                  "Cacna2d3" = sg_box_Cacna2d3,
                  "Kcnk2" = sg_box_Kcnk2,
                  "Kcnq5" = sg_box_Kcnq5,
                  "Kcnt2" = sg_box_Kcnt2,
                  "Scn4b" = sg_box_Scn4b,
                  "Trpc3" = sg_box_Trpc3,
                  "Trpc6" = sg_box_Trpc6)


for (i in rownames(IC_selected_genes)) {
  lista_box[[i]]$names <- factor(lista_box[[i]]$names,
                                 levels = c("MS330","CPU300"))
  ggplot(data = lista_box[[i]],
         aes(x=names, y = points, fill=names, ymin=mean-2*se, ymax=mean+2*se))+
    stat_boxplot(geom = "errorbar",
                 width = 0.15,
                 size=0.1) +
    geom_boxplot(color= "black", show.legend = FALSE, outlier.shape = NA, coef = 0, width = 1, lwd=0.1, position = position_dodge(1))+
    
    # geom_jitter(position = position_dodge2(width = 0.3),
    #             show.legend = FALSE, shape = 21, size = 0.5, stroke = 0.1)+
    theme_classic()+labs(x = NULL, y = NULL) + theme(axis.text = element_text(colour = "black"),
                                                     axis.text.x = element_text(colour = "black"),
                                                     axis.text.y = element_text(colour = "black"),plot.margin = margin(t = 3,  
                                                                                                                       r = 1,  
                                                                                                                       b = 5,  
                                                                                                                       l = 1) ) + theme(panel.grid.major = element_line(colour = NA)) + theme(plot.title = element_text(family = "Helvetica",
                                                                                                                                                                                                                        size = 6, hjust = 0.5)) +labs(title = i)+
    scale_fill_manual(values = c("#FFD662FF","#00539CFF")) + theme(axis.ticks = element_line(colour = "black"),
                                                                   plot.background = element_rect(colour = NA)) + theme(axis.line = element_line(size = 5),
                                                                                                                        axis.title = element_text(size = 5),
                                                                                                                        axis.text = element_text(size = 5)) + theme(axis.line = element_line(size = 0.1),
                                                                                                                                                                    panel.grid.minor = element_line(colour = NA)) + theme(axis.ticks = element_line(size = 0.1),
                                                                                                                                                                                                                          panel.grid.major = element_line(size = 0.1),
                                                                                                                                                                                                                          panel.grid.minor = element_line(size = 0.1)) + theme(axis.text = element_text(angle = 90),
                                                                                                                                                                                                                                                                               plot.title = element_text(face = "italic"))+labs(title = NULL, y = NULL) + theme(axis.text = element_text(angle = 0),
                                                                                                                                                                                                                                                                                                                                                                axis.text.x = element_blank())+ scale_y_continuous(n.breaks = 10)
  ggsave(paste0("Boxplot_IC_", i, '.pdf'), width = 15 , height = 20 , units = "mm")
}

#######receptors#######

IC_selected_genes_up_rec<-receptor_list_szig %>% 
  mutate(MeanCPU = rowMeans(.[8:11])) %>% 
  mutate(MeanMS = rowMeans(.[12:15])) %>% 
  filter(MeanCPU >20 | MeanMS >20) %>%
  arrange(-log2FoldChange) %>% 
  .[1:5,]

IC_selected_genesdown_rec <-receptor_list_szig %>% 
  mutate(MeanCPU = rowMeans(.[8:11])) %>% 
  mutate(MeanMS = rowMeans(.[12:15])) %>% 
  filter(MeanCPU >20 | MeanMS >20) %>% 
  arrange(log2FoldChange) %>% 
  .[1:5,]

IC_selected_genes_rec <- IC_selected_genes_up_rec %>% 
  bind_rows(., IC_selected_genesdown_rec) %>% 
  select(external_gene_name, cpm_CPU300_1:cpm_MS_330_12) %>% 
  data.frame(row.names = 1)

sg_rec <- split(IC_selected_genes_rec, rownames(IC_selected_genes_rec))

for (i in rownames(IC_selected_genes_rec)) {
  assign(paste0("sg_rec_", i), as.data.frame(t(sg_rec[[i]])) %>%
           mutate(names = c((rep(c("CPU300", "MS330"), each=4)))) %>% 
           rename(points = i))
}


for (i in rownames(IC_selected_genes_rec)) {
  assign(paste0("sg_meanse_rec_", i), as.data.frame(t(sg_rec[[i]])) %>%
           mutate(names = c((rep(c("CPU300", "MS330"), each=4)))) %>% 
           rename(points = i) %>% 
           group_by(names) %>% 
           summarise( 
             n=n(),
             mean=mean(points),
             sd=sd(points)
           ) %>%
           mutate(se=sd/sqrt(n)))
}

lista_rec <- list("Adgrf5" = sg_rec_Adgrf5,
                  "Adrb1" = sg_rec_Adrb1, 
                  "Cnr1" = sg_rec_Cnr1, 
                  "Drd1" = sg_rec_Drd1,
                  "Drd2" = sg_rec_Drd2,
                  "Drd5" = sg_rec_Drd5,
                  "Gpr88" = sg_rec_Gpr88,
                  "Htr1a" = sg_rec_Htr1a,
                  "Ngfr" = sg_rec_Ngfr,
                  "Trhr" = sg_rec_Trhr)

meanlista_rec <- list("Adgrf5" = sg_meanse_rec_Adgrf5,
                      "Adrb1" = sg_meanse_rec_Adrb1, 
                      "Cnr1" = sg_meanse_rec_Cnr1,
                      "Drd1" = sg_meanse_rec_Drd1,
                      "Drd2" = sg_meanse_rec_Drd2,
                      "Drd5" = sg_meanse_rec_Drd5,
                      "Gpr88" = sg_meanse_rec_Gpr88,
                      "Htr1a" = sg_meanse_rec_Htr1a,
                      "Ngfr" = sg_meanse_rec_Ngfr,
                      "Trhr" = sg_meanse_rec_Trhr)

boxlista_rec <- purrr::map2(lista_rec, meanlista_rec, dplyr::inner_join, by = "names")


for (i in rownames(IC_selected_genes_rec)) {
  assign(paste0("sg_box_", i), as.data.frame(boxlista_rec[[i]]))
}

lista_box_rec <- list("Adgrf5" = sg_box_Adgrf5,
                      "Adrb1" = sg_box_Adrb1, 
                      "Cnr1" = sg_box_Cnr1, 
                      "Drd1" = sg_box_Drd1,
                      "Drd2" = sg_box_Drd2,
                      "Drd5" = sg_box_Drd5,
                      "Gpr88" = sg_box_Gpr88,
                      "Htr1a" = sg_box_Htr1a,
                      "Ngfr" = sg_box_Ngfr,
                      "Trhr" = sg_box_Trhr)



for (i in rownames(IC_selected_genes_rec)) {
  lista_box_rec[[i]]$names <- factor(lista_box_rec[[i]]$names,
                                     levels = c("MS330","CPU300"))
  ggplot(data = lista_box_rec[[i]],
         aes(x=names, y = points, fill=names, ymin=mean-2*se, ymax=mean+2*se))+
    stat_boxplot(geom = "errorbar",
                 width = 0.15,
                 size=0.1) +
    geom_boxplot(color= "black", show.legend = FALSE, outlier.shape = NA, coef = 0, width = 1, lwd=0.1, position = position_dodge(1))+
    
    # geom_jitter(position = position_dodge2(width = 0.3),
    #             show.legend = FALSE, shape = 21, size = 0.5, stroke = 0.1)+
    theme_classic()+labs(x = NULL, y = NULL) + theme(axis.text = element_text(colour = "black"),
                                                     axis.text.x = element_text(colour = "black"),
                                                     axis.text.y = element_text(colour = "black"),plot.margin = margin(t = 3,  
                                                                                                                       r = 1,  
                                                                                                                       b = 5,  
                                                                                                                       l = 1) ) + theme(panel.grid.major = element_line(colour = NA)) + theme(plot.title = element_text(family = "Helvetica",
                                                                                                                                                                                                                        size = 6, hjust = 0.5)) +labs(title = i)+
    scale_fill_manual(values = c("#FFD662FF","#00539CFF")) + theme(axis.ticks = element_line(colour = "black"),
                                                                   plot.background = element_rect(colour = NA)) + theme(axis.line = element_line(size = 5),
                                                                                                                        axis.title = element_text(size = 5),
                                                                                                                        axis.text = element_text(size = 5)) + theme(axis.line = element_line(size = 0.1),
                                                                                                                                                                    panel.grid.minor = element_line(colour = NA)) + theme(axis.ticks = element_line(size = 0.1),
                                                                                                                                                                                                                          panel.grid.major = element_line(size = 0.1),
                                                                                                                                                                                                                          panel.grid.minor = element_line(size = 0.1)) + theme(axis.text = element_text(angle = 90),
                                                                                                                                                                                                                                                                               plot.title = element_text(face = "italic"))+labs(title = NULL, y = NULL) + theme(axis.text = element_text(angle = 0),
                                                                                                                                                                                                                                                                                                                                                                axis.text.x = element_blank()) + scale_y_continuous(n.breaks = 10)
  ggsave(paste0("Boxplot_receptors_", i, '.pdf'), width = 15 , height = 20 , units = "mm")
}


#####TFs######

IC_selected_genes_up_tf<-TF_met %>%
  filter(padj<0.05) %>% 
  mutate(MeanCPU = rowMeans(.[4:7])) %>% 
  mutate(MeanMS = rowMeans(.[8:11])) %>% 
  filter(MeanCPU >20 | MeanMS >20) %>% 
  .[1:5,]

IC_selected_genesdown_tf <-TF_met %>%
  filter(padj<0.05) %>% 
  mutate(MeanCPU = rowMeans(.[4:7])) %>% 
  mutate(MeanMS = rowMeans(.[8:11])) %>% 
  filter(MeanCPU >20 | MeanMS >20) %>% 
  arrange(log2FoldChange) %>% 
  .[1:5,]

IC_selected_genes_tf <- IC_selected_genes_up_tf %>% 
  bind_rows(., IC_selected_genesdown_tf) %>% 
  select(external_gene_name, cpm_CPU300_1:cpm_MS_330_12) %>% 
  data.frame(row.names = 1)

sg_tf <- split(IC_selected_genes_tf, rownames(IC_selected_genes_tf))

for (i in rownames(IC_selected_genes_tf)) {
  assign(paste0("sg_tf_", i), as.data.frame(t(sg_tf[[i]])) %>%
           mutate(names = c((rep(c("CPU300", "MS330"), each=4)))) %>% 
           rename(points = i))
}


for (i in rownames(IC_selected_genes_tf)) {
  assign(paste0("sg_meanse_tf_", i), as.data.frame(t(sg_tf[[i]])) %>%
           mutate(names = c((rep(c("CPU300", "MS330"), each=4)))) %>% 
           rename(points = i) %>% 
           group_by(names) %>% 
           summarise( 
             n=n(),
             mean=mean(points),
             sd=sd(points)
           ) %>%
           mutate(se=sd/sqrt(n)))
}

lista_tf <- list("Csrp1" = sg_tf_Csrp1,
                 "Hmga1" = sg_tf_Hmga1, 
                 "Meis2" = sg_tf_Meis2, 
                 "Pbx1" = sg_tf_Pbx1,
                 "Pbx3" = sg_tf_Pbx3,
                 "Pou3f3" = sg_tf_Pou3f3,
                 "Rc3h1" = sg_tf_Rc3h1,
                 "Tshz2" = sg_tf_Tshz2,
                 "Tshz3" = sg_tf_Tshz3,
                 "Zbtb18" = sg_tf_Zbtb18)

meanlista_tf <- list("Csrp1" = sg_meanse_tf_Csrp1,
                     "Hmga1" = sg_meanse_tf_Hmga1, 
                     "Meis2" = sg_meanse_tf_Meis2, 
                     "Pbx1" = sg_meanse_tf_Pbx1,
                     "Pbx3" = sg_meanse_tf_Pbx3,
                     "Pou3f3" = sg_meanse_tf_Pou3f3,
                     "Rc3h1" = sg_meanse_tf_Rc3h1,
                     "Tshz2" = sg_meanse_tf_Tshz2,
                     "Tshz3" = sg_meanse_tf_Tshz3,
                     "Zbtb18" = sg_meanse_tf_Zbtb18)

boxlista_tf <- purrr::map2(lista_tf, meanlista_tf, dplyr::inner_join, by = "names")


for (i in rownames(IC_selected_genes_tf)) {
  assign(paste0("sg_box_tf_", i), as.data.frame(boxlista_tf[[i]]))
}

lista_box_tf <-  list("Csrp1" = sg_box_tf_Csrp1,
                      "Hmga1" = sg_box_tf_Hmga1, 
                      "Meis2" = sg_box_tf_Meis2, 
                      "Pbx1" = sg_box_tf_Pbx1,
                      "Pbx3" = sg_box_tf_Pbx3,
                      "Pou3f3" = sg_box_tf_Pou3f3,
                      "Rc3h1" = sg_box_tf_Rc3h1,
                      "Tshz2" = sg_box_tf_Tshz2,
                      "Tshz3" = sg_box_tf_Tshz3,
                      "Zbtb18" = sg_box_tf_Zbtb18)



for (i in rownames(IC_selected_genes_tf)) {
  lista_box_tf[[i]]$names <- factor(lista_box_tf[[i]]$names,
                                    levels = c("MS330","CPU300"))
  ggplot(data = lista_box_tf[[i]],
         aes(x=names, y = points, fill=names, ymin=mean-2*se, ymax=mean+2*se))+
    stat_boxplot(geom = "errorbar",
                 width = 0.15,
                 size=0.1) +
    geom_boxplot(color= "black", show.legend = FALSE, outlier.shape = NA, coef = 0, width = 1, lwd=0.1, position = position_dodge(1))+
    
    # geom_jitter(position = position_dodge2(width = 0.3),
    #             show.legend = FALSE, shape = 21, size = 0.5, stroke = 0.1)+
    theme_classic()+labs(x = NULL, y = NULL) + theme(axis.text = element_text(colour = "black"),
                                                     axis.text.x = element_text(colour = "black"),
                                                     axis.text.y = element_text(colour = "black"),plot.margin = margin(t = 3,  
                                                                                                                       r = 1,  
                                                                                                                       b = 5,  
                                                                                                                       l = 1) ) + theme(panel.grid.major = element_line(colour = NA)) + theme(plot.title = element_text(family = "Helvetica",
                                                                                                                                                                                                                        size = 6, hjust = 0.5)) +labs(title = i)+
    scale_fill_manual(values = c("#FFD662FF","#00539CFF")) + theme(axis.ticks = element_line(colour = "black"),
                                                                   plot.background = element_rect(colour = NA)) + theme(axis.line = element_line(size = 5),
                                                                                                                        axis.title = element_text(size = 5),
                                                                                                                        axis.text = element_text(size = 5)) + theme(axis.line = element_line(size = 0.1),
                                                                                                                                                                    panel.grid.minor = element_line(colour = NA)) + theme(axis.ticks = element_line(size = 0.1),
                                                                                                                                                                                                                          panel.grid.major = element_line(size = 0.1),
                                                                                                                                                                                                                          panel.grid.minor = element_line(size = 0.1)) + theme(axis.text = element_text(angle = 90),
                                                                                                                                                                                                                                                                               plot.title = element_text(face = "italic"))+labs(title = NULL, y = NULL) + theme(axis.text = element_text(angle = 0),
                                                                                                                                                                                                                                                                                                                                                                axis.text.x = element_blank())+ scale_y_continuous(n.breaks = 10)
  ggsave(paste0("Boxplot_tf_", i, '.pdf'), width = 15 , height = 20 , units = "mm")
}

#####Transporters####

IC_selected_genes_up_trans<-transporter_list_szig %>%
  filter(padj<0.05) %>%
  filter(!grepl("Sncg", external_gene_name)) %>%
  filter(!grepl("Fxyd6", external_gene_name)) %>%
  filter(!grepl("Lhfp", external_gene_name)) %>%
  mutate(MeanCPU = rowMeans(.[4:7])) %>% 
  mutate(MeanMS = rowMeans(.[8:11])) %>% 
  filter(MeanCPU >20 | MeanMS >20) %>%
  arrange(-log2FoldChange) %>% 
  .[1:5,]

IC_selected_genesdown_trans <-transporter_list_szig %>%
  filter(padj<0.05) %>%
  filter(!grepl("Penk", external_gene_name)) %>%
  filter(!grepl("Nkain3", external_gene_name)) %>%
  mutate(MeanCPU = rowMeans(.[4:7])) %>% 
  mutate(MeanMS = rowMeans(.[8:11])) %>% 
  filter(MeanCPU >20 | MeanMS >20) %>% 
  arrange(log2FoldChange) %>% 
  .[1:5,]

IC_selected_genes_trans <- IC_selected_genes_up_trans %>% 
  bind_rows(., IC_selected_genesdown_trans) %>% 
  select(external_gene_name, cpm_CPU300_1:cpm_MS_330_12) %>% 
  data.frame(row.names = 1)

sg_trans <- split(IC_selected_genes_trans, rownames(IC_selected_genes_trans))

for (i in rownames(IC_selected_genes_trans)) {
  assign(paste0("sg_trans_", i), as.data.frame(t(sg_trans[[i]])) %>%
           mutate(names = c((rep(c("CPU300", "MS330"), each=4)))) %>% 
           rename(points = i))
}


for (i in rownames(IC_selected_genes_trans)) {
  assign(paste0("sg_meanse_trans_", i), as.data.frame(t(sg_trans[[i]])) %>%
           mutate(names = c((rep(c("CPU300", "MS330"), each=4)))) %>% 
           rename(points = i) %>% 
           group_by(names) %>% 
           summarise( 
             n=n(),
             mean=mean(points),
             sd=sd(points)
           ) %>%
           mutate(se=sd/sqrt(n)))
}

lista_trans <- list("Slc17a8" = sg_trans_Slc17a8,
                    "Slc20a1" = sg_trans_Slc20a1,
                    "Slc24a3" = sg_trans_Slc24a3, 
                    "Slc32a1" = sg_trans_Slc32a1,
                    "Slc6a1" = sg_trans_Slc6a1,
                    "Slc6a11" = sg_trans_Slc6a11,
                    "Slc6a15" = sg_trans_Slc6a15,
                    "Slc7a3" = sg_trans_Slc7a3,
                    "Slc8a1" = sg_trans_Slc8a1,
                    "Slc8a2" = sg_trans_Slc8a2)

meanlista_trans <-  list("Slc17a8" = sg_meanse_trans_Slc17a8,
                         "Slc20a1" = sg_meanse_trans_Slc20a1,
                         "Slc24a3" = sg_meanse_trans_Slc24a3, 
                         "Slc32a1" = sg_meanse_trans_Slc32a1,
                         "Slc6a1" = sg_meanse_trans_Slc6a1,
                         "Slc6a11" = sg_meanse_trans_Slc6a11,
                         "Slc6a15" = sg_meanse_trans_Slc6a15,
                         "Slc7a3" = sg_meanse_trans_Slc7a3,
                         "Slc8a1" = sg_meanse_trans_Slc8a1,
                         "Slc8a2" = sg_meanse_trans_Slc8a2)

boxlista_trans <- purrr::map2(lista_trans, meanlista_trans, dplyr::inner_join, by = "names")


for (i in rownames(IC_selected_genes_trans)) {
  assign(paste0("sg_box_", i), as.data.frame(boxlista_trans[[i]]))
}

lista_box_trans <-  list("Slc17a8" = sg_box_Slc17a8,
                         "Slc20a1" = sg_box_Slc20a1,
                         "Slc24a3" = sg_box_Slc24a3, 
                         "Slc32a1" = sg_box_Slc32a1,
                         "Slc6a1" = sg_box_Slc6a1,
                         "Slc6a11" = sg_box_Slc6a11,
                         "Slc6a15" = sg_box_Slc6a15,
                         "Slc7a3" = sg_box_Slc7a3,
                         "Slc8a1" = sg_box_Slc8a1,
                         "Slc8a2" = sg_box_Slc8a2)


for (i in rownames(IC_selected_genes_trans)) {
  lista_box_trans[[i]]$names <- factor(lista_box_trans[[i]]$names,
                                       levels = c("MS330","CPU300"))
  ggplot(data = lista_box_trans[[i]],
         aes(x=names, y = points, fill=names, ymin=mean-2*se, ymax=mean+2*se))+
    stat_boxplot(geom = "errorbar",
                 width = 0.15,
                 size=0.1) +
    geom_boxplot(color= "black", show.legend = FALSE, outlier.shape = NA, coef = 0, width = 1, lwd=0.1, position = position_dodge(1))+
    
    # geom_jitter(position = position_dodge2(width = 0.3),
    #             show.legend = FALSE, shape = 21, size = 0.5, stroke = 0.1)+
    theme_classic()+labs(x = NULL, y = NULL) + theme(axis.text = element_text(colour = "black"),
                                                     axis.text.x = element_text(colour = "black"),
                                                     axis.text.y = element_text(colour = "black"),plot.margin = margin(t = 3,  
                                                                                                                       r = 1,  
                                                                                                                       b = 5,  
                                                                                                                       l = 1) ) + theme(panel.grid.major = element_line(colour = NA)) + theme(plot.title = element_text(family = "Helvetica",
                                                                                                                                                                                                                        size = 6, hjust = 0.5)) +labs(title = i)+
    scale_fill_manual(values = c("#FFD662FF","#00539CFF")) + theme(axis.ticks = element_line(colour = "black"),
                                                                   plot.background = element_rect(colour = NA)) + theme(axis.line = element_line(size = 5),
                                                                                                                        axis.title = element_text(size = 5),
                                                                                                                        axis.text = element_text(size = 5)) + theme(axis.line = element_line(size = 0.1),
                                                                                                                                                                    panel.grid.minor = element_line(colour = NA)) + theme(axis.ticks = element_line(size = 0.1),
                                                                                                                                                                                                                          panel.grid.major = element_line(size = 0.1),
                                                                                                                                                                                                                          panel.grid.minor = element_line(size = 0.1)) + theme(axis.text = element_text(angle = 90),
                                                                                                                                                                                                                                                                               plot.title = element_text(face = "italic"))+labs(title = NULL, y = NULL) + theme(axis.text = element_text(angle = 0),
                                                                                                                                                                                                                                                                                                                                                                axis.text.x = element_blank())+ scale_y_continuous(n.breaks = 10)
  ggsave(paste0("Boxplot_transporters_", i, '.pdf'), width = 15 , height = 20 , units = "mm")
}



#####Colocalization stack bar graph####

dfstack <- data.frame(bar_kolok_mean)	
dlstack <-data.frame(bar_kolok_points)	


metastack<-full_join(dfstack, dlstack, by = c('Name'='Name')) 

dfstack$Name = factor(dfstack$Name,
                      levels = c("ChAT+ TRH-", "ChAT+ TRH+" ))
ggplot(data = dfstack, 
       aes(x=Title, y = mean,
           fill=Name))+
  geom_bar(stat="identity", width = 0.7, colour="black", size=0.3)+
  labs(x=NULL, y=NULL)+
  theme_classic()+
  scale_y_continuous(expand=c(0, 0), limits = c(0, 100), n.breaks =10)+
  theme(legend.position="right", legend.direction="vertical")+
  geom_errorbar(aes(ymin=Error-SE, ymax=Error+SE), width=0.2, size = 0.3, position = "identity")+
  geom_point(data = metastack, aes(x= Title.x, y=Points), #position = position_dodge2(width = 0.2),
             show.legend = FALSE, shape = 1, size = 1, stroke = 0.5, color= "black",position = position_jitterdodge(jitter.width = 0.1,
                                                                                                                    dodge.width = 0.05))+
  scale_fill_manual(values = c("#0000FE", "forestgreen"))+
  scale_x_discrete(labels="MS") + theme(axis.ticks = element_line(colour = "black"),
                                        axis.text = element_text(colour = "black")) +labs(y = "% of cholinergic neurons")+
  theme(text = element_text(size = 5)) + theme(axis.text = element_text(size = 6),
                                               axis.text.x = element_text(colour = "black"),
                                               axis.text.y = element_text(colour = "black"))+
  theme(legend.key.size = unit(2, 'mm'))+
  theme(legend.text = element_text(size=5))+
  theme(legend.title = element_text(size=5))
ggsave("Rst12.pdf", width = 50 , height = 80 , units = "mm")


#####TRH expression bar plot#####

dftrh_bar <- data.frame(trh_meanse)	
dltrh_bar <-data.frame(trh_points)	

metatrh_bar<-full_join(dftrh_bar, dltrh_bar, by = c('Name'='Name')) 

dftrh_bar$Name = factor(dftrh_bar$Name,
                        levels = c("MS330", "CPU300" ))
ggplot(data = dftrh_bar, 
       aes(x=Name, y = mean,
           fill=Name))+
  geom_bar(stat="identity", width = 0.6, colour="black", show.legend = F)+
  labs(x=NULL, y=NULL)+
  theme_classic()+
  scale_y_continuous(expand=c(0, 0), limits = c(0, 60), n.breaks =10)+
  theme(legend.position="bottom")+
  geom_errorbar(aes(ymin=mean-SEM, ymax=mean+SEM), width=0.05, size = 0.5)+
  geom_point(data = metatrh_bar, aes(x= Name, y=Trh), #position = position_dodge2(width = 0.2),
             show.legend = FALSE, shape = 1, size = 1.2, stroke = 0.5, color= "black",position = position_jitterdodge(jitter.width = 0.1,
                                                                                                                      dodge.width = 0.05))+
  scale_fill_manual(values = c("#0000FE", "forestgreen"))+
  theme(axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black")) +labs(y = "Expression level (cpm)")+
  theme(text = element_text(size = 7)) + theme(axis.text = element_text(size = 7),
                                               axis.text.x = element_text(colour = "black"),
                                               axis.text.y = element_text(colour = "black"))
ggsave("TRH_barplot.pdf", width = 55 , height = 70 , units = "mm")


####create supplementary tables####

###Filter CPU300####

cpm_CPU300 <- Metodikai_antisense_NM_cpm_jav %>% 
  select(cpm_Geneid, cpm_CPU300_1, cpm_CPU300_3, cpm_CPU300_5, cpm_CPU300_7 ) %>% 
  unique() %>% 
  data.frame(row.names = 1)


filter4 <- apply(cpm_CPU300, 1, function(x) length(x[x>1])==4)
cpm_CPU300_4_4 <- cpm_CPU300[filter4,]

filter3 <- apply(cpm_CPU300, 1, function(x) length(x[x>1])==3)
cpm_CPU300_4_3 <- cpm_CPU300[filter3,]

filter2 <- apply(cpm_CPU300, 1, function(x) length(x[x>1])==2)
cpm_CPU300_4_2 <- cpm_CPU300[filter2,]

filter1 <- apply(cpm_CPU300, 1, function(x) length(x[x>1])==1)
cpm_CPU300_4_1 <- cpm_CPU300[filter1,]

cpm_CPU300_4_1 <- rownames_to_column(cpm_CPU300_4_1, var='Ensembl')
cpm_CPU300_4_2 <- rownames_to_column(cpm_CPU300_4_2, var='Ensembl')
cpm_CPU300_4_3 <- rownames_to_column(cpm_CPU300_4_3, var='Ensembl')
cpm_CPU300_4_4 <- rownames_to_column(cpm_CPU300_4_4, var='Ensembl')

cpm_CPU300_4_4 <- cpm_CPU300_4_4 %>% 
  mutate(Mean4 = rowMeans(.[2:5])) %>% 
  arrange(-Mean4) %>% 
  filter(!grepl ("ENSMUSG00000064339", Ensembl),
         !grepl ("ENSMUSG00000064337", Ensembl),
         !grepl ("ENSMUSG00000119584", Ensembl),
         !grepl ("ENSMUSG00000118642", Ensembl),
         !grepl ("ENSMUSG00000035202", Ensembl),
         !grepl ("ENSMUSG00000098973", Ensembl),
         !grepl ("ENSMUSG00000086324", Ensembl)) %>% 
  select(-Mean4)

cpm_CPU300_4_3 <- cpm_CPU300_4_3 %>% 
  mutate(Mean3 = rowMeans(.[2:5])) %>% 
  arrange(-Mean3)%>% 
  select(-Mean3)

cpm_CPU300_4_1 <- cpm_CPU300_4_1 %>% 
  mutate(Mean1 = rowMeans(.[2:5])) %>% 
  arrange(-Mean1)%>% 
  select(-Mean1)

cpm_CPU300_4_2 <- cpm_CPU300_4_2 %>% 
  mutate(Mean2 = rowMeans(.[2:5])) %>% 
  arrange(-Mean2) %>% 
  select(-Mean2)

cpm_CPU300_full <- cpm_CPU300_4_4 %>% 
  bind_rows(., cpm_CPU300_4_3) %>% 
  bind_rows(., cpm_CPU300_4_2) %>%
  bind_rows(., cpm_CPU300_4_1) %>% 
  unique() %>% 
  mutate(Mean = rowMeans(.[2:5]))

cpm_CPU300_full <- cpm_CPU300_full %>% 
  inner_join(Annotiationtable107, by = c("Ensembl" = "ensembl_gene_id")) %>% 
  select(Ensembl, external_gene_name, description, everything() ) %>% 
  unique()

cpm_CPU300_full_jav <- cpm_CPU300_full %>%
  mutate(round(.[, 4:8], digit = 2)) %>% 
  rename("Ensembl ID" = Ensembl) %>% 
  rename("Gene symbol" = external_gene_name) %>% 
  rename("Gene description" = description) %>% 
  rename_at(vars(contains('cpm_')), list( ~ gsub('cpm_', '', .)))


###Filter MS330####

cpm_MS330 <- Metodikai_antisense_NM_cpm_jav %>% 
  select(cpm_Geneid,cpm_MS_330_9:cpm_MS_330_12 ) %>% 
  unique() %>% 
  data.frame(row.names = 1)

filterMS4 <- apply(cpm_MS330, 1, function(x) length(x[x>1])==4)
cpm_MS330_4_4 <- cpm_MS330[filterMS4,]

filterMS3 <- apply(cpm_MS330, 1, function(x) length(x[x>1])==3)
cpm_MS330_4_3 <- cpm_MS330[filterMS3,]

filterMS2 <- apply(cpm_MS330, 1, function(x) length(x[x>1])==2)
cpm_MS330_4_2 <- cpm_MS330[filterMS2,]

filterMS1 <- apply(cpm_MS330, 1, function(x) length(x[x>1])==1)
cpm_MS330_4_1 <- cpm_MS330[filterMS1,]

cpm_MS330_4_1 <- rownames_to_column(cpm_MS330_4_1, var='Ensembl')
cpm_MS330_4_2 <- rownames_to_column(cpm_MS330_4_2, var='Ensembl')
cpm_MS330_4_3 <- rownames_to_column(cpm_MS330_4_3, var='Ensembl')
cpm_MS330_4_4 <- rownames_to_column(cpm_MS330_4_4, var='Ensembl')

cpm_MS330_4_4 <- cpm_MS330_4_4 %>% 
  mutate(Mean4 = rowMeans(.[2:5])) %>% 
  arrange(-Mean4) %>% 
  filter(!grepl ("ENSMUSG00000064339", Ensembl),
         !grepl ("ENSMUSG00000064337", Ensembl),
         !grepl ("ENSMUSG00000119584", Ensembl),
         !grepl ("ENSMUSG00000118642", Ensembl),
         !grepl ("ENSMUSG00000035202", Ensembl),
         !grepl ("ENSMUSG00000098973", Ensembl),
         !grepl ("ENSMUSG00000086324", Ensembl)) %>% 
  select(-Mean4)

cpm_MS330_4_3 <- cpm_MS330_4_3 %>% 
  mutate(Mean3 = rowMeans(.[2:5])) %>% 
  arrange(-Mean3)%>% 
  select(-Mean3)

cpm_MS330_4_1 <- cpm_MS330_4_1 %>% 
  mutate(Mean1 = rowMeans(.[2:5])) %>% 
  arrange(-Mean1)%>% 
  select(-Mean1)

cpm_MS330_4_2 <- cpm_MS330_4_2 %>% 
  mutate(Mean2 = rowMeans(.[2:5])) %>% 
  arrange(-Mean2) %>% 
  select(-Mean2)


cpm_MS330_full <- cpm_MS330_4_4 %>% 
  bind_rows(., cpm_MS330_4_3) %>% 
  bind_rows(., cpm_MS330_4_2) %>%
  bind_rows(., cpm_MS330_4_1) %>% 
  unique() %>% 
  mutate(Mean = rowMeans(.[2:5]))

cpm_MS330_full <- cpm_MS330_full %>% 
  inner_join(Annotiationtable107, by = c("Ensembl" = "ensembl_gene_id")) %>% 
  select(Ensembl, external_gene_name, description, everything() ) %>% 
  unique()

cpm_MS330_full_jav <- cpm_MS330_full %>%
  mutate(round(.[, 4:8], digit = 2)) %>% 
  rename("Ensembl ID" = Ensembl) %>% 
  rename("Gene symbol" = external_gene_name) %>% 
  rename("Gene description" = description )%>% 
  rename_at(vars(contains('cpm_')), list( ~ gsub('cpm_', '', .)))

#####Transcription factors of CPU300####

CPU300 <- unique(CPU300)

CPU300_tf_list <- TF_list %>% 
  inner_join(cpm_CPU300, by = c("ensembl_gene_id" = "cpm_Geneid")) %>% 
  unique() %>% 
  data.frame(row.names = 1)


filter4 <- apply(CPU300_tf_list, 1, function(x) length(x[x>1])==4)
CPU300_tf_list_4_4 <- CPU300_tf_list[filter4,]

filter3 <- apply(CPU300_tf_list, 1, function(x) length(x[x>1])==3)
CPU300_tf_list_4_3 <- CPU300_tf_list[filter3,]

filter2 <- apply(CPU300_tf_list, 1, function(x) length(x[x>1])==2)
CPU300_tf_list_4_2 <- CPU300_tf_list[filter2,]

filter1 <- apply(CPU300_tf_list, 1, function(x) length(x[x>1])==1)
CPU300_tf_list_4_1 <- CPU300_tf_list[filter1,]

CPU300_tf_list_4_1 <- rownames_to_column(CPU300_tf_list_4_1, var='Ensembl')
CPU300_tf_list_4_2 <- rownames_to_column(CPU300_tf_list_4_2, var='Ensembl')
CPU300_tf_list_4_3 <- rownames_to_column(CPU300_tf_list_4_3, var='Ensembl')
CPU300_tf_list_4_4 <- rownames_to_column(CPU300_tf_list_4_4, var='Ensembl')

CPU300_tf_list_4_4 <- CPU300_tf_list_4_4 %>% 
  mutate(Mean4 = rowMeans(.[2:5])) %>% 
  arrange(-Mean4) %>% 
  select(-Mean4)

CPU300_tf_list_4_3 <- CPU300_tf_list_4_3 %>% 
  mutate(Mean3 = rowMeans(.[2:5])) %>% 
  arrange(-Mean3)%>% 
  select(-Mean3)

CPU300_tf_list_4_1 <- CPU300_tf_list_4_1 %>% 
  mutate(Mean1 = rowMeans(.[2:5])) %>% 
  arrange(-Mean1)%>% 
  select(-Mean1)

CPU300_tf_list_4_2 <- CPU300_tf_list_4_2 %>% 
  mutate(Mean2 = rowMeans(.[2:5])) %>% 
  arrange(-Mean2) %>% 
  select(-Mean2)


CPU300_tf_list_full <- CPU300_tf_list_4_4 %>% 
  bind_rows(., CPU300_tf_list_4_3) %>% 
  bind_rows(., CPU300_tf_list_4_2) %>%
  bind_rows(., CPU300_tf_list_4_1) %>% 
  unique() %>% 
  mutate(Mean = rowMeans(.[2:5]))

CPU300_tf_list_full <- CPU300_tf_list_full %>% 
  inner_join(Annotiationtable107, by = c("Ensembl" = "ensembl_gene_id")) %>% 
  select(Ensembl, 'external_gene_name', "description", everything() ) %>% 
  unique()

CPU300_tf_list_full_jav <- CPU300_tf_list_full %>%
  mutate(round(.[, 4:8], digit = 2))


#####receptors of CPU300#####


Annotation_table107 = getBM(attributes = c("ensembl_gene_id","external_gene_name", "description"), 
                            filters = c('ensembl_gene_id') ,
                            values = Annotiation_table104$gene_id,
                            mart = ensembl107)

olfactory_rec <- Annotiationtable107 %>% 
  filter(grepl ("olfactory receptor",  description),
         !grepl ("pseudogene", description)) %>% 
  inner_join(CPU300, by = c("ensembl_gene_id" = "Ensembl"))

olfactory_rec <- olfactory_rec %>% 
  dplyr::select(ensembl_gene_id, `CHAT300-1`, `CHAT300-3`, `CHAT300-5`, `CHAT300-7`) %>% 
  data.frame(row.names = 1)




filter_receptors1 <- apply(olfactory_rec, 1, function(x) length(x[x>1])==1)
olfactory_rec1 <- olfactory_rec[filter_receptors1,]

receptors <- All_receptors %>% 
  na.omit() %>% 
  left_join(Annotiationtable107, by = c("symbol" = "external_gene_name")) %>% 
  select(ensembl_gene_id) %>%
  na.omit()

CPU300receptor_list <- receptors %>% 
  inner_join(cpm_CPU300, by = c("ensembl_gene_id" = "cpm_Geneid")) %>% 
  unique() %>% 
  data.frame(row.names = 1)


filter_receptors4 <- apply(CPU300receptor_list, 1, function(x) length(x[x>1])==4)
CPU300receptor_list_4_4 <- CPU300receptor_list[filter_receptors4,]

filter_receptors3 <- apply(CPU300receptor_list, 1, function(x) length(x[x>1])==3)
CPU300receptor_list_4_3 <- CPU300receptor_list[filter_receptors3,]

filter_receptors2 <- apply(CPU300receptor_list, 1, function(x) length(x[x>1])==2)
CPU300receptor_list_4_2 <- CPU300receptor_list[filter_receptors2,]

filter_receptors1 <- apply(CPU300receptor_list, 1, function(x) length(x[x>1])==1)
CPU300receptor_list_4_1 <- CPU300receptor_list[filter_receptors1,]

CPU300receptor_list_4_1 <- rownames_to_column(CPU300receptor_list_4_1, var='Ensembl')
CPU300receptor_list_4_2 <- rownames_to_column(CPU300receptor_list_4_2, var='Ensembl')
CPU300receptor_list_4_3 <- rownames_to_column(CPU300receptor_list_4_3, var='Ensembl')
CPU300receptor_list_4_4 <- rownames_to_column(CPU300receptor_list_4_4, var='Ensembl')

CPU300receptor_list_4_4 <- CPU300receptor_list_4_4 %>% 
  mutate(Mean4 = rowMeans(.[2:5])) %>% 
  arrange(-Mean4) %>% 
  select(-Mean4)

CPU300receptor_list_4_3 <- CPU300receptor_list_4_3 %>% 
  mutate(Mean3 = rowMeans(.[2:5])) %>% 
  arrange(-Mean3)%>% 
  select(-Mean3)

CPU300receptor_list_4_1 <- CPU300receptor_list_4_1 %>% 
  mutate(Mean1 = rowMeans(.[2:5])) %>% 
  arrange(-Mean1)%>% 
  select(-Mean1)

CPU300receptor_list_4_2 <- CPU300receptor_list_4_2 %>% 
  mutate(Mean2 = rowMeans(.[2:5])) %>% 
  arrange(-Mean2) %>% 
  select(-Mean2)

CPU300receptor_list_full <- CPU300receptor_list_4_4 %>% 
  bind_rows(., CPU300receptor_list_4_3) %>% 
  bind_rows(., CPU300receptor_list_4_2) %>%
  bind_rows(., CPU300receptor_list_4_1) %>% 
  unique() %>% 
  mutate(Mean = rowMeans(.[2:5]))

CPU300receptor_list_full <- CPU300receptor_list_full %>% 
  inner_join(Annotiationtable107, by = c("Ensembl" = "ensembl_gene_id")) %>% 
  select(Ensembl, 'external_gene_name', "description", everything() ) %>% 
  unique()

CPU300receptor_list_full_jav <- CPU300receptor_list_full %>%
  mutate(round(.[, 4:8], digit = 2))

######Transporters of CPU300#####

transporter <- transporter %>% 
  na.omit() %>% 
  inner_join(Annotiationtable107, by = c("Gene symbol" = "external_gene_name")) %>%
  select(ensembl_gene_id)

CPU300transporter_list <- transporter %>% 
  inner_join(cpm_CPU300, by = c("ensembl_gene_id" = "cpm_Geneid")) %>% 
  unique() %>% 
  data.frame(row.names = 1)

filter_transporter4 <- apply(CPU300transporter_list, 1, function(x) length(x[x>1])==4)
CPU300transporter_list_4_4 <- CPU300transporter_list[filter_transporter4,]

filter_transporter3 <- apply(CPU300transporter_list, 1, function(x) length(x[x>1])==3)
CPU300transporter_list_4_3 <- CPU300transporter_list[filter_transporter3,]

filter_transporter2 <- apply(CPU300transporter_list, 1, function(x) length(x[x>1])==2)
CPU300transporter_list_4_2 <- CPU300transporter_list[filter_transporter2,]

filter_transporter1 <- apply(CPU300transporter_list, 1, function(x) length(x[x>1])==1)
CPU300transporter_list_4_1 <- CPU300transporter_list[filter_transporter1,]

CPU300transporter_list_4_1 <- rownames_to_column(CPU300transporter_list_4_1, var='Ensembl')
CPU300transporter_list_4_2 <- rownames_to_column(CPU300transporter_list_4_2, var='Ensembl')
CPU300transporter_list_4_3 <- rownames_to_column(CPU300transporter_list_4_3, var='Ensembl')
CPU300transporter_list_4_4 <- rownames_to_column(CPU300transporter_list_4_4, var='Ensembl')

CPU300transporter_list_4_4 <- CPU300transporter_list_4_4 %>% 
  mutate(Mean4 = rowMeans(.[2:5])) %>% 
  arrange(-Mean4) %>% 
  select(-Mean4)

CPU300transporter_list_4_3 <- CPU300transporter_list_4_3 %>% 
  mutate(Mean3 = rowMeans(.[2:5])) %>% 
  arrange(-Mean3)%>% 
  select(-Mean3)

CPU300transporter_list_4_1 <- CPU300transporter_list_4_1 %>% 
  mutate(Mean1 = rowMeans(.[2:5])) %>% 
  arrange(-Mean1)%>% 
  select(-Mean1)

CPU300transporter_list_4_2 <- CPU300transporter_list_4_2 %>% 
  mutate(Mean2 = rowMeans(.[2:5])) %>% 
  arrange(-Mean2) %>% 
  select(-Mean2)

CPU300transporter_list_full <- CPU300transporter_list_4_4 %>% 
  bind_rows(., CPU300transporter_list_4_3) %>% 
  bind_rows(., CPU300transporter_list_4_2) %>%
  bind_rows(., CPU300transporter_list_4_1) %>% 
  unique() %>% 
  mutate(Mean = rowMeans(.[2:5]))

CPU300transporter_list_full <- CPU300transporter_list_full %>% 
  inner_join(Annotiationtable107, by = c("Ensembl" = "ensembl_gene_id")) %>% 
  select(Ensembl, 'external_gene_name', "description", everything() ) %>% 
  unique()

CPU300transporter_list_full_jav <- CPU300transporter_list_full %>%
  mutate(round(.[, 4:8], digit = 2))

######Ion-channels of CPU300#####

Ion_channel <- Ion_channel %>% 
  na.omit() %>% 
  inner_join(Annotiationtable107, by = c("Gene symbol" = "external_gene_name")) %>%
  select(ensembl_gene_id)

CPU300ion_channel_list <- Ion_channel %>% 
  inner_join(cpm_CPU300, by = c("ensembl_gene_id" = "cpm_Geneid")) %>% 
  unique() %>% 
  data.frame(row.names = 1)


filter_ion_channel4 <- apply(CPU300ion_channel_list, 1, function(x) length(x[x>1])==4)
CPU300ion_channel_list_4_4 <- CPU300ion_channel_list[filter_ion_channel4,]

filter_ion_channel3 <- apply(CPU300ion_channel_list, 1, function(x) length(x[x>1])==3)
CPU300ion_channel_list_4_3 <- CPU300ion_channel_list[filter_ion_channel3,]

filter_ion_channel2 <- apply(CPU300ion_channel_list, 1, function(x) length(x[x>1])==2)
CPU300ion_channel_list_4_2 <- CPU300ion_channel_list[filter_ion_channel2,]

filter_ion_channel1 <- apply(CPU300ion_channel_list, 1, function(x) length(x[x>1])==1)
CPU300ion_channel_list_4_1 <- CPU300ion_channel_list[filter_ion_channel1,]

CPU300ion_channel_list_4_1 <- rownames_to_column(CPU300ion_channel_list_4_1, var='Ensembl')
CPU300ion_channel_list_4_2 <- rownames_to_column(CPU300ion_channel_list_4_2, var='Ensembl')
CPU300ion_channel_list_4_3 <- rownames_to_column(CPU300ion_channel_list_4_3, var='Ensembl')
CPU300ion_channel_list_4_4 <- rownames_to_column(CPU300ion_channel_list_4_4, var='Ensembl')

CPU300ion_channel_list_4_4 <- CPU300ion_channel_list_4_4 %>% 
  mutate(Mean4 = rowMeans(.[2:5])) %>% 
  arrange(-Mean4) %>% 
  select(-Mean4)

CPU300ion_channel_list_4_3 <- CPU300ion_channel_list_4_3 %>% 
  mutate(Mean3 = rowMeans(.[2:5])) %>% 
  arrange(-Mean3)%>% 
  select(-Mean3)

CPU300ion_channel_list_4_1 <- CPU300ion_channel_list_4_1 %>% 
  mutate(Mean1 = rowMeans(.[2:5])) %>% 
  arrange(-Mean1)%>% 
  select(-Mean1)

CPU300ion_channel_list_4_2 <- CPU300ion_channel_list_4_2 %>% 
  mutate(Mean2 = rowMeans(.[2:5])) %>% 
  arrange(-Mean2) %>% 
  select(-Mean2)


CPU300ion_channel_list_full <- CPU300ion_channel_list_4_4 %>% 
  bind_rows(., CPU300ion_channel_list_4_3) %>% 
  bind_rows(., CPU300ion_channel_list_4_2) %>%
  bind_rows(., CPU300ion_channel_list_4_1) %>% 
  unique() %>% 
  mutate(Mean = rowMeans(.[2:5]))

CPU300ion_channel_list_full <- CPU300ion_channel_list_full %>% 
  inner_join(Annotiationtable107, by = c("Ensembl" = "ensembl_gene_id")) %>% 
  select(Ensembl, 'external_gene_name', "description", everything() ) %>% 
  unique()

CPU300ion_channel_list_full_jav <- CPU300ion_channel_list_full %>%
  mutate(round(.[, 4:8], digit = 2))

#####Transcription factors of MS330####

MS330 <- unique(MS330)

MS330_tf_list <- TF_list %>% 
  inner_join(cpm_MS330, by = c("ensembl_gene_id" = "cpm_Geneid")) %>% 
  unique() %>% 
  data.frame(row.names = 1)


filter_tf_list_4 <- apply(MS330_tf_list, 1, function(x) length(x[x>1])==4)
MS330_tf_list_4_4 <- MS330_tf_list[filter_tf_list_4,]

filtertf_list_3 <- apply(MS330_tf_list, 1, function(x) length(x[x>1])==3)
MS330_tf_list_4_3 <- MS330_tf_list[filtertf_list_3,]

filtertf_list_2 <- apply(MS330_tf_list, 1, function(x) length(x[x>1])==2)
MS330_tf_list_4_2 <- MS330_tf_list[filtertf_list_2,]

filtertf_list_1 <- apply(MS330_tf_list, 1, function(x) length(x[x>1])==1)
MS330_tf_list_4_1 <- MS330_tf_list[filtertf_list_1,]

MS330_tf_list_4_1 <- rownames_to_column(MS330_tf_list_4_1, var='Ensembl')
MS330_tf_list_4_2 <- rownames_to_column(MS330_tf_list_4_2, var='Ensembl')
MS330_tf_list_4_3 <- rownames_to_column(MS330_tf_list_4_3, var='Ensembl')
MS330_tf_list_4_4 <- rownames_to_column(MS330_tf_list_4_4, var='Ensembl')

MS330_tf_list_4_4 <- MS330_tf_list_4_4 %>% 
  mutate(Mean4 = rowMeans(.[2:5])) %>% 
  arrange(-Mean4) %>% 
  select(-Mean4)

MS330_tf_list_4_3 <- MS330_tf_list_4_3 %>% 
  mutate(Mean3 = rowMeans(.[2:5])) %>% 
  arrange(-Mean3)%>% 
  select(-Mean3)

MS330_tf_list_4_1 <- MS330_tf_list_4_1 %>% 
  mutate(Mean1 = rowMeans(.[2:5])) %>% 
  arrange(-Mean1)%>% 
  select(-Mean1)

MS330_tf_list_4_2 <- MS330_tf_list_4_2 %>% 
  mutate(Mean2 = rowMeans(.[2:5])) %>% 
  arrange(-Mean2) %>% 
  select(-Mean2)


MS330_tf_list_full <- MS330_tf_list_4_4 %>% 
  bind_rows(., MS330_tf_list_4_3) %>% 
  bind_rows(., MS330_tf_list_4_2) %>%
  bind_rows(., MS330_tf_list_4_1) %>% 
  unique() %>% 
  mutate(Mean = rowMeans(.[2:5]))

MS330_tf_list_full <- MS330_tf_list_full %>% 
  inner_join(Annotiationtable107, by = c("Ensembl" = "ensembl_gene_id")) %>% 
  select(Ensembl, 'external_gene_name', "description", everything() ) %>% 
  unique()

MS330_tf_list_full_jav <- MS330_tf_list_full %>%
  mutate(round(.[, 4:8], digit = 2))


#####Receptors of MS330#####

Annotation_table104 = getBM(attributes = c("ensembl_gene_id","external_gene_name", "description"), 
                            filters = c('ensembl_gene_id') ,
                            values = Annotiation_table104$gene_id,
                            mart = ensembl104)


olfactory_rec <- Annotation_table104 %>% 
  filter(grepl ("olfactory receptor",  description),
         !grepl ("pseudogene", description)) %>% 
  inner_join(MS330, by = c("ensembl_gene_id" = "Ensembl"))

olfactory_rec <- olfactory_rec %>% 
  dplyr::select(ensembl_gene_id, `CHAT300-1`, `CHAT300-3`, `CHAT300-5`, `CHAT300-7`) %>% 
  data.frame(row.names = 1)




filter_receptors1 <- apply(olfactory_rec, 1, function(x) length(x[x>1])==1)
olfactory_rec1 <- olfactory_rec[filter_receptors1,]

receptors <- All_receptors %>% 
  na.omit() %>% 
  left_join(Annotiation_table104, by = c("symbol" = "gene_name")) %>% 
  select(gene_id) %>%
  na.omit()

MS330receptor_list <- receptors %>% 
  inner_join(cpm_MS330, by = c("ensembl_gene_id" = "cpm_Geneid")) %>% 
  unique() %>% 
  data.frame(row.names = 1)


filter_receptors4 <- apply(MS330receptor_list, 1, function(x) length(x[x>1])==4)
MS330receptor_list_4_4 <- MS330receptor_list[filter_receptors4,]

filter_receptors3 <- apply(MS330receptor_list, 1, function(x) length(x[x>1])==3)
MS330receptor_list_4_3 <- MS330receptor_list[filter_receptors3,]

filter_receptors2 <- apply(MS330receptor_list, 1, function(x) length(x[x>1])==2)
MS330receptor_list_4_2 <- MS330receptor_list[filter_receptors2,]

filter_receptors1 <- apply(MS330receptor_list, 1, function(x) length(x[x>1])==1)
MS330receptor_list_4_1 <- MS330receptor_list[filter_receptors1,]

MS330receptor_list_4_1 <- rownames_to_column(MS330receptor_list_4_1, var='Ensembl')
MS330receptor_list_4_2 <- rownames_to_column(MS330receptor_list_4_2, var='Ensembl')
MS330receptor_list_4_3 <- rownames_to_column(MS330receptor_list_4_3, var='Ensembl')
MS330receptor_list_4_4 <- rownames_to_column(MS330receptor_list_4_4, var='Ensembl')

MS330receptor_list_4_4 <- MS330receptor_list_4_4 %>% 
  mutate(Mean4 = rowMeans(.[2:5])) %>% 
  arrange(-Mean4) %>% 
  select(-Mean4)

MS330receptor_list_4_3 <- MS330receptor_list_4_3 %>% 
  mutate(Mean3 = rowMeans(.[2:5])) %>% 
  arrange(-Mean3)%>% 
  select(-Mean3)

MS330receptor_list_4_1 <- MS330receptor_list_4_1 %>% 
  mutate(Mean1 = rowMeans(.[2:5])) %>% 
  arrange(-Mean1)%>% 
  select(-Mean1)

MS330receptor_list_4_2 <- MS330receptor_list_4_2 %>% 
  mutate(Mean2 = rowMeans(.[2:5])) %>% 
  arrange(-Mean2) %>% 
  select(-Mean2)

MS330receptor_list_full <- MS330receptor_list_4_4 %>% 
  bind_rows(., MS330receptor_list_4_3) %>% 
  bind_rows(., MS330receptor_list_4_2) %>%
  bind_rows(., MS330receptor_list_4_1) %>% 
  unique() %>% 
  mutate(Mean = rowMeans(.[2:5]))

MS330receptor_list_full <- MS330receptor_list_full %>% 
  inner_join(Annotiationtable107, by = c("Ensembl" = "ensembl_gene_id")) %>% 
  select(Ensembl, 'external_gene_name', "description", everything() ) %>% 
  unique()

MS330receptor_list_full_jav <- MS330receptor_list_full %>%
  mutate(round(.[, 4:8], digit = 2))

######Transporters of MS330#####

transporter <- transporter %>% 
  na.omit() %>% 
  inner_join(Annotiation_table104, by = c("Gene symbol" = "gene_name")) %>%
  select(gene_id)

MS330transporter_list <- transporter %>% 
  inner_join(cpm_MS330, by = c("ensembl_gene_id" = "cpm_Geneid")) %>% 
  unique() %>% 
  data.frame(row.names = 1)


filter_transporter4 <- apply(MS330transporter_list, 1, function(x) length(x[x>1])==4)
MS330transporter_list_4_4 <- MS330transporter_list[filter_transporter4,]

filter_transporter3 <- apply(MS330transporter_list, 1, function(x) length(x[x>1])==3)
MS330transporter_list_4_3 <- MS330transporter_list[filter_transporter3,]

filter_transporter2 <- apply(MS330transporter_list, 1, function(x) length(x[x>1])==2)
MS330transporter_list_4_2 <- MS330transporter_list[filter_transporter2,]

filter_transporter1 <- apply(MS330transporter_list, 1, function(x) length(x[x>1])==1)
MS330transporter_list_4_1 <- MS330transporter_list[filter_transporter1,]

MS330transporter_list_4_1 <- rownames_to_column(MS330transporter_list_4_1, var='Ensembl')
MS330transporter_list_4_2 <- rownames_to_column(MS330transporter_list_4_2, var='Ensembl')
MS330transporter_list_4_3 <- rownames_to_column(MS330transporter_list_4_3, var='Ensembl')
MS330transporter_list_4_4 <- rownames_to_column(MS330transporter_list_4_4, var='Ensembl')

MS330transporter_list_4_4 <- MS330transporter_list_4_4 %>% 
  mutate(Mean4 = rowMeans(.[2:5])) %>% 
  arrange(-Mean4) %>% 
  select(-Mean4)

MS330transporter_list_4_3 <- MS330transporter_list_4_3 %>% 
  mutate(Mean3 = rowMeans(.[2:5])) %>% 
  arrange(-Mean3)%>% 
  select(-Mean3)

MS330transporter_list_4_1 <- MS330transporter_list_4_1 %>% 
  mutate(Mean1 = rowMeans(.[2:5])) %>% 
  arrange(-Mean1)%>% 
  select(-Mean1)

MS330transporter_list_4_2 <- MS330transporter_list_4_2 %>% 
  mutate(Mean2 = rowMeans(.[2:5])) %>% 
  arrange(-Mean2) %>% 
  select(-Mean2)

MS330transporter_list_full <- MS330transporter_list_4_4 %>% 
  bind_rows(., MS330transporter_list_4_3) %>% 
  bind_rows(., MS330transporter_list_4_2) %>%
  bind_rows(., MS330transporter_list_4_1) %>% 
  unique() %>% 
  mutate(Mean = rowMeans(.[2:5]))

MS330transporter_list_full <- MS330transporter_list_full %>% 
  inner_join(Annotiationtable107, by = c("Ensembl" = "ensembl_gene_id")) %>% 
  select(Ensembl, 'external_gene_name', "description", everything() ) %>% 
  unique()

MS330transporter_list_full_jav <- MS330transporter_list_full %>%
  mutate(round(.[, 4:8], digit = 2))

######Ion-channels of MS330#####

Ion_channel <- Ion_channel1 %>% 
  na.omit() %>% 
  inner_join(Annotiation_table104, by = c("Gene symbol" = "gene_name")) %>%
  select(gene_id)

MS330ion_channel_list <- Ion_channel %>% 
  inner_join(cpm_MS330, by = c("ensembl_gene_id" = "cpm_Geneid")) %>% 
  unique() %>% 
  data.frame(row.names = 1)


filter_ion_channel4 <- apply(MS330ion_channel_list, 1, function(x) length(x[x>1])==4)
MS330ion_channel_list_4_4 <- MS330ion_channel_list[filter_ion_channel4,]

filter_ion_channel3 <- apply(MS330ion_channel_list, 1, function(x) length(x[x>1])==3)
MS330ion_channel_list_4_3 <- MS330ion_channel_list[filter_ion_channel3,]

filter_ion_channel2 <- apply(MS330ion_channel_list, 1, function(x) length(x[x>1])==2)
MS330ion_channel_list_4_2 <- MS330ion_channel_list[filter_ion_channel2,]

filter_ion_channel1 <- apply(MS330ion_channel_list, 1, function(x) length(x[x>1])==1)
MS330ion_channel_list_4_1 <- MS330ion_channel_list[filter_ion_channel1,]

MS330ion_channel_list_4_1 <- rownames_to_column(MS330ion_channel_list_4_1, var='Ensembl')
MS330ion_channel_list_4_2 <- rownames_to_column(MS330ion_channel_list_4_2, var='Ensembl')
MS330ion_channel_list_4_3 <- rownames_to_column(MS330ion_channel_list_4_3, var='Ensembl')
MS330ion_channel_list_4_4 <- rownames_to_column(MS330ion_channel_list_4_4, var='Ensembl')

MS330ion_channel_list_4_4 <- MS330ion_channel_list_4_4 %>% 
  mutate(Mean4 = rowMeans(.[2:5])) %>% 
  arrange(-Mean4) %>% 
  select(-Mean4)

MS330ion_channel_list_4_3 <- MS330ion_channel_list_4_3 %>% 
  mutate(Mean3 = rowMeans(.[2:5])) %>% 
  arrange(-Mean3)%>% 
  select(-Mean3)

MS330ion_channel_list_4_1 <- MS330ion_channel_list_4_1 %>% 
  mutate(Mean1 = rowMeans(.[2:5])) %>% 
  arrange(-Mean1)%>% 
  select(-Mean1)

MS330ion_channel_list_4_2 <- MS330ion_channel_list_4_2 %>% 
  mutate(Mean2 = rowMeans(.[2:5])) %>% 
  arrange(-Mean2) %>% 
  select(-Mean2)

MS330ion_channel_list_full <- MS330ion_channel_list_4_4 %>% 
  bind_rows(., MS330ion_channel_list_4_3) %>% 
  bind_rows(., MS330ion_channel_list_4_2) %>%
  bind_rows(., MS330ion_channel_list_4_1) %>% 
  unique() %>% 
  mutate(Mean = rowMeans(.[2:5]))

MS330ion_channel_list_full <- MS330ion_channel_list_full %>% 
  inner_join(Annotiationtable107, by = c("Ensembl" = "ensembl_gene_id")) %>% 
  select(Ensembl, 'external_gene_name', "description", everything() ) %>% 
  unique()

MS330ion_channel_list_full_jav <- MS330ion_channel_list_full %>%
  mutate(round(.[, 4:8], digit = 2))

for (i in names(sheets1)) {
  
  x <- sheets1[[i]] %>% 
    rename("Ensembl ID" = "Ensembl") %>%
    rename_at(vars(contains('cpm_')), list( ~ gsub('cpm_', '', .))) %>% 
    rename("Gene symbol" = "external_gene_name") %>% 
    rename("Gene description" = "description")
  assign(i,x)  
}