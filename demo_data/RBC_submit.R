# chenghonghan
library(openxlsx)
library(ggplot2)
library(ggsci)
library(ggsignif)
library(reshape2)
library(dplyr)

##### matrix ####
sampleinfo = read.xlsx("Sample_information20241118.xlsx")
rownames(sampleinfo ) = sampleinfo$Raw_name
sampleinfo_RBC = subset(sampleinfo , grepl('RBC', Expriment_ID))

datpro_np23 = read.csv(paste0("RBC_PRP_different_ratio_evaluation_NP23_21samples20241118.pg_matrix.tsv"), sep = '\t', row.names = 1)
rownames(datpro_np23) = paste0(datpro_np23$Protein.Names, '_', datpro_np23$Genes)
datpro_np23 = datpro_np23[5:ncol(datpro_np23)]
dim(datpro_np23)
colnames(datpro_np23)
colnames(datpro_np23) = lapply(colnames(datpro_np23), function(x){
  strsplit(x, '.', fixed = T)[[1]][8]
})
datpre_np23 = read.csv(paste0("RBC_PRP_different_ratio_evaluation_NP23_21samples20241118.pr_matrix.tsv"), sep = '\t' )
colnames(datpre_np23) = lapply(colnames(datpre_np23), function(x){
  x = gsub('D..gaohuanhuan.WAN20241117gaohh_OmniProtV2_24minDIA_platelet.RBC.PPP.NP23.', '',x,  fixed = T)
  x = gsub('.raw', '', x, fixed = T)
})

datpro_neat = read.csv(paste0("RBC_PRP_different_ratio_evaluation_neat_21samples20241118.pg_matrix.tsv"), sep = '\t', row.names = 1)
rownames(datpro_neat) = paste0(datpro_neat$Protein.Names, '_', datpro_neat$Genes)
datpro_neat = datpro_neat[5:ncol(datpro_neat)]
dim(datpro_neat)
colnames(datpro_neat)
colnames(datpro_neat) = lapply(colnames(datpro_neat), function(x){
  strsplit(x, '.', fixed = T)[[1]][8]
})
datpre_neat = read.csv(paste0("RBC_PRP_different_ratio_evaluation_neat_21samples20241118.pr_matrix.tsv"), sep = '\t' )
colnames(datpre_neat) = lapply(colnames(datpre_neat), function(x){
  x = gsub('D..gaohuanhuan.WAN20241117gaohh_OmniProtV2_24minDIA_platelet.RBC.PPP.Neat.', '',x,  fixed = T)
  x = gsub('.raw', '', x, fixed = T)
})
##1 
'HUMAN'
pronum = colSums(!is.na(datpro_np23[grepl('BOVIN', rownames(datpro_np23)), ]))
pronum = data.frame(np23_pro_BOVIN = pronum, row.names = colnames(datpro_np23))
pronum$np23_pro_HUMAN = colSums(!is.na(datpro_np23[grepl('HUMAN', rownames(datpro_np23)), ]))
for(i in colnames(datpro_np23)){
  temp = datpre_np23[c("Protein.Names", "Stripped.Sequence", i)]
  colnames(temp)[3]='cur'
  temp = subset(temp, !is.na(cur) & grepl('BOVIN', Protein.Names))
  #temp = temp %>% distinct(Stripped.Sequence, .keep_all = TRUE)
  pronum[i, 'np23_pre_BOVIN'] = length(unique(temp$Stripped.Sequence))
  
  temp = datpre_np23[c("Protein.Names", "Stripped.Sequence", i)]
  colnames(temp)[3]='cur'
  temp = subset(temp, !is.na(cur) & grepl('HUMAN', Protein.Names))
  pronum[i, 'np23_pre_HUMAN'] = length(unique(temp$Stripped.Sequence))
}
pronum$sample = rownames(pronum)
pronum = melt(pronum, id.vars = 'sample')
pronum$Plasma_type = sampleinfo_RBC[pronum$sample, 'Plasma_type']
pronum$Plasma_ratio = sampleinfo_RBC[pronum$sample, 'Plasma_ratio']
pronum$vol = as.character(sampleinfo_RBC[pronum$sample, 'Vol.%'])
identifyNum = pronum

pronum = colSums(!is.na(datpro_neat[grepl('BOVIN', rownames(datpro_neat)), ]))
pronum = data.frame(neat_pro_BOVIN = pronum, row.names = colnames(datpro_neat))
pronum$neat_pro_HUMAN = colSums(!is.na(datpro_neat[grepl('HUMAN', rownames(datpro_neat)), ]))
for(i in colnames(datpro_neat)){
  temp = datpre_neat[c("Protein.Names", "Stripped.Sequence", i)]
  colnames(temp)[3]='cur'
  temp = subset(temp, !is.na(cur) & grepl('BOVIN', Protein.Names))
  #temp = temp %>% distinct(Stripped.Sequence, .keep_all = TRUE)
  pronum[i, 'neat_pre_BOVIN'] = length(unique(temp$Stripped.Sequence))
  
  temp = datpre_neat[c("Protein.Names", "Stripped.Sequence", i)]
  colnames(temp)[3]='cur'
  temp = subset(temp, !is.na(cur) & grepl('HUMAN', Protein.Names))
  pronum[i, 'neat_pre_HUMAN'] = length(unique(temp$Stripped.Sequence))
}
pronum$sample = rownames(pronum)
pronum = melt(pronum, id.vars = 'sample')
pronum$Plasma_type = sampleinfo_RBC[pronum$sample, 'Plasma_type']
pronum$Plasma_ratio = sampleinfo_RBC[pronum$sample, 'Plasma_ratio']
pronum$vol = as.character(sampleinfo_RBC[pronum$sample, 'Vol.%'])
identifyNum = rbind(identifyNum, pronum)
head(identifyNum)

identifyNum = data.frame(identifyNum)
identifyNum$variable = as.character(identifyNum$variable)
identifyNum$type = as.character(lapply(identifyNum$variable, function(x){ strsplit(x, '_')[[1]][1] } ))
identifyNum$p = as.character(lapply(identifyNum$variable, function(x)strsplit(x, '_')[[1]][2]))
identifyNum$species = as.character(lapply(identifyNum$variable, function(x)strsplit(x, '_')[[1]][3]))

df1.2 = identifyNum
df1.2$x = as.numeric(factor((df1.2$vol)))
unique(df1.2$x)
unique(df1.2$vol)
head(df1.2)

df1.2 = df1.2[c(2,3, 7, 8, 9, 10)]
n = nrow(df1.2)
type='np23'
df1.2[n +1, ] = c(paste0(type,'_pre_HUMAN'), 29019, type,'pre', 'HUMAN' ,7)
df1.2[n +2, ] = c(paste0(type,'_pro_HUMAN'), 3597, type,'pro', 'HUMAN',7 )
df1.2[n +3, ] = c( paste0(type,'_pre_BOVIN'), 0, type,'pre', 'BOVIN' ,7)
df1.2[n +4, ] = c(paste0(type,'_pro_BOVIN'), 0, type,'pro', 'BOVIN' ,7)
type='neat'
n = nrow(df1.2)
df1.2[n+1, ] = c(paste0(type,'_pre_HUMAN'), 11660, type,'pre', 'HUMAN' , 7)
df1.2[n+2, ] = c(paste0(type,'_pro_HUMAN'), 1666, type,'pro', 'HUMAN' , 7)
df1.2[n+3, ] = c(paste0(type,'_pre_BOVIN'), 0, type,'pre', 'BOVIN' , 7)
df1.2[n+4, ] = c(paste0(type,'_pro_BOVIN'), 0, type,'pro', 'BOVIN', 7)
df1.2 = df1.2[df1.2$species== 'HUMAN' , ]
df1.2$x = as.numeric(df1.2$x)
df1.2$value = as.numeric(df1.2$value)
ggplot(df1.2[grepl('pre', df1.2$variable) & df1.2$species == 'HUMAN' , ], aes(x = x, y = value, fill = type)) +    geom_bar(stat = "summary", fun ="mean", position = position_dodge(), alpha=0.7) +
  geom_point()+
  stat_summary(fun.data = 'mean_se', geom = "errorbar", colour = "black", width = 0.3, position = position_dodge( .9))+
  geom_text(stat = "summary", fun = "mean", aes(label = round(..y.., 0)), color = "black", size = 3.5)+
  scale_x_continuous(breaks = c(1, 6, 2, 3, 4, 5, 7 ),
                     labels = c("0", "1e-04", "0.001", "0.01","0.1", "1", 'All_dilution'))+
  scale_y_continuous(n.breaks = 12)+
  theme_classic()+
  labs(y = 'Num', x='')
