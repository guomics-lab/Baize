# chenghonghan
library(openxlsx)
library(ggplot2)
library(ggsci)
library(ggsignif)
library(reshape2)
library(dplyr)

##### matrix ####
sampleinfo = read.xlsx("Sample_information20241118.xlsx")
rownames(sampleinfo) = sampleinfo$Raw_name
sampleinfo_PRP = subset(sampleinfo , grepl('PRP: PPP', Expriment_ID))

datpro_np23 = read.csv(paste0("Platelet_different_ratio_evaluation_np23_21samples20241118.pg_matrix.tsv"), sep = '\t', row.names = 1)
rownames(datpro_np23) = paste0(datpro_np23$Protein.Ids, '_', datpro_np23$Genes)
datpro_np23 = datpro_np23[5:ncol(datpro_np23)]
dim(datpro_np23)
colnames(datpro_np23)
colnames(datpro_np23) = lapply(colnames(datpro_np23), function(x){
  strsplit(x, '.', fixed = T)[[1]][5]
})
datpre_np23 = read.csv(paste0("Platelet_different_ratio_evaluation_np23_21samples20241118.pr_matrix.tsv"), sep = '\t' )
colnames(datpre_np23) = lapply(colnames(datpre_np23), function(x){
  x = gsub('D..gaohuanhuan.WAN20241117gaohh_OmniProtV2_24minDIA_platelet.', '',x,  fixed = T)
  x = gsub('.raw', '', x, fixed = T)
})

datpro_neat = read.csv(paste0("Platelet_different_ratio_evaluation_neat_21samples20241118.pg_matrix.tsv"), sep = '\t', row.names = 1)
rownames(datpro_neat) = paste0(datpro_neat$Protein.Ids, '_', datpro_neat$Genes)
datpro_neat = datpro_neat[5:ncol(datpro_neat)]
dim(datpro_neat)
colnames(datpro_neat)
colnames(datpro_neat) = lapply(colnames(datpro_neat), function(x){
  strsplit(x, '.', fixed = T)[[1]][5]
})
datpre_neat = read.csv(paste0("Platelet_different_ratio_evaluation_neat_21samples20241118.pr_matrix.tsv"), sep = '\t' )
colnames(datpre_neat) = lapply(colnames(datpre_neat), function(x){
  x = gsub('D..gaohuanhuan.WAN20241117gaohh_OmniProtV2_24minDIA_platelet.', '',x,  fixed = T)
  x = gsub('.raw', '', x, fixed = T)
})
##1 
pronum = colSums(!is.na(datpro_np23))
pronum = data.frame(np23_pro = pronum, row.names = colnames(datpro_np23))
for(i in colnames(datpro_np23)){
  temp = datpre_np23[c("Protein.Group", "Stripped.Sequence", i)]
  colnames(temp)[3]='cur'
  temp = subset(temp, !is.na(cur))
  #temp = temp %>% distinct(Stripped.Sequence, .keep_all = TRUE)
  pronum[i, 'np23_pre'] = length(unique(temp$Stripped.Sequence))
}
#pronum = cbind(pronum, sampleinfo_PRP[rownames(pronum), c(1:6)])
#pronum$ps = paste0(pronum$Plasma_ratio, '_', pronum$Sample_ID)
pronum$sample = rownames(pronum)
pronum = melt(pronum, id.vars = 'sample')
pronum$Plasma_type = sampleinfo_PRP[pronum$sample, 'Plasma_type']
#pronum$sample = sampleinfo_PRP[pronum$sample, 'Sample_ID']
pronum$Plasma_ratio = sampleinfo_PRP[pronum$sample, 'Plasma_ratio']
pronum$type = 'np23'
identifyNum = pronum

pronum = colSums(!is.na(datpro_neat))
pronum = data.frame(neat_pro = pronum, row.names = colnames(datpro_neat))
for(i in colnames(datpro_neat)){
  temp = datpre_neat[c("Protein.Group", "Stripped.Sequence", i)]
  colnames(temp)[3]='cur'
  temp = subset(temp, !is.na(cur))
  #temp = temp %>% distinct(Stripped.Sequence, .keep_all = TRUE)
  pronum[i, 'neat_pre'] = length(unique(temp$Stripped.Sequence))
}

pronum$sample = rownames(pronum)
pronum = melt(pronum, id.vars = 'sample')
pronum$Plasma_type = sampleinfo_PRP[pronum$sample, 'Plasma_type']
pronum$Plasma_ratio = sampleinfo_PRP[pronum$sample, 'Plasma_ratio']
pronum$type = 'neat'
identifyNum = rbind(identifyNum, pronum)

head(identifyNum)


df1.2 = identifyNum
df1.2$x = as.numeric(factor((df1.2$Plasma_ratio)))
head(df1.2)
df1.2 %>% na.omit() %>% 
  group_by(variable, x) %>% 
  summarise(mean_value=mean(value), sd_value=sd(value)) -> df1.3
df1.3 = data.frame(df1.3)
df1.3
df1.3$type = rep(c('np23', 'neat' ), each=14)
df1.3$p = rep(c( 'pro', 'pre') , each=7, 2)
df1.3$vpx = paste0(df1.3$variable, df1.3$x)
df1.3$xp = paste0(df1.3$p, df1.3$x)

df1.3[29, ] = c('np23_pre',8, 37792, 0, 'np23', 'pre', 'np23_pre8', 'pre8', 37792)
df1.3[30, ] = c('np23_pro',8, 4909, 0, 'np23', 'pro', 'np23_pro8', 'pre8', 4909)
df1.3[31, ] = c('neat_pre',8, 8706, 0, 'neat', 'pre', 'neat_pre8', 'pre8', 8706+37792)
df1.3[32, ] = c('neat_pro',8, 1333, 0, 'neat', 'pro', 'neat_pre8', 'pre8', 1333+4909)
df1.3$mean_value1 = as.numeric(df1.3$mean_value)

head(df1.2)
ggplot(df1.2[grepl('pro', df1.2$variable), ], aes(x = x, y = value, fill = type)) + 
  geom_bar(stat = "summary", fun ="mean", position = position_dodge(), alpha=0.7) +
  geom_point( )+
  stat_summary(fun.data = 'mean_se', geom = "errorbar", colour = "black", width = 0.3, position = position_dodge( .9))+
  geom_text(stat = "summary", fun = "mean", aes(label = round(..y.., 0)), color = "black", size = 3.5)+
  scale_x_continuous(breaks = c(1, 3, 4, 5, 6, 7, 2, 8 ),
                     labels = c("0:1", "1:1", "1:1.5", "1:2",   "1:9",   "1:99",  "1:0", 'All_dilution'))+
  scale_y_continuous(n.breaks = 10)+
  theme_classic()+
  labs(y = 'Num', x='')