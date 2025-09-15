del_NA_col<-function(mat,ratio=1){
  mat = mat[,apply(mat, 2, function(x) {sum(is.na(x))/length(x) < ratio})]
  return(mat)
}
remove_narow <- function(mat,ratio=1){
  mat = mat[apply(mat, 1, function(x) {sum(is.na(x))/length(x) < ratio}),]
  return(mat)
}

source("D:/Program Files/R/R_function_xuezhangzhi")
library(openxlsx)
samp_info<-read.xlsx("samp_info.xlsx",sheet = 1)
pg_df1<-read.table("DIANN181_OmniProtV2_Heparin_platelet_test_OmniProt_report20240922.pg_matrix.tsv",header = T,sep = "\t")
pg_df1$Protein.Group<-paste(pg_df1$Protein.Group,pg_df1$Genes,sep = "_")
pg_df1<-pg_df1[,-c(2:5)]

pr_df1<-read.table("DIANN181_OmniProtV2_Heparin_platelet_test_OmniProt_report20240922.pr_matrix.tsv",header = T,sep = "\t")
pr_df1<-pr_df1[,10:ncol(pr_df1)]

pg_df2<-read.table("DIANN181_OmniProtV2_Heparin_platelet_test_neat_sample_report20240922.pg_matrix.tsv",header = T,sep = "\t")
pg_df2$Protein.Group<-paste(pg_df2$Protein.Group,pg_df2$Genes,sep = "_")
pg_df2<-pg_df2[,-c(2:5)]

pr_df2<-read.table("DIANN181_OmniProtV2_Heparin_platelet_test_neat_sample_report20240922.pr_matrix.tsv",header = T,sep = "\t")
pr_df2<-pr_df2[,10:ncol(pr_df2)]

names(pg_df1)[2:54]<-sapply(strsplit(names(pg_df1)[2:54],"\\."),function(e){e[9]})
row.names(pg_df1)<-pg_df1$Protein.Group
pg_df1<-pg_df1[,-1]
pg_df1<-data.frame(t(pg_df1))

names(pg_df2)[2:22]<-sapply(strsplit(names(pg_df2)[2:22],"\\."),function(e){e[11]})
row.names(pg_df2)<-pg_df2$Protein.Group
pg_df2<-pg_df2[,-1]
pg_df2<-data.frame(t(pg_df2))

names(pr_df1)[2:54]<-sapply(strsplit(names(pr_df1)[2:54],"\\."),function(e){e[9]})
row.names(pr_df1)<-pr_df1$Precursor.Id
pr_df1<-pr_df1[,-1]
pr_df1<-data.frame(t(pr_df1))

names(pr_df2)[2:22]<-sapply(strsplit(names(pr_df2)[2:22],"\\."),function(e){e[11]})
row.names(pr_df2)<-pr_df2$Precursor.Id
pr_df2<-pr_df2[,-1]
pr_df2<-data.frame(t(pr_df2))

#
library(plyr)
library(dplyr)
group6_1<-pg_df2[row.names(pg_df2)%in%samp_info$Raw.name[samp_info$Group==6],]
group6_2<-pg_df1[row.names(pg_df1)%in%samp_info$Raw.name[samp_info$Group==6],]
group6<-rbind.fill(group6_1,group6_2)
row.names(group6)<-c(row.names(group6_1),row.names(group6_2))

group6$sum<-apply(!is.na(group6),1,sum)
group6$group<-samp_info$Sub.group[match(row.names(group6),samp_info$Raw.name)]


group6$group<-as.factor(group6$group)
group6$Exp.type<-samp_info$Exp.type[match(group6$group,samp_info$Sub.group)]


###########
# Figure 1b/1c
p_data1<-group6[,c(4477,4478,4479)]
p_data1$type<-"prot"
p_data2<-pr_group6[,c(35695,35696,35697)]
p_data2$type<-"pep"
p_data<-rbind(p_data1,p_data2)


# mean and sd
summary_data <- p_data %>%
  group_by(Exp.type, type) %>%
  summarise(
    mean = mean(sum),
    sd = sd(sum),
    .groups = 'drop'
  )

# plot
a <- ggplot() +
  # mean
  geom_col(data = summary_data,
           aes(x = Exp.type, y = mean, fill = type),
           position = position_dodge(0.9),
           width = 0.7) +
  
  geom_errorbar(data = summary_data,
                aes(x = Exp.type, ymin = mean - sd, ymax = mean + sd, group = type),
                position = position_dodge(0.9),
                width = 0.5,
                linewidth = 0.8,
                color = "black") +
  
  # point
  geom_point(data = p_data,
             aes(x = Exp.type, y = sum, group = type),
             position = position_dodge(0.9),
             size = 2.5,
             shape = 21,
             color = "black",
             fill = "black") +
  
  
  geom_text(data = summary_data,
            aes(x = Exp.type, y = mean, group = type, label = round(mean, 1)),
            position = position_dodge(width = 0.9),
            vjust = -0.5,  #
            size = 3.5,    # 
            color = "black") +
  
  # color
  scale_fill_manual(values = c("prot" = "#1f77b4", "pep" = "#ff7f0e"),
                    labels = c("prot" = "Protein", "pep" = "Peptide")) +
  
  # title
  labs(title = "Protein and Peptide Quantification",
       x = "Sample Group",
       y = "Quantity",
       fill = "Measurement") +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.line = element_line(color = "black", linewidth = 0.5),
    legend.position = "top"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

a
ggsave("Fig1D_NP_type_prot_pep_number_column_new_20250704.pdf",a,width = 8,height = 8)
##############


##############
#fig3 
#Figure 6d Figure 6e
library(openxlsx)
library(plyr)
samp_info1<-read.xlsx("lung_samp_info.xlsx",sheet = 1)
file_list<-list.files("E:/Qsync/work/OmniProtV2/OmniProt_V20241012/Figure3_OmniProt/Fig3_v20241014/",pattern = ".pg_matrix.tsv")
temp<-c()
for (i in 1:length(file_list)) {
  indsel<-read.table(paste0("E:/Qsync/work/OmniProtV2/OmniProt_V20241012/Figure3_OmniProt/Fig3_v20241014/",file_list[i]),header = T,sep = "\t",quote="")
  indsel<-indsel[,-c(2,3,4,5)]
  row.names(indsel)<-indsel$Protein.Group
  indsel<-indsel[,-1]
  indsel<-data.frame(t(indsel))
  indsel$ms_name<-row.names(indsel)
  temp<-rbind.fill(temp,indsel) 
}
temp<-temp[,c(9997,1:9996,9998:ncol(temp))]


temp$ms_name<-sapply(strsplit(temp$ms_name,"\\."),function(e){e[5]})

temp$group<-samp_info1$Group_ID[match(temp$ms_name,samp_info1$Raw.name)]
temp$sum<-apply(!is.na(temp[,2:(ncol(temp)-1)]),1,sum)
temp_3A<-temp[temp$group%in%c(5,1,6,2),]




##################
temp_3A<-temp[temp$group%in%c(7,3,8,4),]

# temp_3A<-del_NA_col(temp_3A)

file_list<-list.files("E:/Qsync/work/OmniProtV2/OmniProt_V20241012/Figure3_OmniProt/Fig3_v20241014/",pattern = ".pr_matrix.tsv")
#group 1 2 5 6
file_list<-file_list[c(3,4)]
temp_pep<-c()
for (i in 1:length(file_list)) {
  indsel<-read.table(paste0("E:/Qsync/work/OmniProtV2/OmniProt_V20241012/Figure3_OmniProt/Fig3_v20241014/",file_list[i]),header = T,sep = "\t")
  indsel<-indsel[,10:ncol(indsel)]
  row.names(indsel)<-indsel$Precursor.Id
  indsel<-indsel[,-1]
  indsel<-data.frame(t(indsel))
  indsel$ms_name<-row.names(indsel)
  temp_pep<-rbind.fill(temp_pep,indsel) 
}
temp_pep<-temp_pep[,c(7619,1:7618,7620:ncol(temp_pep))]

temp_pep$ms_name<-sapply(strsplit(temp_pep$ms_name,"\\."),function(e){e[5]})

temp_pep$group<-samp_info1$Group_ID[match(temp_pep$ms_name,samp_info1$Raw.name)]
temp_pep$sum<-apply(!is.na(temp_pep[,2:(ncol(temp_pep)-1)]),1,sum)
temp_pep_3A<-temp_pep[temp_pep$group%in%c(7,3,8,4),]




#########
p_3C_1<-temp_3A[,c(12144,12145)]
p_3C_1$group1<-"prot"
p_3C_2<-temp_pep_3A[,c(25769,25770)]
p_3C_2$group1<-"pep"
p_3C<-rbind(p_3C_1,p_3C_2)
p_3C$group<-factor(p_3C$group,levels = c("7","3","8","4"))
mmm<-data.frame(table(p_3C$group))

summary_data1 <- p_3C %>%
  group_by(group, group1) %>%
  summarise(
    mean = mean(sum),
    sd = sd(sum),
    .groups = 'drop'
  )

a <- ggplot() +
  
  geom_col(data = summary_data1,
           aes(x = group, y = mean, fill = group1),
           position = position_dodge(0.9),
           width = 0.7) +
  
  geom_errorbar(data = summary_data1,
                aes(x = group, ymin = mean - sd, ymax = mean + sd, group = group1),
                position = position_dodge(0.9),
                width = 0.5,
                linewidth = 0.8,
                color = "black") +
  
 
  geom_point(data = p_3C,
             aes(x = group, y = sum, group = group1),
             position = position_dodge(0.9),
             size = 2.5,
             shape = 21,
             color = "black",
             fill = "black") +
  
  
  geom_text(data = summary_data1,
            aes(x = group, y = mean, group = group1, label = round(mean, 1)),
            position = position_dodge(width = 0.9),
            vjust = -0.5,  # 
            size = 3.5,    # 
            color = "black") +
  
  
  scale_fill_manual(values = c("prot" = "#1f77b4", "pep" = "#ff7f0e"),
                    labels = c("prot" = "Protein", "pep" = "Peptide")) +
  
  #
  labs(title = "Protein and Peptide Quantification",
       x = "Sample Group",
       y = "Quantity",
       fill = "Measurement") +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.line = element_line(color = "black", linewidth = 0.5),
    legend.position = "top"
  ) +
  
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

a
ggsave("Fig3C_lung_Neat_NP_plasma_pep_prot_number_column_new_20250709.pdf",a,width = 8,height = 8)



#########################
#neat_plasma, neat_serum, OmniProt_plasma, OmniProt_plasma
fig3H<-temp[temp$group%in%c(3,7,4,8),]

fig3H$type[fig3H$group==3]<-"OmniProt_plasma"
fig3H$type[fig3H$group==7]<-"Neat_plasma"
fig3H$type[fig3H$group==4]<-"OmniProt_serum"
fig3H$type[fig3H$group==8]<-"Neat_serum"
# fig3H<-fig3H[,grepl("\\.",names(fig3H))]
fig3H<-del_NA_col(fig3H)

library(ComplexHeatmap)
library(colorRamp2)
# c("cornflowerblue", "green", "orange","#E7797A","#CDC0DB")
ann_row<-data.frame(fig3H[,3183])
row.names(ann_row)<-fig3H$ms_name
names(ann_row)<-"type"
heatmap_mat<-fig3H[,2:3180]

heatmap_mat[is.na(heatmap_mat)]<-0.5*min(heatmap_mat,na.rm=T)
heatmap_mat<-log2(heatmap_mat)
row.names(heatmap_mat)<-fig3H$ms_name
plot(density(unlist(heatmap_mat)))
col_fun = colorRamp2(c( 0,0.5*max(heatmap_mat,na.rm=T),max(heatmap_mat,na.rm=T)), c( "cornflowerblue",'white', "#E7797A"))
                                                     # lwd=1)))

set.seed(123)
pdf("20241223_omniport_fig3X_complexheatmap.pdf",width = 8,height = 8)

p<-Heatmap(heatmap_mat,
           col = col_fun,
           cluster_columns = T,
           cluster_rows = F,
           
           show_heatmap_legend = T,
           show_column_names = F,
           show_row_names = F,
           border = T,
           
           column_km = 6,
           
           row_names_gp = gpar(fontsize = 8),
           border_gp = gpar(lwd = 2),
           row_split = ann_row,#
           row_gap = unit(2, "mm")#
)
p
dev.off()

set.seed(123)
col_orders2 <- column_order(p)


protein_lists1 <- lapply(col_orders2, function(order) {
  names(heatmap_mat)[order]
})


print(protein_lists)
clus1<-protein_lists1[1]
write.csv(clus1,"20241223_omniport_fig3x_complexheatmap_clust1.csv",row.names = F)
clus2<-protein_lists1[2]
write.csv(clus2,"20241223_omniport_fig3x_complexheatmap_clust2.csv",row.names = F)
clus3<-protein_lists1[3]
write.csv(clus3,"20241223_omniport_fig3x_complexheatmap_clust3.csv",row.names = F)
clus4<-protein_lists1[4]
write.csv(clus4,"20241223_omniport_fig3x_complexheatmap_clust4.csv",row.names = F)
clus5<-protein_lists1[5]
write.csv(clus5,"20241223_omniport_fig3x_complexheatmap_clust5.csv",row.names = F)
clus6<-protein_lists1[6]
write.csv(clus6,"20241223_omniport_fig3x_complexheatmap_clust6.csv",row.names = F)


#########################

####################
#Figure 2g/2h Figure 2i

library(plyr)
file_list<-list.files("E:/Qsync/work/OmniProtV2/Figure3/",pattern="pg_matrix.tsv")
fig_3AA<-c()
for (i in 1:length(file_list)) {
  indsel<-read.table(paste0("E:/Qsync/work/OmniProtV2/Figure3/",file_list[i]),sep="\t",header=T,row.names=1)
  indsel<-indsel[,-c(1:4)]
  indsel<-data.frame(t(indsel))
  row.names(indsel)<-sapply(strsplit(row.names(indsel),"\\."),function(e){e[4]})
  indsel$ms_name<-row.names(indsel)
  indsel$type<-unlist(strsplit(file_list[i],"_"))[1]
  fig_3AA<-rbind.fill(fig_3AA,indsel)
}
fig_3AA<-fig_3AA[,c(738,1:737,740:ncol(fig_3AA),739)]
fig_3AA$sum<-apply(!is.na(fig_3AA[,2:(ncol(fig_3AA)-1)]),1,sum)

file_list<-list.files("E:/Qsync/work/OmniProtV2/Figure3/",pattern="pr_matrix.tsv")
fig_3AA_pep<-c()
for (i in 1:length(file_list)) {
  indsel<-read.table(paste0("E:/Qsync/work/OmniProtV2/Figure3/",file_list[i]),sep="\t",header=T,row.names=10)
  indsel<-indsel[,-c(1:9)]
  indsel<-data.frame(t(indsel))
  row.names(indsel)<-sapply(strsplit(row.names(indsel),"\\."),function(e){e[4]})
  indsel$ms_name<-row.names(indsel)
  indsel$type<-unlist(strsplit(file_list[i],"_"))[1]
  fig_3AA_pep<-rbind.fill(fig_3AA_pep,indsel)
}
which(names(fig_3AA_pep) == "ms_name")

fig_3AA_pep<-fig_3AA_pep[,c(6182,1:6181,6184:ncol(fig_3AA_pep),6183)]
fig_3AA_pep$sum<-apply(!is.na(fig_3AA_pep[,2:(ncol(fig_3AA_pep)-1)]),1,sum)

prot_3AA<-fig_3AA[,3868:3869]
prot_3AA$sub<-"prot"


pep_3AA<-fig_3AA_pep[,28772:28773]
pep_3AA$sub<-"pep"


library(scales)
library(ggh4x)
library(ggplot2)



df_3aa_prot1<-fig_3AA[fig_3AA$type=="Neat",2:(ncol(fig_3AA)-2)]
df_3aa_prot1<-del_NA_col(df_3aa_prot1)

df_3aa_prot2<-fig_3AA[fig_3AA$type=="Top14",2:(ncol(fig_3AA)-2)]
df_3aa_prot2<-del_NA_col(df_3aa_prot2)

df_3aa_prot3<-fig_3AA[fig_3AA$type=="OmniProt",2:(ncol(fig_3AA)-2)]
df_3aa_prot3<-del_NA_col(df_3aa_prot3)



df_3aa_pep1<-fig_3AA_pep[fig_3AA_pep$type=="Neat",2:(ncol(fig_3AA_pep)-2)]
df_3aa_pep1<-del_NA_col(df_3aa_pep1)

df_3aa_pep2<-fig_3AA_pep[fig_3AA_pep$type=="Top14",2:(ncol(fig_3AA_pep)-2)]
df_3aa_pep2<-del_NA_col(df_3aa_pep2)

df_3aa_pep3<-fig_3AA_pep[fig_3AA_pep$type=="OmniProt",2:(ncol(fig_3AA_pep)-2)]
df_3aa_pep3<-del_NA_col(df_3aa_pep3)




#########
p_3A_1<-fig_3AA[,c(3868,3869)]
p_3A_1$group<-"prot"
p_3A_2<-fig_3AA_pep[,c(28772,28773)]
p_3A_2$group<-"pep"
p_3A<-rbind(p_3A_1,p_3A_2)


########

summary_data1 <- p_3A %>%
  group_by(type, group) %>%
  summarise(
    mean = mean(sum),
    sd = sd(sum),
    .groups = 'drop'
  )

a <- ggplot() +
  
  geom_col(data = summary_data1,
           aes(x = type, y = mean, fill = group),
           position = position_dodge(0.9),
           width = 0.7) +
  
  geom_errorbar(data = summary_data1,
                aes(x = type, ymin = mean - sd, ymax = mean + sd, group = group),
                position = position_dodge(0.9),
                width = 0.5,
                linewidth = 0.8,
                color = "black") +
  
  
  geom_point(data = p_3A,
             aes(x = type, y = sum, group = group),
             position = position_dodge(0.9),
             size = 2.5,
             shape = 21,
             color = "black",
             fill = "black") +
  
  
  geom_text(data = summary_data1,
            aes(x = type, y = mean, group = group, label = round(mean, 1)),
            position = position_dodge(width = 0.9),
            vjust = -0.5,  # 向上调整位置
            size = 3.5,    # 字体大小
            color = "black") +
  
  
  scale_fill_manual(values = c("prot" = "#1f77b4", "pep" = "#ff7f0e"),
                    labels = c("prot" = "Protein", "pep" = "Peptide")) +
  
  
  labs(title = "Protein and Peptide Quantification",
       x = "Sample Group",
       y = "Quantity",
       fill = "Measurement") +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.line = element_line(color = "black", linewidth = 0.5),
    legend.position = "top"
  ) +
  
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

a
ggsave("Fig3AA_lung_plasma_Neat_Top14_NP_pep_prot_number_column_20250704.pdf",a,width = 8,height = 8)


##################

#######################

fig_3AA<-fig_3AA[!grepl("\\.",names(fig_3AA))]
fig_line<-fig_3AA[,2:3074]
fig_line1<-aggregate(fig_line[,1:3072],by=list(fig_line$type),mean,na.rm=T)
row.names(fig_line1)<-fig_line1$Group.1

fig_line2<-data.frame(t(fig_line1[,2:ncol(fig_line1)]))
fig_line2$prot<-row.names(fig_line2)
neat_mat<-fig_line2[,c(1,4)]
neat_mat<-neat_mat[!is.na(neat_mat$Neat),]
neat_mat<-neat_mat[order(neat_mat$Neat,decreasing = T),]
neat_mat$order<-seq(1,nrow(neat_mat),1)
neat_mat$log10<-log10(neat_mat$Neat)
neat_mat$group<-"Neat"

omnp_mat<-fig_line2[,c(2,4)]
omnp_mat<-omnp_mat[!is.na(omnp_mat$OmniProt),]
omnp_mat<-omnp_mat[order(omnp_mat$OmniProt,decreasing = T),]
omnp_mat$order<-seq(1,nrow(omnp_mat),1)
omnp_mat$log10<-log10(omnp_mat$OmniProt)
omnp_mat$group<-"OmniProt"

top14_mat<-fig_line2[,c(3,4)]
top14_mat<-top14_mat[!is.na(top14_mat$Top14),]
top14_mat<-top14_mat[order(top14_mat$Top14,decreasing = T),]
top14_mat$order<-seq(1,nrow(top14_mat),1)
top14_mat$log10<-log10(top14_mat$Top14)
top14_mat$group<-"Top14"
names(neat_mat)<-names(omnp_mat)<-names(top14_mat)<-c("intensity","prot","order","log10","group")
df<-rbind(neat_mat,omnp_mat,top14_mat)


library(ggplot2)


a<-ggplot(df, aes(x = order, y = log10, group = group, color = group)) +
  geom_line() +geom_point(shape = 16, size = 3,alpha=0.3) +scale_color_manual(values = c("OmniProt" = "#f9766e", "Neat" = "#5fa664", "Top14" = "#e1c548")) +
  labs(title = "Lung_plasma_neat_Top14_OmniProt_protein_abundance_rank", x = "order", y = "log10 intensity") +
  theme_minimal()+theme_bw()
a
ggsave("Fig3C_Lung_plasma_neat_Top14_OmniProt_protein_abundance_rank20241223.pdf",a,width = 8,height = 8)



###############
#Figure 2c
library(readr)
library(openxlsx)
library(ggplot2)
fig2_mat<-read_tsv("fig2_result.pg_matrix.tsv")

class(fig2_mat)
fig2_mat<-data.frame(fig2_mat)
row.names(fig2_mat)<-fig2_mat$Protein.Group
fig2_mat<-fig2_mat[,-c(1:5)]
names(fig2_mat)<-sapply(strsplit(names(fig2_mat),"\\."),function(e){e[8]})
fig2_mat<-data.frame(t(fig2_mat))


group_fig2<-read.xlsx("nDIA_783files_list_zhangzhi20241121-1349.xlsx",sheet = 1)
fig2_info<-read.xlsx("nDIA_783files_list_zhangzhi20250110-1535.xlsx",sheet = 2)
fig2_info1<-read.xlsx("nDIA_783files_list_zhangzhi20250110-1535.xlsx",sheet = 1)

fig2_info$NP_ID_2<-fig2_info1$NP_ID_1[match(fig2_info$NP_ID,fig2_info1$NP_ID)]
fig2_info$NP_ID_2[is.na(fig2_info$NP_ID_2)]<-fig2_info$NP_ID[is.na(fig2_info$NP_ID_2)]
fig2_info$ms_name<-sapply(strsplit(fig2_info$Raw_name,"\\."),function(e){e[1]})

fig2_mat$NP_ID<-fig2_info$NP_ID_2[match(row.names(fig2_mat),fig2_info$ms_name)]
NP_mat<-fig2_mat[grepl("NP",fig2_mat$NP_ID),]
NP_prot_sum<-data.frame(matrix(NA,nrow = length(unique(NP_mat$NP_ID)),ncol = 2))
names(NP_prot_sum)<-c("NP_type","prot_num")
NP_prot_sum$NP_type<-unique(NP_mat$NP_ID)
for (i in 1:nrow(NP_prot_sum)) {
  indsel<-NP_mat[NP_mat$NP_ID==NP_prot_sum$NP_type[i],-ncol(NP_mat),]
  indsel<-del_NA_col(indsel)
  NP_prot_sum$prot_num[i]<-ncol(indsel)
}

a<-ggplot(NP_prot_sum, aes(x = NP_type, y = prot_num)) +scale_y_continuous(limits =c(0, 8000) ,expand = c(0,0))+
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Protein Numbers",x = "Group",y = "Number", fill = "Type") +theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line())+
  geom_text(aes(label =prot_num), position = position_dodge(width = 0.9), vjust = -0.5)

a
ggsave("Fig2B_NPs_spectral_lib_prot_number_column_new_20250718.pdf",a,width = 8,height = 8)

################
#fig2 heatmap

NP_mat<-NP_mat[order(NP_mat$NP_ID),]
NP_mat<-del_NA_col(NP_mat)
NP_mat_max<-aggregate(NP_mat,by=list(NP_mat$NP_ID),max,na.rm=T)
NP_mat_max<-NP_mat_max[,-1]
write.csv(NP_mat_max,"20241212_test.csv",row.names = F)#check inf 
NP_mat_max<-read.csv("20241212_test.csv",header = T)
row.names(NP_mat_max)<-NP_mat_max$NP_ID
# NP_mat_max[is.infinite(NP_mat_max)]<-NA
annotation_row<-data.frame(NP_mat_max$NP_ID)
row.names(annotation_row)<-row.names(NP_mat_max)
annotation_row$Charge<-fig2_info1$Charge[match(annotation_row$NP_mat_max.NP_ID,fig2_info1$NP_ID_1)]
annotation_row$Hydrophobicity<-fig2_info1$Hydrophobicity[match(annotation_row$NP_mat_max.NP_ID,fig2_info1$NP_ID_1)]
annotation_row$Reaction_class<-fig2_info1$Reaction_class[match(annotation_row$NP_mat_max.NP_ID,fig2_info1$NP_ID_1)]
annotation_row$Functional.group<-fig2_info1$Functional.group[match(annotation_row$NP_mat_max.NP_ID,fig2_info1$NP_ID_1)]
annotation_row$Matrix<-fig2_info1$Matrix[match(annotation_row$NP_mat_max.NP_ID,fig2_info1$NP_ID_1)]
names(annotation_row)[1]<-"group"

heatmap_mat<-data.frame(NP_mat_max[,1:(ncol(NP_mat_max)-1)])

heatmap_mat<-log2(heatmap_mat)
heatmap_mat[is.na(heatmap_mat)]<-0

####
#20241219 complexheatmap and add label
#Figure 2b
library(ComplexHeatmap)
library(colorRamp2)
heatmap_mat<-data.frame(heatmap_mat)

annotation_row <- HeatmapAnnotation(
  df = annotation_row,
  col = list(
    group = c("NP03" = "#5CB85C", "NP10" = "#337AB7", "NP12" = "orange","NP15"= "#CDC0DB", "NP16"= "moccasin", "NP17"= "indianred2", "NP20"= "cornflowerblue",
              "NP24"= "#E7799A", "NP30"= "goldenrod4", "NP34"= "mediumpurple1", "NP40"= "#0072B2", "NP41"= "#009E73", "NP42"= "yellow3", "NP43"= "goldenrod2",
              "NP45"= "olivedrab", "NP48" = "#41b6c4","NP74"= "#FDE724","NP82"= "blue", "NP85" = "red","NP95"= "lightcyan4"),
    Charge = c("Positive" = "#5CB85C", "Negative" = "#337AB7", "Neutral" = "orange"),
    Hydrophobicity = c("Hydrophobic" = "orange", "Hydrophilic" = "#CDC0DB"),
    Reaction_class = c("Methylation" = "#5CB85C", "Free_radical_poly" = "#337AB7", "amino_silane_mod" = "orange","amide_coupling"="#CDC0DB","Sol-Gel"="indianred2"),
    Functional.group = c("-N+(CH3)3" = "#5CB85C", "NO" = "#337AB7", "-NH2" = "orange","-COOH"="#CDC0DB","-SiOH"="indianred2","N-PEA"="goldenrod4"),
    Matrix = c("Polymethacrylate" = "#5CB85C", "Polystyrene" = "#337AB7", "Fe3O4-SiO2 " = "orange","SiO2"="#CDC0DB","SiO2-Al2O3"="indianred","Si-P"="cornflowerblue",
               "Si"="yellow3","Polyamide"= "#009E73","Graphene"="goldenrod4")
  ),
  border = TRUE #
)


row_anno <- rowAnnotation(
  
  group = annotation_row$group,
  
  Charge = annotation_row$Charge,
  
  Hydrophobicity = annotation_row$Hydrophobicity,
  Reaction_class=annotation_row$Reaction_class,
  Functional.group=annotation_row$Functional.group,
  Matrix=annotation_row$Matrix
)
set.seed(123)
heatmap_mat_scaled <- t(apply(heatmap_mat, 1, function(x) {
  (x - mean(x)) / sd(x)
})) 
heatmap_mat_scaled_1 <- scale(heatmap_mat)
col_fun = colorRamp2(c( min(heatmap_mat_scaled_1,na.rm=T),0,max(heatmap_mat_scaled_1,na.rm=T)), c( "cornflowerblue",'white', "#E7797A"))
pdf("Fig2C_NPs_spectral_lib_prot_NPnormalized_heatmap_20250113.pdf",width = 8,height = 8)
p<-Heatmap( heatmap_mat_scaled_1,name = "expression",col = col_fun,
            
         right_annotation = row_anno,
         cluster_columns = T, 
         cluster_rows = F, 
         show_heatmap_legend = T, 
         show_column_names = F, 
         show_row_names = T, 
         border = T, 
         column_km = 5,
         # row_km_repeats = 10,
         column_km_repeats = 10,
         row_names_gp = gpar(fontsize = 2), 
         border_gp = gpar(lwd = 1), 
         row_gap = unit(1, "mm") ,
         legend = list(
           labels_gp = gpar(fontsize = 6) 
         )
         
)
p
dev.off()


##################
#Figure 2d
fig2e_mat<-read.xlsx("fig2_4066prots_concentration_in_human_protein_atlas20241012.xlsx",sheet = 1)
fig2e_lib<-read_tsv("fig2_result_lib_output.tsv")
fig2e_lib<-fig2e_lib[!is.na(fig2e_lib$Genes),]
library(VennDiagram)

venn.plot <- venn.diagram(
  x = list(unique(fig2e_mat$Gene), unique(fig2e_lib$Genes)),
  category.names = c("HPA", "OmniProt"), 
  output = TRUE, 
  filename = NULL, 
  output.type = "ggplot", 
  fontfamily = "serif", 
  fontface = "bold", 
  col = "black", 
  lwd = 2, 
  fill = c("#e5a79a","#abc8e5"), 
  alpha = 0.4 
)

grid.draw(venn.plot)
dev.off()

##################

#Figure 6a ;Figure 6b;  Figure 6c
library(readr)
library(openxlsx)
library(ggplot2)
library(plyr)
library(dplyr)
fig4_info<-read.xlsx("fig4_info.xlsx",sheet = 1)
fig4_mat_neat<-read_tsv("WAN20241021gaohh_OmniprotV2_contamination_neat_treated20241025.pg_matrix.tsv")
fig4_mat_ompt<-read_tsv("ot_treated_report20241024WAN20241021gaohh_OmniprotV2_contamination_OmniPr.gg_matrix.pg_matrix.tsv")
fig4_mat_neat<-fig4_mat_neat[,-c(2:5)]
fig4_mat_ompt<-fig4_mat_ompt[,-c(2:5)]
row.names(fig4_mat_neat)<-fig4_mat_neat$Protein.Group
row.names(fig4_mat_ompt)<-fig4_mat_ompt$Protein.Group
fig4_mat_neat<-data.frame(t(fig4_mat_neat))
fig4_mat_ompt<-data.frame(t(fig4_mat_ompt))
fig4_mat_neat<-fig4_mat_neat[-1,]
fig4_mat_ompt<-fig4_mat_ompt[-1,]
fig4_mat_neat$ms_name<-sapply(strsplit(row.names(fig4_mat_neat),"\\\\"),function(e){e[5]})
fig4_mat_ompt$ms_name<-sapply(strsplit(row.names(fig4_mat_ompt),"\\\\"),function(e){e[5]})

fig4_mat<-rbind.fill(fig4_mat_neat,fig4_mat_ompt)
fig4_mat<-fig4_mat[,c(1857,1:1856,1858:ncol(fig4_mat))]
fig4_mat$ms_name<-gsub(".raw","",fig4_mat$ms_name)
fig4_mat$Project_ID<-fig4_info$Project_ID[match(fig4_mat$ms_name,fig4_info$Raw_name)]
fig4_mat$Group_ID<-fig4_info$Group_ID[match(fig4_mat$ms_name,fig4_info$Raw_name)]
fig4_mat$NP_ID<-fig4_info$NP_ID[match(fig4_mat$ms_name,fig4_info$Raw_name)]

fig4_mat$sum<-apply(!is.na(fig4_mat[,2:4736]),1,sum)


library(scales)
library(ggh4x)




#######
fig4g_mat<-fig4_mat[fig4_mat$Project_ID=="RBC",]

fig4g_mat$samp_id<-sapply(strsplit(fig4g_mat$ms_name,"_"),function(e){e[length(e)]})

fig4g_mat$con<-fig4_info$Con[match(fig4g_mat$ms_name,fig4_info$Raw_name)]
row.names(fig4g_mat)<-fig4g_mat$ms_name

#########################
#BCT_CMP，group 1,4,2,5,3,6 prot indentity bar
fig4i_mat<-fig4_mat[fig4_mat$Project_ID=="BCT_CMP",]

mean_temp_4I <- fig4i_mat %>%
  group_by(Group_ID) %>%
  summarize(
    mean = mean(sum),
    sd = sd(sum),
    n = n(),
    se = sd / sqrt(n),
    lower = mean - 1.96 * se,#
    upper = mean + 1.96 * se
    
  )

fig4I_g1<-fig4i_mat[fig4i_mat$Group_ID==1,2:(ncol(fig4i_mat)-4)]
fig4I_g1<-del_NA_col(fig4I_g1)
fig4I_g2<-fig4i_mat[fig4i_mat$Group_ID==2,2:(ncol(fig4i_mat)-4)]
fig4I_g2<-del_NA_col(fig4I_g2)
fig4I_g3<-fig4i_mat[fig4i_mat$Group_ID==3,2:(ncol(fig4i_mat)-4)]
fig4I_g3<-del_NA_col(fig4I_g3)
fig4I_g4<-fig4i_mat[fig4i_mat$Group_ID==4,2:(ncol(fig4i_mat)-4)]
fig4I_g4<-del_NA_col(fig4I_g4)
fig4I_g5<-fig4i_mat[fig4i_mat$Group_ID==5,2:(ncol(fig4i_mat)-4)]
fig4I_g5<-del_NA_col(fig4I_g5)
fig4I_g6<-fig4i_mat[fig4i_mat$Group_ID==6,2:(ncol(fig4i_mat)-4)]
fig4I_g6<-del_NA_col(fig4I_g6)


library(scales)
library(ggh4x)

##################
#########
p_4I_1<-fig4i_mat[,c(4738,4740)]

p_4I<-p_4I_1
p_4I$Group_ID<-factor(p_4I$Group_ID,levels = c("1","4","2","5","3","6"))


summary_data1 <- p_4I %>%
  group_by(Group_ID) %>%
  summarise(
    mean = mean(sum),
    sd = sd(sum),
    .groups = 'drop'
  )

a <- ggplot() +
 
  geom_col(data = summary_data1,
           aes(x = Group_ID, y = mean),
           position = position_dodge(0.9),
           width = 0.7) +
  
  geom_errorbar(data = summary_data1,
                aes(x = Group_ID, ymin = mean - sd, ymax = mean + sd),
                position = position_dodge(0.9),
                width = 0.5,
                linewidth = 0.8,
                color = "black") +
  
  # 添加单个数据点
  geom_point(data = p_4I,
             aes(x = Group_ID, y = sum),
             position = position_dodge(0.9),
             size = 2.5,
             shape = 21,
             color = "black",
             fill = "black") +
  
  # 
  geom_text(data = summary_data1,
            aes(x = Group_ID, y = mean, label = round(mean, 1)),
            position = position_dodge(width = 0.9),
            vjust = -0.5,  
            size = 3.5,    
            color = "black") +
  
  
  labs(title = "Protein and Peptide Quantification",
       x = "Sample Group",
       y = "Quantity",
       fill = "Measurement") +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.line = element_line(color = "black", linewidth = 0.5),
    legend.position = "top"
  ) +
  # 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

a
ggsave("Fig4I_BCT_neat_omniprot_protein_column_20250709.pdf",a,width = 8,height = 8)





################


listinput <- list("1"=names(fig4I_g1),
                  "4"=names(fig4I_g4),
                  "2"=names(fig4I_g2),
                  "5"=names(fig4I_g5),
                  "3"=names(fig4I_g3),
                  "6"=names(fig4I_g6))
pdf('Fig4J_ BCT_neat_omniprot_protein_upset.pdf',height = 8,width = 8)
upset(fromList(listinput),nsets = 6, order.by = "freq",text.scale =1.5)
dev.off()

fig4i_mat<-del_NA_col(fig4i_mat)
fig4i_mat$Group_ID<-factor(fig4i_mat$Group_ID,levels = c("1","2","3","4","5","6"))
fig4i_mat<-fig4i_mat[order(fig4i_mat$Group_ID),]
fig4i_mat<-fig4i_mat[,!grepl("\\.",names(fig4i_mat))]
fig4i_mat$Plasma_type<-fig4_info$Plasma_type[match(fig4i_mat$ms_name,fig4_info$Raw_name)]
fig4i_mat$group<-paste(fig4i_mat$NP_ID,fig4i_mat$Plasma_type,sep="_")
annotation_row<-data.frame(fig4i_mat$group)
row.names(annotation_row)<-fig4i_mat$ms_name
annotation_row$fig4i_mat.group<-as.factor(annotation_row$fig4i_mat.group)
write.csv(fig4i_mat,"20241128_fig4i_mat.csv",row.names = F)
fig4i_mat<-read.csv("20241128_fig4i_mat.csv",header = T)
heatmap_mat<-log2(fig4i_mat[,2:(ncol(fig4i_mat)-6)])
row.names(heatmap_mat)<-fig4i_mat$ms_name
heatmap_mat[is.na(heatmap_mat)]<-0
heatmap_mat<-as.matrix(heatmap_mat)




library(ComplexHeatmap)
library(colorRamp2)
# c("cornflowerblue", "green", "orange","#E7797A","#CDC0DB")
col_fun = colorRamp2(c( 0,0.5*max(heatmap_mat,na.rm=T),max(heatmap_mat,na.rm=T)), c( "cornflowerblue",'white', "#E7797A"))


set.seed(123)
pdf("20250113_omniport_fig4K_heatmap.pdf",width = 8,height = 8)
p<-Heatmap(heatmap_mat,
           col = col_fun,
           cluster_columns = T,
           cluster_rows = F,
           show_heatmap_legend = T,
           show_column_names = F,
           show_row_names = F,
           border = T,
           column_km  = 5,
           column_km_repeats = 10,
          
           row_names_gp = gpar(fontsize = 8),
           border_gp = gpar(lwd = 2),
           row_split = annotation_row,
           row_gap = unit(2, "mm")
)
p
dev.off()

set.seed(123)
# row_orders <- row_order(p)
column_orders <- column_order(p)


protein_lists <- lapply(column_orders, function(order) {
  colnames(heatmap_mat)[order]
})


print(protein_lists)
clus1<-protein_lists[1]
write.csv(clus1,"20250113_omniport_fig4K_heatmap_clust1.csv",row.names = F)
clus2<-protein_lists[2]
write.csv(clus2,"20250113_omniport_fig4K_heatmap_clust2.csv",row.names = F)
clus3<-protein_lists[3]
write.csv(clus3,"20250113_omniport_fig4K_heatmap_clust3.csv",row.names = F)
clus4<-protein_lists[4]
write.csv(clus4,"20250113_omniport_fig4K_heatmap_clust4.csv",row.names = F)
clus5<-protein_lists[5]
write.csv(clus5,"20250113_omniport_fig4K_heatmap_clust5.csv",row.names = F)

###############
#Figure 2e
library(ggplot2)
library(openxlsx)
df1<-read.xlsx("21NPs_spectral_library_proteins_biological_process.xlsx",sheet = 1)
df1$log10P<- -log10(df1$PValue)

df1<-df1[order(df1$log10P),]
df1$group<-seq(1,nrow(df1),1)


ggp <- ggplot(df1, aes(x = group, y = log10P)) + 
  geom_bar(stat = "identity", fill = "black") +
  geom_line(aes(y = Count / 20), stat = "identity", color = "red", size = 2) +
  labs(title = "Title",
       x = "Term",
       y = "-log10(P)") +
  scale_y_continuous(sec.axis = sec_axis(~ . * 20, name = "ratio")) +
  # scale_x_discrete(breaks = df1$group, labels = df1$Term) +  
  coord_flip() +  
  theme_classic()

print(ggp)
ggsave("21NPs_spectral_library_proteins_biological_process_count_log10P.pdf",ggp,width = 8,height = 8)


################
#Figure 2f
df1<-read.xlsx("DAVID_GO_annotation_KEGG_Top10_pathway.xlsx",sheet = 1)
df1$log10P<- -log10(df1$PValue)

df1<-df1[order(df1$log10P),]
df1$group<-seq(1,nrow(df1),1)

ggp <- ggplot(df1, aes(x = group, y = log10P)) + 
  geom_bar(stat = "identity", fill = "black") +
  geom_line(aes(y = Count / 20), stat = "identity", color = "red", size = 2) +
  labs(title = "Title",
       x = "Term",
       y = "-log10(P)") +
  scale_y_continuous(sec.axis = sec_axis(~ . * 20, name = "ratio")) +
  coord_flip() +  
  theme_classic()

print(ggp)
ggsave("DAVID_GO_annotation_KEGG_Top10_pathway_count_log10P.pdf",ggp,width = 8,height = 8)


































