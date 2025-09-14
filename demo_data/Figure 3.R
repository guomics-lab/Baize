#Code Normalization Outlier


# #Protein Matrix
# 
# #Data
# library(readr)
# Omni <- as.data.frame(read_tsv("OmniProtV2_Bovine_human_mixed_24minDIA_20250110OmniProt-report.tsv"))
# 
# A <- strsplit(Omni$Protein.Names, "_")
# B <- sapply(A, length)
# table(B)
# C <- as.data.frame(Omni$Protein.Names[which(B > 2)])
# Omni <- Omni[-which(B > 2),]
# rm(A, B, C)
# 
# write_tsv(Omni, "Omni_report.tsv")
# 
# library(iq)
# process_long_format("Omni_report.tsv",
#                     sample_id = "File.Name",
#                     intensity_col = "Precursor.Normalised",
#                     output_filename = "Omni_report-pg-global.txt",
#                     annotation_col = c("Protein.Names", "Genes"),
#                     filter_double_less = c("Q.Value" = "0.01", "PG.Q.Value" = "0.05",
#                                            "Lib.Q.Value" = "0.01", "Lib.PG.Q.Value" = "0.01"))
# rm(Omni)
# 
# Neat <- as.data.frame(read_tsv("OmniProtV2_Bovine_human_mixed_24minDIA_20250110-report.tsv"))
# 
# A <- strsplit(Neat$Protein.Names, "_")
# B <- sapply(A, length)
# table(B)
# C <- as.data.frame(Neat$Protein.Names[which(B > 2)])
# Neat <- Neat[-which(B > 2),]
# rm(A, B, C)
# 
# write_tsv(Neat, "Neat_report.tsv")
# 
# process_long_format("Neat_report.tsv",
#                     sample_id = "File.Name",
#                     intensity_col = "Precursor.Normalised",
#                     output_filename = "Neat_report-pg-global.txt",
#                     annotation_col = c("Protein.Names", "Genes"),
#                     filter_double_less = c("Q.Value" = "0.01", "PG.Q.Value" = "0.05",
#                                            "Lib.Q.Value" = "0.01", "Lib.PG.Q.Value" = "0.01"))
# rm(Neat)


#Protein

#Data
library(readr)
Omni <- as.data.frame(read_table("Omni_report-pg-global.txt"))

Species_Omni <- strsplit(Omni$Protein.Names, "_")
Species_Omni <- sapply(Species_Omni, function(x) x[2])
table(Species_Omni)

rownames(Omni) <- Omni$Protein.Group
Omni <- Omni[,-c(1:3)]
colnames(Omni) <- gsub("\\.raw\\.dia$", "", basename(colnames(Omni)))
Omni <- Omni[,-which(colnames(Omni) == "WAN20250110gaohh_OmniProtV2_Bovine_human_mixed_24minDIA_300ng_296")]

Omni <- 2^Omni

Neat <- as.data.frame(read_table("Neat_report-pg-global.txt"))

Species_Neat <- strsplit(Neat$Protein.Names, "_")
Species_Neat <- sapply(Species_Neat, function(x) x[2])
table(Species_Neat)

rownames(Neat) <- Neat$Protein.Group
Neat <- Neat[,-c(1:3)]
colnames(Neat) <- gsub("\\.raw\\.dia$", "", basename(colnames(Neat)))

Neat <- 2^Neat

library(readxl)
Label <- as.data.frame(read_xlsx("Sample_information20250112-1851.xlsx", sheet = 1))
Label <- Label[which(Label$Raw_name %in% c(colnames(Omni), colnames(Neat))),]
rownames(Label) <- Label$Raw_name

Label_Omni <- Label[colnames(Omni),]
Label_Neat <- Label[colnames(Neat),]
rm(Label)

#Data analysis

#Normalization
library(NormalyzerDE)

Prot <- rownames(Omni)
Omni <- as.data.frame(medianNormalization(as.matrix(Omni), noLogTransform = F))
rownames(Omni) <- Prot
Omni <- 2^Omni
rm(Prot)

Prot <- rownames(Neat)
Neat <- as.data.frame(medianNormalization(as.matrix(Neat), noLogTransform = F))
rownames(Neat) <- Prot
Neat <- 2^Neat
rm(Prot)

#Omni
Bovine_omni <- which(Species_Omni == "BOVIN")

omni_mean1 <- apply(Omni[,which(Label_Omni$Plasma_ratio == "1:0")], 1, mean, na.rm = T)
omni_mean2 <- apply(Omni[,which(Label_Omni$Plasma_ratio == "1:1")], 1, mean, na.rm = T)
omni_mean3 <- apply(Omni[,which(Label_Omni$Plasma_ratio == "1:1.5")], 1, mean, na.rm = T)
omni_mean4 <- apply(Omni[,which(Label_Omni$Plasma_ratio == "1:2")], 1, mean, na.rm = T)
omni_mean5 <- apply(Omni[,which(Label_Omni$Plasma_ratio == "1:9")], 1, mean, na.rm = T)
omni_mean6 <- apply(Omni[,which(Label_Omni$Plasma_ratio == "1:99")], 1, mean, na.rm = T)

summary(omni_mean1[Bovine_omni]/omni_mean2[Bovine_omni])
summary(omni_mean2[Bovine_omni]/omni_mean3[Bovine_omni])
summary(omni_mean1[Bovine_omni]/omni_mean4[Bovine_omni])
summary(omni_mean1[Bovine_omni]/omni_mean5[Bovine_omni])
summary(omni_mean1[Bovine_omni]/omni_mean6[Bovine_omni])
summary(omni_mean2[Bovine_omni]/omni_mean5[Bovine_omni])

#Neat
Bovine_neat <- which(Species_Neat == "BOVIN")

neat_mean1 <- apply(Neat[,which(Label_Neat$Plasma_ratio == "1:0")], 1, mean, na.rm = T)
neat_mean2 <- apply(Neat[,which(Label_Neat$Plasma_ratio == "1:1")], 1, mean, na.rm = T)
neat_mean3 <- apply(Neat[,which(Label_Neat$Plasma_ratio == "1:1.5")], 1, mean, na.rm = T)
neat_mean4 <- apply(Neat[,which(Label_Neat$Plasma_ratio == "1:2")], 1, mean, na.rm = T)
neat_mean5 <- apply(Neat[,which(Label_Neat$Plasma_ratio == "1:9")], 1, mean, na.rm = T)
neat_mean6 <- apply(Neat[,which(Label_Neat$Plasma_ratio == "1:99")], 1, mean, na.rm = T)

summary(neat_mean1[Bovine_neat]/neat_mean2[Bovine_neat])
summary(neat_mean2[Bovine_neat]/neat_mean3[Bovine_neat])
summary(neat_mean1[Bovine_neat]/neat_mean4[Bovine_neat])
summary(neat_mean1[Bovine_neat]/neat_mean5[Bovine_neat])
summary(neat_mean1[Bovine_neat]/neat_mean6[Bovine_neat])
summary(neat_mean2[Bovine_neat]/neat_mean5[Bovine_neat])

#Boxplot 1
library(ggplot2)
library(ggsci)
df1 <- data.frame(FC = omni_mean1[Bovine_omni]/omni_mean2[Bovine_omni], Group = rep("OmniProt", length(Bovine_omni)))
df1$FC <- log2(df1$FC)

Omni_Upper <- quantile(df1$FC, 0.75, na.rm = T) + 1.5 *IQR(df1$FC, na.rm = T)
Omni_Lower <- quantile(df1$FC, 0.25, na.rm = T) - 1.5 *IQR(df1$FC, na.rm = T)

df1 <- df1[which(df1$FC < Omni_Upper & df1$FC > Omni_Lower),]
rm(Omni_Upper, Omni_Lower)

df2 <- data.frame(FC = neat_mean1[Bovine_neat]/neat_mean2[Bovine_neat], Group = rep("Neat Plasma", length(Bovine_neat)))
df2$FC <- log2(df2$FC)

Neat_Upper <- quantile(df2$FC, 0.75, na.rm = T) + 1.5 *IQR(df2$FC, na.rm = T)
Neat_Lower <- quantile(df2$FC, 0.25, na.rm = T) - 1.5 *IQR(df2$FC, na.rm = T)

df2 <- df2[which(df2$FC < Neat_Upper & df2$FC > Neat_Lower),]
rm(Neat_Upper, Neat_Lower)

df <- as.data.frame(rbind(df1, df2))
rm(df1, df2)

p1 <- ggplot(df, aes(x = Group, y = FC, fill = Group, color = Group)) +
  geom_boxplot(size = 0.5, color = "black") +
  geom_hline(yintercept = log2(100/50), color = "black", linetype = "dashed", linewidth = 1) +
  theme_classic() +
  scale_fill_jama() +
  scale_color_jama() +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 13, color = "black"), axis.text.y = element_text(size = 13, color = "black")) +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14)) +
  labs(x = "", y = "log2(Fold change)") +
  theme(legend.position = "none") +
  ggtitle("1:0 VS 1:1")
rm(df)

#Boxplot 2
df1 <- data.frame(FC = omni_mean2[Bovine_omni]/omni_mean3[Bovine_omni], Group = rep("OmniProt", length(Bovine_omni)))
df1$FC <- log2(df1$FC)

Omni_Upper <- quantile(df1$FC, 0.75, na.rm = T) + 1.5 *IQR(df1$FC, na.rm = T)
Omni_Lower <- quantile(df1$FC, 0.25, na.rm = T) - 1.5 *IQR(df1$FC, na.rm = T)

df1 <- df1[which(df1$FC < Omni_Upper & df1$FC > Omni_Lower),]
rm(Omni_Upper, Omni_Lower)

df2 <- data.frame(FC = neat_mean2[Bovine_neat]/neat_mean3[Bovine_neat], Group = rep("Neat Plasma", length(Bovine_neat)))
df2$FC <- log2(df2$FC)

Neat_Upper <- quantile(df2$FC, 0.75, na.rm = T) + 1.5 *IQR(df2$FC, na.rm = T)
Neat_Lower <- quantile(df2$FC, 0.25, na.rm = T) - 1.5 *IQR(df2$FC, na.rm = T)

df2 <- df2[which(df2$FC < Neat_Upper & df2$FC > Neat_Lower),]
rm(Neat_Upper, Neat_Lower)

df <- as.data.frame(rbind(df1, df2))
rm(df1, df2)

p2 <- ggplot(df, aes(x = Group, y = FC, fill = Group, color = Group)) +
  geom_boxplot(size = 0.5, color = "black") +
  geom_hline(yintercept = log2(50/40), color = "black", linetype = "dashed", linewidth = 1) +
  theme_classic() +
  scale_fill_jama() +
  scale_color_jama() +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 13, color = "black"), axis.text.y = element_text(size = 13, color = "black")) +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14)) +
  labs(x = "", y = "log2(Fold change)") +
  theme(legend.position = "none") +
  ggtitle("1:1 VS 1:1.5")
rm(df)

#Boxplot 3
df1 <- data.frame(FC = omni_mean1[Bovine_omni]/omni_mean4[Bovine_omni], Group = rep("OmniProt", length(Bovine_omni)))
df1$FC <- log2(df1$FC)

Omni_Upper <- quantile(df1$FC, 0.75, na.rm = T) + 1.5 *IQR(df1$FC, na.rm = T)
Omni_Lower <- quantile(df1$FC, 0.25, na.rm = T) - 1.5 *IQR(df1$FC, na.rm = T)

df1 <- df1[which(df1$FC < Omni_Upper & df1$FC > Omni_Lower),]
rm(Omni_Upper, Omni_Lower)

df2 <- data.frame(FC = neat_mean1[Bovine_neat]/neat_mean4[Bovine_neat], Group = rep("Neat Plasma", length(Bovine_neat)))
df2$FC <- log2(df2$FC)

Neat_Upper <- quantile(df2$FC, 0.75, na.rm = T) + 1.5 *IQR(df2$FC, na.rm = T)
Neat_Lower <- quantile(df2$FC, 0.25, na.rm = T) - 1.5 *IQR(df2$FC, na.rm = T)

df2 <- df2[which(df2$FC < Neat_Upper & df2$FC > Neat_Lower),]
rm(Neat_Upper, Neat_Lower)

df <- as.data.frame(rbind(df1, df2))
rm(df1, df2)

p3 <- ggplot(df, aes(x = Group, y = FC, fill = Group, color = Group)) +
  geom_boxplot(size = 0.5, color = "black") +
  geom_hline(yintercept = log2(100/(100/3)), color = "black", linetype = "dashed", linewidth = 1) +
  theme_classic() +
  scale_fill_jama() +
  scale_color_jama() +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 13, color = "black"), axis.text.y = element_text(size = 13, color = "black")) +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14)) +
  labs(x = "", y = "log2(Fold change)") +
  theme(legend.position = "none") +
  ggtitle("1:0 VS 1:2")
rm(df)

#Boxplot 4
df1 <- data.frame(FC = omni_mean1[Bovine_omni]/omni_mean5[Bovine_omni], Group = rep("OmniProt", length(Bovine_omni)))
df1$FC <- log2(df1$FC)

Omni_Upper <- quantile(df1$FC, 0.75, na.rm = T) + 1.5 *IQR(df1$FC, na.rm = T)
Omni_Lower <- quantile(df1$FC, 0.25, na.rm = T) - 1.5 *IQR(df1$FC, na.rm = T)

df1 <- df1[which(df1$FC < Omni_Upper & df1$FC > Omni_Lower),]
rm(Omni_Upper, Omni_Lower)

df2 <- data.frame(FC = neat_mean1[Bovine_neat]/neat_mean5[Bovine_neat], Group = rep("Neat Plasma", length(Bovine_neat)))
df2$FC <- log2(df2$FC)

Neat_Upper <- quantile(df2$FC, 0.75, na.rm = T) + 1.5 *IQR(df2$FC, na.rm = T)
Neat_Lower <- quantile(df2$FC, 0.25, na.rm = T) - 1.5 *IQR(df2$FC, na.rm = T)

df2 <- df2[which(df2$FC < Neat_Upper & df2$FC > Neat_Lower),]
rm(Neat_Upper, Neat_Lower)

df <- as.data.frame(rbind(df1, df2))
rm(df1, df2)

p4 <- ggplot(df, aes(x = Group, y = FC, fill = Group, color = Group)) +
  geom_boxplot(size = 0.5, color = "black") +
  geom_hline(yintercept = log2(100/10), color = "black", linetype = "dashed", linewidth = 1) +
  theme_classic() +
  scale_fill_jama() +
  scale_color_jama() +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 13, color = "black"), axis.text.y = element_text(size = 13, color = "black")) +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14)) +
  labs(x = "", y = "log2(Fold change)") +
  theme(legend.position = "none") +
  ggtitle("1:0 VS 1:9")
rm(df)

#Boxplot 5
df1 <- data.frame(FC = omni_mean1[Bovine_omni]/omni_mean6[Bovine_omni], Group = rep("OmniProt", length(Bovine_omni)))
df1$FC <- log2(df1$FC)

Omni_Upper <- quantile(df1$FC, 0.75, na.rm = T) + 1.5 *IQR(df1$FC, na.rm = T)
Omni_Lower <- quantile(df1$FC, 0.25, na.rm = T) - 1.5 *IQR(df1$FC, na.rm = T)

df1 <- df1[which(df1$FC < Omni_Upper & df1$FC > Omni_Lower),]
rm(Omni_Upper, Omni_Lower)

df2 <- data.frame(FC = neat_mean1[Bovine_neat]/neat_mean6[Bovine_neat], Group = rep("Neat Plasma", length(Bovine_neat)))
df2$FC <- log2(df2$FC)

Neat_Upper <- quantile(df2$FC, 0.75, na.rm = T) + 1.5 *IQR(df2$FC, na.rm = T)
Neat_Lower <- quantile(df2$FC, 0.25, na.rm = T) - 1.5 *IQR(df2$FC, na.rm = T)

df2 <- df2[which(df2$FC < Neat_Upper & df2$FC > Neat_Lower),]
rm(Neat_Upper, Neat_Lower)

df <- as.data.frame(rbind(df1, df2))
rm(df1, df2)

p5 <- ggplot(df, aes(x = Group, y = FC, fill = Group, color = Group)) +
  geom_boxplot(size = 0.5, color = "black") +
  geom_hline(yintercept = log2(100/1), color = "black", linetype = "dashed", linewidth = 1) +
  theme_classic() +
  scale_fill_jama() +
  scale_color_jama() +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 13, color = "black"), axis.text.y = element_text(size = 13, color = "black")) +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14)) +
  labs(x = "", y = "log2(Fold change)") +
  theme(legend.position = "none") +
  ggtitle("1:0 VS 1:99")
rm(df)

#Boxplot 6
df1 <- data.frame(FC = omni_mean2[Bovine_omni]/omni_mean5[Bovine_omni], Group = rep("OmniProt", length(Bovine_omni)))
df1$FC <- log2(df1$FC)

Omni_Upper <- quantile(df1$FC, 0.75, na.rm = T) + 1.5 *IQR(df1$FC, na.rm = T)
Omni_Lower <- quantile(df1$FC, 0.25, na.rm = T) - 1.5 *IQR(df1$FC, na.rm = T)

df1 <- df1[which(df1$FC < Omni_Upper & df1$FC > Omni_Lower),]
rm(Omni_Upper, Omni_Lower)

df2 <- data.frame(FC = neat_mean2[Bovine_neat]/neat_mean5[Bovine_neat], Group = rep("Neat Plasma", length(Bovine_neat)))
df2$FC <- log2(df2$FC)

Neat_Upper <- quantile(df2$FC, 0.75, na.rm = T) + 1.5 *IQR(df2$FC, na.rm = T)
Neat_Lower <- quantile(df2$FC, 0.25, na.rm = T) - 1.5 *IQR(df2$FC, na.rm = T)

df2 <- df2[which(df2$FC < Neat_Upper & df2$FC > Neat_Lower),]
rm(Neat_Upper, Neat_Lower)

df <- as.data.frame(rbind(df1, df2))
rm(df1, df2)

p6 <- ggplot(df, aes(x = Group, y = FC, fill = Group, color = Group)) +
  geom_boxplot(size = 0.5, color = "black") +
  geom_hline(yintercept = log2(50/10), color = "black", linetype = "dashed", linewidth = 1) +
  theme_classic() +
  scale_fill_jama() +
  scale_color_jama() +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 13, color = "black"), axis.text.y = element_text(size = 13, color = "black")) +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14)) +
  labs(x = "", y = "log2(Fold change)") +
  theme(legend.position = "none") +
  ggtitle("1:1 VS 1:9")
rm(df)

#Patchwork
library(patchwork)
p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(ncol = 3)

rm(list = ls())
gc()

