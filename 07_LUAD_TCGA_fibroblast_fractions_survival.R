# Armon Azizi
#
# Do survival analysis on TCGA fractions deconvolved using fibroblast signature matrix from lambrechts

library(survival)
library(survminer)
library(data.table)

setwd("~/Bioinformatics/sc_prognosis/")

split_fractions<-as.data.frame(fread("cibersort_outputs/deconvolution/luad_rename_fibroblast_split_tcga_deconvolution/CIBERSORTx_Adjusted.txt"))
whole_fractions<-as.data.frame(fread("cibersort_outputs/deconvolution/luad_rename_fibroblast_whole_tcga_deconvolution/CIBERSORTx_Adjusted.txt"))
survival<-as.data.frame(fread("TCGA_LUAD/TCGA-LUAD.survival.tsv"))

samples<-intersect(survival$sample, split_fractions$Mixture)
samples<-intersect(samples, whole_fractions$Mixture)

fractions<-cbind(split_fractions[match(samples,split_fractions$Mixture),c("Mixture","Fibroblast_adv", "Fibroblast_fav", "Immune", "Endothelial", "Tumor")],
                 whole_fractions[match(samples,whole_fractions$Mixture),c("Fibroblasts","Immune", "Endothelial", "Tumor"),drop=FALSE])
colnames(fractions)<-c("Mixture","Fibroblast_adv", "Fibroblast_fav", "Immune", "Endothelial", "Tumor", "Fibroblasts","Immune_whole", "Endothelial_whole", "Tumor_whole")

survival<-survival[match(samples,survival$sample),]
colnames(survival)<-c("sample","censor","patient","os")

fractions$fibroblast_ratio<-fractions$Fibroblast_adv/(fractions$Fibroblast_adv+fractions$Fibroblast_fav)
fractions$fibroblast_difference<-fractions$Fibroblast_adv-fractions$Fibroblast_fav

luad_data<-cbind(fractions,survival[,2:4])

luad_data<-luad_data[!is.na(luad_data$fibroblast_ratio),]

plot_input<-reshape2::melt(as.matrix(luad_data[,2:7]))
ggplot(plot_input,aes(x=Var2,y=value)) +
  geom_violin(fill="grey40") +
  geom_jitter(size=0.5, width = 0.1, alpha=0.3, fill="grey") +
  theme_bw() +
  xlab("celltype") +
  ylab("CSX Fraction") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, color="black", size=10),
        legend.position = "none") 


luad_data$fibroblast_total<-luad_data$Fibroblast_adv+luad_data$Fibroblast_fav
ggplot(luad_data, aes(x=Fibroblasts, y=fibroblast_total)) +
  geom_point() +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  stat_cor(method = "pearson") +
  geom_smooth(method = "lm", se = FALSE, color="red")

ggplot(luad_data, aes(x=Immune, y=Immune_whole)) +
  geom_point() +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  stat_cor(method = "pearson") +
  geom_smooth(method = "lm", se = FALSE, color="red")

ggplot(luad_data, aes(x=Tumor, y=Tumor_whole)) +
  geom_point() +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  stat_cor(method = "pearson") +
  geom_smooth(method = "lm", se = FALSE, color="red")

ggplot(luad_data, aes(x=Endothelial, y=Endothelial_whole)) +
  geom_point() +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  stat_cor(method = "pearson") +
  geom_smooth(method = "lm", se = FALSE, color="red")



hist(luad_data$fibroblast_ratio, breaks=(0:20)/20)

ggplot(luad_data, aes(x=Fibroblast_adv, y=Fibroblast_fav)) +
  geom_point() +
  theme_bw()

m<-coxph(Surv(os, censor) ~ Fibroblasts+Immune_whole+Tumor_whole, luad_data)
m
ggforest(m, main = "Multivariate Cox Analysis")

m<-coxph(Surv(os, censor) ~ fibroblast_ratio+Immune+Tumor, luad_data)
m
ggforest(m, main = "Multivariate Cox Analysis")

m<-coxph(Surv(os, censor) ~ fibroblast_difference+Immune+Tumor, luad_data)
m
ggforest(m, main = "Multivariate Cox Analysis")

m<-coxph(Surv(os, censor) ~ Fibroblast_fav+Immune+Tumor, luad_data)
m
ggforest(m, main = "Multivariate Cox Analysis")

m<-coxph(Surv(os, censor) ~ Fibroblast_adv+Immune+Tumor, luad_data)
m
ggforest(m, main = "Multivariate Cox Hazard Ratio")

luad_data$fibroblast_ratio_split<-rep(1,nrow(luad_data))
luad_data$fibroblast_ratio_split[luad_data$fibroblast_ratio<quantile(luad_data$fibroblast_ratio, probs = c(0.4))]<-0
fit<-survfit(Surv(os, censor) ~ fibroblast_ratio_split, data = luad_data)

luad_data$`Fibroblast Ratio`<-luad_data$fibroblast_ratio_split
m<-coxph(Surv(os, censor) ~ `Fibroblast Ratio`+Immune+Tumor, luad_data)
m
ggforest(m, main = "Multivariate Cox Hazard Ratio", fontsize = 1, )

ggsurvplot(fit, data = luad_data,
           legend.labs = c("Low", "High"),
           pval = TRUE,
           pval.method = TRUE,
           palette = c("goldenrod", "darkblue"),
           ylab="Survival Probability",
           xlab="Time (days)",risk.table = TRUE)


# Calculate optimal split
for(i in (2:10/10)){
  print(i)
  luad_data$split<-rep(1,nrow(luad_data))
  luad_data$split[luad_data$fibroblast_ratio<quantile(luad_data$fibroblast_ratio, probs = c(i))]<-0
  fit<-survfit(Surv(os, censor) ~ split, data = luad_data)
  
  print(coxph(Surv(os, censor) ~ split, luad_data))
  
  print(ggsurvplot(fit, data = luad_data,
                   legend.labs = c("Low", "High"),
                   pval = TRUE,
                   pval.method = TRUE,
                   palette = c("goldenrod", "darkblue"),
                   ylab="Survival Probability",
                   xlab="Time (days)",
                   title=i))
}
