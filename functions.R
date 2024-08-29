
fpkmToTpm <- function(fpkm){
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

plotGene <- function (x, dataset){
  data <- dataset[rownames(dataset)==x,]
  
  datat <- t(data) %>% as.data.frame()
  colnames(datat)<-c("FPKM")
  datat$disease <- factor(datasetmeta$disease_state, levels = c("Normal", "COPD"))

  ggplot(datat, aes(x=disease, y=FPKM, fill=disease))+ geom_jitter() + 
    geom_boxplot(alpha=0.5) +ylab(paste0(x, " Expression"))+ xlab("")+
    ggtitle("GSE57148",subtitle = "Lung homogenates")+ theme_bw() + theme(text = element_text(size = 16))
}

plotGeneTPM <- function (x, tpm){
  data <- tpm[rownames(tpm)==x,]
  
  datat <- t(data) %>% as.data.frame()
  colnames(datat)<-c("TPM")
  datat$disease <- factor(datasetmeta$disease_state, levels = c("Normal", "COPD"))
  
  ggplot(datat, aes(x=disease, y=TPM, fill=disease))+ geom_jitter() + 
    geom_boxplot(alpha=0.5) +ylab(paste0(x, " Expression"))+ xlab("")+
    ggtitle("GSE57148",subtitle = "Lung homogenates")+ theme_bw() + theme(text = element_text(size = 16))
}