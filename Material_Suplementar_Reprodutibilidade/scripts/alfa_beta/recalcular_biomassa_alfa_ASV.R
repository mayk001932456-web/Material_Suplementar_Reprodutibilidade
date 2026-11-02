.script_args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(.script_args)) stop("Run this script with Rscript so its location can be resolved.")
.script_file <- normalizePath(sub("^--file=", "", .script_args[1]), winslash = "/", mustWork = TRUE)
source(file.path(dirname(dirname(.script_file)), "package_config.R"))
valid<-INPUTS_ROOT

outdir<-file.path(RESULTS_ROOT,"biomassa_alfa");dir.create(outdir,showWarnings=FALSE)
alfa<-read.csv(file.path(RESULTS_ROOT,"alfa_bacteriana","alfa_rarefeita_ASV_70554_100iter.csv"),check.names=FALSE,stringsAsFactors=FALSE)
meta<-read.delim(file.path(valid,"metadados_16_raizes.tsv"),stringsAsFactors=FALSE)
alfa<-alfa[,c("SampleID","Cultivar","Batch","observed","shannon")]
alfa$label<-alfa$SampleID
biom<-data.frame(label=c("T1-B2","T1-B4","T2-B1","T2-B4","T3-B2","T3-B4","T4-B1","T4-B4","T5-B2","T5-B4","T6-B1","T6-B4","T7-B1","T7-B3","T8-B1","T8-B4"),
 biomass_g=c(14.0338,3.3530,3.8510,2.8270,3.7740,4.3889,10.8511,26.5423,4.3752,11.5200,1.6725,2.0196,15.6452,1.5354,7.5186,22.7642),stringsAsFactors=FALSE)
# Values above were transcribed from the biomass column of the original root weighing workbook and checked against the existing mapping script.
df<-merge(alfa,meta[,c("label","sample_id","sequencing_batch","cultivar")],by="label",all=TRUE,sort=FALSE)
df<-merge(df,biom,by="label",all=TRUE,sort=FALSE)
stopifnot(nrow(df)==16L,!anyNA(df$sample_id),!anyNA(df$biomass_g),all(df$Cultivar==df$cultivar),all(df$Batch==ifelse(df$sequencing_batch=="Lote 1","L1","L2")))
df<-df[match(meta$label,df$label),]
write.csv(df,file.path(outdir,"biomassa_e_alfa_por_amostra.csv"),row.names=FALSE)
metrics<-c("observed","shannon");res<-list()
for(m in metrics){
  pcor<-cor.test(df$biomass_g,df[[m]],method="pearson")
  scor<-cor.test(df$biomass_g,df[[m]],method="spearman",exact=FALSE)
  fit<-lm(df[[m]]~log(df$biomass_g)+factor(df$Cultivar)+factor(df$sequencing_batch))
  cf<-summary(fit)$coefficients
  term<-"log(df$biomass_g)"
  res[[m]]<-data.frame(metric=m,n=nrow(df),pearson_r=unname(pcor$estimate),pearson_p=pcor$p.value,
    spearman_rho=unname(scor$estimate),spearman_p=scor$p.value,adjusted_log_biomass_slope=cf[term,1],
    adjusted_slope_SE=cf[term,2],adjusted_slope_p=cf[term,4],adjusted_model_R2=summary(fit)$r.squared)
}
out<-do.call(rbind,res);out$pearson_p_holm<-p.adjust(out$pearson_p,"holm");out$spearman_p_holm<-p.adjust(out$spearman_p,"holm");out$adjusted_slope_p_holm<-p.adjust(out$adjusted_slope_p,"holm")
write.csv(out,file.path(outdir,"correlacoes_biomassa_diversidade_alfa.csv"),row.names=FALSE)
summary_cult<-do.call(rbind,lapply(split(df,df$Cultivar),function(x)data.frame(cultivar=x$Cultivar[1],n=nrow(x),mean_biomass_g=mean(x$biomass_g),sd_biomass_g=sd(x$biomass_g),mean_shannon=mean(x$shannon),mean_observed=mean(x$observed))))
write.csv(summary_cult,file.path(outdir,"resumo_por_cultivar.csv"),row.names=FALSE)
print(out,row.names=FALSE);cat("\nBiomass map complete:",nrow(df),"samples; range",range(df$biomass_g),"g\n");print(summary_cult,row.names=FALSE)






