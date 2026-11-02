# Rarefacao sem reposicao e diversidade alfa das 16 amostras radiculares.
# Calcula somente as quatro metricas apresentadas no TCC.
args <- commandArgs(trailingOnly=FALSE)
script_arg <- grep("^--file=", args, value=TRUE)
if (!length(script_arg)) stop("Execute com Rscript para localizar o pacote.")
script_file <- normalizePath(sub("^--file=", "", script_arg[1]), winslash="/", mustWork=TRUE)
source(file.path(dirname(dirname(script_file)), "package_config.R"))

meta <- read.delim(file.path(INPUTS_ROOT,"metadados_16_raizes.tsv"),stringsAsFactors=FALSE,check.names=FALSE)
meta <- meta[meta$role=="raiz",,drop=FALSE]
counts_path <- file.path(INPUTS_ROOT,"matriz_bacterias_ASV_16_raizes.tsv")
dat <- read.delim(counts_path,check.names=FALSE,stringsAsFactors=FALSE)
stopifnot(nrow(meta)==16L,all(meta$sample_id %in% names(dat)))
counts <- as.matrix(dat[,meta$sample_id,drop=FALSE]); storage.mode(counts)<-"integer"
stopifnot(all(is.finite(counts)),all(counts>=0),all(counts==floor(counts)))
depths <- colSums(counts); depth <- min(depths)
if (depth != 70554L) stop("A profundidade minima mudou; confira a matriz e os metadados.")

alpha <- function(v) {
  v <- v[v>0]; richness <- length(v); p <- v/sum(v)
  shannon <- -sum(p*log(p)); f1 <- sum(v==1L); f2 <- sum(v==2L)
  c(observed=richness,chao1=richness+f1*(f1-1)/(2*(f2+1)),
    shannon=shannon,pielou=if(richness>1L) shannon/log(richness) else 0)
}
rarefy_hypergeom <- function(v,size) {
  ix <- which(v>0L); out <- integer(length(v)); remaining_n <- sum(v); remaining_k <- size
  if(length(ix)==1L) { out[ix]<-size; return(out) }
  for(j in seq_len(length(ix)-1L)) {
    i<-ix[j]; draw<-rhyper(1L,v[i],remaining_n-v[i],remaining_k); out[i]<-draw
    remaining_n<-remaining_n-v[i]; remaining_k<-remaining_k-draw
    if(remaining_k==0L) break
  }
  if(remaining_k>0L) out[ix[length(ix)]]<-remaining_k
  out
}
set.seed(20261019); n_iter<-100L
rare <- vector("list",nrow(meta))
for(j in seq_len(nrow(meta))) {
  message(sprintf("Rarefacao %d/%d: %s",j,nrow(meta),meta$label[j]))
  it <- replicate(n_iter,alpha(rarefy_hypergeom(counts[,j],depth)))
  rare[[j]] <- data.frame(SampleID=meta$label[j],Cultivar=meta$cultivar[j],
    Batch=ifelse(meta$sequencing_batch[j]=="Lote 1","L1","L2"),
    as.list(rowMeans(it)),as.list(setNames(apply(it,1,sd),paste0(rownames(it),"_sd"))),
    check.names=FALSE)
}
rare_df<-do.call(rbind,rare)
outdir<-file.path(RESULTS_ROOT,"alfa_bacteriana");dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
write.csv(rare_df,file.path(outdir,"alfa_rarefeita_ASV_70554_100iter.csv"),row.names=FALSE)
summary_mean<-aggregate(rare_df[,c("observed","chao1","shannon","pielou")],list(Cultivar=rare_df$Cultivar),mean)
summary_sd<-aggregate(rare_df[,c("observed","chao1","shannon","pielou")],list(Cultivar=rare_df$Cultivar),sd)
names(summary_mean)[-1]<-paste0(names(summary_mean)[-1],"_mean")
names(summary_sd)[-1]<-paste0(names(summary_sd)[-1],"_sd_between_samples")
summary<-merge(summary_mean,summary_sd,by="Cultivar",sort=FALSE)
summary<-summary[match(c("BL-65","BL-69","BL-74","BL-183"),summary$Cultivar),]
write.csv(summary,file.path(outdir,"resumo_alfa_por_cultivar.csv"),row.names=FALSE)
write.csv(data.frame(sample_id=meta$sample_id,label=meta$label,cultivar=meta$cultivar,
  sequencing_batch=meta$sequencing_batch,library_reads=as.integer(depths),rarefaction_depth=depth,
  iterations=n_iter),file.path(outdir,"mapeamento_amostras_alfa.csv"),row.names=FALSE)
writeLines(c("Rarefacao hipergeometrica sem reposicao; 16 amostras radiculares.",
  paste("Profundidade:",depth),paste("Subamostragens por amostra:",n_iter),
  "Metricas: riqueza observada, Chao1, Shannon e Pielou.",
  "O solo-controle nao foi incluido nas comparacoes entre cultivares.",
  paste("Semente:",20261019),paste("R:",R.version.string)),file.path(outdir,"LEIA-ME.txt"))
