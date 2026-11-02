# Exploratory differential-abundance screen by family and named genus.
# Uses relative abundance of full bacterial libraries and shuffles cultivar labels within sequencing batch.
.script_args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(.script_args)) stop("Run this script with Rscript so its location can be resolved.")
.script_file <- normalizePath(sub("^--file=", "", .script_args[1]), winslash = "/", mustWork = TRUE)
source(file.path(dirname(dirname(.script_file)), "package_config.R"))
valid<-INPUTS_ROOT

outdir<-file.path(RESULTS_ROOT,"triagem_taxons");dir.create(outdir,showWarnings=FALSE)
dat<-read.delim(file.path(valid,"matriz_bacterias_ASV_16_raizes.tsv"),check.names=FALSE,stringsAsFactors=FALSE)
meta<-read.delim(file.path(valid,"metadados_16_raizes.tsv"),stringsAsFactors=FALSE)
samples<-meta$sample_id; counts<-as.matrix(dat[,samples,drop=FALSE]);storage.mode(counts)<-"numeric"
meta$cultivar<-factor(meta$cultivar,levels=c("BL-65","BL-69","BL-74","BL-183"));meta$sequencing_batch<-factor(meta$sequencing_batch,levels=c("Lote 1","Lote 2"))
gcode<-as.integer(meta$cultivar);batch<-meta$sequencing_batch;blocks<-lapply(levels(batch),function(b)which(batch==b));nperm<-9999L
screen_rank<-function(rank,seed){
  tx<-trimws(ifelse(is.na(dat[[rank]]),"",as.character(dat[[rank]])));ok<-nzchar(tx);features<-sort(unique(tx[ok]));
  agg<-rowsum(counts[ok,,drop=FALSE],group=tx[ok],reorder=TRUE);agg<-agg[features,,drop=FALSE]
  Y<-t(sweep(agg,2,colSums(counts),"/"))*100
  # Y: samples by taxon; use all sample reads in the denominator, including unclassified bacteria.
  n<-nrow(Y);grand<-colSums(Y);ss_total<-colSums((Y-matrix(colMeans(Y),n, ncol(Y),byrow=TRUE))^2)
  group_ss<-function(gc){Z<-matrix(0,n,4);Z[cbind(seq_len(n),gc)]<-1;gs<-crossprod(Z,Y);colSums(gs^2/4)-grand^2/n}
  ss_obs<-group_ss(gcode);F_obs<-(ss_obs/3)/((ss_total-ss_obs)/12);eta<-ss_obs/ss_total
  testable<-apply(Y,2,function(v)length(unique(v))>1L);exceed<-integer(ncol(Y));set.seed(seed)
  for(it in seq_len(nperm)){
    gp<-gcode
    for(ii in blocks)gp[ii]<-sample(gp[ii],length(ii),replace=FALSE)
    ss<-group_ss(gp);Fp<-(ss/3)/((ss_total-ss)/12)
    exceed<-exceed+as.integer(Fp>=F_obs-sqrt(.Machine$double.eps))
  }
  p<-(exceed+1)/(nperm+1);p[!testable]<-1
  means<-sapply(1:4,function(g)colMeans(Y[gcode==g,,drop=FALSE]))
  colnames(means)<-paste0("mean_",make.names(levels(meta$cultivar)),"_percent")
  indicator<-levels(meta$cultivar)[max.col(means,ties.method="first")]
  ans<-data.frame(taxon=features,F=F_obs,R2=eta,p=p,q_BH=p.adjust(p,"BH"),prevalence_percent=colMeans(Y>0)*100,
                  max_mean_cultivar=indicator,means,stringsAsFactors=FALSE,check.names=FALSE)
  ans<-ans[order(ans$q_BH,ans$p),]
  write.csv(ans,file.path(outdir,paste0("triagem_",tolower(rank),"_permutacao_lote.csv")),row.names=FALSE)
  data.frame(rank=rank,n_features=nrow(ans),testable=sum(testable),p_raw_lt_0_05=sum(ans$p<.05),q_BH_lt_0_05=sum(ans$q_BH<.05),nperm=nperm,block="within sequencing batch",top_taxon=ans$taxon[1],top_q=ans$q_BH[1],stringsAsFactors=FALSE)
}
summ<-rbind(screen_rank("Familia",20261022L),screen_rank("Genero",20261023L))
write.csv(summ,file.path(outdir,"resumo_triagem_taxons.csv"),row.names=FALSE)
print(summ,row.names=FALSE)
for(r in c("Familia","Genero")){f<-file.path(outdir,paste0("triagem_",tolower(r),"_permutacao_lote.csv"));d<-read.csv(f,check.names=FALSE);cat("\nTOP",r,"\n");print(head(d[,c("taxon","F","R2","p","q_BH","prevalence_percent","max_mean_cultivar")],10),row.names=FALSE)}

