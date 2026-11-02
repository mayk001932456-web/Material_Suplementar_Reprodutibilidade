# Recalculation of PERMDISP using the documented vegan betadisper/permutest algorithm:
# spatial median in PCoA coordinates, Lingoes correction, residual permutation.
# This is a base-R reproduction, not execution of the vegan package.
.script_args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(.script_args)) stop("Run this script with Rscript so its location can be resolved.")
.script_file <- normalizePath(sub("^--file=", "", .script_args[1]), winslash = "/", mustWork = TRUE)
source(file.path(dirname(dirname(.script_file)), "package_config.R"))
valid<-INPUTS_ROOT

nperm<-9999L
meta<-read.delim(file.path(valid,"metadados_16_raizes.tsv"),stringsAsFactors=FALSE)
meta$cultivar<-factor(meta$cultivar,levels=c("BL-65","BL-69","BL-74","BL-183"))
meta$sequencing_batch<-factor(meta$sequencing_batch,levels=c("Lote 1","Lote 2"))
readD<-function(p)as.matrix(read.delim(p,row.names=1,check.names=FALSE))
make_distances<-function(comm){
 rel<-comm/rowSums(comm);n<-nrow(comm);bc<-ja<-matrix(0,n,n,dimnames=list(rownames(comm),rownames(comm)))
 for(i in seq_len(n-1L))for(j in(i+1L):n){bc[i,j]<-bc[j,i]<-sum(abs(rel[i,]-rel[j,]))/sum(rel[i,]+rel[j,]);ai<-comm[i,]>0;aj<-comm[j,]>0;ja[i,j]<-ja[j,i]<-1-sum(ai&aj)/sum(ai|aj)}
 lx<-log(comm+1);clr<-sweep(lx,1L,rowMeans(lx),"-");list(Bray_Curtis=bc,Jaccard=ja,Aitchison=as.matrix(dist(clr)))
}
spatial_median<-function(x,tol=1e-8,maxit=1000L){
 y<-colMeans(x)
 for(i in seq_len(maxit)){dd<-sqrt(rowSums((x-matrix(y,nrow(x),ncol(x),byrow=TRUE))^2));if(any(dd<tol))return(x[which.min(dd),]);w<-1/dd;ynew<-colSums(x*w)/sum(w);if(sqrt(sum((ynew-y)^2))<tol)return(ynew);y<-ynew}
 y
}
restricted_indices<-function(batch){ix<-seq_along(batch);for(b in levels(droplevels(batch))){ii<-which(batch==b);ix[ii]<-sample(ii,length(ii),replace=FALSE)};ix}
permdisp<-function(D,group,batch,seed){
 n<-nrow(D);J<-diag(n)-matrix(1/n,n,n);Dmat<-D
 G0<--0.5*J%*%(Dmat^2)%*%J;ev0<-eigen(G0,symmetric=TRUE,only.values=TRUE)$values;lambda<-max(0,-min(ev0))
 if(lambda>0){D2<-Dmat^2+2*lambda*(1-diag(n));G<--0.5*J%*%D2%*%J}else G<-G0
 ee<-eigen(G,symmetric=TRUE);keep<-ee$values>max(ee$values)*1e-10;coords<-sweep(ee$vectors[,keep,drop=FALSE],2,sqrt(ee$values[keep]),"*")
 z<-numeric(n)
 for(g in levels(group)){ii<-which(group==g);med<-spatial_median(coords[ii,,drop=FALSE]);z[ii]<-sqrt(rowSums((coords[ii,,drop=FALSE]-matrix(med,length(ii),length(med),byrow=TRUE))^2))}
 mod<-lm(z~group);q<-mod$qr;p<-q$rank;rdf<-n-p;Fobs<-unname(summary(mod)$fstatistic[1]);resid<-qr.resid(q,z)
 set.seed(seed);fp<-numeric(nperm)
 for(i in seq_len(nperm)){rp<-resid[restricted_indices(batch)];fit<-qr.fitted(q,rp);err<-qr.resid(q,rp);fp[i]<-(sum((fit-mean(fit))^2)/(p-1))/(sum(err^2)/rdf)}
 c(F=Fobs,p=(1+sum(fp>=Fobs-sqrt(.Machine$double.eps)))/(nperm+1),lingoes=lambda,min_eigenvalue=min(ev0))
}
asvdir<-file.path(valid,"permanova_bloqueada");asvd<-lapply(c("Bray_Curtis","Jaccard","Aitchison"),function(m)readD(file.path(asvdir,paste0("distancia_",m,".tsv"))));names(asvd)<-c("Bray_Curtis","Jaccard","Aitchison")
ids<-rownames(asvd[[1]]);ma<-meta[match(ids,meta$sample_id),];stopifnot(!anyNA(ma$sample_id))
gdat<-read.delim(file.path(valid,"sensibilidade_genero","matriz_genero_com_nao_classificados.tsv"),check.names=FALSE);samples<-grep("^L[12]_Go",names(gdat),value=TRUE);gc<-t(as.matrix(gdat[,samples,drop=FALSE]));storage.mode(gc)<-"numeric";rownames(gc)<-samples;gm<-meta[match(samples,meta$sample_id),];stopifnot(!anyNA(gm$sample_id));gd<-make_distances(gc)
ans<-list();ri<-0L
for(mx in c("ASV","Genero_sensibilidade")){dl<-if(mx=="ASV")asvd else gd;md<-if(mx=="ASV")ma else gm;for(m in names(dl)){cat("PERMDISP",mx,m,"\n");flush.console();z<-permdisp(dl[[m]],md$cultivar,md$sequencing_batch,20261021L+match(m,names(dl)));ri<-ri+1L;ans[[ri]]<-data.frame(matriz=mx,metric=m,F=unname(z["F"]),p=unname(z["p"]),lingoes=unname(z["lingoes"]),permutations=nperm)}}
out<-do.call(rbind,ans);out$p_holm<-ave(out$p,out$matriz,FUN=function(p)p.adjust(p,"holm"));write.csv(out,file.path(RESULTS_ROOT,"permanova","PERMDISP_ASV_e_genero_referencia_vegan.csv"),row.names=FALSE);print(out,row.names=FALSE)
