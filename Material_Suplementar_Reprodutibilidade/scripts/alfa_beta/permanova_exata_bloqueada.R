# Exact enumeration of all distinct cultivar assignments under within-batch restrictions.
# Each batch has two observations per cultivar: 8!/(2!)^4 = 2520 assignments.
# Total exact randomization space: 2520^2 = 6,350,400 assignments.
.script_args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(.script_args)) stop("Run this script with Rscript so its location can be resolved.")
.script_file <- normalizePath(sub("^--file=", "", .script_args[1]), winslash = "/", mustWork = TRUE)
source(file.path(dirname(dirname(.script_file)), "package_config.R"))
valid <- INPUTS_ROOT

meta<-read.delim(file.path(valid,"metadados_16_raizes.tsv"),stringsAsFactors=FALSE)
meta$cultivar<-factor(meta$cultivar,levels=c("BL-65","BL-69","BL-74","BL-183"))
meta$sequencing_batch<-factor(meta$sequencing_batch,levels=c("Lote 1","Lote 2"))
make_assignments<-function() {
  out<-matrix(0L,2520L,8L); r<-0L
  for(a in combn(8,2,simplify=FALSE)) {
    rem1<-setdiff(1:8,a)
    for(b in combn(rem1,2,simplify=FALSE)) {
      rem2<-setdiff(rem1,b)
      for(cc in combn(rem2,2,simplify=FALSE)) {
        rem3<-setdiff(rem2,cc); r<-r+1L
        out[r,a]<-1L; out[r,b]<-2L; out[r,cc]<-3L; out[r,rem3]<-4L
      }
    }
  }
  stopifnot(r==2520L,all(rowSums(out==1L)==2L),all(rowSums(out==2L)==2L),all(rowSums(out==3L)==2L),all(rowSums(out==4L)==2L))
  out
}
assignments<-make_assignments()
readD<-function(p)as.matrix(read.delim(p,row.names=1,check.names=FALSE))
make_distances<-function(comm) {
  rel<-comm/rowSums(comm); n<-nrow(comm)
  bc<-ja<-matrix(0,n,n,dimnames=list(rownames(comm),rownames(comm)))
  for(i in seq_len(n-1L))for(j in (i+1L):n){
    bc[i,j]<-bc[j,i]<-sum(abs(rel[i,]-rel[j,]))/sum(rel[i,]+rel[j,])
    ai<-comm[i,]>0;aj<-comm[j,]>0;ja[i,j]<-ja[j,i]<-1-sum(ai&aj)/sum(ai|aj)
  }
  lx<-log(comm+1);clr<-sweep(lx,1L,rowMeans(lx),"-")
  list(Bray_Curtis=bc,Jaccard=ja,Aitchison=as.matrix(dist(clr)))
}
exact_test<-function(D,group,batch) {
  lev<-levels(group); blev<-levels(droplevels(batch)); i1<-which(batch==blev[1]);i2<-which(batch==blev[2]); n<-nrow(D)
  stopifnot(length(i1)==8L,length(i2)==8L,all(table(group,batch)==2L))
  J<-diag(n)-matrix(1/n,n,n);G<--0.5*J%*%(D^2)%*%J; total<-sum(diag(G))
  G1<-G[i1,i1,drop=FALSE];G2<-G[i2,i2,drop=FALSE];G12<-G[i1,i2,drop=FALSE]
  within_score<-function(Gb,A) {
    v<-numeric(nrow(A))
    for(g in 1:4){M<-1.0*(A==g);v<-v+rowSums((M%*%Gb)*M)/4}
    v
  }
  w1<-within_score(G1,assignments);w2<-within_score(G2,assignments)
  cross<-matrix(0,nrow(assignments),nrow(assignments))
  for(g in 1:4){M1<-1.0*(assignments==g);M2<-M1;cross<-cross+(M1%*%G12)%*%t(M2)}
  ss_all<-outer(w1,w2,"+")+0.5*cross
  obs_codes<-as.integer(group[i1]);obs_codes2<-as.integer(group[i2])
  obs_i<-which(apply(assignments,1,function(x)all(x==obs_codes)))
  obs_j<-which(apply(assignments,1,function(x)all(x==obs_codes2)))
  stopifnot(length(obs_i)==1L,length(obs_j)==1L)
  ss_obs<-ss_all[obs_i,obs_j]
  F_obs<-(ss_obs/3)/((total-ss_obs)/12)
  p_exact<-sum(ss_all>=ss_obs-sqrt(.Machine$double.eps))/length(ss_all)
  c(F=F_obs,R2=ss_obs/total,p_exact=p_exact,exceed=sum(ss_all>=ss_obs-sqrt(.Machine$double.eps)),space=length(ss_all))
}
asvdir<-file.path(valid,"permanova_bloqueada")
asvd<-lapply(c("Bray_Curtis","Jaccard","Aitchison"),function(m)readD(file.path(asvdir,paste0("distancia_",m,".tsv"))))
names(asvd)<-c("Bray_Curtis","Jaccard","Aitchison")
ids<-rownames(asvd[[1]]);ma<-meta[match(ids,meta$sample_id),];stopifnot(!anyNA(ma$sample_id))
gdat<-read.delim(file.path(valid,"sensibilidade_genero","matriz_genero_com_nao_classificados.tsv"),check.names=FALSE)
samples<-grep("^L[12]_Go",names(gdat),value=TRUE);gc<-t(as.matrix(gdat[,samples,drop=FALSE]));storage.mode(gc)<-"numeric";rownames(gc)<-samples
gm<-meta[match(samples,meta$sample_id),];stopifnot(!anyNA(gm$sample_id));gd<-make_distances(gc)
ans<-list();ri<-0L
for(mx in c("ASV","Genero_sensibilidade")){
 dl<-if(mx=="ASV")asvd else gd;md<-if(mx=="ASV")ma else gm
 for(m in names(dl)){
  cat("Enumera??o exata:",mx,m,"\n");flush.console()
  z<-exact_test(dl[[m]],md$cultivar,md$sequencing_batch)
  ri<-ri+1L;ans[[ri]]<-data.frame(matriz=mx,metric=m,F=unname(z["F"]),R2=unname(z["R2"]),
    p_exact=unname(z["p_exact"]),exceed=unname(z["exceed"]),permutation_space=unname(z["space"]))
 }}
out<-do.call(rbind,ans);out$p_holm<-ave(out$p_exact,out$matriz,FUN=function(p)p.adjust(p,"holm"))
write.csv(out,file.path(RESULTS_ROOT,"permanova","PERMANOVA_exata_bloqueada.csv"),row.names=FALSE)
print(out,row.names=FALSE)
