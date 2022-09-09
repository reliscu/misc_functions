kme2fidelity <- function(kme, gene_cols=1, kme_col_grep, sample_col_grep) {
  
  rtoz <- function(x){0.5*log((1+x)/(1-x))}
  ztor <- function(x){(exp(2*x)-1)/(exp(2*x)+1)}
  
  z <- apply(kme[,grep(kme_col_grep,colnames(kme))], 2, rtoz)
  n_samples_adj <- apply(kme[,grep(sample_col_grep, colnames(kme))], 2, function(x) x-3)
  mean_z <- rowSums(mapply(function(idx) z[,idx]*n_samples_adj[,idx], 1:ncol(z)), na.rm=T)/rowSums(n_samples_adj, na.rm=T)
  fid <- mean_z/sqrt(1/rowSums(n_samples_adj, na.rm=T))
  n_datasets <- apply(z, 1, function(x) sum(!is.na(x)))
  
  df <- data.frame(kme[,gene_cols], Fidelity=fid, No.datasets=n_datasets)
  return(df)
  
}
