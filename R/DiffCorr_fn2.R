#edited comp.2.cc.fdr function  for one dataframe
comp.2.cc.fdr_fn2<-function ( 
  data,
  method = "pearson", 
  p.adjust.methods = "BH", 
  threshold = 0.05 
) 
{
  #Compute metabolite-wise correlation
  cc1 <- cor(t(data), method = method) 

  #Extract just the lower.tri
  ccc1 <- as.vector(cc1[lower.tri(cc1)])

  n1 <- ncol(data)
  n <- nrow(data)
  N <- n * (n - 1)/2
  
  #######################
  p1 <- rep(1, N)

  #######################
  
  mol.names <- rownames(cc1)
  p1 <- cor2.test(n1, ccc1)


  #   if (p.adjust.methods == "local") {
  #   p1.lfdr <- get.lfdr(p1)$lfdr
  # 
  #     }
  # #  else (p.adjust.methods == "BH" | p.adjust.methods == "bh") {
  # else {  
    p1.lfdr <- p.adjust(p1, method = p.adjust.methods)
    
#  }

  myindex <- which((lower.tri(cc1)) == TRUE, arr.ind = TRUE)
  mol.names1 <- mol.names[myindex[, 2]]
  mol.names2 <- mol.names[myindex[, 1]]
  
  
  fin.ind <- p1.lfdr < threshold
  res <- data.frame(mol.names1[fin.ind], mol.names2[fin.ind], ccc1[fin.ind], 
                    p1[fin.ind],p1.lfdr[fin.ind])
  
  head <- c("metabolite X", "metabolite Y", "correlation", "p-value", "p-adjusted")
  colnames(res) <- head
  # write.table(res, file = output.file, row.names = FALSE, col.names = FALSE, 
  #             sep = "\t", quote = FALSE)
  return(res)
}