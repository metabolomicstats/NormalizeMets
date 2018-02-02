#shiftCor function in the statTarget package, edited as commented

shiftCor_fn<-function (samPeno, samFile, Frule = 0.8, QCspan = 0.75, degree = 2, 
                       imputeM = "KNN", saveoutput, outputname,lg) 
{
  cat("\nData File Checking Start..., Time: ", date(), "\n")
  samPeno <- read.csv(samPeno, header = TRUE, check.names = FALSE, 
                      stringsAsFactors = FALSE)
  samPeno <- as.data.frame(samPeno)
  samFile <- read.csv(samFile, header = FALSE, check.names = FALSE, 
                      stringsAsFactors = FALSE)
  samFile <- t(samFile)
  colnames(samFile) <- samFile[1, ]
  
  colnames(samFile)[1]<-"name" #ADDED
  
  samFile <- as.data.frame(samFile[-1, ])
  rownames(samFile) <- samFile$name
  # message("\n", dim(samPeno)[1], " Pheno Samples vs ", dim(samFile)[1], 
  #         " Profile samples", sep = "")
  message("\n", dim(samPeno)[1], " samples in sampledata and ", dim(samFile)[1], 
          "samples in featuredata", sep = "")
  
    message("\nThe samples list (*NA, missing data from featuredata)")
  mcdat <- samFile[, 1][match(samPeno[, 1], samFile[, 1])]
  print(as.vector(mcdat))
  if (any(is.na(mcdat))) {
    stop("Missing data from featuredata! Check your data please!!")
  }
  if (dim(samFile)[1] - dim(samPeno)[1] > 0) {
    message("\nWarning: The sample size in featuredata is larger than sampledata! ")
  }
  else if (dim(samFile)[1] - dim(samPeno)[1] < 0) {
    stop("The sample size in featuredata should be no less than sampledata!\n      Check your data please!!")
  }
  samFP <- samFile[samPeno$sample, ]
  if (sum(is.na(samPeno$class)) <= 0) {
    stop("There were no QC samples in your data!")
  }
  samPeno_stat <- samPeno
  rownames(samPeno_stat) <- samPeno[, 1]
  qc_seq <- rownames(samPeno_stat)
  qc_seq_tmp <- grep("QC", qc_seq)
  sam_seq_tmp <- samPeno_stat[-c(qc_seq_tmp), ]
  if (sum(is.na(samPeno_stat$batch)) > 0) {
    stop("There were missing values in batch! Check your data please!\n")
  }
  if (sum(is.na(sam_seq_tmp$class)) > 0) {
    stop("There were missing values (ungrouped data) in sample class! \n         Check your data please!\n")
  }
  samPeno_stat$class[is.na(samPeno_stat$class)] <- "QC"
  data_stat = aggregate(samPeno_stat$class, by = list(Category = samPeno_stat$class), 
                        FUN = length)
  batch_stat = aggregate(samPeno_stat$class, by = list(Category = samPeno_stat$batch), 
                         FUN = length)
  colnames(data_stat) <- c("Class", "No.")
  colnames(batch_stat) <- c("Batch", "No.")
  message("\nsampledata information:")
  print(data_stat)
  print(batch_stat)
  num_sam <- dim(samFile[, 2:ncol(samFile)])
  num_sam <- data.frame(num_sam)
  colnames(num_sam) <- c("no.")
  rownames(num_sam) <- c("QC and samples", "Metabolites")
  message("\nfeaturedata information:")
  print(num_sam)
  samFP <- as.matrix(samFP)
  samFP[samFP == 0] <- NA
  message("\nWarning: The zeros are assigned missing") #Added shiftCor_fn
  
  cat("\nstatTarget: shiftCor start...Time: ", date(), "\n\nStep 1: Evaluation of missing value...")
  message("\nThe number of NA value in Data Profile before QC-RLSC: ", 
          sum(is.na(samFP)))
  imsamFP <- samFP
  FilterMV = function(m, degree) {
    dx <- c()
    for (i in 1:ncol(m)) {
      freq <- as.vector(tapply(m[, i], degree, function(x) {
        sum(is.na(x))/length(x)
      }))
      if (sum(freq > Frule) > 0) 
        dx <- c(dx, i)
    }
    if (length(dx) > 0) 
      m <- m[, -dx]
    return(m)
  }
  classF <- as.factor(samPeno$class)
  classF = addNA(classF)
  imsamFPF = FilterMV(imsamFP, classF)
  Frule_warning = paste("\nThe number of variables including", 
                        Frule * 100, "% of missing value :", sep = " ")
  message(Frule_warning, " ", dim(imsamFP)[2] - dim(imsamFPF)[2])
  imsamFP = as.matrix(imsamFPF)
  cat("\nStep 2: Imputation start...\n")
  if (imputeM == "KNN") {
    mvd <- impute::impute.knn(imsamFP[, 2:ncol(imsamFP)], 
                              rowmax = 0.99, colmax = 0.99, maxp = 15000)
    inputedData <- mvd$data
  }
  else if (imputeM == "min") {
    minValue <- function(x, group) {
      group = as.factor(as.numeric(group))
      for (i in 1:dim(x)[1]) {
        for (j in 3:dim(x)[2]) {
          if (is.na(x[i, j]) == TRUE) {
            x[i, j] <- tapply(as.numeric(x[, j]), group, 
                              min, na.rm = TRUE)[group[i]]
          }
        }
      }
      return(x)
    }
    inputedData = missvalue(imsamFP, classF)
    inputedData = inputedData[, -1]
  }
  else if (imputeM == "median") {
    missvalue <- function(x, group) {
      group = as.factor(as.numeric(group))
      for (i in 1:dim(x)[1]) {
        for (j in 3:dim(x)[2]) {
          if (is.na(x[i, j]) == TRUE) {
            x[i, j] <- tapply(as.numeric(x[, j]), group, 
                              median, na.rm = TRUE)[group[i]]
          }
        }
      }
      return(x)
    }
    cat("\n", "The imputation method was set at 'median'")
    inputedData = missvalue(imsamFP, classF)
    inputedData = inputedData[, -1]
  }
  message("\nThe number of NA value in Data Profile after the initial imputation: ", 
          sum(is.na(inputedData)))
  if (sum(is.na(inputedData)) > 0) {
    mvd2 <- impute::impute.knn(inputedData[, 1:ncol(inputedData)], 
                               rowmax = 0.99, colmax = 0.99, maxp = 15000)
    inputedData <- mvd2$data
    message("\nThe number of NA value in Data Profile after the second imputation (KNN): ", 
            sum(is.na(inputedData) | as.matrix(inputedData) == 
                  0))
  }
  message("\nImputation Finished!")
  cat("\nStep 3: QC-RLSC Start... Time: ", date())
  dat <- as.matrix(t(inputedData))
  numX <- 1:dim(dat)[2]
  if (QCspan > 0) {
    message("\nWarning: The QCspan was set at ", QCspan, 
            "\n")
    loessFit = function(x, y, QCspan, degree) {
      cn <- colnames(x)
      st_QC <- grep("QC", cn[1])
      ed_QC <- grep("QC", cn[length(cn)])
      if (length(st_QC) == 0) {
        stop("the first sample must be QC sample; please check ......")
      }
      if (length(ed_QC) == 0) {
        stop("the sample at the end of sequence must be QC sample; \n          please check ......")
      }
      qcid <- grep("QC", cn)
      pb <- txtProgressBar(min = 1, max = dim(x)[1], style = 3)
      for (i in 1:dim(x)[1]) {
        loe <- stats::loess(x[i, qcid] ~ qcid, span = QCspan, 
                            degree = degree)
        yf <- stats::predict(loe, y)
        x[i, ] <- as.numeric(x[i, ])/yf
        setTxtProgressBar(pb, i)
      }
      close(pb)
      loessDat = x
    }
    loessDat <- loessFit(x = dat, y = numX, QCspan = QCspan, 
                         degree = degree)
  }
  else if (QCspan <= 0) {
    message("\nWarning: The QCspan was set at '0'.\n", "\nThe GCV was used to", 
            " avoid overfitting the observed data\n")
    autoFit <- function(xl, y) {
      cn <- colnames(xl)
      st_QC <- grep("QC", cn[1])
      ed_QC <- grep("QC", cn[length(cn)])
      if (length(st_QC) == 0) {
        stop("the first sample must be QC sample; please check ......")
      }
      if (length(ed_QC) == 0) {
        stop("the sample at the end of sequence must be QC sample; \n          please check ......")
      }
      qcid <- grep("QC", cn)
      pb <- txtProgressBar(min = 1, max = dim(xl)[1], style = 3)
      for (i in 1:dim(xl)[1]) {
        Sys.sleep(1e-06)
        loe1 <- loess(xl[i, qcid] ~ qcid)
        env <- environment()
        sploe <- function(sp) {
          loe2 <- get("loe1", envir = env)
          mod <- stats::update(loe2, span = sp)
          CVspan = loessGCV_fn(mod)[["gcv"]]#statTarget:::loessGCV(mod)[["gcv"]] #loessGCV(mod)[["gcv"]] GO 5/6
        }
        sp <- c(seq(0.2, 0.75, 0.01))
        CVspan = as.matrix(lapply(sp, sploe))
        CVspan[!is.finite(as.numeric(CVspan))] <- NA
        minG <- data.frame(sp, CVspan)
        minspan <- minG[which.min(minG[, 2]), 1]
        minspan
        loeN <- stats::update(loe1, span = minspan)
        yf <- predict(loeN, y)
        xl[i, ] <- as.numeric(xl[i, ])/yf
        setTxtProgressBar(pb, i)
      }
      close(pb)
      loessDat = xl
    }
    loessDat <- autoFit(xl = dat, y = numX)
  }
  loessDat <- as.matrix(loessDat)
  if (length(loessDat[loessDat < 0L]) > 0) {
    loessDat[loessDat < 0L] <- 0
  }
  loessDat[loessDat == 0L] <- NA
  if (sum(is.na(loessDat)) > 0) {
    mvd2 <- impute::impute.knn(loessDat[, 1:ncol(loessDat)], 
                               rowmax = 0.99, colmax = 0.99, maxp = 15000)
    loessDat <- mvd2$data
    message("\nThe number of NA value in Data Profile after QC-RLS Correction (KNN): ", 
            sum(is.na(loessDat)))
  }
  
  
  #For saving the output
  if (saveoutput){ 
    dirout.uni = paste(getwd(), "/", outputname, "/", sep = "")
    dirsc.ID = getwd()
    dir.create(dirout.uni)
    # dirout.w = paste(getwd(), "/statTarget/shiftCor", sep = "")
    # dir.create(dirout.w)
    # dirout.Bs = paste(getwd(), "/statTarget/shiftCor/Before_shiftCor", 
    #                   sep = "")
    # dir.create(dirout.Bs)
    # dirout.As = paste(getwd(), "/statTarget/shiftCor/After_shiftCor", 
    #                   sep = "")
    # dir.create(dirout.As)
    dirout.loplot <- paste(getwd(), "/", outputname, "/", "figures", 
                           sep = "")
    dir.create(dirout.loplot)
    
    cat("\nHigh-resolution images output...")
    for (i in 1:dim(dat)[1]) {
      loplot_fn(dat, loessDat, i,dirout.loplot=dirout.loplot)
    }
  }
  #Calculation of CV and RSD  
  #   raw_temp <- cbind(samPeno, inputedData)
  #   nam_qc <- rownames(raw_temp)
  #   QC_temp_raw <- grep("QC", nam_qc)
  #   QC_temp_raw <- raw_temp[c(QC_temp_raw), ]
  #   raw_temp_qc <- QC_temp_raw[, -c(3, 4)]
  #   rownames(raw_temp_qc) <- NULL
  # #  RSD30_CV = paste("shift_QC_raw", ".csv", sep = "")
  # #  write.csv(raw_temp_qc, paste(dirout.Bs, RSD30_CV, sep = "/"))
  #   
  #   # cat("\n\nCalculation of CV distribution of raw peaks (QC)...\n\n")
  #   # Rsdist_QC_raw = statTarget:::RsdCal(raw_temp_qc, batch = TRUE, DistPattern = TRUE, 
  #   #                        output = FALSE)
  #   sam_temp_raw <- grep("QC", nam_qc)
  #   sam_temp_raw <- raw_temp[-c(sam_temp_raw), ]
  #   raw_temp_sam <- sam_temp_raw[, -c(3, 4)]
  #   rownames(raw_temp_sam) <- NULL
  # #  RSD30_CV = paste("shift_sam_raw", ".csv", sep = "")
  # #  write.csv(raw_temp_sam, paste(dirout.Bs, RSD30_CV, sep = "/"))
  #   # Rsdist_sam_raw = statTarget:::RsdCal(raw_temp_sam, batch = FALSE, DistPattern = FALSE, 
  #   #                         output = FALSE)
  lo_temp <- cbind(samPeno, t(loessDat))
  
  #Added outdata
  outdata<-list()
  outdata$sampledata<-samPeno
  outdata$metabolitedata<-NULL
  if (lg){
    outdata$featuredata<-LogTransform(t(loessDat))$featuredata
#    message("\nFeaturedata is logtransformed after rlsc normalization.")
  }  else
    outdata$featuredata<-t(loessDat)
  rownames(outdata$featuredata)<-rownames(outdata$sampledata)
  outdata$featuredata<-editcolnames(outdata$featuredata)
  
  #   nam_qc <- rownames(lo_temp)
  #   QC_temp <- grep("QC", nam_qc)
  #   QC_temp <- lo_temp[-c(QC_temp), ]
  #   lo_temp_sam <- QC_temp[, -c(2, 4)]
  #   rownames(lo_temp_sam) <- NULL
  # #  RSD30_CV = paste("shift_sample_cor", ".csv", sep = "")
  # #  write.csv(lo_temp_sam, paste(dirout.uni, RSD30_CV, sep = "/"))
  # #  Rsdist_sam_cor = statTarget:::RsdCal(lo_temp_sam, batch = FALSE, DistPattern = FALSE, 
  # #                          output = FALSE)
  #   QC_temp <- grep("QC", nam_qc)
  #   QC_temp <- lo_temp[c(QC_temp), ]
  #   lo_temp_qc <- QC_temp[, -c(3, 4)]
  #   rownames(lo_temp_qc) <- NULL
  # #  RSD30_CV = paste("shift_QC_cor", ".csv", sep = "")
  # #  write.csv(lo_temp_qc, paste(dirout.uni, RSD30_CV, sep = "/"))
  #   # cat("\n\nCalculation of CV distribution of corrected peaks (QC)...\n\n")
  #   # Rsdist_QC_cor = statTarget:::RsdCal(lo_temp_qc, batch = TRUE, DistPattern = TRUE, 
  #   #                        output = TRUE)
  #   # statTarget:::RSDdist(Rsdist_sam_raw, Rsdist_sam_cor, Rsdist_QC_raw, Rsdist_QC_cor)
  lo_temp_all <- lo_temp[, -c(2, 4)]
  rownames(lo_temp_all) <- NULL
  if (saveoutput){
    RSD30_CV = paste(outputname, ".csv", sep = "")
    write.csv(lo_temp_all, paste(dirout.uni, RSD30_CV, sep = "/"))
    cat("\n\nCorrection Finished! Time: ", date())
    setwd(dirsc.ID)
    tmpfilesc = paste(getwd(), "/tmp", sep = "")
    unlink(tmpfilesc, recursive = TRUE)
  }
  
  return(outdata)
}
