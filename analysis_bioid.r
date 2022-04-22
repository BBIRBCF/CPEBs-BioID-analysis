###################################
######### BioID analysis ##########
###################################

### required R packages
library(ggplot2)
library(impute)
library(hwriter)
library(limma)
library(openxlsx)

### source functions4BioIDanalysis.R
source(pathfunctions4BioIDanalysis)

###########################################################################
## 1. Import data - iBAQ columns
###########################################################################

## dataDir <- setPath2Data  #(ProteinGroups with iBAQ abundances and SampleNames with metadata)

DAT <- read.table(paste0(dataDir, "proteinGroups.txt"), sep="\t", header=T, stringsAsFactors=F)
sel <- grepl("iBAQ", colnames(DAT))
iBACmat <- DAT[,c(sel)]
rownames(iBACmat) <- DAT[,2]
names(iBACmat)[2:length(names(iBACmat))] <- sub("iBAQ.","2",names(iBACmat)[2:length(names(iBACmat))])
sampleinfo <- read.table(paste0(dataDir, "SampleNames.txt"), sep="\t", header=T, stringsAsFactors=F)
sampleinfo$Raw.file <- gsub("-",".",sampleinfo$Raw.file,fixed=TRUE)
rownames(sampleinfo) <- sampleinfo$Raw.file
iBACmat <- iBACmat[,sampleinfo$Raw.file]
iBACmat <- iBACmat[which(apply(is.na(iBACmat),1,sum)==0),]
gene.names <- DAT$Gene.names
names(gene.names) <- DAT[,2]

## res <- setPath2Results; dir.create(res)

#################################################################################
## 2. Missing values treatment and differential expression between conditions - case study 2591
#################################################################################

controltests <- TRUE
Contrasts2591 <- cbind(c("BirA-RBD-ZZ", "RBD-ZZ-BirA", "xCPEB1-BirA", "BirA-RBD-ZZ"),
                       c("BirA-xCPEB1", "xCPEB1-BirA", "BirA-xCPEB1", "RBD-ZZ-BirA"))

study <- 2591

for(ite in seq_len(dim(eval(parse(text = paste0("Contrasts",study))))[1])){
    Cond1 <- eval(parse(text = paste0("Contrasts",study)))[ite,1]
    Cond2 <- eval(parse(text = paste0("Contrasts",study)))[ite,2]
    nameCont <- paste0(gsub("-",".",Cond1), "_", gsub("-",".",Cond2))

    ## create folder for results
    resp2 <- paste0(res,nameCont,"/")
    dir.create(resp2)

    x <- as.matrix(iBACmat[,sampleinfo$sample == Cond1])
    y <- as.matrix(iBACmat[,sampleinfo$sample == Cond2])
    C <- as.matrix(iBACmat[,sampleinfo$sample == "BirA" & sampleinfo$study == study])


    ## understanding missings values
    t1 <- table(apply(x==0,1,sum),apply(C ==0,1,sum))
    t1a <- array(as.numeric(t1),dim = dim(t1))
    rownames(t1a) <- paste0(Cond1, " = ", rownames(t1))
    colnames(t1a) <- paste0("Control", " = ", colnames(t1))

    t2 <- table(apply(y==0,1,sum),apply(C ==0,1,sum))
    t2a <- array(as.numeric(t2),dim = dim(t2))
    rownames(t2a) <- paste0(Cond2, " = ", rownames(t2))
    colnames(t2a) <- paste0("Control", " = ", colnames(t2))

    t3 <- table(apply(x==0,1,sum),apply(y ==0,1,sum))
    t3a <- array(as.numeric(t3),dim = dim(t3))
    rownames(t3a) <- paste0(Cond1, " = ", rownames(t3))
    colnames(t3a) <- paste0(Cond2, " = ", colnames(t3))

    p <- openPage(paste(resp2, 'missings_table.html', sep=''));
    hwrite(paste('<br>Number of missing values by protein<br><br>'), p);
    for (j in 1:3)
        {
            ts <- eval(parse(text = c(paste0("t",j,"a"))))
            hwrite(paste('<br><br><br>',"table ", j,'<br><br>'), p);
            hwrite(as.matrix(ts), page=p, center=F, row.names=T, col.names=T,
                   col.style=c("text-align:center", rep("text-align:center", ncol(ts)-1)),
                   col.width=c(100, rep(125, ncol(ts)-1)));
        }
    closePage(p);

    mns <- apply(log10(x),1, function(x) mean(x[is.finite(x)]))
    mss <- apply(x==0, 1, sum)
    mss <- mss[!is.na(mns)]
    mns <- mns[!is.na(mns)]
    png(paste0(resp2, "missings_randomness_",Cond1,".png"))
    boxplot(mns ~ mss, lwd = 2, ylab = 'avg log10 iBAC', xlab="n missings", main = Cond1)
    stripchart(mns ~ mss, vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = 'blue')
    dev.off()

    mns <- apply(log10(y),1, function(y) mean(y[is.finite(y)]))
    mss <- apply(y==0, 1, sum)
    mss <- mss[!is.na(mns)]
    mns <- mns[!is.na(mns)]
    png(paste0(resp2, "missings_randomness_",Cond2,".png"))
    boxplot(mns ~ mss, lwd = 2, ylab = 'avg log10 iBAC', xlab="n missings", main = Cond2)
    stripchart(mns ~ mss, vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = 'blue')
    dev.off()

    ## frog names
    frog.x <- vapply(strsplit(colnames(x),"_"),function(x) x[5], character(1))
    frog.y <- vapply(strsplit(colnames(y),"_"),function(x) x[5], character(1))
    frog.C <- vapply(strsplit(colnames(C),"_"),function(x) x[5], character(1))
    colnames(x) <- frog.x
    colnames(y) <- frog.y
    colnames(C) <- frog.C
    nam <- sort(frog.x)

    ## log 10
    x <- log10(x[,nam])
    y <- log10(y[,nam])
    C <- log10(C[,nam])
    x[!is.finite(x)] <- NA
    y[!is.finite(y)] <- NA
    C[!is.finite(C)] <- NA

    ## normalization
    N <- length(nam)
    jd <- cbind(x,y,C)
    colnames(jd) <- paste0(colnames(jd), rep(paste0("_",c(Cond1,Cond2,"C")),each = N))
    jd.norm <- mdQuantileNormalization(jd)

    png(paste0(resp2,"boxplot_intensities_",nameCont, ".png"), width=1200, height = 800)
    par(mfrow=c(1,2))
    boxplot(jd, xlab = "", ylab = "Intensity", main = "Before normalization", las=2)
    boxplot(jd.norm, xlab = "", ylab = "", main = "After normalization", las=2)
    dev.off()

    x <- jd.norm[,1:N]
    y <- jd.norm[,1:N + N]
    C <- jd.norm[,1:N + N*2]

    ############
    ### A. subdata with not many missings: impute missing values with KNN when there is at least 2 non-missing values
    ############
    cond <- apply(!is.na(x),1,sum) > 1 & apply(!is.na(y),1,sum) > 1
    x2 <- x[cond,]
    y2 <- y[cond,]
    C2 <- C[cond,]

    x.imp <- impute.knn(x2 ,k = 10, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=23)
    y.imp <- impute.knn(y2 ,k = 10, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=23)
    C.imp <- impute.knn(C2 ,k = 10, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=23)

    C.imp$data[apply(is.na(C2),1,sum)>3,] <- NA

    ## Joint data
    ed <- cbind(x.imp$data, y.imp$data,C.imp$data )

    pd <- data.frame(Condition = as.factor(c(rep(Cond1, 4), rep(Cond2, 4), rep("control",4))),
                     BioReplicate = as.factor(rep(colnames(x.imp$data),3)))
    rownames(pd) <- colnames(ed)
    pd$sample.id <- rownames(pd)

    ## PCA
    pc <- prcomp(t(ed[apply(is.na(ed),1,sum)==0,]))
    png(paste0(resp2,"pca_",nameCont, ".png"), width=1400, height = 800)
    m <- matrix(c(1,2,3,3),nrow = 2,ncol = 2,byrow = TRUE)
    layout(mat = m,heights = c(0.8,.2))
    s <- 2
    plot(pc$x[,1],pc$x[,s], pch = " ", xlab = "pcomp 1", ylab = paste0("pcomp ", s))
    text(pc$x[,1],pc$x[,s], rownames(pd), col = as.numeric(pd$Condition))
    s <- 3
    plot(pc$x[,1],pc$x[,s], pch = " ", xlab = "pcomp 1", ylab = paste0("pcomp ", s))
    text(pc$x[,1],pc$x[,s], rownames(pd), col = as.numeric(pd$Condition))
    plot(1, type = "n", axes=FALSE, xlab="", ylab="")
    legend("top", inset = 0, legend=c(Cond1,Cond2,"Control"), col = c(1,2,3), pch= c(21,21,21), horiz = TRUE,cex =2)
    dev.off()

    ## Differential expression
    Nall <- N*2
    pd2 <- pd
    pd <- pd[1:Nall,]
    pd$Condition <- droplevels(pd$Condition)
    form <-  paste0("~ -1 + BioReplicate + Condition ", sep='')
    design <- model.matrix(as.formula(form), data=pd)
    lmf <- lmFit(ed[,1:Nall], design)
    betas <- coef(lmf)

    lgroups <- list(
        Cond1_Cond2 = list(list(Condition = Cond1), list( Condition = Cond2))
    );
    lsigns <- list(c(-1, 1))
    names(lsigns) <- names(lgroups);

    conts <- sapply(1:length(lsigns), function(j)
        make.contrasts(d=pd, form=as.formula(form),
                       lsel=lgroups[[j]],
                       signs=lsigns[[j]]), simplify=F)
    cont.mat <- sapply(conts, function(o) o$cont)
    colnames(cont.mat) <- names(lgroups)

    lmc <- contrasts.fit(lmf, cont.mat)
    eb <- eBayes(lmc)
    ds <- topTable(eb, adjust="BH", number=nrow(ed), coef=colnames(cont.mat)[1], confint=T)

    ## Prepare results tables (imputed values marked with *)
    Cond1.miss <- apply(is.na(x[rownames(ds),]),1,sum)
    Cond2.miss <- apply(is.na(y[rownames(ds),]),1,sum)
    Control.miss <- apply(is.na(C[rownames(ds),]),1,sum)
    ad <- array("", dim = dim(ed[,1:Nall]))
    ad[is.na(cbind(x[rownames(ds),],y[rownames(ds),]))] <- "*"
    rownames(ad) <- rownames(ds)

    imp <- array(paste0(round(ed[rownames(ds),1:Nall],3),ad), dim = dim(ed[,1:Nall]))
    colnames(imp) <- colnames(ed)[1:Nall]
    rownames(imp) <- rownames(ds)

    ## Prepare results tables (average values)
    avg.vals <- data.frame(Cond1.mean = apply(ed[rownames(ds),1:N],1,mean),Cond2.mean = apply(ed[rownames(ds),1:N+N],1,mean),
                           Control.mean = apply(ed[rownames(ds),1:N+N*2],1,mean))

    ## Prepare results tables (add all columns)
    rt <- data.frame(Protein = rownames(ds),gene.names = gene.names[rownames(ds)], ds[,1:4], avg.vals, ds[,5:ncol(ds)], imp, Cond1.miss, Cond2.miss, Control.miss)
    rt[, c(3:5)] <- 10^abs(rt[, c(3:5)])*sign(rt[, c(3:5)])
    colnames(rt)[3] <- "FC"
    colnames(rt)[regexpr("Cond1",colnames(rt))>0] <- sub("Cond1", Cond1,colnames(rt)[regexpr("Cond1",colnames(rt))>0])
    colnames(rt)[regexpr("Cond2",colnames(rt))>0] <- sub("Cond2", Cond2,colnames(rt)[regexpr("Cond2",colnames(rt))>0])

    ## Save Rdata and csv
    save(rt,file = paste(resp2, "d", nameCont, ".Rdata", sep=''));
    write.xlsx(rt,file = paste0(resp2, "d", nameCont, ".xlsx"), sep="\t",  row.names=FALSE, quote=FALSE)


    ############
    ### B. subdata with many missings in only one of the two conditions: no imputation.
    ############
    cond2 <- (apply(!is.na(x),1,sum) >= 2 & apply(!is.na(y),1,sum) < 2) | (apply(!is.na(y),1,sum) >= 2 & apply(!is.na(x),1,sum) < 2)
    xe <- (x[cond2,])
    colnames(xe) <-  paste0(colnames(xe), rep(paste0("_",c(Cond1)),each = N))
    ye <- (y[cond2,])
    colnames(ye) <-  paste0(colnames(ye), rep(paste0("_",c(Cond2)),each = N))
    Ce <- (C[cond2,])
    colnames(Ce) <-  paste0(colnames(Ce), rep(paste0("_",c("control")),each = N))
    rt2 <- data.frame(Protein = rownames(xe),gene.names = gene.names[rownames(xe)], xe, ye,
                      Cond1.miss = apply(is.na(xe),1,sum), Cond2.miss = apply(is.na(ye),1,sum), Cond1.avg = mean(ed[,1:N]), Cond2.avg = mean(ed[,1:N+N]))
    colnames(rt2)[regexpr("Cond1",colnames(rt2))>0] <- sub("Cond1", Cond1, colnames(rt2)[regexpr("Cond1",colnames(rt2))>0])
    colnames(rt2)[regexpr("Cond2",colnames(rt2))>0] <- sub("Cond2", Cond2, colnames(rt2)[regexpr("Cond2",colnames(rt2))>0])

    write.xlsx(rt2,file = paste0(resp2, "d", nameCont, "largemiss.xlsx"), sep="\t",  row.names=FALSE, quote=FALSE)

}



#################################################################################
## 3. Missing values treatment and differential expression between condition and control - case study 2591
#################################################################################
Contrasts2591 <- cbind(c("Control", "Control", "Control", "Control"),
                       c("BirA-xCPEB1", "xCPEB1-BirA", "BirA-RBD-ZZ", "RBD-ZZ-BirA"))
study <- 2591

for(ite in seq_len(dim(eval(parse(text = paste0("Contrasts",study))))[1])){
    Cond1 <- eval(parse(text = paste0("Contrasts",study)))[ite,1]
    Cond2 <- eval(parse(text = paste0("Contrasts",study)))[ite,2]
    nameCont <- paste0(gsub("-",".",Cond1), "_", gsub("-",".",Cond2))

    ### create folder for results
    resp2 <- paste0(res,nameCont,"/")
    dir.create(resp2)

    y <- as.matrix(iBACmat[,sampleinfo$sample == Cond2])
    C <- as.matrix(iBACmat[,sampleinfo$sample == "BirA" & sampleinfo$study == study])

    ### understanding missings values
    t2 <- table(apply(y==0,1,sum),apply(C ==0,1,sum))
    t2a <- array(as.numeric(t2),dim = dim(t2))
    rownames(t2a) <- paste0(Cond2, " = ", rownames(t2))
    colnames(t2a) <- paste0("Control", " = ", colnames(t2))

    p <- openPage(paste(resp2, 'missings_table.html', sep=''));
    hwrite(paste('<br>Number of missing values by protein<br><br>'), p);
    for (j in 2)
        {
            ts <- eval(parse(text = c(paste0("t",j,"a"))))
            hwrite(paste('<br><br><br>',"table ", j,'<br><br>'), p);
            hwrite(as.matrix(ts), page=p, center=F, row.names=T, col.names=T,
                   col.style=c("text-align:center", rep("text-align:center", ncol(ts)-1)),
                   col.width=c(100, rep(125, ncol(ts)-1)));
        }
    closePage(p);

    mns <- apply(log10(C),1, function(x) mean(x[is.finite(x)]))
    mss <- apply(C==0, 1, sum)
    mss <- mss[!is.na(mns)]
    mns <- mns[!is.na(mns)]
    png(paste0(resp2, "missings_randomness_control.png"))
    boxplot(mns ~ mss, lwd = 2, ylab = 'avg log10 iBAC', xlab="n missings", main = Cond1)
    stripchart(mns ~ mss, vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = 'blue')
    dev.off()

    mns <- apply(log10(y),1, function(y) mean(y[is.finite(y)]))
    mss <- apply(y==0, 1, sum)
    mss <- mss[!is.na(mns)]
    mns <- mns[!is.na(mns)]
    png(paste0(resp2, "missings_randomness_",Cond2,".png"))
    boxplot(mns ~ mss, lwd = 2, ylab = 'avg log10 iBAC', xlab="n missings", main = Cond2)
    stripchart(mns ~ mss, vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = 'blue')
    dev.off()

    ## frog names
    frog.y <- vapply(strsplit(colnames(y),"_"),function(x) x[5], character(1))
    frog.C <- vapply(strsplit(colnames(C),"_"),function(x) x[5], character(1))
    colnames(y) <- frog.y
    colnames(C) <- frog.C
    nam <- sort(frog.y)

    ## log 10
    y <- log10(y[,nam])
    C <- log10(C[,nam])
    y[!is.finite(y)] <- NA
    C[!is.finite(C)] <- NA

    ## normalization
    N <- length(nam)
    jd <- cbind(y,C)
    colnames(jd) <- paste0(colnames(jd), rep(paste0("_",c(Cond2,"C")),each = N))
    jd.norm <- mdQuantileNormalization(jd)

    png(paste0(resp2,"boxplot_intensities_",nameCont, ".png"), width=1200, height = 800)
    par(mfrow=c(1,2))
    boxplot(jd, xlab = "", ylab = "Intensity", main = "Before normalization", las=2)
    boxplot(jd.norm, xlab = "", ylab = "", main = "After normalization", las=2)
    dev.off()

    y <- jd.norm[,1:N ]
    C <- jd.norm[,1:N+N]


    ############
    ### A. subdata with not many missings: impute missing values with KNN when there is at least 2 non-missing values
    ############

    cond <- apply(!is.na(C),1,sum) > 1 & apply(!is.na(y),1,sum) > 1
    y2 <- y[cond,]
    C2 <- C[cond,]

    y.imp <- impute.knn(y2 ,k = 10, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=23)
    C.imp <- impute.knn(C2 ,k = 10, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=23)

    ## Joint data
    ed <- cbind(y.imp$data, C.imp$data)

    pd <- data.frame(Condition = as.factor(c(rep(Cond2, N), rep("control",N))),
                     BioReplicate = as.factor(rep(colnames(y.imp$data),2)))
    rownames(pd) <- colnames(ed)
    pd$sample.id <- rownames(pd)

    ## PCA
    pc <- prcomp(t(ed))
    png(paste0(resp2,"pca_",nameCont, ".png"), width=1400, height = 800)
    m <- matrix(c(1,2,3,3),nrow = 2,ncol = 2,byrow = TRUE)
    layout(mat = m,heights = c(0.8,.2))
    s <- 2
    plot(pc$x[,1],pc$x[,s], pch = " ", xlab = "pcomp 1", ylab = paste0("pcomp ", s))
    text(pc$x[,1],pc$x[,s], rownames(pd), col = as.numeric(pd$Condition))
    s <- 3
    plot(pc$x[,1],pc$x[,s], pch = " ", xlab = "pcomp 1", ylab = paste0("pcomp ", s))
    text(pc$x[,1],pc$x[,s], rownames(pd), col = as.numeric(pd$Condition))
    plot(1, type = "n", axes=FALSE, xlab="", ylab="")
    legend("top", inset = 0, legend=c(Cond2,"Control"), col = c(1,2), pch= c(21,21), horiz = TRUE,cex =2)
    dev.off()

    ## Differential expression
    form <-  paste0("~ -1 + BioReplicate + Condition ", sep='')
    design <- model.matrix(as.formula(form), data=pd)
    lmf <- lmFit(ed, design)
    betas <- coef(lmf)

    lgroups <- list(
        Cond1_Cond2 = list(list(Condition = "control"), list( Condition = Cond2))
    );
    lsigns <- list(c(-1, 1))
    names(lsigns) <- names(lgroups);

    conts <- sapply(1:length(lsigns), function(j)
        make.contrasts(d=pd, form=as.formula(form),
                       lsel=lgroups[[j]],
                       signs=lsigns[[j]]), simplify=F)
    cont.mat <- sapply(conts, function(o) o$cont)
    colnames(cont.mat) <- names(lgroups)

    lmc <- contrasts.fit(lmf, cont.mat)
    eb <- eBayes(lmc)
    ds <- topTable(eb, adjust="BH", number=nrow(ed), coef=colnames(cont.mat)[1], confint=T)


    ## Prepare results tables (imputed values marked with *)
    Nall <- N * 2
    Control.miss <- apply(is.na(C[rownames(ds),]),1,sum)
    Cond2.miss <- apply(is.na(y[rownames(ds),]),1,sum)
    ad <- array("", dim = dim(ed[,1:Nall]))
    ad[is.na(cbind(y[rownames(ds),],C[rownames(ds),]))] <- "*"
    rownames(ad) <- rownames(ds)

    imp <- array(paste0(round(ed[rownames(ds),1:Nall],3),ad), dim = dim(ed[,1:Nall]))
    colnames(imp) <- colnames(ed)[1:Nall]
    rownames(imp) <- rownames(ds)

    ## Prepare results tables (average values)
    avg.vals <- data.frame(Control.mean = apply(ed[rownames(ds),1:N+N],1,mean),Cond2.mean = apply(ed[rownames(ds),1:N],1,mean))

    rt <- data.frame(Protein = rownames(ds),gene.names = gene.names[rownames(ds)], ds[,1:4], avg.vals, ds[,5:ncol(ds)],imp[,1:N+N], imp[,1:N], Control.miss, Cond2.miss)
    rt[, c(3:5)] <- 10^abs(rt[, c(3:5)])*sign(rt[, c(3:5)])
    colnames(rt)[3] <- "FC"
    colnames(rt)[regexpr("Cond2",colnames(rt))>0] <- sub("Cond2", Cond2,colnames(rt)[regexpr("Cond2",colnames(rt))>0])

    ## Save Rdata and csv
    save(rt,file = paste(resp2, "d", nameCont, ".Rdata", sep=''));
    write.xlsx(rt,file = paste0(resp2, "d", nameCont, ".xlsx"), sep="\t",  row.names=FALSE, quote=FALSE)

    ############
    ### B. subdata with many missings in only one of the two conditions: no imputation.
    ############

    cond2 <- (apply(!is.na(C),1,sum) >= 2 & apply(!is.na(y),1,sum) < 2) | ( apply(!is.na(y),1,sum) >= 2 & apply(!is.na(C),1,sum) < 2 )
    Ce <- (C[cond2,])
    colnames(Ce) <-  paste0(colnames(Ce), rep(paste0("_",c("control")),each = N))
    ye <- (y[cond2,])
    colnames(ye) <-  paste0(colnames(ye), rep(paste0("_",c(Cond2)),each = N))
    rt2 <- data.frame(Protein = rownames(Ce),gene.names = gene.names[rownames(Ce)], Ce, ye, control.miss = apply(is.na(Ce),1,sum),Cond2.miss = apply(is.na(ye),1,sum),
                      control.avg = mean(ed[,1:N+N]), Cond2.avg = mean(ed[,1:N]))
    colnames(rt2)[regexpr("Cond2",colnames(rt2))>0] <- sub("Cond2", Cond2, colnames(rt2)[regexpr("Cond2",colnames(rt2))>0])

    write.xlsx(rt2,file = paste0(resp2, "d", nameCont, "largemiss.xlsx"), sep="\t",  row.names=FALSE, quote=FALSE)
}


#################################################################################
## 4. Missing values treatment and differential expression between condition and control - case study 2591 vs 2833
#################################################################################
studies <- cbind(c(2591,2591,2591,2591,2591,2591),
                 c(2833,2833,2833,2833,2833,2833))

Conds <- cbind(c("BirA-xCPEB1","BirA-xCPEB1","BirA-xCPEB1","xCPEB1-BirA", "xCPEB1-BirA","xCPEB1-BirA"),
               c("BirA-xCPEB2","BirA-xCPEB3","BirA-xCPEB4","xCPEB2-BirA", "xCPEB3-BirA", "xCPEB4-BirA"))

for(LA in 1:dim(studies)[1]){
    Cond1 <- Conds[LA,1]; study1 <- studies[LA,1]
    Cond2 <- Conds[LA,2]; study2 <- studies[LA,2]

    ## create folder for results
    nameCont <- paste0(gsub("-",".",Cond1), "_", gsub("-",".",Cond2))
    resp2 <- paste0(res,nameCont,"/")
    dir.create(resp2)

    x <- as.matrix(iBACmat[,sampleinfo$sample == Cond1])
    y <- as.matrix(iBACmat[,sampleinfo$sample == Cond2])
    Cx <- as.matrix(iBACmat[,sampleinfo$sample == "BirA" & sampleinfo$study == study1])
    Cy <- as.matrix(iBACmat[,sampleinfo$sample == "BirA" & sampleinfo$study == study2])

    t1 <- table(apply(x==0,1,sum), apply(Cx ==0,1,sum))
    t1a <- array(as.numeric(t1), dim = dim(t1))
    rownames(t1a) <- paste0(Cond1, " = ", rownames(t1))
    colnames(t1a) <- paste0("Control", " = ", colnames(t1))

    t2 <- table(apply(y==0,1,sum),apply(Cy ==0,1,sum))
    t2a <- array(as.numeric(t2),dim = dim(t2))
    rownames(t2a) <- paste0(Cond2, " = ", rownames(t2))
    colnames(t2a) <- paste0("Control", " = ", colnames(t2))

    ## frog names
    frog.x <- vapply(strsplit(colnames(x),"_"),function(x) x[5], character(1))
    frog.y <- vapply(strsplit(colnames(y),"_"),function(x) x[5], character(1))
    frog.Cx <- vapply(strsplit(colnames(Cx),"_"),function(x) x[5], character(1))
    frog.Cy <- vapply(strsplit(colnames(Cy),"_"),function(x) x[5], character(1))
    colnames(x) <- frog.x
    colnames(y) <- frog.y
    colnames(Cx) <- frog.Cx
    colnames(Cy) <- frog.Cy
    nam1 <- sort(frog.x)
    nam2 <- sort(frog.y)

    ## log 10
    x <- log10(x[,nam1])
    y <- log10(y[,nam2])
    Cx <- log10(Cx[,nam1])
    Cy <- log10(Cy[,nam2])

    x[!is.finite(x)] <- NA
    y[!is.finite(y)] <- NA
    Cx[!is.finite(Cx)] <- NA
    Cy[!is.finite(Cy)] <- NA

    ## normalization
    N1 <- length(nam1)
    N2 <- length(nam2)
    jd <- cbind(x,y, Cx, Cy)
    jd.norm <- mdQuantileNormalization(jd)

    png(paste0(resp2,"boxplot_intensities_",nameCont, ".png"), width=1200, height = 800)
    par(mfrow=c(1,2))
    boxplot(jd, xlab = "", ylab = "Intensity", main = "Before normalization", las=2)
    boxplot(jd.norm, xlab = "", ylab = "", main = "After normalization", las=2)
    dev.off()

    N <- N1
    x <- jd.norm[,1:N]
    y <- jd.norm[,1:N + N]
    Cx <- jd.norm[,1:N + N*2]
    Cy <- jd.norm[,1:N + N*3]

    ############
    ### A. subdata with not many missings: impute missing values with KNN when there is at least 2 non-missing values
    ############
    cond <- apply(!is.na(x),1,sum) > 1 & apply(!is.na(y),1,sum) > 1 & apply(!is.na(Cx),1,sum) > 1 & apply(!is.na(Cy),1,sum) > 1
    x2 <- x[cond,]
    y2 <- y[cond,]
    Cx2 <- Cx[cond,]
    Cy2 <- Cy[cond,]

    ## Imputing missing values through KNN
    x.imp <- impute.knn(x2 ,k = 10, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=23)
    y.imp <- impute.knn(y2 ,k = 10, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=23)
    Cx.imp <- impute.knn(Cx2 ,k = 10, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=23)
    Cy.imp <- impute.knn(Cy2 ,k = 10, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=23)

    ed <- cbind(x.imp$data, y.imp$data, Cx.imp$data, Cy.imp$data)
    colnames(ed) <- paste0(colnames(ed),c(rep(Cond1, N1), rep(Cond2, N2), rep("Cont1",N1), rep("Cont2",N2)))

    pd <- data.frame(Condition = as.factor(c(rep(Cond1, N1), rep(Cond2, N2), rep("Cont1",N1), rep("Cont2",N2))),
                     BioReplicate = as.factor(rep(c(colnames(x.imp$data),colnames(y.imp$data)),2)))
    rownames(pd) <- colnames(ed)
    pd$sample.id <- rownames(pd)

    ## PCA
    pc <- prcomp(t(ed))

    png(paste0(resp2,"pca_",nameCont, ".png"), width=1400, height = 800)
    m <- matrix(c(1,2,3,3),nrow = 2,ncol = 2,byrow = TRUE)
    layout(mat = m,heights = c(0.8,.2))
    s <- 2
    plot(pc$x[,1],pc$x[,s], pch = " ", xlab = "pcomp 1", ylab = paste0("pcomp ", s))
    text(pc$x[,1],pc$x[,s], rownames(pd), col = as.numeric(pd$Condition))
    s <- 3
    plot(pc$x[,1],pc$x[,s], pch = " ", xlab = "pcomp 1", ylab = paste0("pcomp ", s))
    text(pc$x[,1],pc$x[,s], rownames(pd), col = as.numeric(pd$Condition))
    plot(1, type = "n", axes=FALSE, xlab="", ylab="")
    legend("top", inset = 0, legend=c(Cond1,Cond2, "Control1", "Control2"), col = c(1,2,3,4), pch= c(21,21,21,21), horiz = TRUE,cex =2)
    dev.off()


    ## center to control
    Nall <- N*2
    edn <- cbind(ed[,1:N]-ed[,1:N +Nall],ed[,1:N+N]-ed[,1:N +Nall+N])
    edn <- edn - min(edn)+0.1
    pdn <- pd[1:Nall,]
    pdn$Condition <- droplevels(pdn$Condition)

    ## Differential expression
    form <-  paste0("~ -1  + Condition ", sep='')
    design <- model.matrix(as.formula(form), data=pdn)
    lmf <- lmFit(edn, design)
    betas <- coef(lmf)

    lgroups <- list(
        Cond1_Cont1 = list(list(Condition = Cond1), list( Condition = Cond2))      );
    lsigns <- list(c(-1, 1))
    names(lsigns) <- names(lgroups);

    conts <- sapply(1:length(lsigns), function(j)
        make.contrasts(d=pdn, form=as.formula(form),
                       lsel=lgroups[[j]],
                       signs=lsigns[[j]]), simplify=F)
    cont.mat <- sapply(conts, function(o) o$cont)
    colnames(cont.mat) <- names(lgroups)

    lmc <- contrasts.fit(lmf, cont.mat)
    eb <- eBayes(lmc)
    ds <- topTable(eb, adjust="BH", number=nrow(edn), coef=colnames(cont.mat)[1], confint=T)

    ## Prepare results tables (imputed values marked with *)
    Cond1.miss <- apply(is.na(x[rownames(ds),]),1,sum)
    Cond2.miss <- apply(is.na(y[rownames(ds),]),1,sum)
    Cont1.miss <- apply(is.na(Cx[rownames(ds),]),1,sum)
    Cont2.miss <- apply(is.na(Cy[rownames(ds),]),1,sum)
    ad <- array("", dim = dim(ed[,1:(Nall*2)]))
    ad[is.na(cbind(x[rownames(ds),],y[rownames(ds),],Cx[rownames(ds),],Cy[rownames(ds),]))] <- "*"
    rownames(ad) <- rownames(ds)

    imp <- array(paste0(round(ed[rownames(ds),1:(Nall*2)],3),ad), dim = dim(ed[,1:(Nall*2)]))
    colnames(imp) <- colnames(ed)[1:(Nall*2)]
    rownames(imp) <- rownames(ds)

    ## Prepare results tables (average values)
    avg.vals <- data.frame(Cond1.mean = apply(ed[rownames(ds),1:N],1,mean), Cond2.mean = apply(ed[rownames(ds),1:N + N],1,mean),
                           Cont1.mean = apply(ed[rownames(ds),1:N*2],1,mean), Cont2.mean = apply(ed[rownames(ds),1:N + N*3],1,mean))

    rt <- data.frame(Protein = rownames(ds), gene.names = gene.names[rownames(ds)], ds[,1:4],avg.vals, ds[,5:ncol(ds)], imp, Cond1.miss, Cond2.miss, Cont1.miss, Cont2.miss)
    rt[, c(3:5)] <- 10^abs(rt[, c(3:5)])*sign(rt[, c(3:5)])
    colnames(rt)[3] <- "FC"
    colnames(rt)[regexpr("Cond1",colnames(rt))>0] <- sub("Cond1", Cond1,colnames(rt)[regexpr("Cond1",colnames(rt))>0])
    colnames(rt)[regexpr("Cond2",colnames(rt))>0] <- sub("Cond2", Cond2,colnames(rt)[regexpr("Cond2",colnames(rt))>0])

    ## Save Rdata and csv
    save(rt,file = paste(resp2, "d", nameCont, ".Rdata", sep=''));
    write.xlsx(rt,file = paste0(resp2, "d", nameCont, ".xlsx"), sep="\t",  row.names=FALSE, quote=FALSE)

    ############
    ### B. subdata with many missings in only one of the two conditions: no imputation.
    ############
    conde <- !cond
    xe <-  cbind(x,y,Cx, Cy)[conde,]
    colnames(xe) <- paste0(c(rep(Cond1,N),rep(Cond2,4), rep("Cont1",4), rep("Cont2",4)), "_", colnames(xe))
    rt2 <- data.frame(Protein = rownames(xe),gene.names = gene.names[rownames(xe)], xe, Cond1.miss = apply(is.na(xe[,1:N]),1,sum), Cond2.miss = apply(is.na(xe[,1:N+N]),1,sum),
                      Cont1.miss = apply(is.na(xe[,1:N + 2*N]),1,sum), Cont2.miss = apply(is.na(xe[,1:N + 3*N]),1,sum),
                      Cond1.avg = rep(mean(xe[,1:N],na.rm=TRUE),dim(xe)[1]), Cond2.avg = rep(mean(xe[,1:N+N],na.rm=TRUE),dim(xe)[1]),
                      Cont1.avg = rep(mean(xe[,1:N +2*N],na.rm=TRUE),dim(xe)[1]),
                      Cont2.avg = rep(mean(xe[,1:N+3*N],na.rm=TRUE),dim(xe)[1]))
    colnames(rt2)[regexpr("Cond1",colnames(rt2))>0] <- sub("Cond1", Cond1,colnames(rt2)[regexpr("Cond1",colnames(rt2))>0])
    colnames(rt2)[regexpr("Cond2",colnames(rt2))>0] <- sub("Cond2", Cond2,colnames(rt2)[regexpr("Cond2",colnames(rt2))>0])

    write.xlsx(rt2,file = paste0(resp2, "d", nameCont, "largemiss.xlsx"), sep="\t",  row.names=FALSE, quote=FALSE)

}

