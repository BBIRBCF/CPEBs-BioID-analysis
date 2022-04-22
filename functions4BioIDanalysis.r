#### Quantile normalization
mdQuantileNormalization <- function(x){
    x.nm <- apply(x,2,function(x) x[!is.na(x)])
    percentiles <- sapply(x.nm, function(s) quantile(s,c(0,seq_len(100)/100)))
    avg.percentiles <- apply(percentiles,1,mean)
    y.nm <- sapply(1:length(x.nm), function(o)  approx(percentiles[,o], avg.percentiles, x.nm[[o]])$y)
    for(o in seq_len(length(x.nm))) x[names(x.nm[[o]]),o] <- y.nm[[o]]
    return(x)
}

formatApval <- function(x){
    st = "ns"
    st <- ifelse(x < 0.25, "+", st)
    st <- ifelse(x < 0.10, "*", st)
    st <- ifelse(x < 0.05, "**", st)
    st <- ifelse(x < 0.01, "***", st)
    if(is.na(st)) st <- "nc"
    return(st)
}


#### make contrasts for limma modelling
make.contrasts <- function(d, form, lsel, signs, type='balanced', ngroups=NULL)
{

    nv <- strsplit(paste(form)[2], split=' \\+ ')[[1]];
    nv <- nv[regexpr("\\:", nv) < 0];
    nv <- nv[nv!='-1']
    for (n in nv) if (is.character(d[, n])) d[, n] <- factor(d[, n])

    nvf <- nv[sapply(nv, function(v, d) is.factor(d[, v]), d)];
    nvn <- nv[!nv%in%nvf];

    daf <- expand.grid(sapply(nvf, function(o) levels(d[, o]), simplify=F));
    dan <- matrix(rep(apply(d[, nvn, drop=F], 2, mean), each=nrow(daf)), nrow=nrow(daf))
    colnames(dan) <- nvn;
    colnames(daf) <- nvf;
    dal <- cbind(daf, dan)[nv];

    dsg <- model.matrix(form, dal);

    lss <- sapply(lsel, function(sel, dal)
              {
                  ss <- sapply(1:length(sel), function(j, sel, dal)
                           {
                               nm <- names(sel)[[j]];
                               if (is.factor(dal[, nm])) paste("'", sel[[j]], "'", sep="");
                           }, sel, dal, simplify=F);
                  names(ss) <- names(sel);
                  ss;
              }, dal, simplify=F);

    ds0 <- sapply(lss, function(ss, d)
             {
                 eval(parse(text=paste(sapply(1:length(ss), function(j, ss, d)
                        {
                            o <- ss[[j]];
                            names(o) <- names(ss)[j];
                            paste("d$", names(ss)[j], "%in%c(",
                                  paste(o, collapse=', ', sep=''), ")", sep='');
                        },  ss, dal), collapse=' & ')))
             }, d)

    if (any(apply(ds0, 1, sum) > 1)) stop("Error: Overlapping groups!");

    if (type=='balanced')
    {
        ds <- sapply(lss, function(ss, dal)
                 {
                     eval(parse(text=paste(sapply(1:length(ss), function(j, ss, dal)
                            {
                                o <- ss[[j]];
                                names(o) <- names(ss)[j];
                                paste("dal$", names(ss)[j], "%in%c(",
                                      paste(o, collapse=', ', sep=''), ")", sep='');
                            },  ss, dal), collapse=' & ')))
                 }, dal)

        grs <- t(apply(ds, 2, function(s, dsg) apply(dsg[s, , drop=F], 2, mean), dsg));
        cont <- (t(grs)%*%signs)[, 1]

    }else
    {
        dsg0 <- model.matrix(form, d);

        nvm <- colnames(dsg)[regexpr("\\:", colnames(dsg)) < 0];
        lnve <- sapply(ldsg, function(ds, nvm)
                   {
                       nvm[apply(ds[, nvm, drop=F], 2, function(o) length(unique(o)) == 1)];
                   }, nvm);
        lnvd <- sapply(ldsg, function(ds, nvm)
                   {
                       nvm[apply(ds[, nvm, drop=F], 2, function(o) length(unique(o)) > 1)];
                   }, nvm);
        lmed.nvd <- sapply(lnvd, function(v, dsg0) apply(dsg0[, v, drop=F], 2, mean), dsg0);

        lp <- sapply(1:length(med.nvd), function(j, nvm)
                 {
                     p <- data.frame(t(c(unique(ldsg[[j]][, lnve[[j]]]), lmed.nvd[[j]])));
                     colnames(p) <- c(lnve[[j]], lnvd[[j]]);
                     p[, nvm, drop=F];
                 }, nvm);

        grs <- t(apply(ds, 2, function(s, dsg) apply(dsg[s, , drop=F], 2, mean), dsg));
        cont <- (t(grs)%*%signs)[, 1]

    }

    if (!is.null(ngroups)) rownames(grs) <- colnames(ds0) <- ngroups;
    list(grs=grs, cont=cont, sels=ds0);

}

