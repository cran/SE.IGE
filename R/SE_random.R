SE_random <-
  function(N,n,vpd,vps,h2d,h2s,rgds,reds,nw,r) {
    vad <- h2d*vpd
    vas <- h2s*vps
    ved <- vpd-vad
    ves <- vps-vas
    cads <- rgds*sqrt(vad*vas)
    ceds <- reds*sqrt(ved*ves)
    cpds <- cads+ceds
    vp <- vpd+(nw-1.0)*vps
    vtbv <- vad+2*(nw-1.0)*cads+(nw-1.0)**2*vas
    covw <- 2*cpds+(nw-2.0)*vps
    ## DGE
    vary <- vp
    varf <- r*vad
    vvarf <- vvarb(vary,varf,N,n)
    SEvad <- sqrt(vvarf)/r
    ## IGE
    vary <- (vp+(nw-2)*covw)/(nw-1)
    varf <- r*vas
    vvarf <- vvarb(vary,varf,N,n)
    SEvas <- sqrt(vvarf)/r
    ## TBV
    vary <- nw*(vp + (nw-1)*covw)
    varf <- r*vtbv
    vvarf <- vvarb(vary,varf,N,n)
    SEvtbv <- sqrt(vvarf)/r
    ## T2
    SET2 <- SEvtbv/vp
    ## COVARIANCE
    SEcads <- sqrt(0.5*(SEvad*SEvas)+cads^2/(N-1))
    ## CORRELATION
    SErgds <- SE_corr(vad,vas,cads,SEvad,SEvas,SEcads)
    x <- c(SEvad,SEvas,SEvtbv,SET2,SEcads,SErgds)
    names(x) <- c("SEvad","SEvas","SEvtbv","SET2","SEcads","SErgds")
    return(x)
  }

