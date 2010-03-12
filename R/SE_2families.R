SE_2families <-
  function(N,n,vpd,vps,h2d,h2s,rgds,reds,nw,r) {
    vad <- h2d*vpd
    vas <- h2s*vps
    ved <- vpd-vad
    ves <- vps-vas
    cads <- rgds*sqrt(vad*vas)
    ceds <- reds*sqrt(ved*ves)
    cpds <- cads+ceds
    vp <- vpd+(nw-1)*vps+2*r*(0.5*nw-1)*(cads+(0.5*nw-1)*vas)
    vtbv <- vad+2*(nw-1)*cads+(nw-1.0)**2*vas
    covwf <- 2*cpds+(nw-2)*vps + 
      r*( vad + 2*(0.5*nw-2)*cads + (0.5*nw**2-2*nw+3)*vas )
    covwnf <- 2*cpds+(nw-2)*vps + 
      2*r*(0.5*nw-1)*(cads+(0.5*nw-1)*vas)
    vpavg <- (vp+(0.5*nw-1)*covwf)/(0.5*nw)
    ngf <- n/(0.5*nw)		
    ## DGE
    phi <- (0.5*nw-1)/(0.5*nw)
    vary <- (1+phi**2)*vpavg - 2*phi*covwnf
    varf <- r*vad
    vvarf <- vvarb(vary,varf,N,ngf)
    SEvad <- sqrt(vvarf)/r
    ## IGE
    vary <- vpavg/(0.25*nw**2)
    varf <- r*vas
    vvarf <- vvarb(vary,varf,N,ngf)
    SEvas <- sqrt(vvarf)/r
    ## TBV
    vary <- (4/nw)*(vp + (0.5*nw-1)*covwf + 0.5*nw*covwnf)
    varf <- r*vtbv
    vvarf <- vvarb(vary,varf,N,ngf)
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

