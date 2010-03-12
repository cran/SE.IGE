nopt_2families <-
  function(T,vpd,vps,h2d,h2s,rgds,reds,nw,r) {
                                        #basic parameters
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
    ## Direct effect
    phi <- (0.5*nw-1)/(0.5*nw)
    vary <- (1+phi**2)*vpavg - 2*phi*covwnf
    varf <- r*vad
    t=varf/vary
    nopt_d=1/t					#note: this is optimum ngf
    md <- T/(nopt_d*0.5*nw)
    vvarf <- vvarb(vary,varf,md,nopt_d)
    SEvad <- sqrt(vvarf)/r
    ## SOCIAL EFFECT
    vary <- vpavg/(0.25*nw**2)
    varf <- r*vas
    t=varf/vary
    nopt_s <- 1/t					#note: this is optimum ngf
    ms <- T/(nopt_s*0.5*nw)
    vvarf <- vvarb(vary,varf,ms,nopt_s)
    SEvas <- sqrt(vvarf)/r
    ## TBV
    vary <- (4/nw)*(vp + (0.5*nw-1)*covwf + 0.5*nw*covwnf)
    varf <- r*vtbv
    t <- varf/vary
    nopt_t <- 1/t					#note: this is optimum ngf
    mt=T/(nopt_t*0.5*nw)
    vvarf <- vvarb(vary,varf,mt,nopt_t)
    SEvtbv <- sqrt(vvarf)/r
    ## T2
    SET2 <- SEvtbv/vp
    return(matrix(c(md,nopt_d*0.5*nw,SEvad,ms,nopt_s*0.5*nw,SEvas,
                    mt,nopt_t*0.5*nw,SEvtbv,mt,nopt_t*0.5*nw,SET2),
                    ncol=3,nrow=4,byrow=TRUE,
                    dimnames = list(c("direct","indirect","total","T2"),
                                    c("#families", "fam size", "SE")) ))
  }

