nopt_random <-
  function(T,vpd,vps,h2d,h2s,rgds,reds,nw,r) {
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
    t <- varf/vary
    nopt_d <- 1/t
    md <- T/nopt_d
    vvarf <- vvarb(vary,varf,md,nopt_d)
    SEvad <- sqrt(vvarf)/r
    ## IGE
    vary <- (vp+(nw-2)*covw)/(nw-1)
    varf <- r*vas
    t <- varf/vary
    nopt_s <- 1/t
    ms <- T/nopt_s
    vvarf <- vvarb(vary,varf,ms,nopt_s)
    SEvas <- sqrt(vvarf)/r
    ## TBV
    vary <- nw*(vp + (nw-1)*covw)
    varf <- r*vtbv
    t <- varf/vary
    nopt_t <- 1/t
    mt <- T/nopt_t
    vvarf <- vvarb(vary,varf,mt,nopt_t)
    SEvtbv <- sqrt(vvarf)/r
    ## T2
    SET2 <- SEvtbv/vp
    return(matrix(c(md,nopt_d,SEvad,ms,nopt_s,SEvas,mt,nopt_t,SEvtbv,
                    mt,nopt_t,SET2),ncol=3,nrow=4,byrow=TRUE, 
                    dimnames = list(c("direct","indirect","total","T2"),
                                    c("#families", "fam size", "SE"))   ))
}

