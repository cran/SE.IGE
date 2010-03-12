SE_corr <-
function(v1,v2,c12,SEv1,SEv2,SEc12) {
         vsd1sd2 <- (SEv1^2/(4*v1)+v1)*(SEv2^2/(4*v2)+v2) - v1*v2
         return(sqrt( SEc12^2/(v1*v2) + c12^2*vsd1sd2/(v1^2*v2^2) ))
 }

