vvarb <-
  function(var,varb,N,n) {
    vare <- var-varb
    return((2/(N-1))*(varb^2 + 2*varb*vare/n + vare^2/(n*(n-1))))
  }

