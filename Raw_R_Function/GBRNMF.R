###NMF work
#x: the non-negative data (if any negatives exist, the whole data is shifted)

#q: the number of factors

#eps: convergence criteria (lack of progress on sum squared error)

#maxit: convergence criteria (max iteration when lack of progress not met)

#w: n by q loadings on factors (amount of chemical/molecule/etc)

#a: q by q scaling auxiliary matrix

#s: q by p observed factors (basis chemical/molecule/etc spectra)

#By default (when u, h NULL) both u and h are initialized randomly

library(zoo)

cnmf <- function(x, q, eps=0.0000001, maxit=2000, w=NULL, a=NULL, s=matrix(nrow = 0,ncol = 0), xax=1:ncol(x)){
  
  n <- nrow(x)
  p <- ncol(x)
  
  known_factors <- nrow(as.matrix(s))

  if(any(x<0)){
    x <- as.matrix(x)+abs(min(x))
  }
  else{
    x <- as.matrix(x)
  }
  
  if(any(s<0)){
    s <- as.matrix(s) + abs(min(s))
  }
  else{
    s <- as.matrix(s)
  }
  
  if(is.null(w)){
    known_groupings <- 0
    w <- matrix(runif(n*q, min(x), max(x)), n, q)
    
    if(is.null(a)){
      a <- diag(nrow = q)
    }
    
    if(q > known_factors){
      extra <- matrix(runif((q-nrow(s))*p, min(s), max(s)), q-nrow(s), p)
      s <- s_save <- rbind(s,extra)
    }
    else{
      q = known_factors
      s_save <- s
    }
  }
  else{
    known_groupings <- ncol(w)
    addition <- matrix(runif((q-ncol(w))*n, min(x), max(x)), n, q-ncol(w))
    w <- cbind(w,addition)
    
    if(is.null(a)){
      a <- diag(nrow = q)
      #a <- matrix(1, nrow = q, ncol = q)
      # for (i in 1:q){
      #   a[i,i] <- 1
      # }
    }
    
    if(q > known_factors){
      if(known_factors>0){
        extra1 <- matrix(runif(known_groupings*p, min(s), max(s)), known_groupings, p)
        extra2 <- matrix(runif((q-known_groupings-nrow(s))*p, min(s), max(s)), q-known_groupings-nrow(s), p)
        s <- s_save <- rbind(extra1,s,extra2)
      }
      else{
        s <- s_save <- matrix(runif(q*p, min(x), max(x)), q, p)
      }
    }
    else{q = known_factors+known_groupings}
  }
  
  ed <- vector()
  ed[1] <- sum((x-w%*%a%*%s)^2)
  conv <- FALSE
  ctr <- 1
  
  while(!conv){
    ctr <- ctr+1
    
    a <- a * (t(w) %*% x %*% t(s)) / (t(w) %*% w %*% t(a) %*% s %*% t(s))
    
    
    if(known_groupings > 0){
      temp_w <- w * (x %*% t(s) %*% t(a)) / (w %*% a %*% s %*% t(s) %*% t(a))
      for(i in (known_groupings+1):q){
        w[,i] <- temp_w[,i]
      }
      if(q > known_factors + known_groupings){
        temp <- s * (t(a) %*% t(w) %*% x) / (t(a) %*% t(w) %*% w %*% a %*% s)
        for(i in 1:known_groupings){
          s[i,] <- temp[i,]
        }
        for(i in (known_groupings+known_factors+1):q){
          s[i,] <- temp[i,]
        }
      }
      else if(q== known_factors + known_groupings){
        temp <- s * (t(a) %*% t(w) %*% x) / (t(a) %*% t(w) %*% w %*% a %*% s)
        for(i in 1:known_groupings){
          s[i,] <- temp[i,]
        }
      }
    }
    else{
      w <- w * (x %*% t(s) %*% t(a)) / (w %*% a %*% s %*% t(s) %*% t(a))
      if(q > known_factors){
        temp <- s * (t(a) %*% t(w) %*% x) / (t(a) %*% t(w) %*% w %*% a %*% s)
        for(i in (known_factors+1):q){
          s[i,] <- temp[i,]
        }
      }
    }
  
    wh <- w%*%a%*%s
    ed[ctr] <- sum((x-wh)^2)
    
    if(ctr >= 3 & ((ed[ctr-1]-ed[ctr] < eps)|(ctr==maxit))){
      conv <- TRUE
    }
  }
  
  column_sum_w <- colSums(w)
  scaled_w_col <- diag(nrow(w)/column_sum_w)
  scaled_w_col_inv <- diag(column_sum_w/nrow(w))
  
  w <- w%*%scaled_w_col
  a <- scaled_w_col_inv %*% a
  
  for(id in 1:nrow(s)){
    AUC <- sum(diff((xax))*rollmean(s[id,],2))
    s[id,] <- s[id,] / AUC
    a[id,id] <- a[id,id]*AUC
  }
  
  list(x = x, ed = ed, w = w, a = a, s = s, reconstructedx = wh, s_init = s_save)
  
}
