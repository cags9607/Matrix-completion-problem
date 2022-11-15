fill.correlations = function(A, replace=FALSE, eps=1e-10) {
  library("CVXR")   # A library for Disciplined Convex Programming
  
  # Parameters:
  
  # A -------- A partially defined symmetric matrix with NA's in the missed entries
  # replace --- A logic value that indicates whether the output matrix should be corrected
  # eps -------- An epsilon value to consider as 0 all entries whose absolute value is less than eps
  
  # Outputs:
  
  # X0 ------- A correlation matrix with maximum determinant and entries x0ij = aij for specified entries aij of A.
  
  n = dim(A)[1]
  X = Variable(n, n, PSD = T) # this is the matrix that we are looking for
  constraints = list(diag(X) == 1)
  for (i in 1:n){
    for (j in 1:n){
      if (j < i) { # We update the constraints. We have one constraint for each specified correlation
        if (!is.na(A[i,j])) constraints = c(constraints,X[i,j] == A[i,j])
      }
    }
  }
  objective = Maximize(log_det(X)) # This is a monotone increasing function of det and a concave function (so -log_det is convex)
  # This ensures that the problem is a SDP, as det is not a convex nor concave function of the entries of X
  # However, maximizing log_det is the same as maximing det. 
  
  program = Problem(objective,constraints)
  result = solve(program)
  status = result$status
  X0 = as.matrix(result$getValue(X))
  diag(X0) = 1
  if (replace==FALSE) {   # To recover the original values of the specified entries
    X0[!is.na(A)] = A[!is.na(A)]  
    X0[abs(X0) < eps] = 0  
  }
  if (status == "optimal") print("Solver succesfully converged!!!") else print(paste("Warning!!! (The solution is ",status,")"))
  return (list(Sigma_maxdet = X0, Sigma = A, status = status))
}
is.PD = function(X, eps = 1e-12) {
  X[is.na(X)] = 0
  try(min(eigen(X, symmetric = T, only.values = T)$values) > eps, F)
}
is.PSD =  function(X, eps = 1e-12) {
  X[is.na(X)] = 0
  try(min(eigen(X, symmetric = T, only.values = T)$values) > -eps, F)
}

find.lambda = function(eps, Sigma_maxdet, Sigma_guess){
  d = ncol(Sigma_maxdet)
  if(!is.PD(Sigma_guess)){
    alpha  = uniroot(function(alpha,R,Rfixed,tol=eps){
      if(anyNA(R)){
        R[is.na(R)] = 0
      }
      f = alpha*R + (1-alpha)*Rfixed
      min(eigen(f)$values) - tol
    },c(0,1),R=Sigma_guess,Rfixed=Sigma_maxdet,tol = eps)$root
  } else {alpha = 1}
  alpha
}

find.lambda.max = function(lambdas, Sigma_maxdet, Sigma_guess){
  pds = sapply(lambdas, function(lambda) is.PD(lambda*Sigma_guess + (1-lambda)*Sigma_maxdet))
  max(lambdas[pds==1])
}

randcorr_sylv = function(N = 1000, Sigma, G = seq(-1,1,by = 0.00001), eps = 1e-12, verbose = T){
  matrices = list()
  for (i in 1:N) {
    Sigma_c = Sigma
    num_elem = dim(Sigma_c)[1]^2
    size = row_num = col_num_ind = 3; total_size = dim(Sigma_c)[1]
    here = Sigma_c
    
    while (size <= total_size) {
      here = Sigma_c[(1:size), (1:size)] # take the upper left size x size corner
      here[is.na(here)] = sample(x = G, size = length(here[is.na(here)]), replace = TRUE)
      here[upper.tri(here)] = t(here)[upper.tri(here)]
      
      # Let's perform rejection sampling
      
      while (det(here) < 0){
        here = Sigma_c[(1:size), (1:size)]  
        here[is.na(here)] = sample(x = G, size = length(here[is.na(here)]), replace = TRUE)
        here[upper.tri(here)] = t(here)[upper.tri(here)]
      }
      
      Sigma_c[1:row_num, 1:col_num_ind] = here 
      row_num = row_num + 1
      col_num_ind = col_num_ind + 1
      size = size + 1
    }
    
    Sigma_c = here
    matrices[[i]] = Sigma_c
    if (verbose) cat(paste(i*100/N, "% of progress", "\n"))
    
  }
  return(matrices)
}

randcorr_maxdet_maxlambda = function(N, Sigma, G = seq(-1,1,by = 0.00001), eps = 1e-12,
                                verbose = T){
  Sigma_maxdet = fill.correlations(Sigma)
  status = Sigma_maxdet$status
  if (status != "optimal") return(warning("Cannot sample from Lebesgue 0 measure sets!!!"))
  Sigma_maxdet = Sigma_maxdet$Sigma_maxdet
  n = dim(Sigma)[1]
  matrices = list()
  lambdas = list()
  for (k in 1:N) {  
    Sigma_rand = Sigma_new = Sigma
    Sigma_rand[is.na(Sigma_rand)] = sample(x = G, size = length(Sigma_rand[is.na(Sigma_rand)]), replace = TRUE)
    Sigma_rand[upper.tri(Sigma_rand)] = t(Sigma_rand)[upper.tri(Sigma_rand)]
    Sigma_new = Sigma_rand
    
    if(!is.PD(Sigma_new)) {lambda_max = find.lambda(0.00001, Sigma_maxdet, Sigma_new)} else {lambda_max = 1}
    lambdas[[k]] = lambda_max
    lambda = runif(1,0,lambda_max)
    entries = lambda*Sigma_rand + (1-lambda)*Sigma_maxdet
    Sigma_new[is.na(Sigma)] = entries[is.na(Sigma)]
    matrices[[k]] = Sigma_new
    if (verbose) cat(paste(k*100/N, "% of progress", "\n"))
  }
  return(list(matrices = matrices, lambdas = lambdas))
}
randcorr_maxdet_unif = function(N, Sigma, G = seq(-1,1,by = 0.00001), eps = 1e-12,
                                     verbose = T){
  Sigma_maxdet = fill.correlations(Sigma)
  status = Sigma_maxdet$status
  if (status != "optimal") return(warning("Cannot sample from Lebesgue 0 measure sets!!!"))
  Sigma_maxdet = Sigma_maxdet$Sigma_maxdet
  n = dim(Sigma)[1]
  matrices = list()
  lambdas = list()
  for (k in 1:N) {  
    lambda = 1
    Sigma_rand = Sigma_new = Sigma
    Sigma_rand[is.na(Sigma_rand)] = sample(x = G, size = length(Sigma_rand[is.na(Sigma_rand)]), replace = TRUE)
    Sigma_rand[upper.tri(Sigma_rand)] = t(Sigma_rand)[upper.tri(Sigma_rand)]
    Sigma_new = Sigma_rand
    while(!is.PD(Sigma_new)){
      lambda = runif(1,0,lambda)
      entries = lambda*Sigma_rand + (1-lambda)*Sigma_maxdet
      Sigma_new[is.na(Sigma)] = entries[is.na(Sigma)]
  }
    matrices[[k]] = Sigma_new
    lambdas[[k]] = lambda
    if (verbose) cat(paste(k*100/N, "% of progress", "\n"))
  }
  return(list(matrices = matrices, lambdas = lambdas))
}

randcorr = function(N, Sigma, G = seq(-1,1,by = 0.00001), verbose = T){
  matrices = list()
  for (i in 1:N){
    Sigma_new = Sigma
    repeat{
      Sigma_new[is.na(Sigma)] = sample(x = G, size = length(Sigma[is.na(Sigma)]), replace = TRUE)
      Sigma_new[upper.tri(Sigma_new)] = t(Sigma_new)[upper.tri(Sigma_new)]
      if(is.PSD(Sigma_new)){break}
    }
    matrices[[i]] = Sigma_new
    if (verbose) cat(paste(i*100/N, "% of progress", "\n"))
    
  }
  matrices
}
