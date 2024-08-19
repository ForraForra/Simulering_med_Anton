###################################BYGGA KERNEL

# Funktion för att välja kernelfunktion.
kernel_typ <- function(kernel_type) {
  switch(
    kernel_type,
    "gaussian" = function(u) dnorm(u, mean = 0, sd = 1),
    "kvadratisk" = function(u) 0.75 * (1 - u^2),
    "uniform" = function(u)  0.5,
    stop("uniform,gaussian o epa funkar ba")
  )
}

# Gasser-Muller estimatorn.
kernel_smooth <- function(x, y, kernel = "gaussian", bandwidth = 0.5) {
  # x är domänet och y funktionsvärdet, n är antalet punkter och kernel_typ 
  # specificerar den valda kernel:en utifrån funktionen kernel_typ ovan. 
  n <- length(x)
  kernel_function <- kernel_typ(kernel)
  
  kernel_vikt_S <- 0 #tom vektor för de slutgiltiga resultaten. 
  
  
  # I loopen nedan är "  for (i in 1:n) " för varje punkt där funktionen ska skattas. 
  # "  for (j in 1:n) " är observationerna som ska användas för att skatta i punkterna i.
  
  for (i in 1:n) {
    resultat <- 0 # Sparar vikterna från varje observation för skattningen i punkten i.
    for (j in 1:n) {
      # Övre undre integrationsgräns 
      ti_litet <- ifelse(j == 1, x[1], (x[j] + x[j - 1]) / 2) 
      ti_stort <- ifelse(j == n, x[n], (x[j + 1] + x[j]) / 2) 
      # göra Gasser-Muller för vikterna
      if (kernel == "gaussian"){
        resultat[j] <- integrate(function(u) kernel_function((u - x[i]) / bandwidth), ti_litet, ti_stort)$value / bandwidth
      } else {
        if(abs(x[j]-x[i])<=bandwidth){
          if (kernel == "uniform"){
            resultat[j] <- 0.5*(ti_stort-ti_litet) / bandwidth
          } else {
            resultat[j] <- integrate(function(u) kernel_function((u - x[i]) / bandwidth), ti_litet, ti_stort)$value / bandwidth
          }
        } else {
          resultat[j] <- 0
        }
      }
      
    }
    
    # Vikterna multipliceras med de observerade värdena. 
    kernel_vikt_S[i] <- sum(resultat * y)
  }
  # returnerar de skattade funktionsvärdena över x.
  return(kernel_vikt_S)
}
# Estimator för kernel över 300 punkter i intervallet.
kernel_smooth_MISE <- function(x, x_out, y, kernel = "gaussian", bandwidth = 0.5) {
  # Precis som innan fast nu är det 300 punkter (n1) av intresse x_out och datasetet 
  # x har bara n2 antal punkter.
  n1 <- length(x_out)
  n2 <- length(x)
  kernel_function <- kernel_typ(kernel)
  
  kernel_vikt_S <- 0
  
  
  # Skillnaden i den här loopen från tidigare är att n2 observationer används för 
  # att skatta funktionsvärdet i 300 punkter (n1).
  for (i in 1:n1) {
    resultat <- 0
    for (j in 1:n2) {
      
      ti_litet <- ifelse(j == 1, x[1], (x[j] + x[j - 1]) / 2) 
      ti_stort <- ifelse(j == n2, x[n2], (x[j + 1] + x[j]) / 2) 
      
      if (kernel == "gaussian"){
        resultat[j] <- integrate(function(u) kernel_function((u - x_out[i]) / bandwidth), ti_litet, ti_stort)$value / bandwidth
      } else{
        if (abs(x[j]-x_out[i])<=bandwidth){
          if (kernel == "uniform"){
            resultat[j] <- 0.5*(ti_stort-ti_litet) / bandwidth
          } else {
            resultat[j] <- integrate(function(u) kernel_function((u - x_out[i]) / bandwidth), ti_litet, ti_stort)$value / bandwidth
          }
        } else {
          resultat[j] <- 0
        }
      }
    }
    
    kernel_vikt_S[i] <- sum(resultat * y)
  }
  
  return(kernel_vikt_S)
}

# funktion för att skatta värdet i en enskild punkt som inte är en del av x-vektorn.
# Används för att genomföra LOOCV
kernel_pred <- function(x_eval, x, y, kernel = "gaussian", bandwidth){
  # Allt görs som tidigare bara att den enda punkten av intresse är x_eval så 
  # ingen loop för i behöver göras utan det blir bara en skattad punkt som returneras
  # av funktionen.
  n <- length(x)
  kernel_function <- kernel_typ(kernel)
  
  resultat <- 0
  for (j in 1:n) {
    ti_litet <- ifelse(j == 1, x[1], (x[j] + x[j - 1]) / 2) 
    ti_stort <- ifelse(j == n, x[n], (x[j + 1] + x[j]) / 2) 
    # 
    if (kernel == "gaussian"){
      resultat[j] <- integrate(function(u) kernel_function((u - x_eval) / bandwidth), ti_litet, ti_stort)$value / bandwidth
    } else{
      if (abs(x[j]-x_eval)<=bandwidth){
        if (kernel == "uniform"){
          resultat[j] <- 0.5*(ti_stort-ti_litet) / bandwidth
        } else {
          resultat[j] <- integrate(function(u) kernel_function((u - x_eval) / bandwidth), ti_litet, ti_stort)$value / bandwidth
        }
      } else {
        resultat[j] <- 0
      }
    }
  }
  return(sum(resultat*y))
}

# Skatta GCV för kernels
kernel_GCV <- function(kernel_smooth, x, y, kernel = "gaussian", bandwidth = 0.5){
  # SSE beräknas för kernel_smooth objektet av intresse
  SSE <- sum((kernel_smooth - y)^2)
  N <- length(kernel_smooth)
  
  # Smoothingmatrix skapas genom att kolumn för kolumn beräkna hur stor påverkan 
  # varje observation har på skattningen i punkten x[j].
  
  L <- matrix(nrow = length(X), ncol = length(X))
  for (j in 1:length(x)) {
    yi <- rep_len(0, length(x))
    yi[j] <- 1
    L[, j] <- kernel_smooth(x, yi, kernel = kernel, bandwidth = bandwidth)
  }
  # De effektiva frihetsgraderna är summan av diagonalen av smoothing matrix.
  vv <- sum(diag(L))
  
  # Beräkna GCV med Wassermans formel.
  return((SSE / (1 - vv / N)^2) / N)
}

# Skatta GCV för basfunktionsmetoderna.
basis_GCV <- function(x, y, nbasis, type.basis = "fourier"){
  # Samma sak som innan förutom att andra funktioner används för att skatta
  # funktionsvärden för SSE och när smoothing matrix ska skapas.
  
  xyfdata1 <- fdata(y,x)
  basis_obj <- optim.basis(xyfdata1, numbasis = nbasis, type.basis = type.basis)
  SSE <- sum((basis_obj$fdata.est$data - y)^2)
  N <- length(x)
  
  
  L <- matrix(nrow = length(x), ncol = length(x))
  for (j in 1:length(x)) {
    yi <- rep_len(0, length(x))
    yi[j] <- 1
    xyfdata2 <- fdata(yi, x)
    L[, j] <- optim.basis(xyfdata2, numbasis = nbasis, type.basis = type.basis)$fdata.est$data
  }
  vv <- sum(diag(L))
  
  
  
  return((SSE / (1 - vv / N)^2) / N)
}

# Skatta GCV för smoothing spline och fourier
smooth_GCV <- function(x, y, nbasis, lambda, type.basis = "fourier"){
  if (type.basis == "fourier"){
    basisobjekt <- create.fourier.basis(c(0,100), nbasis)
    fdParobjekt <- fdPar(basisobjekt, 2, lambda)
    fdobjekt <- smooth.basis(x, y, fdParobjekt)$fd$coefs
    
    SSE <- sum((y - eval.basis(x, basisobjekt)%*%fdobjekt)^2)
    N <- length(x)
    
    L <- matrix(nrow = length(x), ncol = length(x))
    for (j in 1:length(x)) {
      yi <- rep_len(0, length(x))
      yi[j] <- 1
      L[, j] <- eval.basis(x, basisobjekt)%*%smooth.basis(X, yi, fdParobjekt)$fd$coefs
    }
  }else {
    basisobjekt <- create.bspline.basis(c(0,100), nbasis, norder = 4, x)
    fdParobjekt <- fdPar(basisobjekt, 2, lambda)
    fdobjekt <- smooth.basis(x, y, fdParobjekt)$fd$coefs
    
    SSE <- sum((y - eval.basis(x, basisobjekt)%*%fdobjekt)^2)
    N <- length(x)
    
    L <- matrix(nrow = length(x), ncol = length(x))
    for (j in 1:length(x)) {
      yi <- rep_len(0, length(x))
      yi[j] <- 1
      L[, j] <- eval.basis(x, basisobjekt)%*%smooth.basis(x, yi, fdParobjekt)$fd$coefs
    }
  }
  vv <- sum(diag(L))
  
  return((SSE / (1 - vv / N)^2) / N)
}