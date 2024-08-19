library(fda)
library(fda.usc)
#Skapa alla vektorer
case_1_Fourier_GCV <- 0 # Vektor för GCV beroende på antal basfunktioner.
case_1_Fourier_GCV_MSE <- 0 # Fourier MSE-mätning med GCV.
case_1_Fourier_GCV_MISE_hat <- 0 # Fourier MISE_hat-mätning med GCV.
case_1_Fourier_GCV_time <- 0 # Tidvektor för fourier med GCV.
case_1_Fourier_GCV_basis <- 0 # Fourier, antal basfunktioner med GCV.

case_1_Fourier_CV_pred <- 0 # Fourier, För att hjälpa med LOOCV hyperparameterskattning.
case_1_Fourier_CV_hyper_MSE <- 0 # Fourier, LOOCV-MSE skattning beroende på antalet basfunktioner.
case_1_Fourier_CV_minbasis <- 0 # Fourier, antalet basfunktioner med LOOCV.
case_1_Fourier_CV_MSE <- 0 # Fourier MSE-mätning med LOOCV.
case_1_Fourier_CV_MISE_hat <- 0 # Fourier MISE_hat-mätning med LOOCV.
case_1_Fourier_CV_time <- 0 # Tidvektor för fourier med LOOCV.

case_1_spline_GCV <- 0 # Vektor för GCV beroende på antal basfunktioner
case_1_Spline_GCV_MSE <- 0 # Spline MSE-mätning med GCV.
case_1_Spline_GCV_MISE_hat <- 0 # Spline MISE_hat-mätning med GCV.
case_1_Spline_GCV_time <- 0 # Tidvektor för Spline med GCV.
case_1_Spline_GCV_basis <- 0 # Spline, antal basfunktioner med GCV.

case_1_Spline_CV_pred <- 0 # Spline, För att hjälpa med LOOCV hyperparameterskattning.
case_1_Spline_CV_hyper_MSE <- 0 # Spline, LOOCV-MSE skattning beroende på antalet basfunktioner.
case_1_Spline_CV_minbasis <- 0 # Spline, antalet basfunktioner med LOOCV.
case_1_Spline_CV_MSE <- 0 # Spline MSE-mätning med LOOCV.
case_1_Spline_CV_MISE_hat <- 0 # Spline MISE_hat-mätning med LOOCV.
case_1_Spline_CV_time <- 0 # Tidvektor för Spline med LOOCV.

case_1_Fourier_smooth_GCV <- 0 # hjälpvektor för val av lambda för GCV.
case_1_Fourier_smooth_GCV_lambda <- 0 # Det värde på lambda som ger lägst GCV.
case_1_Fourier_smooth_GCV_MSE <- 0 # MSE-skattning för GCV.
case_1_Fourier_smooth_GCV_MISE_hat <- 0 # MISE_hat skattning med GCV.
case_1_Fourier_smooth_GCV_time <- 0 # tidvektor för GCV.

case_1_Fourier_smooth_CV_pred <- 0 # hjälpvektor för val av lambda med LOOCV
case_1_Fourier_smooth_CV_hyper_MSE <- 0 # samma som ovan
case_1_Fourier_smooth_CV_minlambda <- 0 # Det värde på lambda som ger lägst LOOCV.
case_1_Fourier_smooth_CV_MSE <- 0 # MSE-skattning för LOOCV.
case_1_Fourier_smooth_CV_MISE_hat <- 0 # MISE_hat skattning med LOOCV.
case_1_Fourier_smooth_CV_time <- 0 # tidvektor för LOOCV.

case_1_Spline_smooth_GCV <- 0 # hjälpvektor för val av lambda för GCV.
case_1_Spline_smooth_GCV_lambda <- 0 # Det värde på lambda som ger lägst GCV.
case_1_Spline_smooth_GCV_MSE <- 0 # MSE-skattning för GCV.
case_1_Spline_smooth_GCV_MISE_hat <- 0 # MISE_hat skattning med GCV.
case_1_Spline_smooth_GCV_time <- 0 # tidvektor för GCV.

case_1_Spline_smooth_CV_pred <- 0 # hjälpvektor för val av lambda med LOOCV
case_1_Spline_smooth_CV_hyper_MSE <- 0 # samma som ovan
case_1_Spline_smooth_CV_minlambda <- 0 # Det värde på lambda som ger lägst LOOCV.
case_1_Spline_smooth_CV_MSE <- 0 # MSE-skattning för LOOCV.
case_1_Spline_smooth_CV_MISE_hat <- 0 # MISE_hat skattning med LOOCV.
case_1_Spline_smooth_CV_time <- 0 # tidvektor för LOOCV.

case_1_Unif_Kernel_GCV_GCV <- 0 # Hjälpvektor när bästa h ska väljas med GCV
case_1_Unif_Kernel_GCV_minh <- 0 # Det h som ger lägst GCV
case_1_Unif_Kernel_GCV_MSE <- 0 # Unifkernel MSE-mätning med GCV
case_1_Unif_Kernel_GCV_MISE_hat <- 0 # Unifkernel MISE_hat-mätning med GCV
case_1_Unif_Kernel_GCV_time <- 0 # tidvektor för Unifkernel med GCV

case_1_Unif_Kernel_CV_pred <- 0 # Hjälpvektor när bästa h ska väljas med LOOCV
case_1_Unif_Kernel_CV_hyper_MSE <- 0 # Vektor med LOOCV-MSE beroende på h
case_1_Unif_Kernel_CV_minh <- 0 # Det h som ger lägst LOOCV-MSE
case_1_Unif_Kernel_CV_MSE <- 0 # Unifkernel MSE-mätning med LOOCV
case_1_Unif_Kernel_CV_MISE_hat <- 0 # Unifkernel MISE_hat-mätning med LOOCV
case_1_Unif_Kernel_CV_time <- 0 # Tidverktor för Unifkernel med LOOCV

case_1_Quad_Kernel_GCV_GCV <- 0 # Hjälpvektor när bästa h ska väljas med GCV
case_1_Quad_Kernel_GCV_minh <- 0 # Det h som ger lägst GCV
case_1_Quad_Kernel_GCV_MSE <- 0 # Quad-kernel MSE-mätning med GCV
case_1_Quad_Kernel_GCV_MISE_hat <- 0 # Quad-kernel MISE_hat-mätning med GCV
case_1_Quad_Kernel_GCV_time <- 0 # tidvektor för Quad-kernel med GCV

case_1_Quad_Kernel_CV_pred <- 0 # Hjälpvektor när bästa h ska väljas med LOOCV
case_1_Quad_Kernel_CV_hyper_MSE <- 0 # Vektor med LOOCV-MSE beroende på h
case_1_Quad_Kernel_CV_minh <- 0 # Det h som ger lägst LOOCV-MSE
case_1_Quad_Kernel_CV_MSE <- 0 # Quad-kernel MSE-mätning med LOOCV
case_1_Quad_Kernel_CV_MISE_hat <- 0 # Quad-kernel MISE_hat-mätning med LOOCV
case_1_Quad_Kernel_CV_time <- 0 # Tidverktor för Quad-kernel med LOOCV

case_1_Normal_Kernel_GCV_GCV <- 0 # Hjälpvektor när bästa h ska väljas med GCV
case_1_Normal_Kernel_GCV_minh <- 0 # Det h som ger lägst GCV
case_1_Normal_Kernel_GCV_MSE <- 0 # Normalkernel MSE-mätning med GCV
case_1_Normal_Kernel_GCV_MISE_hat <- 0 # Normalkernel MISE_hat-mätning med GCV
case_1_Normal_Kernel_GCV_time <- 0 # tidvektor för Normalkernel med GCV

case_1_Normal_Kernel_CV_pred <- 0 # Hjälpvektor när bästa h ska väljas med LOOCV
case_1_Normal_Kernel_CV_hyper_MSE <- 0 # Vektor med LOOCV-MSE beroende på h
case_1_Normal_Kernel_CV_minh <- 0 # Det h som ger lägst LOOCV-MSE
case_1_Normal_Kernel_CV_MSE <- 0 # Normalkernel MSE-mätning med LOOCV
case_1_Normal_Kernel_CV_MISE_hat <- 0 # Normalkernel MISE_hat-mät
case_1_Normal_Kernel_CV_time <- 0 #tidvektor för normalkernel med LOOCV


X_exact <- seq(0, 100, length.out = 300)

eps_funktion <- function(X, sd_value) {
  n <- length(X)
  eps <- rnorm(n, 0, sd_value)
  return(eps)
}
X_5 <- seq(0,100,length.out=5)
X_25 <- seq(0,100, length.out=25)
fx1_5 <- sin(0.02*pi*X_5)
fx2_5 <- sin(0.02*pi*X_5) + cos(0.02*pi*X_5) + sin(0.08*pi*X_5) + cos(0.08*pi*X_5) +
  sin(0.14*pi*X_5) + cos(0.14*pi*X_5) + sin(0.2*pi*X_5) + cos(0.2*pi*X_5)
fx1_25 <- sin(0.02*pi*X_25)
fx2_25 <- sin(0.02*pi*X_25) + cos(0.02*pi*X_25) + sin(0.08*pi*X_25) + cos(0.08*pi*X_25) +
  sin(0.14*pi*X_25) + cos(0.14*pi*X_25) + sin(0.2*pi*X_25) + cos(0.2*pi*X_25)

cases <- list(
  case1 <- list(X = X_5, fx = fx1_5, fx_exact = sin(0.02*pi*X_exact), 
                eps_funk = function() eps_funktion(X_5, 0.5)),
  case2 <- list(X = X_5, fx = fx1_5, fx_exact = sin(0.02*pi*X_exact),
                eps_funk = function() eps_funktion(X_5, 0.2)),
  case3 <- list(X = X_5, fx = fx2_5, fx_exact = sin(0.02*pi*X_exact) + cos(0.02*pi*X_exact) + 
                  sin(0.08*pi*X_exact) + cos(0.08*pi*X_exact) + sin(0.14*pi*X_exact) + 
                  cos(0.14*pi*X_exact) + sin(0.2*pi*X_exact) + cos(0.2*pi*X_exact),
                eps_funk = function() eps_funktion(X_5, 0.5)),
  case4 <- list(X = X_5, fx = fx2_5, fx_exact = sin(0.02*pi*X_exact) + cos(0.02*pi*X_exact) + 
                  sin(0.08*pi*X_exact) + cos(0.08*pi*X_exact) + sin(0.14*pi*X_exact) + 
                  cos(0.14*pi*X_exact) + sin(0.2*pi*X_exact) + cos(0.2*pi*X_exact),
                eps_funk = function() eps_funktion(X_5, 0.2)),
  case5 <- list(X = X_25, fx = fx1_25, fx_exact = sin(0.02*pi*X_exact),
                eps_funk = function() eps_funktion(X_25, 0.5)),
  case6 <- list(X = X_25, fx = fx1_25, fx_exact = sin(0.02*pi*X_exact),
                eps_funk = function() eps_funktion(X_25, 0.2)),
  case7 <- list(X = X_25, fx = fx2_25, fx_exact = sin(0.02*pi*X_exact) + cos(0.02*pi*X_exact) + 
                  sin(0.08*pi*X_exact) + cos(0.08*pi*X_exact) + sin(0.14*pi*X_exact) + 
                  cos(0.14*pi*X_exact) + sin(0.2*pi*X_exact) + cos(0.2*pi*X_exact),
                eps_funk = function() eps_funktion(X_25, 0.5)),
  case8 <- list(X = X_25, fx = fx2_25, fx_exact = sin(0.02*pi*X_exact) + cos(0.02*pi*X_exact) + 
                  sin(0.08*pi*X_exact) + cos(0.08*pi*X_exact) + sin(0.14*pi*X_exact) + 
                  cos(0.14*pi*X_exact) + sin(0.2*pi*X_exact) + cos(0.2*pi*X_exact),
                eps_funk = function() eps_funktion(X_25, 0.2))
)
resultat3 <- list() # Denna ska spara de vi vill "för varje case"
for (case_index in seq_along(cases)) {
  nuvarande_case <- cases[[case_index]]
  fx <- nuvarande_case$fx
  X <- nuvarande_case$X
  fx_exact <- nuvarande_case$fx_exact
  
  
  if (length(X) == 5){
    K_fourier <- c(1, 3)
    K_spline <- c(4)
    h_kernel <- c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
  } else{
    K_fourier <- c(1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23)
    K_spline <- c(4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)
    h_kernel <- c(1,3,5,7,9,11,13,15,17,19,21,23)
  }
  
  lambda_values_spline <- seq(0.01, 0.15, 0.01) 
  lambda_values_fourier <- seq(0.1, 1.5, 0.1)
  
  
  for (i in 1:500) {
    eps <- nuvarande_case$eps_funk()
    Y <- fx + eps
    XY <- fdata(Y, X)
    
    ########## FOURIER ##########
    
    # Välj antal basfunktioner med hjälp av GCV.
    start.time <- Sys.time()
    
    for (K in K_fourier) {
      case_1_Fourier_GCV[K] <- basis_GCV(X, Y, nbasis = K, type.basis = "fourier")
    }
    case_1_Fourier_GCV_basis[i] <- which.min(case_1_Fourier_GCV)
    fourier_GCV <- optim.basis(XY, numbasis = case_1_Fourier_GCV_basis[i], type.basis = "fourier")
    case_1_Fourier_GCV_MSE[i] <- mean((Y-fourier_GCV$fdata.est$data)^2)
    
    end.time <- Sys.time()
    
    # Göra 'exakt' skattning av MISE för GCV.
    fourier_GCV_MISE_basis <- create.fourier.basis(c(0, 100), nbasis = fourier_GCV$base.opt$nbasis)
    fourier_GCV_MISE <- smooth.basis(X, Y, fourier_GCV_MISE_basis)
    
    case_1_Fourier_GCV_MISE_pred <- eval.basis(X_exact, fourier_GCV_MISE_basis)%*%fourier_GCV_MISE$fd$coefs
    case_1_Fourier_GCV_MISE_hat[i] <- mean((fx_exact-case_1_Fourier_GCV_MISE_pred)^2)
    
    
    case_1_Fourier_GCV_time[i] <- round(end.time - start.time, 10)
    
    
    
    
    # Välja hyperparametrar med hjälp av LOOCV.
    start.time <- Sys.time()
    
    for (nbasis in K_fourier) {
      fourier_CV_basis <- create.fourier.basis(c(0, 100), nbasis = nbasis)
      for (l in 1:length(X)) {
        X_CV <- X[-l]
        Y_CV <- Y[-l]
        fourier_CV <- smooth.basis(X_CV, Y_CV, fourier_CV_basis)
        case_1_Fourier_CV_pred[l] <- eval.basis(X[l], fourier_CV_basis)%*%fourier_CV$fd$coefs
      }
      case_1_Fourier_CV_hyper_MSE[nbasis] <- mean((Y-case_1_Fourier_CV_pred)^2)
    }
    case_1_Fourier_CV_minbasis[i] <- which.min(case_1_Fourier_CV_hyper_MSE)
    fourier_CV <- optim.basis(XY, numbasis = case_1_Fourier_CV_minbasis[i], type.basis = "fourier")
    case_1_Fourier_CV_MSE[i] <- mean((Y-fourier_CV$fdata.est$data)^2)
    
    end.time <- Sys.time()
    
    # Göra 'exakt' skattning av MISE för LOOCV.
    fourier_CV_MISE_basis <- create.fourier.basis(c(0, 100), nbasis = case_1_Fourier_CV_minbasis[i])
    fourier_CV_MISE <- smooth.basis(X, Y, fourier_CV_MISE_basis)
    
    case_1_Fourier_CV_MISE_pred <- eval.basis(X_exact, fourier_CV_MISE_basis)%*%fourier_CV_MISE$fd$coefs
    case_1_Fourier_CV_MISE_hat[i] <- mean((fx_exact-case_1_Fourier_CV_MISE_pred)^2)
    
    
    
    case_1_Fourier_CV_time[i] <- round(end.time - start.time, 10)
    
    
    
    
    
    
    
    
    ########## SPLINE ##########
    # Välj antal basfunktioner med hjälp av GCV.
    start.time <- Sys.time()
    
    for (K in K_spline) {
      case_1_spline_GCV[K] <- basis_GCV(X, Y, nbasis = K, type.basis = "bspline")
    }
    case_1_spline_GCV[1] <- NA # Väl medveten om att det finns bättre sätt att hantera
    # detta på men det här är lättast.
    case_1_Spline_GCV_basis[i] <- which.min(case_1_spline_GCV)
    Spline_GCV <- optim.basis(XY, numbasis = case_1_Spline_GCV_basis[i], type.basis = "bspline")
    case_1_Spline_GCV_MSE[i] <- mean((Y-Spline_GCV$fdata.est$data)^2)
    
    end.time <- Sys.time()
    
    # Göra 'exakt' skattning av MISE för GCV.
    Spline_GCV_MISE_basis <- create.bspline.basis(c(0, 100), nbasis = Spline_GCV$base.opt$nbasis)
    Spline_GCV_MISE <- smooth.basis(X, Y, Spline_GCV_MISE_basis)
    
    case_1_Spline_GCV_MISE_pred <- eval.basis(X_exact, Spline_GCV_MISE_basis)%*%Spline_GCV_MISE$fd$coefs
    case_1_Spline_GCV_MISE_hat[i] <- mean((fx_exact-case_1_Spline_GCV_MISE_pred)^2)
    
    
    case_1_Spline_GCV_time[i] <- round(end.time - start.time, 10)
    
    
    
    # Välja hyperparametrar med hjälp av LOOCV.
    start.time <- Sys.time()
    
    for (nbasis in K_spline) {
      Spline_CV_basis <- create.bspline.basis(c(0, 100), nbasis = nbasis)
      for (l in 1:length(X)) {
        X_CV <- X[-l]
        Y_CV <- Y[-l]
        Spline_CV <- smooth.basis(X_CV, Y_CV, Spline_CV_basis)
        case_1_Spline_CV_pred[l] <- eval.basis(X[l], Spline_CV_basis)%*%Spline_CV$fd$coefs
      }
      case_1_Spline_CV_hyper_MSE[nbasis] <- mean((Y-case_1_Spline_CV_pred)^2)
    }
    case_1_Spline_CV_hyper_MSE[1] <- NA # Återigen dålig hantering
    case_1_Spline_CV_minbasis[i] <- which.min(case_1_Spline_CV_hyper_MSE)
    Spline_CV <- optim.basis(XY, numbasis = case_1_Spline_CV_minbasis[i], type.basis = "bspline")
    case_1_Spline_CV_MSE[i] <- mean((Y-Spline_CV$fdata.est$data)^2)
    
    end.time <- Sys.time()
    
    # Göra 'exakt' skattning av MISE för LOOCV.
    Spline_CV_MISE_basis <- create.bspline.basis(c(0, 100), nbasis = case_1_Spline_CV_minbasis[i])
    Spline_CV_MISE <- smooth.basis(X, Y, Spline_CV_MISE_basis)
    
    case_1_Spline_CV_MISE_pred <- eval.basis(X_exact, Spline_CV_MISE_basis)%*%Spline_CV_MISE$fd$coefs
    case_1_Spline_CV_MISE_hat[i] <- mean((fx_exact-case_1_Spline_CV_MISE_pred)^2)
    
    
    
    case_1_Spline_CV_time[i] <- round(end.time - start.time, 10)
    
    
    
    
    
    
    
    
    
    
    ########## Smooth Fourier ##########
    
    # Välj antal basfunktioner med hjälp av GCV.
    start.time <- Sys.time()
    
    nbasis <- length(X) - 1
    
    
    for (lambdaindex in 1:length(lambda_values_fourier)) {
      lambda <- lambda_values_fourier[lambdaindex]
      case_1_Fourier_smooth_GCV[lambdaindex] <- smooth_GCV(X, Y, nbasis, lambda, type.basis = "fourier")
    }
    case_1_lambdaindex <- which.min(case_1_Fourier_smooth_GCV)
    case_1_Fourier_smooth_GCV_lambda[i] <- lambda_values_fourier[case_1_lambdaindex]
    basisobjekt <- create.fourier.basis(c(0,100), nbasis)
    fdParobjekt <- fdPar(basisobjekt, 2, case_1_Fourier_smooth_GCV_lambda[i])
    fdobjekt <- smooth.basis(X, Y, fdParobjekt)$fd$coefs
    
    case_1_Fourier_smooth_GCV_MSE[i] <- mean((Y-eval.basis(X, basisobjekt)%*%fdobjekt)^2)
    
    end.time <- Sys.time()
    
    # Göra 'exakt' skattning av MISE för GCV.
    
    case_1_Fourier_smooth_GCV_MISE_pred <- eval.basis(X_exact, basisobjekt)%*%fdobjekt
    
    case_1_Fourier_smooth_GCV_MISE_hat[i] <- mean((fx_exact-case_1_Fourier_smooth_GCV_MISE_pred)^2)
    
    case_1_Fourier_smooth_GCV_time[i] <- round(end.time - start.time, 10)
    
    
    
    
    # Välja hyperparametrar med hjälp av LOOCV.
    start.time <- Sys.time()
    
    for (lambdaindex in 1:length(lambda_values_fourier)) {
      lambda <- lambda_values_fourier[lambdaindex]
      fdParobjekt <- fdPar(basisobjekt, 2, lambda)
      for (l in 1:length(X)) {
        X_CV <- X[-l]
        Y_CV <- Y[-l]
        fdobjekt <- smooth.basis(X_CV, Y_CV, fdParobjekt)$fd
        case_1_Fourier_smooth_CV_pred[l] <- eval.basis(X[l], basisobjekt)%*%fdobjekt$coefs
      }
      case_1_Fourier_smooth_CV_hyper_MSE[lambdaindex] <- mean((Y-case_1_Fourier_smooth_CV_pred)^2)
    }
    case_1_lambdaindex <- which.min(case_1_Fourier_smooth_CV_hyper_MSE)
    case_1_Fourier_smooth_CV_minlambda[i] <- lambda_values_fourier[case_1_lambdaindex]
    
    fdParobjekt <- fdPar(basisobjekt, 2, case_1_Fourier_smooth_CV_minlambda[i])
    fdobjekt <- smooth.basis(X, Y, fdParobjekt)$fd$coefs
    
    case_1_Fourier_smooth_CV_MSE[i] <- mean((Y-eval.basis(X, basisobjekt)%*%fdobjekt)^2)
    
    end.time <- Sys.time()
    
    # Göra 'exakt' skattning av MISE för LOOCV.
    
    case_1_Fourier_smooth_CV_MISE_hat[i] <- mean((fx_exact-eval.basis(X_exact, basisobjekt)%*%fdobjekt)^2)
    
    
    case_1_Fourier_smooth_CV_time[i] <- round(end.time - start.time, 10)
    
    
    
    
    
    
    
    
    
    ########## Smooth Spline ##########
    
    # Välj antal basfunktioner med hjälp av GCV.
    start.time <- Sys.time()
    
    nbasis <- length(X) + 2
    
    
    for (lambdaindex in 1:length(lambda_values_spline)) {
      lambda <- lambda_values_spline[lambdaindex]
      case_1_Spline_smooth_GCV[lambdaindex] <- smooth_GCV(X, Y, nbasis, lambda, type.basis = "bspline")
    }
    case_1_lambdaindex <- which.min(case_1_Spline_smooth_GCV)
    case_1_Spline_smooth_GCV_lambda[i] <- lambda_values_spline[case_1_lambdaindex]
    basisobjekt <- create.bspline.basis(c(0,100), nbasis, norder = 4, X)
    fdParobjekt <- fdPar(basisobjekt, 2, case_1_Spline_smooth_GCV_lambda[i])
    fdobjekt <- smooth.basis(X, Y, fdParobjekt)$fd$coefs
    
    case_1_Spline_smooth_GCV_MSE[i] <- mean((Y-eval.basis(X, basisobjekt)%*%fdobjekt)^2)
    
    end.time <- Sys.time()
    
    # Göra 'exakt' skattning av MISE för GCV.
    
    case_1_Spline_smooth_GCV_MISE_pred <- eval.basis(X_exact, basisobjekt)%*%fdobjekt
    
    case_1_Spline_smooth_GCV_MISE_hat[i] <- mean((fx_exact-case_1_Spline_smooth_GCV_MISE_pred)^2)
    
    case_1_Spline_smooth_GCV_time[i] <- round(end.time - start.time, 10)
    
    
    
    
    # Välja hyperparametrar med hjälp av LOOCV.
    start.time <- Sys.time()
    
    for (lambdaindex in 1:length(lambda_values_spline)) {
      lambda <- lambda_values_spline[lambdaindex]
      fdParobjekt <- fdPar(basisobjekt, 2, lambda)
      for (l in 1:length(X)) {
        X_CV <- X[-l]
        Y_CV <- Y[-l]
        fdobjekt <- smooth.basis(X_CV, Y_CV, fdParobjekt)$fd$coefs
        case_1_Fourier_smooth_CV_pred[l] <- eval.basis(X[l], basisobjekt)%*%fdobjekt
      }
      case_1_Fourier_smooth_CV_hyper_MSE[lambdaindex] <- mean((Y-case_1_Fourier_smooth_CV_pred)^2)
    }
    case_1_lambdaindex <- which.min(case_1_Fourier_smooth_CV_hyper_MSE)
    case_1_Fourier_smooth_CV_minlambda[i] <- lambda_values_spline[case_1_lambdaindex]
    
    fdParobjekt <- fdPar(basisobjekt, 2, case_1_Fourier_smooth_CV_minlambda[i])
    fdobjekt <- smooth.basis(X, Y, fdParobjekt)$fd$coefs
    
    case_1_Spline_smooth_CV_MSE[i] <- mean((Y-eval.basis(X, basisobjekt)%*%fdobjekt)^2)
    
    end.time <- Sys.time()
    
    # Göra 'exakt' skattning av MISE för LOOCV.
    
    case_1_Spline_smooth_CV_MISE_hat[i] <- mean((fx_exact-eval.basis(X_exact, basisobjekt)%*%fdobjekt)^2)
    
    
    case_1_Spline_smooth_CV_time[i] <- round(end.time - start.time, 10)
    
    
    
    
    
    
    
    
    
    
    ############################### KERNEL UNIF 
    
    
    for (hindex in 1:length(h_kernel)) {
      h <- h_kernel[hindex]
      Unif_Kernel_GCV <- kernel_smooth(X, Y, kernel = "uniform", bandwidth = h)
      case_1_Unif_Kernel_GCV_GCV[hindex] <- kernel_GCV(Unif_Kernel_GCV, X, Y, kernel = "uniform", bandwidth = h)
    }
    case_1_unif_Kernel_GCV_index <- which.min(case_1_Unif_Kernel_GCV_GCV)
    case_1_Unif_Kernel_GCV_minh[i] <- h_kernel[case_1_unif_Kernel_GCV_index]
    Unif_Kernel_GCV <- kernel_smooth(X, Y, kernel = "uniform", bandwidth = case_1_Unif_Kernel_GCV_minh[i])
    
    
    case_1_Unif_Kernel_GCV_MSE[i] <- mean((Y-Unif_Kernel_GCV)^2)
    
    end.time <- Sys.time()
    
    # Göra 'exakt' skattning av MISE för GCV.
    Unif_Kernel_GCV_MISE <- kernel_smooth_MISE(X, X_exact, Y, kernel = "uniform", bandwidth = case_1_Unif_Kernel_GCV_minh[i])
    
    case_1_Unif_Kernel_GCV_MISE_hat[i] <- mean((fx_exact-Unif_Kernel_GCV_MISE)^2)
    
    
    case_1_Unif_Kernel_GCV_time[i] <- round(end.time - start.time, 10)
    
    
    
    # Välja hyperparametrar med hjälp av LOOCV.
    start.time <- Sys.time()
    
    for (hindex in 1:length(h_kernel)) {
      h <- h_kernel[hindex]
      for (l in 1:length(X)) {
        X_CV <- X[-l]
        Y_CV <- Y[-l]
        case_1_Unif_Kernel_CV_pred[l] <- kernel_pred(X[l], X_CV, Y_CV, kernel = "uniform", bandwidth = h)
      }
      case_1_Unif_Kernel_CV_hyper_MSE[hindex] <- mean((Y-case_1_Unif_Kernel_CV_pred)^2)
    }
    case_1_Unif_Kernel_CV_index <- which.min(case_1_Unif_Kernel_CV_hyper_MSE)
    case_1_Unif_Kernel_CV_minh[i] <- h_kernel[case_1_Unif_Kernel_CV_index]
    case_1_Unif_Kernel_CV <- kernel_smooth(X, Y, kernel = "uniform", bandwidth = case_1_Unif_Kernel_CV_minh[i])
    case_1_Unif_Kernel_CV_MSE[i] <- mean((Y-case_1_Unif_Kernel_CV)^2)
    
    end.time <- Sys.time()
    
    # Göra 'exakt' skattning av MISE för LOOCV.
    Unif_Kernel_CV_MISE <- kernel_smooth_MISE(X, X_exact, Y, kernel = "uniform", bandwidth = case_1_Unif_Kernel_CV_minh[i])
    
    case_1_Unif_Kernel_CV_MISE_hat[i] <- mean((fx_exact-Unif_Kernel_CV_MISE)^2)
    
    case_1_Unif_Kernel_CV_time[i] <- round(end.time - start.time, 10)
    
    
    
    
    
    
    
    
    ########## KERNEL QUAD ########## 
    
    
    
    ## GCV:s MSE & MAE & MISE
    start.time <- Sys.time()
    
    for (hindex in 1:length(h_kernel)) {
      h <- h_kernel[hindex]
      Quad_Kernel_GCV <- kernel_smooth(X, Y, kernel = "kvadratisk", bandwidth = h)
      case_1_Quad_Kernel_GCV_GCV[hindex] <- kernel_GCV(Quad_Kernel_GCV, X, Y, kernel = "kvadratisk", bandwidth = h) 
    }
    case_1_Quad_Kernel_GCV_index <- which.min(case_1_Quad_Kernel_GCV_GCV)
    case_1_Quad_Kernel_GCV_minh[i] <- h_kernel[case_1_Quad_Kernel_GCV_index]
    Quad_Kernel_GCV <- kernel_smooth(X, Y, kernel = "kvadratisk", bandwidth = case_1_Quad_Kernel_GCV_minh[i])
    
    case_1_Quad_Kernel_GCV_MSE[i] <- mean((Y-Quad_Kernel_GCV)^2)
    
    end.time <- Sys.time()
    
    
    # MISE för GCV.
    Quad_Kernel_GCV_MISE <- kernel_smooth_MISE(X, X_exact, Y, kernel = "kvadratisk", bandwidth = case_1_Quad_Kernel_GCV_minh[i])
    
    case_1_Quad_Kernel_GCV_MISE_hat[i] <- mean((fx_exact-Quad_Kernel_GCV_MISE)^2)
    
    case_1_Quad_Kernel_GCV_time[i] <- round(end.time - start.time, 10)
    
    
    
    # LOOCV Quad
    
    start.time <- Sys.time()
    for (hindex in 1:length(h_kernel)) {
      h <- h_kernel[hindex]
      for (l in 1:length(X)) {
        X_CV <- X[-l]
        Y_CV <- Y[-l]
        case_1_Quad_Kernel_CV_pred[l] <- kernel_pred(X[l], X_CV, Y_CV, kernel = "kvadratisk", bandwidth = h)
      }
      case_1_Quad_Kernel_CV_hyper_MSE[hindex] <- mean((Y-case_1_Quad_Kernel_CV_pred)^2)
    }
    case_1_Quad_Kernel_CV_index <- which.min(case_1_Quad_Kernel_CV_hyper_MSE)
    case_1_Quad_Kernel_CV_minh[i] <- h_kernel[case_1_Quad_Kernel_CV_index]
    Quad_Kernel_CV <- kernel_smooth(X, Y, kernel = "kvadratisk", bandwidth = case_1_Quad_Kernel_CV_minh[i])
    case_1_Quad_Kernel_CV_MSE[i] <- mean((Y-Quad_Kernel_CV)^2)
    
    end.time <- Sys.time()
    
    # Göra 'exakt' skattning av MISE för LOOCV.
    Quad_Kernel_CV_MISE <- kernel_smooth_MISE(X, X_exact, Y, kernel = "kvadratisk", bandwidth = case_1_Quad_Kernel_CV_minh[i])
    
    case_1_Quad_Kernel_CV_MISE_hat[i] <- mean((fx_exact-Quad_Kernel_CV_MISE)^2)
    
    
    case_1_Quad_Kernel_CV_time[i] <- round(end.time - start.time, 10)
    
    
    
    
    
    ########## KERNEL NORMAL ##########
    
    
    # GCV NORMAL
    
    start.time <- Sys.time()
    for (hindex in 1:length(h_kernel)) {
      h <- h_kernel[hindex]
      Normal_Kernel_GCV <- kernel_smooth(X, Y, kernel = "gaussian", bandwidth = h)
      case_1_Normal_Kernel_GCV_GCV[hindex] <- kernel_GCV(Normal_Kernel_GCV, X, Y, kernel = "gaussian", bandwidth = h)
    }
    case_1_Normal_Kernel_GCV_index <- which.min(case_1_Normal_Kernel_GCV_GCV)
    case_1_Normal_Kernel_GCV_minh[i] <- h_kernel[case_1_Normal_Kernel_GCV_index]
    Normal_Kernel_GCV <- kernel_smooth(X, Y, kernel = "gaussian", bandwidth = case_1_Normal_Kernel_GCV_minh[i])
    
    
    case_1_Normal_Kernel_GCV_MSE[i] <- mean((Y-Normal_Kernel_GCV)^2)
    
    end.time <- Sys.time()
    
    # Göra 'exakt' skattning av MISE för GCV.
    Normal_Kernel_GCV_MISE <- kernel_smooth_MISE(X, X_exact, Y, kernel = "gaussian", bandwidth = case_1_Normal_Kernel_GCV_minh[i])
    
    case_1_Normal_Kernel_GCV_MISE_hat[i] <- mean((fx_exact-Normal_Kernel_GCV_MISE)^2)
    
    
    case_1_Normal_Kernel_GCV_time[i] <- round(end.time - start.time, 10)
    
    
    
    
    
    # LOOCV Normal
    
    start.time <- Sys.time()
    
    for (hindex in 1:length(h_kernel)) {
      h <- h_kernel[hindex]
      for (l in 1:length(X)) {
        X_CV <- X[-l]
        Y_CV <- Y[-l]
        Normal_Kernel_CV <- kernel_smooth(X_CV, Y_CV, kernel = "gaussian", bandwidth = h)
        case_1_Normal_Kernel_CV_pred[l] <- kernel_pred(X[l], X_CV, Y_CV, kernel = "gaussian", bandwidth = h)
      }
      case_1_Normal_Kernel_CV_hyper_MSE[hindex] <- mean((Y-case_1_Normal_Kernel_CV_pred)^2)
    }
    case_1_Normal_Kernel_CV_index <- which.min(case_1_Normal_Kernel_CV_hyper_MSE)
    case_1_Normal_Kernel_CV_minh[i] <- h_kernel[case_1_Normal_Kernel_CV_index]
    Normal_Kernel_CV <- kernel_smooth(X, Y, kernel = "gaussian", bandwidth = case_1_Normal_Kernel_CV_minh[i])
    case_1_Normal_Kernel_CV_MSE[i] <- mean((Y-Normal_Kernel_CV)^2)
    
    end.time <- Sys.time()
    
    # MISE LOOCV Normal
    Normal_Kernel_CV_MISE <- kernel_smooth_MISE(X, X_exact, Y, kernel = "gaussian", bandwidth = case_1_Normal_Kernel_CV_minh[i])
    
    case_1_Normal_Kernel_CV_MISE_hat[i] <- mean((fx_exact-Normal_Kernel_CV_MISE)^2)
    
    case_1_Normal_Kernel_CV_time[i] <- round(end.time - start.time, 10)
  }
  
  
  resultat3[[case_index]] <- list(Fourier_GCV_MSE <- mean(case_1_Fourier_GCV_MSE),
                                  Fourier_GCV_MISE <- mean(case_1_Fourier_GCV_MISE_hat),
                                  Fourier_GCV_tid <- mean(case_1_Fourier_GCV_time),
                                  Fourier_GCV_baser <- mean(case_1_Fourier_GCV_basis),
                                  Fourier_CV_MSE <- mean(case_1_Fourier_CV_MSE),
                                  Fourier_CV_MISE <- mean(case_1_Fourier_CV_MISE_hat),
                                  Fourier_CV_tid <- mean(case_1_Fourier_CV_time),
                                  Fourier_CV_baser <- mean(case_1_Fourier_CV_minbasis),
                                  Spline_GCV_MSE <- mean(case_1_Spline_GCV_MSE),
                                  Spline_GCV_MISE <- mean(case_1_Spline_GCV_MISE_hat),
                                  Spline_GCV_tid <- mean(case_1_Spline_GCV_time),
                                  Spline_GCV_baser <- mean(case_1_Spline_GCV_basis),
                                  Spline_CV_MSE <- mean(case_1_Spline_CV_MSE),
                                  Spline_CV_MISE <- mean(case_1_Spline_CV_MISE_hat),
                                  Spline_CV_tid <- mean(case_1_Spline_CV_time),
                                  Spline_CV_baser <- mean(case_1_Spline_CV_minbasis),
                                  Smoothing_Fourier_GCV_MSE <- mean(case_1_Fourier_smooth_GCV_MSE),
                                  Smoothing_Fourier_GCV_MISE <- mean(case_1_Fourier_smooth_GCV_MISE_hat),
                                  Smoothing_Fourier_GCV_tid <- mean(case_1_Fourier_smooth_GCV_time),
                                  Smoothing_Fourier_GCV_lambda <- mean(case_1_Fourier_smooth_GCV_lambda),
                                  Smoothing_Fourier_CV_MSE <- mean(case_1_Fourier_smooth_CV_MSE),
                                  Smoothing_Fourier_CV_MISE <- mean(case_1_Fourier_smooth_CV_MISE_hat),
                                  Smoothing_Fourier_CV_tid <- mean(case_1_Fourier_smooth_CV_time),
                                  Smoothing_Fourier_CV_lambda <- mean(case_1_Fourier_smooth_CV_minlambda),
                                  Smoothing_Spline_GCV_MSE <- mean(case_1_Spline_smooth_GCV_MSE),
                                  Smoothing_Spline_GCV_MISE <- mean(case_1_Spline_smooth_GCV_MISE_hat),
                                  Smoothing_Spline_GCV_tid <- mean(case_1_Spline_smooth_GCV_time),
                                  Smoothing_Spline_GCV_lambda <- mean(case_1_Spline_smooth_GCV_lambda),
                                  Smoothing_Spline_CV_MSE <- mean(case_1_Spline_smooth_CV_MSE),
                                  Smoothing_Spline_CV_MISE <- mean(case_1_Spline_smooth_CV_MISE_hat),
                                  Smoothing_Spline_CV_tid <- mean(case_1_Spline_smooth_CV_time),
                                  Smoothing_Spline_CV_lambda <- mean(case_1_Spline_smooth_CV_minlambda),
                                  Unif_Kernel_GCV_MSE <- mean(case_1_Unif_Kernel_GCV_MSE),
                                  Unif_Kernel_GCV_MISE <- mean(case_1_Unif_Kernel_GCV_MISE_hat),
                                  Unif_Kernel_GCV_tid <- mean(case_1_Unif_Kernel_GCV_time),
                                  Unif_Kernel_GCV_h <- mean(case_1_Unif_Kernel_GCV_minh),
                                  Unif_Kernel_CV_MSE <- mean(case_1_Unif_Kernel_CV_MSE),
                                  Unif_Kernel_CV_MISE <- mean(case_1_Unif_Kernel_CV_MISE_hat),
                                  Unif_Kernel_CV_tid <- mean(case_1_Unif_Kernel_CV_time),
                                  Unif_Kernel_CV_h <- mean(case_1_Unif_Kernel_CV_minh),
                                  Quad_Kernel_GCV_MSE <- mean(case_1_Quad_Kernel_GCV_MSE),
                                  Quad_Kernel_GCV_MISE <- mean(case_1_Quad_Kernel_GCV_MISE_hat),
                                  Quad_Kernel_GCV_tid <- mean(case_1_Quad_Kernel_GCV_time),
                                  Quad_Kernel_GCV_h <- mean(case_1_Quad_Kernel_GCV_minh),
                                  Quad_Kernel_CV_MSE <- mean(case_1_Quad_Kernel_CV_MSE),
                                  Quad_Kernel_CV_MISE <- mean(case_1_Quad_Kernel_CV_MISE_hat),
                                  Quad_Kernel_CV_tid <- mean(case_1_Quad_Kernel_CV_time),
                                  Quad_Kernel_CV_h <- mean(case_1_Quad_Kernel_CV_minh),
                                  Normal_Kernel_GCV_MSE <- mean(case_1_Normal_Kernel_GCV_MSE),
                                  Normal_Kernel_GCV_MISE <- mean(case_1_Normal_Kernel_GCV_MISE_hat),
                                  Normal_Kernel_GCV_tid <- mean(case_1_Normal_Kernel_GCV_time),
                                  Normal_Kernel_GCV_h <- mean(case_1_Normal_Kernel_GCV_minh),
                                  Normal_Kernel_CV_MSE <- mean(case_1_Normal_Kernel_CV_MSE),
                                  Normal_Kernel_CV_MISE <- mean(case_1_Normal_Kernel_CV_MISE_hat),
                                  Normal_Kernel_CV_tid <- mean(case_1_Normal_Kernel_CV_time),
                                  Normal_Kernel_CV_h <- mean(case_1_Normal_Kernel_CV_minh))
}