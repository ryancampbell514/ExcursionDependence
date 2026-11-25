# A clone of "proportion.range" found in https://github.com/pony-helping-treat/Local-Excursion-Set-Coefficients
proportion.range <- function(ep, proportion = 1) {
  left <- 1
  right <- length(ep$y)
  result <- NULL
  
  while (left <= right) {
    mid <- floor((left + right) / 2)
    
    if (ep$y[mid] < proportion) {
      result <- mid
      right <- mid - 1  # Keep searching to the left for the first occurrence
    } else {
      left <- mid + 1  # Search to the right
    }
  }
  if (is.null(result)) return(Inf)
  return(ep$x[result])
}

# Metropolis Algorithm to simulate from Gaussian location and scale mixutres,
# adapted from Algorithm 1 of https://arxiv.org/abs/2602.09512
MetropolisAlg = function(lambda0,target.dens,nburn=5000,
                            prop.fn = function(loc) rnorm(1,mean=loc),
                            plot=FALSE){
  
  lambda <- lambda0
  dlam <- target.dens(lambda)
  accepts <- 0
  
  i=1
  
  lambdas=c(lambda)
  
  while (i<nburn) {
    
    # Propose (log-normal random walk)
    lambda_prop <- prop.fn(loc=lambda)
    
    dlam_prop <- target.dens(lambda_prop)
    
    # Metropolis–Hastings acceptance ratio
    # Hastings correction: q(x | x') / q(x' | x) = x_prop / x
    alpha = dlam_prop/dlam
    alpha=min(0,alpha)
    
    unif.obs = runif(1)
    
    if (unif.obs < alpha) {
      lambda <- lambda_prop
      dlam <- target.dens(lambda)
      print(lambda)
    }
    i = i+1
    
    lambdas=c(lambdas,lambda)
    
  }
  if(plot){
    plot(lambdas,type="l")
  }
  return(lambda)
}

# Same as above, but on the log-scale. May improve performance in some instances.
logMetropolisAlg = function(lambda0,log.target.dens,nburn=5000,
                          prop.fn = function(loc) rnorm(1,mean=loc),
                          plot=FALSE){
  
  lambda <- lambda0
  dlam <- log.target.dens(lambda)
  accepts <- 0
  
  i=1
  
  lambdas=c(lambda)
  
  while (i<nburn) {
    
    # Propose (log-normal random walk)
    lambda_prop <- prop.fn(loc=lambda)
    
    dlam_prop <- log.target.dens(lambda_prop)
    
    # Metropolis–Hastings acceptance ratio
    # Hastings correction: q(x | x') / q(x' | x) = x_prop / x
    alpha = dlam_prop - dlam
    # alpha=min(0,alpha)

    unif.obs = log(runif(1)) 
    
    if (unif.obs < alpha) {
      lambda <- lambda_prop
      dlam <- log.target.dens(lambda)
      # print(lambda)
    }
    i = i+1
    
    lambdas=c(lambdas,lambda)
    
  }
  if(plot){
    plot(lambdas,type="l")
  }
  return(lambda)
}


# A copy of the generateGaussian function in https://github.com/pony-helping-treat/Local-Excursion-Set-Coefficients,
# but using the truncnorm package is faster for very large conditioning thresholds.
generateGaussian.truncnorm = function (generator, cond_val = NULL, exceedance = FALSE)
{
  require(truncnorm)
  if (!is.null(cond_val) & exceedance) {
    stopifnot(nrow(generator$cond_sites) == 1)
    cond_val = rtruncnorm(n=1,a=cond_val)
  }
  Z <- rnorm(ncol(generator$L))
  if (!is.null(generator$cond_sites)) {
    if (is.null(cond_val)) {
      if (nrow(generator$cond_sites) == 1) {
        cond_val = rnorm(1)
      }
      else {
        stop("generateRF: Please provide values for the conditioning sites.")
      }
    }
    mu = generator$reg_coeff %*% cond_val
  }
  else {
    mu = 0
  }
  values <- as.vector(generator$L %*% Z + mu)
  return(random.field(values, generator$grid))
}

# A copy of the plotRandomField function in https://github.com/pony-helping-treat/Local-Excursion-Set-Coefficients,
# but with improved visual features.
plotRandomField.v2 = function (rf, x, y, thresh = NULL, draw.boundary=FALSE, main=NA) {
  require(ExcursionSets)
  if (!is.random.field(rf)) {
    rf = random.field(as.matrix(rf))
  }
  z = rf$values
  dim(z) = rf$dimensions
  if (is.null(thresh)) {
    # colors
    cols <- rev(heat.colors(100))
    fields::image.plot(x, y, z,col=cols,xlab=NA,ylab=NA,main=main)#, asp = 1)
    
  }
  else {
    par(pty="s")
    image(x, y, z > thresh, col = gray.colors(2, start = 1, 
                                              end = 0.5),
          xlab=NA,ylab=NA,main=main)
    
    if(draw.boundary){
      # add lines (indicating perimeter)
      contour(x, y, z,
              levels = thresh,
              add = TRUE,
              drawlabels = FALSE,
              col = "black",
              lwd = 2)
    }
  }
}

# Matern covariance function
matern_cov <- function(r, nu, l, sigma2 = 1) {
  r[r == 0] <- 1e-12
  kappa <- sqrt(2 * nu) / l
  val <- sigma2 * (2^(1 - nu) / gamma(nu)) *
    (kappa * r)^nu * besselK(kappa * r, nu)
  val[r < 1e-10] <- sigma2
  val
}

# An R translation of the Matlab code used to simulate Gaussian random fields on large grids,
# with original Matlab code accomanying the manuscript https://projecteuclid.org/journals/electronic-journal-of-statistics/volume-13/issue-1/Lipschitz-Killing-curvatures-of-excursion-sets-for-two-dimensional-random/10.1214/19-EJS1530.full
stationary_Gaussian_process_EDiBcode <- function(
    nx, ny,
    xlim = c(0, 1),
    ylim = c(0, 1),
    rho,
    anis=NULL,
    check.lam=TRUE
) {
  # Reference:
    # Kroese, D. P., & Botev, Z. I. (2015). Spatial Process Simulation.
    # In Stochastic Geometry, Spatial Statistics and Random Fields (pp. 369-404)
    # Springer International Publishing, DOI: 10.1007/978-3-319-10064-7_12
  
  # Create spatial grid
  tx <- seq(xlim[1], xlim[2], length.out = nx)
  ty <- seq(ylim[1], ylim[2], length.out = ny)
  
  dx <- tx[2] - tx[1]
  dy <- ty[2] - ty[1]
  
  # Build covariance samples (IMPORTANT: scale by dx, dy)
  if(is.null(anis)){
    Rows <- outer(ty - ty[1], tx - tx[1],
                  Vectorize(function(y, x) rho(c(x, y))))
    Cols <- outer(ty - ty[1], tx - tx[1],
                  Vectorize(function(y, x) rho(c(-x, y))))
  } else {
    Rows <- outer(ty - ty[1], tx - tx[1],
                  Vectorize(function(y, x) {
                      inp.vec = as.numeric(anis %*% matrix(c(x,y),ncol=1))
                      rho(inp.vec)
                    }))
    Cols <- outer(ty - ty[1], tx - tx[1],
                  Vectorize(function(y, x) {
                      inp.vec = as.numeric(anis %*% matrix(c(-x,y),ncol=1))
                      rho(inp.vec)
                    }))
  }
  
  
  # Circulant embedding
  BlkCirc_row <- rbind(
    cbind(Rows, Cols[, ny:2]),
    cbind(Cols[nx:2, ], Rows[nx:2, ny:2])
  )
  
  # Eigenvalues
  lam <- Re(fft(BlkCirc_row)) / ((2*nx - 1) * (2*ny - 1))
  
  if(check.lam){
    if (min(lam) < -1e-15) {
      stop("Embedding not positive definite — try larger domain or smoother covariance")
    }
  }
  
  lam[lam < 0] <- 0
  lam <- sqrt(lam)
  
  # Complex Gaussian noise
  Z <- complex(
    real = matrix(rnorm((2*nx-1)*(2*ny-1)), 2*nx-1),
    imaginary = matrix(rnorm((2*nx-1)*(2*ny-1)), 2*nx-1)
  )
  
  # Simulation
  F <- fft(lam * Z)
  field <- Re(F[1:nx, 1:ny])
  
  return(list(field = field, x = tx, y = ty))
}

# A conditional version of the above, using krigging.
stationary_Gaussian_process_EDiBcode_cond <- function(nx,ny,xlim,ylim,rho,
                                                      obs_loc=matrix(c(0,0),nrow=1), 
                                                      obs_val=qnorm(0.99)) {
  
  field = stationary_Gaussian_process_EDiBcode(nx,ny,xlim,ylim,rho,check.lam=F)
  field = field$field
  
  x=seq(xlim[1],xlim[2],length.out=nx)
  y=seq(ylim[1],ylim[2],length.out=ny)
  grid <- expand.grid(x = x, y = y)
  
  # simulated values at observation locations
  sim_vals <- numeric(nrow(obs_loc))
  
  for (i in 1:nrow(obs_loc)) {
    ix <- which.min(abs(x - obs_loc[i,1]))
    iy <- which.min(abs(x - obs_loc[i,2]))
    sim_vals[i] <- field[ix, iy]  # the value of the unconditional random field at s0
  }
  
  # residuals
  r <- obs_val - sim_vals  # target value - cirrent value
  
  # covariance among observations
  # D_oo <- as.matrix(dist(obs_loc))
  # C_oo <- rho(D_oo)#matern_cov(D_oo, rho)
  C_oo = rho(obs_loc)
  Cinv <- solve(C_oo)
  
  # covariance grid → observations
  C_go <- matrix(0, nrow(grid), nrow(obs_loc))
  
  for (i in 1:nrow(obs_loc)) {
    
    # d <- sqrt((grid$x - obs_loc[i,1])^2 +
    #             (grid$y - obs_loc[i,2])^2)
    # 
    # C_go[,i] <- matern_cov(d, nu, l, sigma2)
    C_go[,i] = apply(cbind(grid$x - obs_loc[i,1], grid$y - obs_loc[i,2]),1,rho)
  }
  
  # kriging correction
  correction <- C_go %*% Cinv %*% r
  
  n = length(x)
  field_cond <- field + matrix(correction, n, n)
  
  field_cond
}

# define the Matern covariance function
make_rho_matern <- function(nu, l, sigma2 = 1) {
  function(h) {
    r <- sqrt(h[1]^2 + h[2]^2)  # Euclidean distance
    matern_cov(r, nu, l, sigma2)
  }
}

# define the squared-exponential covariance function
make_rho_squexp <- function(l, sigma2 = 1) {
  function(h) {
    r <- sqrt(h[1]^2 + h[2]^2)  # Euclidean distance
    sigma2 * exp(-r^2 / (2*(l^2)))
  }
}



simulate_BrownResnick <- function(n, L, ell, K = 1000) {
  
  stopifnot(n >= 2, L > 0, ell > 0, K >= 1)
  
  # 1. Grid
  x <- seq(-L / 2, L / 2, length.out = n)
  dx <- L / (n - 1)
  
  i0 <- which.min(abs(x))
  
  # 2. Circulant embedding
  # Embedding dimension (power of 2)
  m <- 2^ceiling(log2(2 * n))
  
  j <- 0:(n - 1)
  
  # Squared-exponential covariance
  c1 <- exp(-0.5 * (j * dx / ell)^2)
  
  # Symmetric circulant embedding
  cemb <- numeric(m)
  cemb[1:n] <- c1
  cemb[m - j[-1] + 1] <- c1[-1]
  
  # Eigenvalues
  lambda1 <- Re(fft(cemb))
  
  if (min(lambda1) < -1e-8) {
    stop("Circulant embedding is not positive definite.")
  }
  
  lambda1 <- pmax(lambda1, 0)
  
  # Separable 2D covariance
  sqrt_lambda2 <- sqrt(outer(lambda1, lambda1))
  
  # 3. Half variogram relative to anchor
  half_vario <- 1 - exp(
    -outer(
      (x - x[i0])^2,
      (x - x[i0])^2,
      "+"
    ) / (2 * ell^2)
  )
  
  # 4. Poisson points
  Gamma <- cumsum(rexp(K))
  
  # 5. Running maximum
  logZ <- matrix(-Inf, n, n)
  
  for (i in seq_len(K)) {
    
    eps <- matrix(rnorm(m * m), m, m)
    
    # FFT Gaussian simulation
    Gfull <- Re(
      fft(
        sqrt_lambda2 * fft(eps),
        inverse = TRUE
      )
    ) / m^2
    
    # Crop to original grid
    G <- Gfull[1:n, 1:n]
    
    # Stationary increments
    W <- G - G[i0, i0]
    
    # Brown-Resnick spectral function
    logY <- -log(Gamma[i]) + W - half_vario
    
    logZ <- pmax(logZ, logY)
  }
  
  list(
    Z = exp(logZ),
    x = x,
    y = x,
    n = n,
    L = L,
    ell = ell,
    K = K
  )
}


# setup.GRF.inla = function(nu=2.5,range=3,alpha=2,sigma0=1,
#                           xlim=c(-15,15),ylim=c(-15,15),
#                           nx=60,ny=60,thresh=NULL,aniso=NULL){
#   # get the INLA precision matrix
#   
#   # x <- seq(xlim[1], xlim[2], length.out = nx)
#   # y <- seq(ylim[1], ylim[2], length.out = ny)
#   
#   # ---------------------------
#   # 2. Build mesh (coarse!)
#   # ---------------------------
#   coords = expand.grid(
#     x = seq(xlim[1], xlim[2], length.out = nx),
#     y = seq(ylim[1], ylim[2], length.out = ny)
#   )
#   if(!is.null(aniso)){
#     aniso.mat = matrix(c(aniso[1],1,1,aniso[2]),2,2)
#     coords = t(apply(coords,1,function(vec){
#       vec = matrix(vec,ncol=1)
#       return(as.numeric(aniso.mat %*% vec))
#     }))
#   }
# 
#   # mesh <- inla.mesh.2d(
#   #   loc = coords,
#   #   max.edge = c(1, 3),
#   #   cutoff = 0.1
#   # )
#   mesh = fmesher::fm_mesh_2d_inla(loc=coords,max.edge = c(1, 3),
#                            cutoff = 0.1)
# 
#   n.mesh = mesh$n
#   
#   # ---------------------------
#   # 3. Matérn parameters
#   # ---------------------------
#   # SPDE uses:
#   # alpha <- 3   # <-- approximation (nu ≈ 2)
#   
#   kappa <- sqrt(2 * nu) / range   # from your definition
#   tau <- 1/(2*sqrt(pi)*kappa*sigma0)
#   # ---------------------------
#   # 4. Define SPDE model
#   # (parameters effectively fixed)
#   # ---------------------------
#   kappa.vec <- rep(kappa, n.mesh)
#   tau.vec <- rep(tau, n.mesh)
#   B.kappa <- matrix(log(kappa.vec),n.mesh, 1)
#   B.tau <- matrix(log(tau.vec), n.mesh, 1)
#   spde = inla.spde2.matern(mesh, alpha = alpha, B.tau=B.tau, B.kappa=B.kappa)
#   # range0=2*range
#   # spde = inla.spde2.pcmatern(mesh,alpha=alpha,prior.range = c(range0, 0.5), prior.sigma = c(sigma0, 0.01) )
#   
#   # ---------------------------
#   # 5. Sample from GRF
#   # ---------------------------
# 
#   Q <- inla.spde2.precision(spde)
#   return(list(prec.mat=Q,
#               mesh=mesh))
#   
# }
# 
# 
# sim.GRF.inla = function(n=1,inla.prec.mat,thresh=NULL){
#   field_sample <- inla.qsample(n=n, inla.prec.mat)
# }
# 
# proj.GRF.inla = function(grf.sims){
#   stop("DO NOT RUN THIS. TOO SLOW.")
#   
#   # ---------------------------
#   # 6. Project to 1024x1024 grid
#   # ---------------------------
#   proj <- inla.mesh.projector(
#     mesh,
#     xlim = c(-15, 15),
#     ylim = c(-15, 15),
#     dims = c(nx, ny)
#   )
#   
#   field <- inla.mesh.project(proj, field_sample[,1])
#   
#   return(field)
# }
# 
# 
# setup.cond.GRF.inla = function(nu=2.5,range=3,alpha=2,sigma0=1,
#                           xlim=c(-15,15),ylim=c(-15,15),
#                           nx=60,ny=60,s0=matrix(c(0, 0), ncol = 2),thresh=1.0){
#   # get the INLA precision matrix
#   
#   # x <- seq(xlim[1], xlim[2], length.out = nx)
#   # y <- seq(ylim[1], ylim[2], length.out = ny)
#   
#   # ---------------------------
#   # 2. Build mesh (coarse!)
#   # ---------------------------
#   mesh <- inla.mesh.2d(
#     loc = expand.grid(
#       x = seq(xlim[1], xlim[2], length.out = 60),
#       y = seq(ylim[1], ylim[2], length.out = 60)
#     ),
#     max.edge = c(1, 3),
#     cutoff = 0.1
#   )
#   n.mesh = mesh$n
#   
#   # ---------------------------
#   # 3. Matérn parameters
#   # ---------------------------
#   # SPDE uses:
#   # alpha <- 3   # <-- approximation (nu ≈ 2)
#   
#   kappa <- sqrt(2 * nu) / range   # from your definition
#   tau <- 1/(2*sqrt(pi)*kappa*sigma0)
#   # ---------------------------
#   # 4. Define SPDE model
#   # (parameters effectively fixed)
#   # ---------------------------
#   kappa.vec <- rep(kappa, n.mesh)
#   tau.vec <- rep(tau, n.mesh)
#   B.kappa <- matrix(log(kappa.vec),n.mesh, 1)
#   B.tau <- matrix(log(tau.vec), n.mesh, 1)
#   spde = inla.spde2.matern(mesh, alpha = alpha, B.tau=B.tau, B.kappa=B.kappa)
#   # range0=2*range
#   # spde = inla.spde2.pcmatern(mesh,alpha=alpha,prior.range = c(range0, 0.5), prior.sigma = c(sigma0, 0.01) )
#   
#   # # ---------------------------
#   # # 5. Sample from GRF
#   # # ---------------------------
#   # 
#   # Q <- inla.spde2.precision(spde)
#   # return(list(prec.mat=Q,
#   #             mesh=mesh))
#   
#   # ---------------------------
#   # 4. Conditioning point
#   # ---------------------------
#   # s0 <- matrix(c(0, 0), ncol = 2)  # location
#   # thresh <- 1.0                        # value to enforce
#   
#   # Projector matrix at observation
#   A_obs <- inla.spde.make.A(mesh, loc = s0)
#   
#   # ---------------------------
#   # 5. Build stack
#   # ---------------------------
#   stack <- inla.stack(
#     data = list(y = thresh),
#     A = list(A_obs, 1),
#     effects = list(
#       field = 1:spde$n.spde,
#       intercept = 1
#     )
#   )
#   
#   # ---------------------------
#   # 6. Fit model with tiny noise
#   # ---------------------------
#   result <- inla(
#     y ~ 0 + intercept + f(field, model = spde),
#     data = inla.stack.data(stack),
#     control.predictor = list(
#       A = inla.stack.A(stack),
#       compute = TRUE
#     ),
#     family = "gaussian",
#     control.family = list(
#       hyper = list(
#         prec = list(initial = log(1e12), fixed = TRUE)
#       )
#     ),
#     control.compute = list(config = TRUE)
#     # control.compute = list(config = TRUE, openmp.strategy = "safe")
#   )
#   
#   # return the model fit.
#   return(list(model.fit=result,mesh=mesh,spde=spde))
#   
# }
# 
# # update_conditioning <- function(z0, mesh, spde, s0) {
# #   
# #   A_obs <- inla.spde.make.A(mesh, loc = s0)
# #   
# #   stack <- inla.stack(
# #     data = list(y = z0),
# #     A = list(A_obs, 1),
# #     effects = list(
# #       field = 1:spde$n.spde,
# #       intercept = 1
# #     )
# #   )
# #   
# #   result <- inla(
# #     y ~ 0 + intercept + f(field, model = spde),
# #     data = inla.stack.data(stack),
# #     control.predictor = list(
# #       A = inla.stack.A(stack),
# #       compute = TRUE
# #     ),
# #     family = "gaussian",
# #     control.family = list(
# #       hyper = list(
# #         prec = list(initial = log(1e12), fixed = TRUE)
# #       )
# #     ),
# #     control.compute = list(config = TRUE)
# #   )
# #   
# #   return(result)
# # }
# 
# gauss.filter = function(X,sigma=1.5,nx=2^10,ny=2^10,delta){
#   # sigma <- 1.5  # controls how much \nu increases (not exact)
#   # nx = 2^10; ny=2^10
#   kx <- if (nx %% 2 == 0) {
#     c(0:(nx/2 - 1), -(nx/2): -1)
#   } else {
#     c(0:((nx-1)/2), -((nx-1)/2): -1)
#   }
#   
#   ky <- if (ny %% 2 == 0) {
#     c(0:(ny/2 - 1), -(ny/2): -1)
#   } else {
#     c(0:((ny-1)/2), -((ny-1)/2): -1)
#   }
#   wx <- (2 * pi / (nx * delta)) * kx
#   wy <- (2 * pi / (ny * delta)) * ky
#   wx2 <- outer(wx^2, rep(1, nx))
#   wy2 <- outer(rep(1, ny), wy^2)
#   omega2 <- wx2 + wy2
#   # Fourier transform of Gaussian kernel
#   K_hat <- exp(-0.5 * sigma^2 * omega2)
#   # FFT convolution
#   field_fft <- fft(X)
#   smoothed_fft <- field_fft * K_hat
#   X.smooth <- Re(fft(smoothed_fft, inverse = TRUE)) / (nx * ny)
#   # X.smooth <- X.smooth * sqrt(4 * pi * sigma^2)
#   X.smooth = (X.smooth - mean(X.smooth,na.rm=T))/sd(X.smooth,na.rm=T)
#   return(X.smooth)
# }
# 
# 
# est.spec2 = function(X,n,dx){
#   # Compute gradients
#   dx_forward <- function(M) {
#     (M[c(2:n,1), ] - M[c(n,1:(n-1)), ]) / (2*dx)
#   }
#   
#   dy_forward <- function(M) {
#     (M[, c(2:n,1)] - M[, c(n,1:(n-1))]) / (2*dx)
#   }
#   
#   grad_x <- dx_forward(X)
#   grad_y <- dy_forward(X)
#   grad_sq <- grad_x^2 + grad_y^2
#   m2_grad <- mean(grad_sq)
#   
#   return(m2_grad)
#   
# }

# use the terra package, obtain LER for rasters. Much faster than the 
# code in https://github.com/pony-helping-treat/Local-Excursion-Set-Coefficients
terra.ler = function(rf=NULL,bf=NULL,thresh,half.width,s0){
  if(is.null(bf)){
    if(is.null(dim(rf))){
      stop("rf must be a matrix or a raster")
    }
    # if(!("matrix" %in% class(s0))){
    #   stop("s0 must be a matrix with 1 row and 2 columns")
    # }
    is.not.rast = !("SpatRaster" %in% class(rf))
    if(is.not.rast){
      bf = terra::rast(x=rf>thresh)
    } else {
      bf = rf>thresh  # the binary field, this takes time
    }
  }
  
  if(is.null(half.width)){
    # terra::ext(bf) <- c(0,ncol(bf),0,nrow(bf))  # set the dimensions
  } else {
    terra::ext(bf) <- c(-half.width, half.width, -half.width, half.width)  # set the dimensions
  }
  # plot(bf)
  # points(s0,col="pink",pch=20)
  d = terra::gridDist(bf, target = FALSE)
  # project into raster CRS

  # stop("TO DO: fix the s0 projection")
  
  # extract distance
  ler = as.numeric(terra::extract(d, s0))[-1]
  if(!is.null(half.width)){
    ler = ifelse(is.na(ler),sqrt(2*(half.width^2)),ler)
  }
  
  return(ler) 
}

raster.crop = function(rf.obs,crop.len,lon,lat){
  half.crop.len=ceiling(crop.len/2)

  loc <- terra::vect(data.frame(lon=lon, lat=lat),crs="EPSG:4326")
  s0_proj <- terra::project(loc, crs(rf.obs))
  cell <- terra::cellFromXY(rf.obs, terra::crds(s0_proj))
  xy   <- terra::xyFromCell(rf.obs, cell)
  res_r <- terra::res(rf.obs)
  win <- terra::ext(
    xy[1] - half.crop.len*res_r[1], xy[1] + half.crop.len*res_r[1],
    xy[2] - half.crop.len*res_r[2], xy[2] + half.crop.len*res_r[2]
  )

  # this takes the most time...
  rf.crop <- terra::resample(terra::crop(rf.obs, win),
                             terra::rast(win,
                                         resolution=res_r,
                                         crs=terra::crs(rf.obs)))

  return(rf.crop)
}

raster.cc = function(rf.obs,thresh,lon,lat){
  # this function takes in a random field ans threshold, and returns the
  # binary field with connected component containing lon and lat.

  # warning("this can cause problems")
  exc <- (rf.obs > thresh)#*1

  # plot(exc)
  # points(s0_proj, col="red", pch=16, cex=1.2)
  cl <- terra::patches(exc,values=T)

  # cl
  loc <- terra::vect(data.frame(lon=lon, lat=lat),crs="EPSG:4326")
  s0_proj <- terra::project(loc, terra::crs(rf.obs))
  loc_cell <- terra::cellFromXY(cl, terra::crds(s0_proj))

  loc_id <- terra::values(cl)[loc_cell]

  # loc_id
  exc_loc <- cl
  exc_loc[exc_loc != loc_id] <- 0
  exc_loc[exc_loc == loc_id] <- 1

  bf = as.vector(exc_loc)
  bf <- matrix(bf, nrow = sqrt(length(bf)), ncol = sqrt(length(bf)), byrow = TRUE)
  bf = t(bf)
  bf <- bf[, nrow(bf):1]

  return(bf)
}
# 
# terra.perim = function(rf.obs,thresh,lon,lat){
#   # excursion set
#   exc <- rf.obs > thresh
#   
#   # par(mfrow=c(1,2),pty="s")
#   # plot(exc)
#   
#   # keep only the component containing lon and lat
#   # plot(exc)
#   # points(s0_proj, col="red", pch=16, cex=1.2)
#   cl <- terra::patches(exc,values=T)
#   
#   # cl
#   loc <- terra::vect(data.frame(lon=lon, lat=lat),crs="EPSG:4326")
#   s0_proj <- terra::project(loc, terra::crs(rf.obs))
#   loc_cell <- terra::cellFromXY(cl, terra::crds(s0_proj))
#   
#   loc_id <- terra::values(cl)[loc_cell]
#   
#   # loc_id
#   exc_loc <- cl
#   exc_loc[exc_loc != loc_id] <- 0
#   exc_loc[exc_loc == loc_id] <- 1  
#   
#   # plot(exc_loc)
#   # points(s0_proj,col="red")
#   
#   # p <- terra::as.polygons(exc_loc, dissolve=TRUE)
#   # if (is.null(p) || nrow(p) == 0) return(0)
#   # 
#   # b <- terra::as.lines(p)
#   # 
#   # # window boundary
#   # w <- terra::ext(rf.obs)
#   # s_poly <- terra::as.polygons(w)
#   # s_boundary <- terra::as.lines(s_poly)
#   # 
#   # # remove boundary induced by window clipping
#   # b_clean <- terra::erase(b, s_boundary)
#   # 
#   # return((terra::perim(b_clean))[1])
#   
#   # remove FALSE cells
#   exc_loc[exc_loc == 0] <- NA
# 
#   # convert connected patches to polygons
#   p <- as.polygons(exc_loc, dissolve = FALSE)
# 
#   # if no excursion region exists
#   if (is.null(p) || nrow(p) == 0) {
#     return(0)
#   }
# 
#   # perimeter of polygons
#   perim <- perim(p)
# 
#   # total perimeter
#   perim.full = sum(perim)
#   
#   
#   S <- as.polygons(ext(rf.obs))
#   bS <- as.lines(S)
#   bE <- as.lines(p)
#   b_touch <- intersect(bE, bS)
#   L_touch <- perim(b_touch)#sum(terra::lengths(b_touch, geodesic=FALSE), na.rm=TRUE)
#   
#   return(perim.full-L_touch)
# }
# 
# perimeter_minkowski_bc <- function(rf.obs, thresh, dx=1, dy=1, lon=NULL, lat=NULL) {
#   
#   exc <- rf.obs > thresh
#   
#   # par(mfrow=c(1,2),pty="s")
#   # plot(exc)
#   
#   # keep only the component containing lon and lat
#   # plot(exc)
#   # points(s0_proj, col="red", pch=16, cex=1.2)
#   cl <- terra::patches(exc,values=T)
#   
#   # cl
#   loc <- terra::vect(data.frame(lon=lon, lat=lat),crs="EPSG:4326")
#   s0_proj <- terra::project(loc, terra::crs(rf.obs))
#   loc_cell <- terra::cellFromXY(cl, terra::crds(s0_proj))
#   
#   loc_id <- terra::values(cl)[loc_cell]
#   
#   # loc_id
#   exc_loc <- cl
#   exc_loc[exc_loc != loc_id] <- 0
#   exc_loc[exc_loc == loc_id] <- 1 
#   
#   x = as.matrix(exc_loc)
#   nr <- nrow(x)
#   nc <- ncol(x)
#   
#   # remove boundary rows/cols (this is the window correction)
#   x2 <- x[2:(nr-1), 2:(nc-1)]
#   
#   nr <- nrow(x2)
#   nc <- ncol(x2)
#   
#   Eh <- sum(abs(x2[, -nc] - x2[, -1]), na.rm=TRUE)
#   Ev <- sum(abs(x2[-nr, ] - x2[-1, ]), na.rm=TRUE)
#   
#   return(dx * Ev + dy * Eh)
# }
# 
# perim.L1 <- function(bf, dx = 1, dy = 1) {
#   # bf is an observed binary field
#   # dx and dy are the pixel dimensions
#   
#   nr <- nrow(bf)
#   nc <- ncol(bf)
#   
#   # horizontal edges
#   p_vert <-
#     subf(bf[, -1] != bf[, -nc]) * dy +
#     subf(bf[, 1]) * dy +
#     subf(bf[, nc]) * dy
#   
#   # vertical edges
#   p_horz <-
#     subf(bf[-1, ] != bf[-nr, ]) * dx +
#     subf(bf[1, ]) * dx +
#     subf(bf[nr, ]) * dx
#   
#   p_vert + p_horz
# }

# L1-based perimeter estimation, ignoring boundaries
perim.L1.na.rm <- function(bf, dx = 1, dy = 1, incl_boundary=FALSE) {
  
  nr <- nrow(bf)
  nc <- ncol(bf)
  
  if(incl_boundary){
    bf[is.na(bf)] <- 0
  }
  
  # left-right neighbors
  A <- bf[, -1, drop = FALSE]
  B <- bf[, -nc, drop = FALSE]
  
  vertical_transitions = sum((A == 1 & B == 0) | (A == 0 & B == 1), na.rm = !incl_boundary)
  
  # up-down neighbors
  A <- bf[-1, , drop = FALSE]
  B <- bf[-nr, , drop = FALSE]
  
  horizontal_transitions = sum((A == 1 & B == 0) | (A == 0 & B == 1), na.rm = !incl_boundary)
  
  dy * vertical_transitions +
    dx * horizontal_transitions
}
