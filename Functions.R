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

# define the squared-exponential covariance function
make_rho_squexp <- function(l, sigma2 = 1) {
  function(h) {
    r <- sqrt(h[1]^2 + h[2]^2)  # Euclidean distance
    sigma2 * exp(-r^2 / (2*(l^2)))
  }
}

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
