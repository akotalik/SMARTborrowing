

########################################################################################
########### functions define for the SMARTborrow.R file #######################################################
########################################################################################

#' @importFrom foreach getDoParWorkers

num_workers <- function(factor = 1) {
  
  getDoParWorkers() * factor
  
}



boa_hpd <- function(x, alpha) {
  
  n <- length(x)
  
  m <- max(1, ceiling(alpha * n))
  
  y <- sort(x)
  
  a <- y[seq_len(m)]
  
  b <- y[(n - m + 1):n]
  
  i <- order(b - a)[1]
  
  structure(c(a[i], b[i]), names = c("Lower Bound", "Upper Bound"))
  
}



mem_mat <- function(indices, mod_mat, h) {
  
  m <- matrix(NA, h, h)
  
  diag(m) <- rep(1, dim(m)[1])
  
  if (h == 2) {
    
    m[1, 2] <- m[2, 1] <- mod_mat[[1]][indices[1]]
    
    return(m)
    
  }
  
  for (i in seq_len(length(indices) - 1)) {
    
    m[(i + 1):dim(m)[2], i] <- m[i, (i + 1):dim(m)[2]] <-
      
      mod_mat[[i]][indices[i], ]
    
  }
  
  m[dim(m)[2], i + 1] <- m[i + 1, dim(m)[2]] <-
    
    c(0, 1)[indices[length(indices)]]
  
  m
  
}



log_marg_dens <- function(i_indices, mod_mat, xvec, nvec, avec, bvec) {
  
  m <- mem_mat(i_indices, mod_mat, length(xvec))
  
  marg_vec <- rep(NA, dim(m)[1])
  
  # calculate the product portion of integrated marginal likelihood
  
  prod_vec <- beta(xvec + avec, nvec + bvec - xvec) / beta(avec, bvec)
  
  for (i in seq_len(dim(m)[1])) {
    
    p_vec <- prod(prod_vec^(1 - m[i, ]))
    
    marg_vec[i] <- (beta(avec[i] + m[i, ] %*% xvec, bvec[i] + m[i, ] %*%
                           
                           (nvec - xvec)) / beta(avec[i], bvec[i])) *
      
      p_vec
    
  }
  
  sum(log(marg_vec))
  
}



ess <- function(x, n, omega, a, b) {
  
  alpha <- a + omega %*% x
  
  beta <- b + (omega %*% n - omega %*% x)
  
  alpha + beta
  
}



mem_prior <- function(i_indices, mod_mat, prior_inclusion) {
  
  m <- mem_mat(i_indices, mod_mat, nrow(prior_inclusion))
  
  mem <- m[upper.tri(m)]
  
  source_vec <- prior_inclusion[upper.tri(prior_inclusion)]
  
  s_in <- source_vec # source inclusion probability
  
  s_ex <- 1 - source_vec # source exclusion probability
  
  prod(s_in^(mem) * s_ex^(1 - mem))
  
}



mem_ij <- function(i_indices, mod_mat, prior_inclusion, i, j) {
  
  mem_mat(i_indices, mod_mat, nrow(prior_inclusion))[i, j]
  
}



ij_pep <- function(j, mod_i, mod_mat, prior_inclusion, i, log_marg, prior) {
  
  i_1 <- apply(mod_i, MARGIN = 1, FUN = mem_ij, mod_mat, prior_inclusion, i, j)
  
  ((exp(log_marg) * prior) %*% i_1) / sum(exp(log_marg) * prior)
  
}



comp_pep <- function(i, j_indices, mod_i, mod_mat, prior_inclusion, 
                     
                     log_marg, prior) {
  
  vapply((i + 1):j_indices,
         
         FUN = ij_pep, FUN.VALUE = NA_real_,
         
         mod_i, mod_mat, prior_inclusion, i, log_marg, prior
         
  )
  
}



mem_j <- function(i_indices, mod_mat, prior_inclusion, j, m) {
  
  mm <- mem_mat(i_indices, mod_mat, nrow(prior_inclusion))
  
  as.numeric(sum(as.numeric(mm[j, ] == m)) == length(m))
  
}



j_weight_1 <- function(i, mod_i, mod_mat, prior_inclusion, j, log_marg, prior) {
  
  if (length(log_marg) == 2) {
    
    i_m <- apply(mod_i,
                 
                 MARGIN = 1, FUN = mem_j, mod_mat, prior_inclusion, j,
                 
                 c(1, mod_mat[[1]][i])
                 
    )
    
  } else {
    
    i_m <- apply(mod_i,
                 
                 MARGIN = 1, FUN = mem_j, mod_mat, prior_inclusion, j,
                 
                 c(1, mod_mat[[1]][i, ])
                 
    )
    
  }
  
  ((exp(log_marg) * prior) %*% i_m) / sum(exp(log_marg) * prior)
  
}



j_weight_j <- function(i, mod_i, mod_mat, prior_inclusion, j, log_marg, prior) {
  
  u <- mod_mat[[1]][i, ]
  
  i_m <- apply(mod_i,
               
               MARGIN = 1, FUN = mem_j, mod_mat, prior_inclusion, j,
               
               c(u[1:(j - 1)], 1, u[j:length(u)])
               
  )
  
  ((exp(log_marg) * prior) %*% i_m) / sum(exp(log_marg) * prior)
  
}



j_weight_j_indices <- function(i, mod_i, mod_mat, prior_inclusion, j, 
                               
                               log_marg, prior) {
  
  if (length(log_marg) == 2) {
    
    i_m <- apply(mod_i,
                 
                 MARGIN = 1, FUN = mem_j, mod_mat, prior_inclusion, j,
                 
                 c(mod_mat[[1]][i], 1)
                 
    )
    
  } else {
    
    i_m <- apply(mod_i,
                 
                 MARGIN = 1, FUN = mem_j, mod_mat, prior_inclusion, j,
                 
                 c(mod_mat[[1]][i, ], 1)
                 
    )
    
  }
  
  ((exp(log_marg) * prior) %*% i_m) / sum(exp(log_marg) * prior)
  
}



post_weights <- function(j, j_indices, mod_i, mod_mat, prior_inclusion,
                         
                         log_marg, prior) {
  
  
  
  if (length(log_marg) == 2) {
    
    if (j == 1) {
      
      out <- vapply(1:2,
                    
                    FUN = j_weight_1, FUN.VALUE = NA_real_,
                    
                    mod_i, mod_mat, prior_inclusion, j, log_marg, prior
                    
      )
      
    } else if (j == j_indices) {
      
      out <- vapply(1:2,
                    
                    FUN = j_weight_j_indices, FUN.VALUE = NA_real_,
                    
                    mod_i, mod_mat, prior_inclusion, j, log_marg, prior
                    
      )
      
    }
    
    return(out)
    
  }
  
  
  
  if (j == 1) {
    
    out <- vapply(seq_len(nrow(mod_mat[[1]])),
                  
                  FUN = j_weight_1, FUN.VALUE = NA_real_,
                  
                  mod_i, mod_mat, prior_inclusion, j, log_marg, prior
                  
    )
    
  } else if (j == j_indices) {
    
    out <- vapply(seq_len(nrow(mod_mat[[1]])),
                  
                  FUN = j_weight_j_indices, FUN.VALUE = NA_real_,
                  
                  mod_i, mod_mat, prior_inclusion, j, log_marg, prior
                  
    )
    
  } else {
    
    out <- vapply(seq_len(nrow(mod_mat[[1]])),
                  
                  FUN = j_weight_j, FUN.VALUE = NA_real_,
                  
                  mod_i, mod_mat, prior_inclusion, j, log_marg, prior
                  
    )
    
  }
  
  out
  
}



euc_dist <- function(x1, x2, w = c(1, 1)) {
  
  if (sum(is.na(x1)) > 1) {
    
    Inf
    
  } else {
    
    sqrt(sum(w * ((x1 - x2)^2)))
    
  }
  
}



#' @importFrom stats qbeta

dist_beta_hpd <- function(ess, fit, alpha, jj) {
  
  al <- fit$mean_est[jj] * ess
  
  al <- max(1e-2, al)
  
  be <- ess - al
  
  be <- max(1e-2, be)
  
  euc_dist(fit$hpd[, jj], qbeta(
    
    c(alpha / 2, 1 - alpha / 2),
    
    al, be
    
  ))
  
}



#' @importFrom GenSA GenSA

ess_from_hpd_i <- function(jj, fit, alpha) {
  
  opt <-
    
    GenSA(
      
      par = 1,
      
      fn = dist_beta_hpd,
      
      lower = 0,
      
      upper = 10000000,
      
      fit = fit,
      
      alpha = alpha,
      
      jj = jj
      
    )
  
  opt$par
  
}



#' @importFrom foreach %dopar%

ess_from_hpd <- function(fit, alpha) {
  
  ## fit is list with median vec and HPD vec ##
  
  i <- NULL
  
  foreach(i = seq_along(fit$mean_est), .combine = c) %dopar% {
    
    ess_from_hpd_i(i, fit, alpha)
    
  }
  
}



####################################################################

########## MCMC for Bayes with Metropolis-Hastings #################

####################################################################



flip_mem <- function(v, m_old, m) {
  
  out <- m_old
  
  if (m_old[which(m == v)[1]] == 1) {
    
    out[which(m == v)] <- 0
    
  } else {
    
    out[which(m == v)] <- 1
    
  }
  
  out
  
}



mem_prior_mat <- function(m, mod_mat, prior_inclusion) {
  
  mem <- m[upper.tri(m)]
  
  source_vec <- prior_inclusion[upper.tri(prior_inclusion)]
  
  s_in <- source_vec # source inclusion probability
  
  s_ex <- 1 - source_vec # source exclusion probability
  
  prod(s_in^(mem) * s_ex^(1 - mem))
  
}



i_models <- function(hh, models, samp) {
  
  k <- length(models[1, ])
  
  if (hh == 1) {
    
    ii <- seq_len(k)
    
  } else if (hh == k) {
    
    ii <- c(k, seq_len(k - 1))
    
  } else {
    
    ii <- c(hh, seq_len(hh - 1), (hh + 1):k)
    
  }
  
  which(apply(models,
              
              MARGIN = 1,
              
              FUN = function(x, t) {
                
                sum(x == t)
                
              },
              
              t = samp[hh, ii]
              
  ) == k)
  
}



models_count <- function(samp, models) {
  
  out <- matrix(0, nrow(models), ncol(models))
  
  u <- vapply(seq_len(ncol(models)),
              
              FUN = i_models,
              
              FUN.VALUE = NA_integer_,
              
              models = models,
              
              samp = samp
              
  )
  
  for (i in seq_along(u)) {
    
    out[u[i], i] <- 1
    
  }
  
  out
  
}



#' @importFrom crayon red

mem_post_prob <- function(model, fit) {
  
  if (model$alternative == "greater") {
    
    out <-
      
      vapply(
        
        seq_len(ncol(fit$samples)),
        
        FUN = function(j, x, t) {
          
          return(sum(x[, j] > t[j]) / length(x[, j]))
          
        },
        
        FUN.VALUE = NA_real_,
        
        x = fit$samples,
        
        t = model$p0
        
      )
    
  } else if (model$alternative == "less") {
    
    out <-
      
      vapply(
        
        seq_len(ncol(fit$samples)),
        
        FUN = function(j, x, t) {
          
          return(sum(x[, j] < t[j]) / length(x[, j]))
          
        },
        
        FUN.VALUE = NA_real_,
        
        x = fit$samples,
        
        t = model$p0
        
      )
    
  }
  
  else {
    
    stop(red("Alternative must be either \"greater\" or \"less\"."))
    
  }
  
  
  
  names(out) <- model$name
  
  return(out)
  
}





#############################################################

### Added Shape and Direction to Posterior Prob #############

#############################################################



#' @importFrom stats pbeta

eval.Post <- function(p0, x, n, omega, w, a, b, alternative = "greater") {
  
  alpha <- a + omega %*% x
  
  beta <- b + (omega %*% n - omega %*% x)
  
  if (alternative == "greater") {
    
    out <- sum((1 - pbeta(p0, alpha, beta)) * w)
    
  } else {
    
    out <- sum((pbeta(p0, alpha, beta)) * w)
    
  }
  
  out
  
}



#############################################################

### Fixed shape paramter specification in gen_post() ########

#############################################################



#' @importFrom stats rbeta

gen_post <- function(x, n, omega, a, b) {
  
  alpha <- a + omega %*% x
  
  beta <- b + (omega %*% n - omega %*% x)
  
  rbeta(1, alpha, beta)
  
}



#' @importFrom stats rmultinom

samp_post <- function(x, n, omega, w, a, b) {
  
  return(gen_post(x, n, omega[which(rmultinom(1, 1, w) == 1), ], a, b))
  
}





sample_posterior_model <- function(model, num_samples = 100000) {
  
  ret <- replicate(
    
    num_samples,
    
    samp_post(
      
      model$responses,
      
      model$size,
      
      model$models,
      
      model$pweights[[1]],
      
      model$shape1[1],
      
      model$shape2[1]
      
    )
    
  )
  
  k <- length(model$responses)
  
  if (k > 2) {
    
    ret <- rbind(
      
      ret,
      
      foreach(j = seq_len(k - 1)[-1], .combine = rbind) %do% {
        
        ii <- c(j, seq_len(j - 1), (j + 1):k)
        
        replicate(
          
          num_samples,
          
          samp_post(
            
            model$responses[ii],
            
            model$size[ii],
            
            model$models,
            
            model$pweights[[j]],
            
            model$shape1[j],
            
            model$shape2[j]
            
          )
          
        )
        
      }
      
    )
    
  }
  
  j <- k
  
  ii <- c(j, 1:(j - 1))
  
  ret <- rbind(
    
    ret,
    
    replicate(
      
      num_samples,
      
      samp_post(
        
        model$responses[ii],
        
        model$size[ii],
        
        model$models,
        
        model$pweights[[j]],
        
        model$shape1[j],
        
        model$shape2[j]
        
      )
      
    )
    
  )
  
  ret <- t(ret)
  
  dimnames(ret) <- list(NULL, model$name)
  
  ret
  
}



#' @importFrom stats rbeta

samp_one_group <- function(x, n, a, b, num_samples = 100000) {
  
  return(rbeta(num_samples, x + a, n - x + b))
  
}



#' @importFrom stats pbeta

eval_post_one_group <- function(p0, x, n, a, b, alternative = "greater") {
  
  alpha <- a + x
  
  beta <- b + n - x
  
  if (alternative == "greater") {
    
    out <- 1 - pbeta(p0, alpha, beta)
    
  } else {
    
    out <- pbeta(p0, alpha, beta)
    
  }
  
  out
  
}



#' Cluster Baskets Based on the Posterior Exchangeabilities

#'

#' This is the default function used to cluster cohorts in the

#' \code{basket}, \code{mem_mcmc}, and \code{mem_exact} functions.

#' The approach creates a graph where each vertex is a cohort and the

#' weight between two cohorts is determined by their posterior exchangeability

#' probability. The graph is then clustered using \pkg{igraph}'s

#' \code{louvain} function, which determines the number of clusters and

#' the cluster memberships, and has been shown to perform well with

#' real clinical data.

#' @param m the adjacency matrix.

#' @return A factor variable with cluster memberships for each cohort in

#' the study.

#' @importFrom igraph graph_from_adjacency_matrix E cluster_louvain

#' @seealso basket mem_mcmc mem_exact

#' @export

cluster_membership <- function(m) {
  
  graph <-
    
    graph_from_adjacency_matrix(m,
                                
                                mode = "undirected",
                                
                                weighted = TRUE,
                                
                                diag = FALSE
                                
    )
  
  factor(cluster_louvain(graph, weights = E(graph)$weight)$membership)
  
}



cluster_comp <- function(basket, cluster_function) {
  
  n_vec <- length(basket$size)
  
  map <- basket$map
  
  name <- basket$name
  
  all_samples <- basket$samples
  
  result <- cluster_function(map)
  
  num_clusters <- length(levels(result))
  
  
  
  cluster_sample <- list()
  
  cluster_name <- c()
  
  cluster_element <- list()
  
  for (k in 1:num_clusters) {
    
    rank <- which(result == levels(result)[k])
    
    
    
    cluster_element[[k]] <- name[rank]
    
    cluster_name <- c(cluster_name, paste0("Cluster ", k))
    
    
    
    if (n_vec == 1) {
      
      sample_vector <- as.vector(all_samples)
      
    } else {
      
      sample_vector <- as.vector(all_samples[, rank])
      
    }
    
    cluster_sample[[k]] <- sample_vector
    
  }
  
  names(cluster_sample) <- cluster_name
  
  
  
  ret <- basket
  
  ret$p0 <- unique(ret$p0)
  
  ret$name <- cluster_name
  
  ret$cluster <- cluster_element
  
  ret$samples <- cluster_sample
  
  
  
  p0_test <- unique(ret$p0)
  
  all_cdf <- matrix(0, 0, num_clusters)
  
  for (kk in seq_along(p0_test)) {
    
    if (ret$alternative == "greater") {
      
      res <- unlist(lapply(
        
        seq_len(num_clusters),
        
        FUN = function(j, x, t) {
          
          return(sum(x[[j]] > t) / length(x[[j]]))
          
        },
        
        x = cluster_sample,
        
        t = p0_test[kk]
        
      ))
      
    } else {
      
      res <- unlist(lapply(
        
        seq_len(num_clusters),
        
        FUN = function(j, x, t) {
          
          return(sum(x[[j]] > t) / length(x[[j]]))
          
        },
        
        x = cluster_sample,
        
        t = p0_test[kk]
        
      ))
      
    }
    
    all_cdf <- rbind(all_cdf, res)
    
  }
  
  colnames(all_cdf) <- cluster_name
  
  rownames(all_cdf) <- p0_test
  
  ret$post_prob <- all_cdf
  
  
  
  ret$mean_est <- unlist(lapply(ret$samples, mean))
  
  ret$median_est <- unlist(lapply(ret$samples, median))
  
  ret$hpd <- matrix(unlist(
    
    lapply(
      
      ret$samples,
      
      FUN = boa_hpd,
      
      alpha = ret$alpha
      
    )
    
  ), nrow = 2, byrow = FALSE)
  
  colnames(ret$hpd) <- cluster_name
  
  row.names(ret$hpd) <- c("Lower Bound", "Upper Bound")
  
  ret$ess <- ess_from_hpd(fit = ret, alpha = ret$alpha)
  
  names(ret$ess) <- cluster_name
  
  class(ret) <- c("mem_exact", "exchangeability_model")
  
  ret
  
}





update_mh <- function(m_old,
                      
                      m,
                      
                      xvec,
                      
                      nvec,
                      
                      avec,
                      
                      bvec,
                      
                      mod_mat,
                      
                      prior_ep,
                      
                      beta_vec, old_dens, prod_vec) {
  
  k <- max(m, na.rm = TRUE) + 1
  
  v <- sample(
    
    seq_len(k - 1),
    
    seq_len(k - 1)[which(rmultinom(
      
      1, 1,
      
      (rev(seq_len(k - 1))^3) / sum(rev(seq_len(k - 1))^3)
      
    ) == 1)]
    
  )
  
  
  
  m_prop <- flip_mem(v[1], m_old, m)
  
  if (length(v) > 1) {
    
    for (ii in seq_along(v)[-1]) {
      
      m_prop <- flip_mem(v[ii], m_prop, m)
      
    }
    
  }
  
  
  
  if (is.na(old_dens)) {
    
    old_dens <-
      
      log_marg_dens_mcmc(m_old, mod_mat, xvec, nvec, avec, bvec, beta_vec,
                         
                         prod_vec) +
      
      log(mem_prior_mat(m_old, mod_mat, prior_ep))
    
  }
  
  new_dens <-
    
    log_marg_dens_mcmc(m_prop, mod_mat, xvec, nvec, avec, bvec, beta_vec, 
                       
                       prod_vec) +
    
    log(mem_prior_mat(m_prop, mod_mat, prior_ep))
  
  
  
  rho <- exp(new_dens - old_dens)
  
  
  
  if (rho >= 1) {
    
    old_dens <- new_dens
    
    out <- m_prop
    
  } else {
    
    if (rbinom(1, 1, rho) == 1) {
      
      old_dens <- new_dens
      
      out <- m_prop
      
    } else {
      
      out <- m_old
      
    }
    
  }
  
  list(out, old_dens)
  
}



log_marg_dens_mcmc <- function(m, mod_mat, xvec, nvec, avec, bvec, beta_vec,
                               
                               prod_vec) {
  
  marg_vec <- rep(NA, dim(m)[1])
  
  
  
  # calculate the product portion of integrated marginal likelihood
  
  for (i in seq_len(dim(m)[1])) {
    
    p_vec <- prod(prod_vec^(1 - m[i, ]))
    
    marg_vec[i] <-
      
      (beta(avec[i] + m[i, ] %*% xvec, bvec[i] + m[i, ] %*% (nvec - xvec)) /
         
         beta_vec[i]) * p_vec
    
  }
  
  sum(log(marg_vec))
  
}

#priormat <- matrix(c(1,0.1, 0.1, 1), ncol=2)
#responses=c(50, 50); size=c(100, 100); name=c("basket1", "basket2"); prior = priormat; #mcmc_burnin = 1000; mcmc_iter = 10000; shape1 = 1; shape2=1

#pr <- 0.5
#priormat <- matrix(c(1,pr, pr, pr, 1, pr, pr, pr, 1), ncol=3)
#responses=c(50, 55, 70); size=c(100, 100, 100); name=c("basket1", "basket2", "basket3"); prior #= priormat; mcmc_burnin = 1000; mcmc_iter = 10000; shape1 = 1;shape2=1

#pr <- 0.1
#priormat <- matrix(pr, nrow=4, ncol=4)
#diag(priormat) <- 1
#responses=c(50, 55, 50, 55); size=c(100, 100, 100, 100); name=c("basket1", "basket2", "basket3","basket4"); prior = priormat; mcmc_burnin = 1000; mcmc_iter = 10000; shape1 = 1;shape2=1

########################################################################################
########### function MAIN #######################################################
########################################################################################


mem_mcmc2 <- function (responses, size, name, p0 = 0.15, shape1 = 0.5, shape2 = 0.5, 
                       prior = diag(length(responses))/2 + matrix(0.5, nrow = length(responses), 
                                                                  ncol = length(responses)), hpd_alpha = 0.05, alternative = "greater", 
                       mcmc_iter = 2e+05, mcmc_burnin = 50000, initial_mem = round(prior - 
                                                                                     0.001), seed = 1000, cluster_analysis = FALSE, call = NULL, 
                       cluster_function = cluster_membership) 
{
  
  set.seed(seed)
  k <- NULL
  if (is.null(getDoParName())) {
    registerDoSEQ()
  }
  if (length(responses) != length(size)) {
    stop(red("The length of the responses and size parameters", 
             "must be equal."))
  }
  if (length(shape1) == 1) {
    shape1 <- rep(shape1, length(responses))
  }
  if (length(shape2) == 1) {
    shape2 <- rep(shape2, length(responses))
  }
  if (length(p0) == 1) {
    p0 <- rep(p0, length(responses))
  }
  size1 <- size[size != 0]
  alp <- hpd_alpha
  if (length(size1) < 1) {
    stop(red("The length of the responses must be equal or greater than 1"))
  }
  if (length(size1) == 1) {
    ind <- which(size != 0)
    n_vec <- length(size)
    prior_inclusion <- prior
    pep <- matrix(1, n_vec, n_vec)
    colnames(pep) <- rownames(pep) <- name
    if (n_vec > 1) {
      pep[, ind] <- 0
      pep[ind, ] <- 0
      pep[ind, ind] <- 1
    }
    maximizer <- map <- pep
    pweights <- rep(0, n_vec)
    pweights[ind] <- 1
    hpd <- matrix(NA, 2, n_vec)
    p_ess <- post_prob <- rep(NA, n_vec)
    names(p_ess) <- names(post_prob) <- name
    rownames(hpd) <- c("lower", "upper")
    colnames(hpd) <- name
    a <- shape1
    b <- shape2
    samp <- samp_one_group(responses[1], size[1], a[1], b[1])
    hpd[, 1] <- boa_hpd(samp, alp)
    p_ess[1] <- a[1] + b[1] + size[1]
    t <- eval_post_one_group(p0[1], responses[1], size[1], 
                             a[1], b[1], alternative)
    post_prob[1] <- t
    if (n_vec > 1) {
      for (i in 2:n_vec) {
        group_sample <- samp_one_group(responses[i], 
                                       size[i], a[i], b[i])
        samp <- cbind(samp, group_sample)
        hpd[, i] <- boa_hpd(group_sample, alp)
        p_ess[i] <- a[i] + b[i] + size[i]
        post_prob[i] <- eval_post_one_group(p0[i], responses[i], 
                                            size[i], a[i], b[i], alternative)
      }
      colnames(samp) <- name
    }
    if (is.null(call)) {
      call <- match.call()
    }
    ret <- list(maximizer = maximizer, prior = prior_inclusion, 
                map = map, pep = pep, post_prob = post_prob, hpd = hpd, 
                responses = responses, size = size, name = name, 
                p0 = p0, alpha = hpd_alpha, alternative = alternative, 
                pweights = pweights, shape1 = shape1, shape2 = shape2, 
                call = call)
    ret$samples <- samp
    if (n_vec > 1) {
      ret$mean_est <- colMeans(ret$samples)
      ret$median_est <- apply(ret$samples, 2, median)
    }
    else {
      ret$mean_est <- mean(ret$samples)
      ret$median_est <- median(ret$samples)
    }
    ret$ess <- p_ess
    class(ret) <- c("mem_basket", "mem")
    if (cluster_analysis) {
      cluster_ret <- cluster_comp(ret, cluster_function)
      class(cluster_ret) <- c("mem_cluster", "mem")
      result <- list(call = call, basket = ret, cluster = cluster_ret)
    }
    else {
      result <- list(call = call, basket = ret, cluster = NA)
    }
    class(result) <- c("mem_exact", "exchangeability_model")
    return(result)
  }
  if (!isTRUE(all.equal(diag(prior), rep(1, ncol(prior))))) {
    stop(red("Elements on the main diagonal of `prior` must be 1."))
  }
  mod_mat <- foreach(k = rev(seq_len(length(responses) - 1))) %do% 
    {
      mem_sample_space <- as.matrix(expand.grid(rep(list(c(0, 
                                                           1)), k)))
      mem_sample_space[order(rowSums(mem_sample_space)), 
                       ]
    }
  if (!isTRUE(all.equal(initial_mem, t(initial_mem)))) {
    stop(red("The `initial_mem` matrix must be symmetric."))
  }
  if (!isTRUE(all(diag(initial_mem) == 1))) {
    stop(red("The main diagonal of the `initial_mem` matrix must be 1's."))
  }
  m_init <- initial_mem
  m_old <- m_init
  m <- diag(NA, nrow(m_old))
  k <- 1
  for (ii in seq_len(nrow(m_old) - 1)) {
    for (jj in (ii + 1):ncol(m_old)) {
      m[ii, jj] <- m[jj, ii] <- k
      k <- k + 1
    }
  }
  n_chg <- 0
  mod_mat[[1]] <- as.matrix(mod_mat[[1]])
  models <- cbind(rep(1, dim(mod_mat[[1]])[1]), mod_mat[[1]])
  mweights <- matrix(0, nrow(models), length(responses))
  if (missing(name)) {
    name <- paste("basket", seq_along(size))
  }
  if (is.factor(name)) {
    name <- as.character(name)
  }
  colnames(mweights) <- name
  mem_samp <- list(m_old)
  mweights <- mweights + models_count(samp = mem_samp[[1]], 
                                      models = models)
  map_list <- list(mem_samp[[1]])
  map_count <- c(1)
  map_hash <- new.env()
  map_hash[[toString(m_old)]] <- 1
  old_dens <- NA
  xvec <- responses
  nvec <- size
  beta_vec <- beta(shape1, shape2)
  prod.vec <- beta(xvec + shape1, nvec + shape2 - xvec)/beta(shape1, 
                                                             shape2)
  for (i in 1:mcmc_burnin) {
    t <- update_mh(m_old, m, responses, size, shape1, shape2, 
                   mod_mat, prior, beta_vec, old_dens, prod.vec)
    m_old <- t[[1]]
    old_dens <- t[[2]]
  }
  t <- update_mh(m_old, m, responses, size, shape1, shape2, 
                 mod_mat, prior, beta_vec, old_dens, prod.vec)
  mem_samp[[2]] <- t[[1]]
  old_dens <- t[[2]]
  mweights <- mweights + models_count(samp = mem_samp[[2]], 
                                      models = models)
  samp_sum <- mem_samp[[1]] + mem_samp[[2]]
  if (sum(mem_samp[[2]] == mem_samp[[1]]) < length(mem_samp[[2]])) {
    n_chg <- n_chg + 1
  }
  new <- mem_samp[[2]]
  key <- toString(new)
  if (!is.null(map_hash[[key]])) {
    index <- map_hash[[key]]
    map_count[index] <- map_count[index] + 1
  } else {
    map_list[[length(map_list) + 1]] <- mem_samp[[2]]
    map_count <- c(map_count, 1)
    map_hash[[key]] <- length(map_list)
  }
  for (kk in seq_len(mcmc_iter)[-(1:2)]) {
    t <- update_mh(mem_samp[[kk - 1]], m, responses, size, 
                   shape1, shape2, mod_mat, prior, beta_vec, old_dens, 
                   prod.vec)
    mem_samp[[kk]] <- t[[1]]
    old_dens <- t[[2]]
  }
  it <- NULL
  models_count <- foreach(it = isplitVector(seq_len(mcmc_iter)[-(1:2)], 
                                            chunks = num_workers()), .combine = c) %dopar% {
                                              foreach(k = it) %do% {
                                                models_count(samp = mem_samp[[k]], models = models)
                                              }
                                            }
  for (kk in seq_len(mcmc_iter)[-(1:2)]) {
    mweights <- mweights + models_count[[kk - 2]]
    samp_sum <- samp_sum + mem_samp[[kk]]
    if (sum(mem_samp[[kk]] == mem_samp[[kk - 1]]) < length(mem_samp[[kk - 
                                                                     1]])) {
      n_chg <- n_chg + 1
    }
    new <- mem_samp[[kk]]
    key <- toString(new)
    if (!is.null(map_hash[[key]])) {
      index <- map_hash[[key]]
      map_count[index] <- map_count[index] + 1
    } else {
      map_list[[length(map_list) + 1]] <- mem_samp[[kk]]
      map_count <- c(map_count, 1)
      map_hash[[key]] <- length(map_list)
    }
  }
  pweights <- list()
  for (kk in seq_len(ncol(mweights))) {
    pweights[[kk]] <- mweights[, kk]/mcmc_iter
  }
  model <- list(responses = responses, size = size, name = name, 
                shape1 = shape1, shape2 = shape2, models = models, pweights = pweights, 
                p0 = p0, alpha = hpd_alpha, alternative = alternative)
  pep <- samp_sum/mcmc_iter
  rownames(pep) <- colnames(pep) <- model$name
  map <- map_list[[order(map_count, decreasing = TRUE)[1]]]
  rownames(map) <- colnames(map) <- model$name
  if (is.null(call)) {
    call <- match.call()
  }
  ret <- list(responses = responses, size = size, name = name, 
              p0 = p0, alpha = hpd_alpha, alternative = alternative, 
              shape1 = shape1, shape2 = shape2, prior = prior, call = call)
  ret$mod_mat <- mod_mat
  ret$initial <- m_init
  rownames(ret$initial) <- colnames(ret$initial) <- model$name
  ret$models <- models
  ret$pweights <- pweights
  ret$samples <- sample_posterior_model(model)
  ret$accept_rate <- (n_chg)/mcmc_iter
  ret$mean_est <- colMeans(ret$samples)
  ret$median_est <- apply(ret$samples, 2, median)
  ret$pep <- pep
  ret$map <- map
  ret$hpd <- apply(ret$samples, MARGIN = 2, FUN = boa_hpd, 
                   alpha = model$alpha)
  ret$post_prob <- mem_post_prob(model, fit = ret)
  ret$ess <- ess_from_hpd(fit = ret, alpha = model$alpha)
  names(ret$ess) <- model$name
  class(ret) <- c("mem_basket", "mem")
  if (cluster_analysis) {
    cluster_ret <- cluster_comp(ret, cluster_function)
    class(cluster_ret) <- c("mem_cluster", "mem")
    result <- list(call = call, basket = ret, cluster = cluster_ret, 
                   seed = seed, map_list=map_list, mem_samp=mem_samp)
  }
  else {
    result <- list(call = call, basket = ret, cluster = NA, 
                   seed = seed, map_list=map_list, mem_samp=mem_samp)
  }
  class(result) <- c("mem_mcmc", "exchangeability_model")
  result
}