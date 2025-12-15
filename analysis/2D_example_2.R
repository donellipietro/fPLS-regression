
graphics.off()
rm(list = ls())

library(MASS)
source("src/utils/pls.R")

# classification
library(pROC)

set.seed(1)
N <- 300
# X1 <- mvrnorm(N, c(50, 160), matrix(c(10, 4, 4, 2), ncol = 2, byrow = TRUE))
# X2 <- mvrnorm(N, c(75, 180), matrix(c(2, 4, 4, 10), ncol = 2, byrow = TRUE))

X1 <- mvrnorm(N, c(50, 160), matrix(c(1, 2, 2, 10), ncol = 2, byrow = TRUE))
X2 <- mvrnorm(N, c(52, 160), matrix(c(1, 2, 2, 10), ncol = 2, byrow = TRUE))

X <- rbind(X1, X2)
Y <- c(rep(1, N), rep(-1, N))
Ycl <- c(rep("Gruppo_1", N), rep("Gruppo_2", N))
colors <- c(rep("darkgreen", N), rep("purple", N))

# Data exploration
plot(NaN, xlim = range(X[,1]), ylim = range(X[,2]),
     xlab = "Weight", ylab = "Height", main = "Data", asp = 1)
grid()
points(X, col = colors)

# Data centering
Xc <- X - rep(1, N*2) %*% t(colMeans(X))
# Xc[1:N, ] <- X[1:N, ] - rep(1, N) %*% t(colMeans(X[1:N, ]))
# Xc[(N+1):(2*N), ] <- X[(N+1):(2*N), ] - rep(1, N) %*% t(colMeans(X[(N+1):(2*N), ]))
plot(NaN, xlim = range(Xc[,1]), ylim = range(Xc[,2]),
     xlab = "Weight", ylab = "Height", main = "Centered Data", asp = 1)
grid()
points(Xc, col = colors)

# Fit model
model_PLS <- PLS(Y, Xc, center = TRUE, n_comp = 2)

# Classifier
Y_hat <- model_PLS$Y_hat[[2]]
roc_curve <- roc(Ycl, as.vector(Y_hat))
optimal_coords <- coords(roc_curve, "best", ret = c("threshold", "sensitivity", "specificity"))
optimal_treshold <- optimal_coords$threshold[1]

# Make sure x-axis can accommodate both sets of data
df <- data.frame(
  s = Y_hat[, 1],
  group = Ycl
)
ggplot(df, aes(x = s, color = group, fill = group)) +
  geom_density(alpha = 0.5) +
  # match your colors
  scale_color_manual(values = c("darkgreen", "purple")) +
  scale_fill_manual(values = c("darkgreen", "purple")) +
  # same x-range, same y-limit as in your base R plot
  # coord_cartesian(xlim = range_all, ylim = c(0, 1)) +
  labs(x = "", y = "") +
  theme_minimal()  + 
  geom_vline(xintercept = optimal_treshold, linetype = "dashed")

# Latent space
C <- model_PLS$Y_loadings
TT <- model_PLS$X_latent_scores
plot(NaN, xlim = range(TT[, 1]), ylim = range(TT[, 2]),
     xlab = "t1", ylab = "t2", main = "Latent space", asp = 1)
grid()
points(TT, col = colors)
abline(a = optimal_treshold/C[2], b = -C[1]/C[2], lty = 2)

# Transormation
generate_basis <- function(M) {
  # Input: A vector v from R^d
  # Output: An orthonormal basis for R^d containing v
  
  v <- M[, 1]
  d <- length(v) # Dimension of the space
  basis <- matrix(0, nrow = d, ncol = d) # To store the basis vectors
  
  # Normalize the initial vector
  basis[, 1] <- v / sqrt(sum(v^2))
  
  # Generate additional vectors orthogonal to the initial vector
  for (i in 2:d) {
    # Create a random vector
    new_vec <- M[, i]
    
    # Orthogonalize against all previous basis vectors
    for (j in 1:(i - 1)) {
      new_vec <- new_vec - sum(basis[, j] * new_vec) * basis[, j]
    }
    
    # Normalize the new vector
    basis[, i] <- new_vec / sqrt(sum(new_vec^2))
  }
  
  return(basis)
}
Q <- generate_basis(cbind(t(C), c(0, 1)))
SS <- as.matrix(TT[,1:2]) %*% Q + rep(1, nrow(X)) %*% t(c(- optimal_treshold / (Q %*% t(C))[1], 0))
colnames(SS) <- c("s1", "s2")

plot(NaN, xlim = range(SS[, 1]), ylim = range(SS[, 2]),
     xlab = "s1", ylab = "s2", main = "Transformed Latent space", asp = 1)
grid()
# abline(h = 0, lty = 2)
abline(v = 0, lty = 2)
points(SS, col = colors)


# Centered data and loadings
R <- model_PLS$X_loadings
W <- model_PLS$X_space_directions
r1 <- R[, 1]
r2 <- R[, 2]
w1 <- W[, 1]
w2 <- W[, 2]
plot(NaN, xlim = range(Xc[,1]), ylim = range(Xc[,2]),
     xlab = "Weight", ylab = "Height", main = "Centered data", asp = 1)
grid()
points(Xc, col = colors)
abline(a = 0, b = r1[2]/r1[1], lty = 2)
abline(a = 0, b = r2[2]/r2[1], lty = 3)
#abline(a = 0, b = w1[2]/w1[1], lty = 2)
#abline(a = 0, b = w2[2]/w2[1], lty = 3)

# Comparison PCA
model_PCA <- princomp(Xc)
V <- model_PCA$loadings
v1 <- V[, 1]
v2 <- V[, 2]
plot(NaN, xlim = range(Xc[,1]), ylim = range(Xc[,2]),
     xlab = "Weight", ylab = "Height", main = "Centered data", asp = 1)
grid()
points(Xc, col = colors)
abline(a = 0, b = r1[2]/r1[1], lty = 2)
abline(a = 0, b = r2[2]/r2[1], lty = 3)
abline(a = 0, b = v1[2]/v1[1], lty = 2, col = "blue")
abline(a = 0, b = v2[2]/v2[1], lty = 3, col = "blue")

svd(t(Y)%*%Xc)


alphas <- seq(0,pi, length = 100)
TT_alpha <- NULL
cov <- NULL
for(alpha in alphas){
  v <- c(cos(alpha), sin(alpha))
  t <- Xc %*% v
  TT_alpha <- cbind(TT_alpha, t)
  cov <- c(cov, (t(t)%*%Y)^2)
}

par(mfrow = c(2, 1), mai = c(0.5, 1, 0.5, 1))
matplot(NaN, xlim = c(0,180), ylim = range(TT_alpha),  xlab = "", ylab = "scores")
matlines(alphas/(2*pi)*360, t(TT_alpha), col = colors, type = "l")
abline(v = atan(w1[2]/w1[1])/(2*pi)*360)
grid()
matplot(alphas/(2*pi)*360, cov, type = "l", xlab = "alpha", ylab = "cov")
abline(v = atan2(w1[2], w1[1])/(2*pi)*360)
grid()
