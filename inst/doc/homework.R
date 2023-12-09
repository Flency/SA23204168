## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
n <- 1000
u <- runif(n)
index <- which(u >= 0.5)
x <- c(-log(2-2 * u[index]), log(2 * u[-index]))
x_neg<-seq(-10,0,0.1)
x_pos<-seq(0,10,0.1)
hist(x,freq = FALSE)
lines(x_neg, 0.5 * exp(x_neg))
lines(x_pos, 0.5 * exp(-x_pos))

## -----------------------------------------------------------------------------
rbeta <- function(n, a, b) {
n <- 1000
k <- 0
y <- numeric(n)
while (k < n) {
  u <- runif(1)
  x <- runif(1)
    if (x^(a - 1) * (1 - x)^(b - 1) > u) {
      k <- k + 1
      y[k] <- x
    }
 }
return(y)
} 

## -----------------------------------------------------------------------------
beta_sample<-rbeta(1000,3,2)
hist(beta_sample,freq = F)
y <- seq(0, 1, 0.01)
fy <- 12 * y^2 * (1 - y)
lines(y, fy)

## -----------------------------------------------------------------------------
data_gen<-function(n){
  x<-rep(0,n)
  u1<-runif(n,-1,1)
  u2<-runif(n,-1,1)
  u3<-runif(n,-1,1)
  u2index<-which(abs(u3)>=abs(u2)&abs(u3)>=abs(u1))
  u3index<-setdiff(1:n,u2index)
  x[u3index]<-u3[u3index]
  x[u2index]<-u2[u2index]
  return(x)
}

x<-data_gen(10000)
hist(x,freq = F)
a <- seq(-1, 1, 0.001)
fe <- 0.75 * (1 - a^2)
lines(a, fe)

## -----------------------------------------------------------------------------
mysample<-function(data,sample_size,prob){
  
  # 放回
  # 使用逆变换法生成有放回抽样的索引
  sampled_indices <- integer(sample_size)
  # 生成[0, 1)之间的均匀随机数
  U <- runif(sample_size)
  # 计算累积概率
  cumulative_prob <- cumsum(prob)
  # 找到所在区间
  index <- findInterval(U,cumulative_prob)+1
  # 根据索引获取抽样的数据
  sampled_data <- data[index]
  return(sampled_data)
}

# 测试
# 准备数据集
data <- c("A", "B", "C", "D", "E")
# 指定每个元素的抽样概率
probabilities <- c(0.2, 0.3, 0.1, 0.2, 0.2)


# 放回

# 指定抽样大小
sample_size <- 3
mysample(data=data,sample_size = sample_size,prob = probabilities)


## -----------------------------------------------------------------------------
# Calculate the variance of the estimator
pihat.var <- function(rho){
  d <- 1
  l <- d*rho
  n <- 1e6
  M=100 #repeat 100 times
  pihat <- numeric(M)
  for (i in 1:M) {
    X <- runif(n,0,d/2)
    Y <- runif(n,0,pi/2)
    pihat[i] <- 2*l/d/mean(l/2*sin(Y)>X)
  }
  var(pihat)
}

## -----------------------------------------------------------------------------
# rho=0.5
pihat.var(0.5)
# rho=0.8
pihat.var(0.8)
# rho=1
pihat.var(1)

## -----------------------------------------------------------------------------
m <- 10000
simple_mc<-function(){
  m<-1e5
  times<-1e3
  simple_mc<-numeric(times)
  for (i in 1:times) {
    U<-runif(m)
    simple_mc[i]<-mean(exp(U))
  }
  simple_mc
}

anti_variate<-function(){
  m<-1e5
  times<-1e3
  anti_variate<-numeric(times)
  for (i in 1:times) {
    U<-runif(m/2)
    anti_variate[i]<-mean((exp(U)+exp(1-U))/2)
  }
  anti_variate
}

simple_mc_result<-simple_mc()
anti_variate_result<-anti_variate()

# Estimates from both methods
c(mean(simple_mc_result), mean(anti_variate_result))
# Variance of the two methods
c(var(simple_mc_result), var(anti_variate_result))
# Proportion of variance reduction
(var(simple_mc_result) - var(anti_variate_result))/var(simple_mc_result)

## -----------------------------------------------------------------------------
x <- seq(1, 10, 0.01)
y <- x^2 * exp(-x^2/2)/sqrt(2 * pi)
plot(x, y, type = "l",ylim=c(0,1))
lines(x, 2*dt(x-1,df = 1), lty = 2)
lines(x, 2*dnorm(x,mean=1), lty = 3)
legend("topright", legend = c("g(x)", "f1","f2"), lty = 1:3)

## -----------------------------------------------------------------------------
plot(x, y/(2*dt(x-1,df = 1)), type = "l", lty = 2,ylab = "")
lines(x, y/(2*dnorm(x, mean=1)), lty = 3)
legend("topright",  legend = c("f1", "f2"),lty = 2:3)

## -----------------------------------------------------------------------------
m <- 10000

Importance_sampling1 <- replicate(1000, expr = {
x <- abs(rt(m,df=1))+1
f <- 2*dt(x-1,df=1)
g <- x^2 * exp(-x^2/2)/sqrt(2 * pi)
mean(g/f)
})
Importance_sampling2 <- replicate(1000, expr = {
x <- abs(rnorm(m)) + 1
f <- 2*dnorm(x, mean=1)
g <- x^2 * exp(-x^2/2)/sqrt(2 * pi)
mean(g/f)
})

c(mean(Importance_sampling1), mean(Importance_sampling2))
c(var(Importance_sampling1), var(Importance_sampling2))
var(Importance_sampling2)/var(Importance_sampling1)


## -----------------------------------------------------------------------------
M <- 10000
k <- 5 #层数
m <- M/k #每层样本
theta.hat <- numeric(k)
se <- numeric(k)
g <- function(x) exp(-x)/(1 + x^2)
f <- function(x) (k/(1 - exp(-1))) * exp(-x)
for (j in 1:k) {
u <- runif(m, (j - 1)/k, j/k)
x <- -log(1 - (1 - exp(-1)) * u)
fg <- g(x)/f(x)
theta.hat[j] <- mean(fg)
se[j] <- sd(fg)
}
sum(theta.hat)
mean(se)

## -----------------------------------------------------------------------------
k <- 1
m <- M/k
theta.hat <- numeric(k)
se <- numeric(k)
for (j in 1:k) {
u <- runif(m, (j - 1)/k, j/k)
x <- -log(1 - (1 - exp(-1)) * u)
fg <- g(x)/f(x)
theta.hat[j] <- mean(fg)
se[j] <- sd(fg)
}
sum(theta.hat)
mean(se)

## -----------------------------------------------------------------------------
n <- 20
t0 <- qt(c(0.025, 0.975), df = n - 1)
CI <- replicate(10000, expr = {
x <- rchisq(n, df = 2)
ci <- mean(x) + t0 * sd(x)/sqrt(n)
})
LCL <- CI[1, ]
UCL <- CI[2, ]
mean(LCL < 2 & UCL > 2)


## -----------------------------------------------------------------------------
num_simulations <- 10000
alpha <- 0.05
population_mean <- 1

# 初始化
rejection_count <- rep(0,3)

# Monte Carlo模拟
for (i in 1:num_simulations) {
  # (i) χ²(1)分布
  sample_data_chisq <- rchisq(30, df = 1)
  # (ii) Uniform(0,2)分布
  sample_data_unif <- runif(30, min = 0, max = 2)
  # (iii) 指数分布
  sample_data_exp <- rexp(30, rate = 1)
  
  # t检验
  t_test_result_chisq <- t.test(sample_data_chisq, mu = population_mean, alternative = "two.sided")
  t_test_result_unif <- t.test(sample_data_unif, mu = population_mean, alternative = "two.sided")
  t_test_result_exp <- t.test(sample_data_exp, mu = population_mean, alternative = "two.sided")
  # 判断是否拒绝零假设
  if (t_test_result_chisq$p.value < alpha) {
    rejection_count[1] <- rejection_count[1] + 1
  }
  if (t_test_result_unif$p.value < alpha) {
    rejection_count[2] <- rejection_count[2] + 1
  }
  if (t_test_result_exp$p.value < alpha) {
    rejection_count[3] <- rejection_count[3] + 1
  }
}

# 计算经验I型错误率
empirical_type_1_error_rate <- rejection_count / num_simulations
cat("Empirical Type I Error Rate:", empirical_type_1_error_rate, "\n")

## -----------------------------------------------------------------------------
# 设置模拟参数
m <- 1000
alpha <- .1 # 显著性水平
# 模拟一千次
rep<- replicate(1000,{
  H0 <- runif(.95*m,0,1)# H0假设的p值，占总数的95%
  H1 <- rbeta(.05*m,.1,1)# H1假设的p值，占总数的5%
  p <- c(H0,H1)
  p.adj1 <- p.adjust(p,method='BH')# Benjamini-Hochberg 校正
  p.adj2 <- p.adjust(p,method='bonferroni')# Bonferroni 校正

  FWER <- FDR <- TPR <- numeric(2)
  FDR[1] <- sum(p.adj1[1:950]<alpha)/sum(p.adj1<alpha)
  FDR[2] <- sum(p.adj2[1:950]<alpha)/sum(p.adj2<alpha)
  FWER[1] <- (sum(p.adj1[1:950]<alpha)>0)
  FWER[2] <- (sum(p.adj2[1:950]<alpha)>0)
  TPR[1] <- sum(p.adj1[-(1:950)]<alpha)/50
  TPR[2] <- sum(p.adj2[-(1:950)]<alpha)/50
  c(FWER,FDR,TPR)
})
result <- apply(rep,1,FUN=mean)

A <- matrix(round(result,3),ncol=3,nrow = 2)
colnames(A)<- c("FWER","FDR","TPR")
rownames(A)<- c("B-H","Bonferroni")
print(A)

## -----------------------------------------------------------------------------
# 设置参数
true_lambda <- 2          # 真实参数 lambda
sample_sizes <- c(5, 10, 20)  # 样本大小
B <- 1000                 # Bootstrap重抽样次数
m <- 1000                 # 用于Bootstrap的样本数量

# 创建空数组以存储结果
mean_bias_results <- numeric(length(sample_sizes))      # 均值偏差结果
standard_error_results <- numeric(length(sample_sizes))  # 标准误差结果

# 模拟循环
for (i in 1:length(sample_sizes)) {
  n <- sample_sizes[i]   # 当前样本大小
  biases <- numeric(m)   # 存储每次模拟的均值偏差
  standard_errors <- numeric(m)  # 存储每次模拟的标准误差
  
  for (j in 1:m) {
    # 从指数分布生成一个随机样本
    data <- rexp(n, rate = 1/true_lambda)
    sample_mean <- mean(data)
    
    # 计算 lambda 的极大似然估计
    lambda_hat <- 1 / sample_mean
    
    # Bootstrap重抽样
    bootstrap_estimates <- replicate(B, {
      bootstrap_sample <- sample(data, replace = TRUE)
      lambda_b <- 1 / mean(bootstrap_sample)
      lambda_b
    })
    
    # 计算偏差和标准误差
    bias <- mean(bootstrap_estimates) - lambda_hat
    std_error <- sd(bootstrap_estimates)
    
    biases[j] <- bias
    standard_errors[j] <- std_error
  }
  
  mean_bias_results[i] <- mean(biases)
  standard_error_results[i] <- mean(standard_errors)
}

# 计算理论值
theoretical_biases <- true_lambda / (sample_sizes - 1)
theoretical_standard_errors <- true_lambda * sqrt(sample_sizes) / ((sample_sizes - 1) * sqrt(sample_sizes - 2))

# 创建一个数据框以存储结果
results_df <- data.frame(
  Sample_Size = sample_sizes,
  Simulated_Bias = mean_bias_results,
  Theoretical_Bias = theoretical_biases,
  Simulated_Standard_Error = standard_error_results,
  Theoretical_Standard_Error = theoretical_standard_errors
)
print(results_df)

## -----------------------------------------------------------------------------
# 加载必要的包
library(boot)
library(bootstrap)
attach(law)  # 将数据附加到工作环境中

# 自定义函数，用于计算每个Bootstrap样本的相关性统计量和方差估计
cor.stat <- function(x, i = 1:NROW(x)) {
  cor(x[i, 1], x[i, 2])  # 计算相关性统计量
}

# 自定义函数，用于计算相关性统计量的方差估计
cor.stat2 <- function(x, i = 1:NROW(x)) {
  # 在每个Bootstrap样本上执行自定义函数 cor.stat
  o <- boot(x[i, ], cor.stat, R = 25)
  n <- length(i)  # 样本大小
  # 返回相关性统计量和方差估计
  c(o$t0, var(o$t) * (n - 1) / n^2)
}

# 使用boot函数执行Bootstrap分析
b <- boot(law, statistic = cor.stat2, R = 1000)

# 利用boot.ci函数计算Bootstrap置信区间，使用"stud"类型，方差默认在返回的Bootstrap对象的第二个位置
result <- boot.ci(b, type = "stud")
result
detach(law)

## -----------------------------------------------------------------------------

attach(aircondit)
data <- aircondit[1]

# 定义一个函数meant，用于计算Bootstrap样本的均值
mean.boot <- function(x, i) return(mean(as.matrix(x[i, ])))

# 运行Bootstrap方法，R参数指定Bootstrap采样的次数
result <- boot(data, statistic = mean.boot, R = 2000)

# Bootstrap估计结果
result

# 计算不同类型的Bootstrap置信区间，这里包括"norm"、"perc"、"basic"和"bca"
boot.ci(result, type = c("norm", "perc", "basic", "bca"))
detach(aircondit)

## -----------------------------------------------------------------------------
attach(scor)
data <- scor
n <- nrow(data)

# 创建一个空的数值向量theta.jack，用于存储Jackknife估计的结果
theta.jack <- numeric(n)

# 计算数据集的特征值（eigenvalues）
eigenvalues <- eigen(cov(data))$values

# 计算全样本数据的估计值theta.hat
theta.hat <- eigenvalues[1] / sum(eigenvalues)

# 开始进行Jackknife估计的循环，逐一排除每个样本点
for (i in 1:n) {
  # 计算删除一个观测后的数据集的特征值
  eigenvalues <- eigen(cov(data[-i, ]))$values

  # 计算Jackknife估计的结果，并将其存储在theta.jack向量中
  theta.jack[i] <- eigenvalues[1] / sum(eigenvalues)
}

# 计算Jackknife估计的偏差（bias.jack）
bias.jack <- (n - 1) * (mean(theta.jack) - theta.hat)

# 计算Jackknife估计的标准误差（se.jack）
se.jack <- sqrt((n - 1) / n * sum((theta.jack - mean(theta.jack))^2))

# 估计结果：估计值(theta.hat)，偏差(bias.jack)，和标准误差(se.jack)
cat("theta hat:", theta.hat, "\n bias:", bias.jack, "\n se:", se.jack)

detach(scor)

## -----------------------------------------------------------------------------
library(DAAG)
attach(ironslag)
n <- length(magnetic)
e1 <- e2 <- e3 <- e4 <- numeric(choose(n, 2))

# for n-fold cross validation
# fit models on leave-one-out samples
t <- 1
for (i in 1:(n - 1)) for (j in (i + 1):n) {
k <- c(i, j)
y <- magnetic[-k]
x <- chemical[-k]
J1 <- lm(y ~ x)
yhat1 <- J1$coef[1] + J1$coef[2] * chemical[k]
e1[t] <- sum((magnetic[k] - yhat1)^2)
J2 <- lm(y ~ x + I(x^2))
yhat2 <- J2$coef[1] + J2$coef[2] * chemical[k] +
J2$coef[3] * chemical[k]^2
e2[t] <- sum((magnetic[k] - yhat2)^2)
J3 <- lm(log(y) ~ x)
logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[k]
yhat3 <- exp(logyhat3)
e3[t] <- sum((magnetic[k] - yhat3)^2)
J4 <- lm(log(y) ~ log(x))
logyhat4 <- J4$coef[1] + J4$coef[2] * log(chemical[k])
yhat4 <- exp(logyhat4)
e4[t] <- sum((magnetic[k] - yhat4)^2)
t <- t + 1
}
c(mean(e1), mean(e2), mean(e3), mean(e4))
detach(ironslag)


## -----------------------------------------------------------------------------
attach(chickwts)
x <- sort(as.vector(weight[feed == "soybean"]))  # x 包含 "soybean" 饲料类型的体重数据
y <- sort(as.vector(weight[feed == "linseed"]))  # y 包含 "linseed" 饲料类型的体重数据
detach(chickwts)

# 初始化参数
R <- 999  # 模拟次数
n <- length(x)  # 样本 x 的大小
m <- length(y)  # 样本 y 的大小
z <- c(x, y)  # 将两个样本合并为一个样本
N <- n + m  # 合并后样本的大小
K <- Fn <- Gm <- numeric(N)  

# 计算原始样本的统计量
for (i in 1:N) {
  Fn[i] <- mean((z[i] <= x))  # 计算 Fn
  Gm[i] <- mean((z[i] <= y))  # 计算 Gm
}
C0 <- ((n * m) / N) * sum((Fn - Gm)^2)  # 计算原始样本的统计量

# 执行 R 次模拟
C <- replicate(R, {
  k <- sample(N, size = n, replace = FALSE)  # 随机抽取样本 x 中的数据点
  x1 <- z[k]  # 根据抽取的索引构建新的 x1
  y1 <- z[-k]  # 剩余的数据构成 y1
  for (i in 1:N) {
    Fn[i] <- mean((z[i] <= x1))  # 计算 Fn
    Gm[i] <- mean((z[i] <= y1))  # 计算 Gm
  }
  ((n * m) / N) * sum((Fn - Gm)^2)  # 计算模拟样本的统计量
})

# 输出原始样本的统计量和 p-value
cat("statistic", C0, "\n")
cat("p.value", mean(c(C0, C) >= C0))  # 计算 p-value


## -----------------------------------------------------------------------------
maxout <- function(x,y){
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  return(max(c(outx,outy)))
}

maxout_revised <- function(x, y, R = 999) {
  n <- length(x)  # 获取第一个样本的大小
  m <- length(y)
  z <- c(x, y)  # 合并两个样本的数据
  N <- n+m  # 合并后的总样本大小
  M0 <- maxout(x, y)  # 计算原始数据的统计量

  # 使用 replicate 函数执行 R 次模拟
  M <- replicate(R, {
    k <- sample(N, size = n, replace = FALSE)  # 随机抽取索引，保持原始样本大小
    x1 <- z[k]  # 基于抽取的索引构建新的 x1 样本
    y1 <- z[-k]  # 剩余的数据构成新的 y1 样本
    maxout(x1, y1)  # 计算新样本的统计量
  })

  # 计算 p-value，表示原始数据的统计量在模拟中的分布情况
  p.value <- mean(c(M0, M) >= M0)

  # 返回估计值和 p-value
  return(list(estimate = M0, p.value = p.value))
}



## -----------------------------------------------------------------------------
# 方差相同
set.seed(0)
n1 <- 100
n2 <- 200
mu1 <- mu2 <- 0
sigma1 <- sigma2 <- 1
x <- rnorm(n1, mu1, sigma1)
y <- rnorm(n2, mu2, sigma2)
maxout_revised(x, y)

## -----------------------------------------------------------------------------
# 方差不同
set.seed(0)
n1 <- 100
n2 <- 200
mu1 <- mu2 <- 0
sigma1 <- 1
sigma2 <- 1.5
x <- rnorm(n1, mu1, sigma1)
y <- rnorm(n2, mu2, sigma2)
maxout_revised(x, y)

## -----------------------------------------------------------------------------
# 定义函数计算 alpha
alpha <- function(N, b1, b2, b3, f0) {
  # 生成符合分布的随机样本
  x1 <- rpois(N, 1)       # X1 服从泊松分布
  x2 <- rexp(N)           # X2 服从指数分布
  x3 <- rbinom(N, 1, 0.5)  # X3 服从二项分布

  # 定义函数 g
  g <- function(alpha) {
    tmp <- exp(-alpha - b1 * x1 - b2 * x2 - b3 * x3)
    p <- 1 / (1 + tmp)
    mean(p) - f0
  }

  # 使用 uniroot 求解方程
  solution <- uniroot(g, c(-20, 10))

  # 返回解的 alpha
  return(alpha = solution$root)
}

## -----------------------------------------------------------------------------
N <- 10^6
b1 <-0
b2 <- 1
b3 <- -1
f0 <- c(.1,.01,.001,.0001)
alpha(N,b1,b2,b3,f0[1])
alpha(N,b1,b2,b3,f0[2])
alpha(N,b1,b2,b3,f0[3])
alpha(N,b1,b2,b3,f0[4])

## -----------------------------------------------------------------------------
#生成 f0 的序列
f0 <- seq(0.0001, 0.99, length = 100)

# 初始化存储 alphas 的向量
alphas <- numeric(length(f0))

# 循环计算 alpha
for (i in 1:length(f0)) {
  alphas[i] <- alpha(N, b1, b2, b3, f0[i])
}

# 绘制图形
plot(-log(f0), alphas, xlab = "-log(f0)", ylab = "alpha", type = "l", main = "Plot of -log(f0) vs alpha")

## -----------------------------------------------------------------------------
set.seed(1)
# Metropolis采样函数
rw.Metropolis <- function(sigma, x0, N) {
    # 初始化向量存储采样结果
    x <- numeric(N)
    x[1] <- x0
    # 初始化接受次数
    k <- 0
    # 生成均匀分布的随机数
    u <- runif(N)
    # Metropolis采样过程
    for (i in 2:N) {
        # 从正态分布生成候选样本
        y <- rnorm(1, x[i-1], sigma)
        # 判断是否接受候选样本
        if (u[i] <= exp(abs(x[i-1]) - abs(y))) {
            x[i] <- y  
            k <- k + 1
        } else {
            x[i] <- x[i-1]
        }
    }
    # 返回结果
    return(list(x = x, k = k))
}

N <- 2000
sigma <- c(.05, .5, 2.5,  16)

x0 <- 25
rw1 <- rw.Metropolis(sigma[1], x0, N)
rw2 <- rw.Metropolis(sigma[2], x0, N)
rw3 <- rw.Metropolis(sigma[3], x0, N)
rw4 <- rw.Metropolis(sigma[4], x0, N)

#计算每个链条的接受率
data.frame(sigma=sigma,accept.rate=c(rw1$k, rw2$k, rw3$k, rw4$k)/N)

## -----------------------------------------------------------------------------
rline <- c(-qexp(.975),qexp(.975))
rw <- cbind(rw1$x, rw2$x, rw3$x,rw4$x)
for (j in 1:4) {
    plot(rw[, j], type = "l", xlab = bquote(sigma == .(round(sigma[j], 3))), ylab = "X", ylim = range(rw[, j]))
    abline(h = rline)
}

## -----------------------------------------------------------------------------
# 生成概率点
p <- ppoints(200)

# 生成Laplace分布的分位数
y <- qexp(p, 1)
z <- c(-rev(y), y)

# 从每个Metropolis采样结果中提取一部分样本
y1 <- rw1$x[501:N]
y2 <- rw2$x[501:N]
y3 <- rw3$x[501:N]
y4 <- rw4$x[501:N]


# 绘制Q-Q图
Q1 <- quantile(y1, p)
qqplot(z, Q1, main = bquote(sigma == .(round(sigma[1], 3))), xlab = "Laplace Quantiles", ylab = "Sample Quantiles")
abline(0, 1)

Q2 <- quantile(y2, p)
qqplot(z, Q2, main = bquote(sigma == .(round(sigma[2], 3))), xlab = "Laplace Quantiles", ylab = "Sample Quantiles")
abline(0, 1)

Q3 <- quantile(y3, p)
qqplot(z, Q3, main = bquote(sigma == .(round(sigma[3], 3))), xlab = "Laplace Quantiles", ylab = "Sample Quantiles")
abline(0, 1)

Q4 <- quantile(y4, p)
qqplot(z, Q4, main = bquote(sigma == .(round(sigma[4], 3))), xlab = "Laplace Quantiles", ylab = "Sample Quantiles")
abline(0, 1)


## -----------------------------------------------------------------------------
# 初始化常量和参数
N <- 5000    # 链的长度
burn <- 1000  # 预烧期长度
X <- matrix(0, N, 2)  
rho <- 0.9   # 相关性
mu1 <- 0
mu2 <- 0
sigma1 <- 1
sigma2 <- 1
s1 <- sqrt(1 - rho^2) * sigma1
s2 <- sqrt(1 - rho^2) * sigma2

# 生成链
X[1, ] <- c(mu1, mu2)  # 初始化
for (i in 2:N) {
  x2 <- X[i - 1, 2]
  m1 <- mu1 + rho * (x2 - mu2) * sigma1 / sigma2
  X[i, 1] <- rnorm(1, m1, s1)
  x1 <- X[i, 1]
  m2 <- mu2 + rho * (x1 - mu1) * sigma2 / sigma1
  X[i, 2] <- rnorm(1, m2, s2)
}
b <- burn + 1
x <- X[b:N, ]

# 绘制图形
plot(x[, 1], type = 'l', col = 1, lwd = 2, xlab = 'Index', ylab = 'Random numbers')
lines(x[, 2], col = 2, lwd = 2)
legend('bottomright', c(expression(X), expression(Y)), col = 1:2, lwd = 2)
X <- x[, 1]
Y <- x[, 2]
fit <- lm(Y ~ X)
summary(fit)

## -----------------------------------------------------------------------------
plot(fit)

## -----------------------------------------------------------------------------
Gelman.Rubin <- function(psi) {
# psi[i,j] is the statistic psi(X[i,1:j])
# for chain in i-th row of X
psi <- as.matrix(psi)
n <- ncol(psi)
k <- nrow(psi)
psi.means <- rowMeans(psi) #row means
B <- n * var(psi.means) #between variance est.
psi.w <- apply(psi, 1, "var") #within variances
W <- mean(psi.w) #within est.
v.hat <- W*(n-1)/n + (B/n) #upper variance est.
r.hat <- v.hat / W #G-R statistic
return(r.hat)
}

# 定义概率密度函数
f <- function(x, sigma) {
  if (x < 0)
    return(0)
  stopifnot(sigma > 0)
  return((x/sigma^2) * exp(-x^2/(2 * sigma^2)))
}

# Metropolis-Hastings链生成函数
Rayleigh.MH.chain <- function(sigma, m, x0) {
  x <- numeric(m)
  x[1] <- x0
  u <- runif(m)
  for (i in 2:m) {
    xt <- x[i - 1]
    y <- rchisq(1, df = xt)
    num <- f(y, sigma) * dchisq(xt, df = y)
    den <- f(xt, sigma) * dchisq(y, df = xt)
    if (u[i] <= num/den)
      x[i] <- y
    else x[i] <- xt
  }
  return(x)
}

# 设定参数
sigma <- 4
x0 <- c(1/sigma^2, 1/sigma, sigma^2, sigma^3)
k <- 4
m <- 2000
X <- matrix(0, nrow = k, ncol = m)

# 生成Metropolis-Hastings链
for (i in 1:k) X[i, ] <- Rayleigh.MH.chain(sigma, m, x0[i])

# 计算路径的累积均值
psi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(psi)) psi[i, ] <- psi[i, ] / (1:ncol(psi))

# 计算 Gelman-Rubin 统计量
rhat <- Gelman.Rubin(psi)

# 打印结果
cat("Gelman-Rubin 统计量：", rhat, "\n")


## ----warning=FALSE------------------------------------------------------------
# 导入 coda 包
library(coda)

# 将每个链转换为 coda::mcmc 对象
X1 <- as.mcmc(X[1, ])
X2 <- as.mcmc(X[2, ])
X3 <- as.mcmc(X[3, ])
X4 <- as.mcmc(X[4, ])

# 将所有链组成 mcmc.list
Y <- mcmc.list(X1, X2, X3, X4)

# 打印 Gelman-Rubin 统计量
print(gelman.diag(Y))

# 绘制 Gelman-Rubin 统计量的图形
gelman.plot(Y, col = c(1, 1))


## -----------------------------------------------------------------------------
# 区间观测数据
data <- matrix(c(11, 8, 27, 13, 16, 0, 23, 10, 24, 2,
                 12, 9, 28, 14, 17, 1, 24, 11, 25, 3), ncol = 2, byrow = FALSE)
logL <- function(lambda) {
  log_likelihood<-sum((data[,2]*exp(-lambda*data[,2])-data[,1]*exp(-lambda*data[,1]))/(exp(-lambda*data[,1])-exp(-lambda*data[,2])))
  return(log_likelihood)
}

solution <- uniroot(logL, c(0, 20))

print(paste("Lambda 的 MLE 数值解为:", round(solution$root,4)))

## -----------------------------------------------------------------------------
# E-M 算法
n <- nrow(data)
E_step <- function(lambda){
  sum(lambda*(data[,1]*exp(-lambda*data[,1])-data[,2]*exp(-lambda*data[,2]))/(exp(-lambda*data[,1])-exp(-lambda*data[,2])))
}
M_step <- function(lambda){
  n*lambda/(E_step(lambda)+n)
}
lambda <- 1 #初始值
iter <- 1
epsilon<-1e-6
criteria<-1e10
while (criteria>epsilon && iter < 1000) {
  lambda_new <- M_step(lambda)
  iter <- iter+1
  criteria<-abs(lambda_new-lambda)
  lambda <- lambda_new
}
print(paste("E-M算法得到的Lambda的MLE解为:", round(lambda_new,3)))


## -----------------------------------------------------------------------------
solve.game <- function(A) {
  min.A <- min(A)
  A <- A - min.A
  max.A <- max(A)
  A <- A/max(A)
  m <- nrow(A)
  n <- ncol(A)
  it <- n^3
  a <- c(rep(0, m), 1)
  A1 <- -cbind(t(A), rep(-1, n))
  b1 <- rep(0, n)
  A3 <- t(as.matrix(c(rep(1, m), 0)))
  b3 <- 1
  sx <- simplex(a = a, A1 = A1, b1 = b1, A3 = A3, b3 = b3,
                maxi = TRUE, n.iter = it)
  a <- c(rep(0, n), 1)
  A1 <- cbind(A, rep(-1, m))
  b1 <- rep(0, m)
  A3 <- t(as.matrix(c(rep(1, n), 0)))
  b3 <- 1
  sy <- simplex(a = a, A1 = A1, b1 = b1, A3 = A3, b3 = b3,
                maxi = FALSE, n.iter = it)
  soln <- list(A = A * max.A + min.A, x = sx$soln[1:m],
               y = sy$soln[1:n], v = sx$soln[m + 1] * max.A + min.A)
  soln
}

# 示例输入数据
A <- matrix(c(0, -2, -2, 3, 0, 0, 4, 0, 0,
              2, 0, 0, 0, -3, -3, 4, 0, 0,
              2, 0, 0, 3, 0, 0, 0, -4, -4,
              -3, 0, -3, 0, 4, 0, 0, 5, 0, 0,
              3, 0, -4, 0, -4, 0, 5, 0, 0, 3,
              0, 0, 4, 0, -5, 0, -5, -4, -4, 0,
              0, 0, 5, 0, 0, 6, 0, 0, 4, -5, -5,
              0, 0, 0, 6, 0, 0, 4, 0, 0, 5, -6, -6, 0), 9, 9)


library(boot) 
B <- A + 2
s <- solve.game(B)
s$v  # 输出解决方案的值
round(cbind(s$x, s$y), 7)
round(s$x * 61, 7)

## -----------------------------------------------------------------------------
# 可以
a<-data.frame(x=0.1,y=0.2)
a[FALSE,]
a[,FALSE]
a[FALSE,FALSE]

## -----------------------------------------------------------------------------
scale01 <- function(x) {
rng <- range(x, na.rm = TRUE)
(x - rng[1]) / (rng[2] - rng[1])
}

## -----------------------------------------------------------------------------
# 先构造两个dataframe 
# 数值数据框
x<-1:5
y<-6:10
z<-c("A","B","C","D","E")
df<-data.frame(x=x,y=y)
# 混合数据框
mix.df<-data.frame(x=x,y=y,z=z)

# 对于数据框的所有列应用函数
apply(df, 2, scale01)

# 对于数据框的数值列应用函数
data.frame(lapply(mix.df, function(x) if (is.numeric(x)) scale01(x) else x))

## -----------------------------------------------------------------------------
# 仍然利用Exercises 2 (page 204, Advanced R)构造的数据框
vapply(df, sd, numeric(1))

## -----------------------------------------------------------------------------
# 仍然利用Exercises 2 (page 204, Advanced R)构造的数据框
vapply(mix.df[vapply(mix.df, is.numeric, logical(1))],sd, numeric(1))

## ----eval=FALSE---------------------------------------------------------------
#  gibbsR <- function(N,a,b,n) {
#    mat <- matrix(nrow = N, ncol = 2)
#    x <- y <- 0
#    for (i in 1:N) {
#      x <- rbinom(1, n, y)
#      y <- rbeta(1, x+a, n-x+b)
#      mat[i, ] <- c(x, y)
#    }
#    mat
#  }

## ----eval=FALSE---------------------------------------------------------------
#  NumericMatrix gibbsC(int N,  int a, int b, int n) {
#    NumericMatrix mat(N, 2);
#    double x = 0, y = 0;
#    for(int i = 0; i < N; i++) {
#      x = rbinom(1, n, y)[0];
#      y = rbeta(1, x+a, n-x+b)[0];
#      mat(i, 0) = x;
#      mat(i, 1) = y;
#    }
#    return(mat);
#  }

## -----------------------------------------------------------------------------
library(Rcpp)
library(microbenchmark)
library(SA23204168)
ts <- microbenchmark(gibbsR=gibbsR(100,2,3,10),gibbsC=gibbsC(100,2,3,10))
#summary(ts)[,c(1,3,5,6)]

