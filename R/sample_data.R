set.seed(2020)
n = 100
p = 10
k = 1:p
beta = c(1, -1, rep(x = 0,
                    times = p-2))
x = matrix(data = rnorm(n*p),
           ncol = p)
colnames(x) = paste0("X", 1:p)
y = rbinom(n = n,
           size = 1,
           prob = expit(x %*% beta))

sampleData = data.frame(x, y)
