N <- 10
R <- 3
n <- 4

x <- 0:4

com1 <- choose(R, x)
com2 <- choose(N - R, n - x)
com3 <- choose(N,n)

com1*com2/com3

dhyper(x, R, N - R, n)



