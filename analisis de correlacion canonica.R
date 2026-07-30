library(CCA) # CCA: An R Package to Extend Canonical Correlation Analysis
data(mtcars)

X = mtcars[,c(2,3,5,10,8,11)] # Design features variables
Y = mtcars[,c(1,4,6,7,9)] # Driver features variables

ccs <- cc(X, Y)

ufvcc1 <- ccs$scores$xscores[ , 1]
dfvcc1 <- ccs$scores$yscores[ , 1]

cat ("Correlation of the first canonical pair: ",cor(ufvcc1, dfvcc1))

ufvcc2 <- ccs$scores$xscores[ , 2]
dfvcc2 <- ccs$scores$yscores[ , 2]

cat ("\nCorrelation of the second canonical pair: ",cor(ufvcc2, dfvcc2))

sdr <- sort(ufvcc1)
sdr <- sdr[c(1, length(sdr) - 1)] # first and next-to-last
ext <- match(sdr, ufvcc1)

plot( ufvcc1, dfvcc1, cex.lab = 1.5,
      xlab = "User features canonical scores", 
      ylab = "Driver features canonical scores",
      pch = 16, cex = 1, col = "red", 
      xlim = sdr * c(1.3, 2), 
      frame.plot=FALSE)

text(ufvcc1[ext], dfvcc1[ext], 
     labels = rownames(mtcars)[ext],
     pos = c(1, 2), 
     cex = 1, 
     col = "blue")


library(CCP)
plt.cc(ccs, d1 = 1, d2 = 2, type = "v", var.label = TRUE)



p.perm(X, Y, nboot = 999, rhostart = 3, type = "Wilks")

#Test Ho : ρ2* = ρ3* = · · · = rhop = 0 vs Ha: at least one is not zero. Interpret results
n = 5
for(i in 2:n){
  cat("rhostart = ", i,"\n")
  p.perm(X, Y, nboot = 999, rhostart = i, type = "Wilks")
  cat("\n")
}