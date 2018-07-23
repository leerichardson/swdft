devtools::load_all(".")

lfp <- read.table(file = "http://www.stat.cmu.edu/~kass/KEB/data/lfp-ryan.dat")
lfp$t <- 0:(nrow(lfp) - 1)
names(lfp) <- c("val", "time")

inds <- 1:5000
window_size <- 2^9

lfp <- lfp[inds, ]

a_lfp <- swdft::swdft(x=lfp$val, n=window_size)

par(mfrow=c(3,1))
plot(lfp$time, lfp$val, type="l")
abline(v=0)
abline(v=0 + window_size)
plot_mvts(a=a_lfp, legend=FALSE, only_unique=TRUE, take_log=TRUE)
plot_swdft(a=a_lfp, take_log=TRUE, only_unique=TRUE)

amod_lfp <- log( Mod(a_lfp)^2 )
plot(amod_lfp[2, ], ylim = c(0, 20))
lines(amod_lfp[3, ])
lines(amod_lfp[4, ])
lines(amod_lfp[5, ])
lines(amod_lfp[6, ])
