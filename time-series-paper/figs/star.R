devtools::load_all(pkg = "/home/lee/Dropbox/swdft/r/swdft")
data(star)

window_size <- 2^7
a_star <- swdft::swdft(x = star, n = window_size, normalize = (1 / sqrt(window_size)))
n <- nrow(a_star)
P <- ncol(a_star)

png("/home/lee/Dropbox/thesis/writing/swft_timeseries_paper/doc/images/stars.png")

  par(mfrow = c(3, 1))

  plot(0:(N - 1), star, type="l", xlab="Time", ylab="Brightness",
       main="Time-Series of Variable Star Data")

  pad_array <- matrix(data = NA_complex_, nrow = n, ncol = window_size - 1)
  a_pad <- cbind(pad_array, a_star)
  amod_star <- Mod(a_pad)^2
  amod_ylim <- c(0, max(amod_star[2:nrow(amod_star), ], na.rm = TRUE))

  plot(0:(N - 1), amod_star[5, ], type = "l", ylim = amod_ylim,
       lty=2, col="black", xlab="Time", ylab="Energy") # 32 Days
  lines(0:(N - 1), amod_star[6, ], col = "grey", lty=1) # 25.6 Days

  legend("topleft", c("25.6 Days", "32 Days"), col = c("black", "grey"),
         lty=c(2, 1), cex = .8)

  plot_swdft(a=a_star, take_log=TRUE, only_unique=TRUE, use_fields=FALSE,
             title= "SWDFT of Variable Star Data")

dev.off()


# From chapter three of bloomfields
N <- length(star)
star_dft <- fft(star) * (1 / N)
star_mod <- Mod(star_dft)^2
freqs <- (0:(N - 1)) / N
M <- floor(N / 2)
star_df <- data.frame(freqs, star_mod)
star_inds <- 2:M

plot(star_df$freqs[star_inds], log(star_df$star_mod[star_inds]), pch=19, cex=.5)

