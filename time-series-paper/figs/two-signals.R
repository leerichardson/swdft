devtools::load_all("/home/lee/Dropbox/swdft/r/swdft")

# R Local Periodic Signals ---
png("/home/lee/Dropbox/thesis/writing/swft_timeseries_paper/doc/images/two_signals.png", width = 720)

    N <- 2^8
    n <- 2^4
    m <- floor(n / 2)
    A_vec <- c(1, 1, 1, 1)
    Fr_vec <- c(18.3, 48, 18.3, 48)
    phi_vec <- c(0, 0, 0, 0)
    S_vec <- c(30, 150, 70, 100)
    L_vec <- c(90, 50, 90, 50)

    x1 <- swdft::local_signal(N=N, A=A_vec[1], Fr=Fr_vec[1], phase=phi_vec[1], S=S_vec[1], L=L_vec[1])
    x2 <- swdft::local_signal(N=N, A=A_vec[2], Fr=Fr_vec[2], phase=phi_vec[2], S=S_vec[2], L=L_vec[2])
    y1 <- x1 + x2
    b1 <- swdft::swdft(x=x1, n=n, normalize=(1 / sqrt(n)))
    b2 <- swdft::swdft(x=x2, n=n, normalize=(1 / sqrt(n)))
    a1 <- swdft::swdft(x=y1, n=n, normalize=(1 / sqrt(n)))

    x3 <- swdft::local_signal(N=N, A=A_vec[3], Fr=Fr_vec[3], phase=phi_vec[3], S=S_vec[3], L=L_vec[3])
    x4 <- swdft::local_signal(N=N, A=A_vec[4], Fr=Fr_vec[4], phase=phi_vec[4], S=S_vec[4], L=L_vec[4])
    y2 <- x3 + x4
    b3 <- swdft::swdft(x=x3, n=n, normalize=(1 / sqrt(n)))
    b4 <- swdft::swdft(x=x4, n=n, normalize=(1 / sqrt(n)))
    a2 <- swdft::swdft(x=y2, n=n, normalize=(1 / sqrt(n)))

    layout(mat = matrix(data = 1:8, nrow = 4, ncol = 2,  byrow = FALSE))

    plot(y1, pch = 19, cex = 1, ylim = c(-2, 2),
       main = "Two Local Periodic Signals: No Overlap", xlab = "", ylab = "")
    lines(x1, col = "grey", lwd = 1.5)
    lines(x2, col = "black", lwd = 1.5)

    swdft::plot_swdft(a = b1, title = "SWDFT of Signal 1", zlim= c(0, max(Mod(a1)^2)), use_fields = FALSE, ylab = "")
    swdft::plot_swdft(a = b2, title = "SWDFT of Signal 2", zlim= c(0, max(Mod(a1)^2)), use_fields = FALSE, ylab = "")
    swdft::plot_swdft(a = a1, title = "SWDFT of Signal 1 + Signal 2", zlim= c(0, max(Mod(a1)^2)), use_fields = FALSE, ylab = "")

    plot(y2, pch = 19, cex = 1, ylim = c(-2, 2),
       main = "Two Local Periodic Signals: Overlap", xlab = "", ylab = "")
    lines(x3, col = "grey", lwd = 1.5)
    lines(x4, col = "black", lwd = 1.5)

    swdft::plot_swdft(a = b3, title = "SWDFT of Signal 1", zlim= c(0, max(Mod(a1)^2)), use_fields = FALSE, ylab = "")
    swdft::plot_swdft(a = b4, title = "SWDFT of Signal 2", zlim= c(0, max(Mod(a1)^2)), use_fields = FALSE, ylab = "")
    swdft::plot_swdft(a = a2, title = "SWDFT of Signal 1 + Signal 2", zlim= c(0, max(Mod(a1)^2)), use_fields = FALSE, ylab = "")

dev.off()
