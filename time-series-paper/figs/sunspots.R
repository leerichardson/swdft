devtools::load_all(pkg = "/home/lee/Dropbox/swdft/r/swdft")

data("sunspot.year")

sunspots <- sunspot.year
sunspot_years <- 1700:1988
window_size_sunspots <- 2^6

a_sunspots <- swdft::swdft(x = sunspots,
                           n = window_size_sunspots,
                           normalize = (1 / sqrt(window_size_sunspots)))

data("lynx")

lynx <- lynx
lynx_years <- 1821:1934
window_size_lynx <- 2^5
a_lynx <- swdft::swdft(x = lynx, n = window_size_lynx)

png("/home/lee/Dropbox/thesis/writing/swft_timeseries_paper/doc/images/sunspots_lynx.png")

  layout(mat = matrix(data = 1:4, nrow = 2, ncol = 2,  byrow = FALSE))
  par(mar = c(4,4,4,4))
  plot(lynx, main = "Canadian Lynx Time-Series", ylab = "Annual Lynx Trappings", xlab="")
  abline(v=lynx_years[1], lty=2)
  abline(v=lynx_years[1] + window_size_lynx, lty=2)

  custom_x_axis_lynx <- (lynx_years[1] + window_size_lynx - 1):lynx_years[length(lynx_years)]

  plot_swdft(a = a_lynx, only_unique=TRUE, take_log=FALSE, use_fields = FALSE, pad_array=TRUE,
           title = "SWDFT of Canadian Lynx", xlab = "",
           custom_xaxis=lynx_years)

  custom_x_axis_sunspots <- (sunspot_years[1] + window_size_sunspots - 1):sunspot_years[length(sunspot_years)]
  plot(sunspots, main = "Yearly Sunspots Time-Series", ylab = "Number of Sunspots", xlab="")
  abline(v=sunspot_years[1], lty=2)
  abline(v=sunspot_years[1] + window_size_sunspots, lty=2)
  plot_swdft(a=a_sunspots, take_log=FALSE, only_unique=TRUE, use_fields=FALSE, pad_array=TRUE,
           title = "SWDFT of Sunspots", xlab="", custom_xaxis=sunspot_years)

dev.off()
