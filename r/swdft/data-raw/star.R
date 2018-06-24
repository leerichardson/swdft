# From the Following Website:
# https://www.york.ac.uk/depts/maths/data/ts/welcome.htm
star_raw <- read.table("https://www.york.ac.uk/depts/maths/data/ts/ts26.dat")

n <- nrow(star_raw)
m <- ncol(star_raw)
star <- vector(mode = "numeric", length = n * m)

iter <- 0
for (i in 1:n) {
  for (j in 1:m) {
    iter <- iter + 1
    cat("i: ", i, " j: ", j, " \n")
    star[iter] <- star_raw[i, j]
  }
}

devtools::use_data(star)
