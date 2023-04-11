library(terra)


## NOTE 1 is water, 17 is barren, 16 snow-ice. It is shifted. snow-ice is counted as land here
mcd12_files <- list.files("G:/MCD12C1/2020/reprocessed/percent", full.names = TRUE)

mcd12 <- rast(mcd12_files)

test_map <- sum(mcd12)

plot(test_map)

plot(mcd12[[1]])

test <- mcd12[[1]]

m <- c(0, 100, 0)
m <- matrix(m, ncol=3, byrow=TRUE)
mcd12[[1]] <- classify(mcd12[[1]], m, include.lowest=TRUE)


land_cover_percent <- sum(mcd12)

writeCDF(land_cover_percent, "G:/MCD12C1/2020/reprocessed/land_cover_percent.nc", overwrite = TRUE)


mcd12[[17]]
