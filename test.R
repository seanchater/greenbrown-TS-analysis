library(raster)
library(reshape2)
library(rgdal)
library(scales)
library(dplyr)


r1 <- raster(ncols=5, nrows=10)
values(r1) <- 1#:ncell(r1)
plot(r1)
dev.print(png, file = "outputs/myplot.png", width = 1024, height = 768)
# png("test/test.png")

# r2 <- raster(ncols=5, nrows=10)
# r2[] <- 2#NA
# plot(r2)
# 
# r3 <- raster(ncols=5, nrows=10)
# values(r3) <- 3#(1:ncell(r2))*3
# plot(r3)
# 
# 
# stack <- stack(r1,r2,r3)
# stack
# stack[]
# plot(stack)
# 
# for (i in 1:nlayers(stack)) {
# 	stack[[i]] <- stack[[i]] + 1
# }
# 
# # replace NA values by row mean in raster stack
# for( rl in 1:nrow(stack) ) { 
# 	v <- getValues(stack, rl, 1)
# 	v[is.na(v)] <- rowMeans(v, na.rm=TRUE)
# 	stack[rl,] <- v
# }
# stack
# stack[]
# plot(stack)
