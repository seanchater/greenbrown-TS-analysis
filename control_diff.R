library(sf)
library(ggplot2)
library(reshape2)


# data_path <- 'farm5_merged_2005-2022'
data_path <- 'farm5_actual_ndvi_2005-2022'
raster_files = list.files(path=data_path, pattern='.tif$', full.names = TRUE)
raster_files

ndvi_stack <- stack(raster_files)
crs(ndvi_stack)


# fill gaps 
impute.loess <- function(y, x.length = NULL, s = 0.80, 
						 smooth.data = FALSE, ...) {
	if(is.null(x.length)) { x.length = length(y) }
	options(warn=-1)
	x <- 1:x.length
	p <- loess(y ~ x, span = s, data.frame(x=x, y=y))
	if(smooth.data == TRUE) {
		y <- predict(p, x)
	} else {
		na.idx <- which( is.na(y) )
		if( length(na.idx) > 1 ) {
			y[na.idx] <- predict(p, data.frame(x=na.idx))
		}
	}   
	return(y)
}

ndvi.new <- ndvi_stack
ndvi.new[] <- NA
for( rl in 1:nrow(ndvi_stack) ) { 
	v <- getValues(ndvi_stack, rl, 1)
	ndvi.new[rl,] <- as.matrix(t(apply(v, MARGIN=1, FUN=impute.loess)))
}

# set the points where you want to extract values from the raster stack
df <- data.frame(point = c("5.1.2", "5.2.7"), lat = c(-30.3021059112123, -30.3031838895532), long = c(21.7661344184518, 21.7672123967927))
df_sf <- st_as_sf(x = df,
				  coords = c("long", "lat"),
				  crs = 4326)

# extract the point values
values <- raster::extract(ndvi.new, df_sf)

complete <- cbind(df_sf, values)

final <- as.data.frame(t(complete[,-1])) # transpose to have the observations for each point as rows instead of cols
row.names(final) <- NULL # remove the filenames as rows from the df
final <- cbind(rownames(final), data.frame(final, row.names=NULL))
colnames(final) <- c("obs", "cleared512", "control527")
final <- final[-217,] # remove the row containing coords from the df

final <- data.frame(lapply(final,as.numeric)) # convert all cols to numeric format and return as a dataframe
final$diff <- (final$control527 - final$cleared512)

final2 <- reshape::melt(final, id.var="obs") # reshape to long format
ggplot(data=final2, mapping=aes(x=obs, y=value, color=variable)) + 
	geom_line() +
	scale_color_manual(values=c('red', 'green', "blue")) +
	geom_point()

boxplot(final$cleared512)
boxplot(final$control527)

# gap filling



# SG filtering
library(signal) # for sgolayfilt function for smoothing and removing outliers
sg <- sgolayfilt(final$cleared512) # default: sgolayfilt(x, p = 3, n = p + 3 - p%%2, m = 0, ts = 1)
sg_n7 <- sgolayfilt(final$cleared512, n=9)

plot(final$obs, final$cleared512, type="l", main="unfiltered signal", xaxt='n', xlab='month', ylab='NDVI')
axis(side=1, at = seq(0,216,12))
plot(sg, type="l", main="sg default", xaxt='n', xlab='month', ylab='NDVI')
axis(side=2, at=c(-0.15, -0.1, -0.05, 0, 0.05, 0.1, 0.15))
axis(side=1, at = seq(0,216,12))
plot(sg_n7, type="l", main="sg n=7", xaxt='n', xlab='month', ylab='NDVI')
axis(side=1, at = seq(0,216,12))




ggplot(data=final, mapping=aes(x=obs, y=diff)) + geom_line() + geom_point()

#temporal filling
# remove outliers
# gee evi - savitzky-golay EVI filter

