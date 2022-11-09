library(raster)  
library(reshape2)
library(rgdal)
library(scales)
library(dplyr)
library(greenbrown)
library(ggplot2)
library(signal) # for sgolayfilt function for gap filling
library(sf)

####################################################################################
###########      Load in Data      ###########
####################################################################################

load_data <- function(data_path, control_path) {
	data_files = list.files(path=data_path, pattern='.tif$', full.names=TRUE)
	control_files = list.files(path=control_path, pattern='.tif$', full.names=TRUE)
	data_files
	
	# create a stack from the file list
	data_stack <<- stack(data_files)
	control_stack <<- stack(control_files)
	crs(data_stack)
	
	# Add 1 to all values in each raster stack
	# This is to make all NDVI values positive to simplify calculations further on such as with differencing section
	for (i in 1:nlayers(data_stack)) {
		data_stack[[i]] <<- data_stack[[i]] + 1
	}
	for (i in 1:nlayers(control_stack)) {
		control_stack[[i]] <<- control_stack[[i]] + 1
	}
}

####################################################################################
###########      Other thingy      ###########
####################################################################################

gb_gapfill_filter <- function(x) {
	v=as.vector(x)
	ndvi_ts = ts(v, start=c(2000,1), end=c(2022,12), frequency=12)
	x=TSGFlinear(ndvi_ts, interpolate = TRUE)
}


original_filter <- function() {
	# https://matinbrandt.wordpress.com/2014/12/02/smoothingfiltering-a-ndvi-time-series-using-a-savitzky-golay-filter-and-r/
	# check for gaps and fill - SG filtering https://www.sciencedirect.com/science/article/pii/S0034425716304692#s0120
	# data_stack_sg <- sgolayfilt(data_gapfill[[2:nlayers(data_gapfill)]])
	
	
	data_filtered <<- data_stack
	data_filtered[] <<- NA
	for( rl in 1:nrow(data_stack) ) { 
		v <- getValues(data_stack, rl, 1)
		data_filtered[rl,] <<- as.matrix(t(apply(v, MARGIN=1, FUN=gb_gapfill_filter)))
	}
	
	control_filtered <<- control_stack
	control_filtered[] <<- NA
	for( rl in 1:nrow(control_stack) ) { 
		v <- getValues(control_stack, rl, 1)
		control_filtered[rl,] <<- as.matrix(t(apply(v, MARGIN=1, FUN=gb_gapfill_filter)))
	}
}
####################################################################################
###########      Gap filling      ###########
####################################################################################

# imput NA values
#### Loess impute function - https://gis.stackexchange.com/questions/279354/ndvi-time-series-with-missing-values
# y            Vector with NA values or to be smoothed
# x.length     length of resulting vector
# s            Smoothing parameter (see loess span argument)
# smooth.data  (FALSE/TRUE) smooth all of the data
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

fill_gaps <- function() {
	data_gapfill <<- data_stack
	data_gapfill[] <<- NA
	for( rl in 1:nrow(data_stack) ) { 
		v <- getValues(data_stack, rl, 1)
		data_gapfill[rl,] <<- as.matrix(t(apply(v, MARGIN=1, FUN=impute.loess)))
	}
	
	control_gapfill <<- control_stack
	control_gapfill[] <<- NA
	for( rl in 1:nrow(control_stack) ) { 
		v <- getValues(control_stack, rl, 1)
		control_gapfill[rl,] <<- as.matrix(t(apply(v, MARGIN=1, FUN=impute.loess)))
	}
	
	############ ISSUE WITH THIS AS THEN GRAPHS PLOTTING A SINGLE PIXEL ARE NOT GOING TO BE ACCURATE :
	############ NEED TO INFILL WITH MEAN OR SOMETHING SO THAT TIME SERIES IS STILL 276 (2000-2022) MONTHS LONG
	# for (i in 1:nlayers(data_gapfill)) {
	# 	s <- summary(data_gapfill[[i]])
	# 	n <- s[[6]]
	# 	if (n > 0) {
	# 		print(paste("layer " , i , " has " , n , " NA's ... DROPPING LAYER"))
	# 		data_gapfill <- dropLayer(data_gapfill, i)
	# 	}
	# }
	
	# gap fill function may have trouble with imputing beginning and end of time series so
	# just replace the leftover NA values with the mean of each pixels TS
	for( rl in 1:nrow(data_gapfill) ) { 
		v <- getValues(data_gapfill, rl, 1)
		v[is.na(v)] <- rowMeans(v, na.rm=TRUE)
		data_gapfill[rl,] <<- v
	}
	
	for( rl in 1:nrow(control_gapfill) ) { 
		v <- getValues(control_gapfill, rl, 1)
		v[is.na(v)] <- rowMeans(v, na.rm=TRUE)
		control_gapfill[rl,] <<- v
	}
}

####################################################################################
###########     Filtering      ###########
####################################################################################

sg_filter <- function(x) {
	v=as.vector(x)
	ndvi_ts = ts(v, start=c(2000,1), end=c(2022,12), frequency=12)
	x=sgolayfilt(ndvi_ts, p=1, n=3, ts=30)
}

filter <- function() {
	# https://matinbrandt.wordpress.com/2014/12/02/smoothingfiltering-a-ndvi-time-series-using-a-savitzky-golay-filter-and-r/
	# check for gaps and fill - SG filtering https://www.sciencedirect.com/science/article/pii/S0034425716304692#s0120
	# data_stack_sg <- sgolayfilt(data_gapfill[[2:nlayers(data_gapfill)]])
	
	
	data_filtered <<- data_gapfill
	data_filtered[] <<- NA
	for( rl in 1:nrow(data_gapfill) ) { 
		v <- getValues(data_gapfill, rl, 1)
		data_filtered[rl,] <<- as.matrix(t(apply(v, MARGIN=1, FUN=sg_filter)))
	}
	
	control_filtered <<- control_gapfill
	control_filtered[] <<- NA
	for( rl in 1:nrow(control_gapfill) ) { 
		v <- getValues(control_gapfill, rl, 1)
		control_filtered[rl,] <<- as.matrix(t(apply(v, MARGIN=1, FUN=sg_filter)))
	}
}

# TODO: rename layers in rasterbricks to appropriate names

####################################################################################
###########     Control difference     ###########
####################################################################################

difference <- function(output_dir, farm) {
	# set the points where you want to extract values from the raster stack
	df <- NULL
	points <- NULL
	control <- NULL
	if (farm == 5) {
		# farm 5
		control <- 5
		points <- c("5.2.1", "5.2.2", "5.2.3", "5.2.4", "5.2.5")
		df <- data.frame(point = points,
						 lat = c( -30.3039923733089, -30.3037228787237, -30.3037228787237, -30.3034533841385, -30.3034533841385),
						 lon = c(21.7666734076222, 21.7666734076222, 21.7669429022075, 21.7666734076222, 21.7669429022075))
	} else if (farm == 19) {
		# farm 19
		control <- 17
		points <- c("17.1.1", "17.1.2", "17.1.3", "17.1.4", "17.1.5", "17.1.6")
		df <- data.frame(point = points,
						 lat = c(-30.89768894, -30.89768894, -30.89795844, -30.89795844, -30.89795844, -30.89795844),
						 lon = c(22.08683297, 22.08656348, 22.08656348, 22.08683297, 22.08710247, 22.08737196))
	} else if (farm == 11) {
		control <- 11
		points <- c("11.4.1", "11.4.2", "11.4.3", "11.4.4", "11.4.5", "11.4.6", "11.4.7", "11.4.8")
		df <- data.frame(point = points,
						 lat = c(-30.6909866, -30.6907171, -30.6907171, -30.69044761, -30.69044761, -30.69017811, -30.69017811, -30.69017811),
						 lon = c(22.0555716, 22.0555716, 22.0558411, 22.0558411, 22.0555716, 22.0555716, 22.0558411, 22.05611059))
	} else if (farm == 14) {
		control <- 18
		points <- c("18.1.1", "18.1.2", "18.1.3", "18.1.4", "18.1.9", "18.1.10", "18.1.11", "18.1.12", "18.1.13")
		df <- data.frame(point = points,
						 lat = c(-30.6610727, -30.6610727, -30.6608032, -30.66134219, -30.66053371, -30.66026421, -30.65999472, -30.65999472, -30.66026421),
						 lon = c(22.22643117, 22.22616168, 22.22616168, 22.22643117, 22.22643117, 22.22643117, 22.22643117, 22.22616168, 22.22616168))
	} else if (farm == 1) {
		control <- 12
		points <- c("12.1.1", "12.1.2", "12.1.3", "12.1.4", "12.1.5", "12.1.6", "12.1.7", "12.1.8", "12.1.9")
		df <- data.frame(point = points,
						 lat = c(-30.47835537,-30.47862486,-30.47862486,-30.47862486,-30.47862486,
						 		-30.47862486,-30.47862486,-30.47862486,-30.47835537),
						 lon = c(22.32048478,22.32075427,22.32102377,22.32129326,22.32156276,
						 		22.32183225,22.32210175,22.32237124,22.32183225))
	}
	
	df_sf <- st_as_sf(x = df,
					  coords = c("lon", "lat"),
					  crs = 4326)
	
	# extract the point values
	values <- raster::extract(control_filtered, df_sf)
	
	# TODO: combine with coordinates: not sure if this is neccessary?
	# complete <- cbind(df_sf, values)
	
	final <- as.data.frame(t(values)) # transpose to have the observations for each point as rows instead of cols
	# need to average the pixel values of (minimum) 5 different control points per layer per month 
	row_av <- apply(final, 1, mean)
	row_sd <- apply(final, 1, sd)
	final <- cbind(final, mean=row_av, sd=row_sd)
	row.names(final) <- NULL # remove the filenames as rows from the df
	final <- cbind(rownames(final), data.frame(final, row.names=NULL)) # add the row numbers as it's own col
	colnames(final) <- c("months", points, "mean", "sd")
	final <- data.frame(lapply(final,as.numeric)) # convert all cols to numeric format and return as a dataframe
	# final$diff <- (final$cleared512 - final$control521)
	var <- cbind(obs=final$obs, mean=final$mean, sd = final$sd) # create list with variables for plotting
	
	# TODO: plot mean and SD for averaged control sites
	print(ggplot(final, aes(x=months, y=mean, colour="red")) + 
		geom_line()+
		geom_linerange(aes(ymin=mean-sd, ymax=mean+sd, colour="grey"), width=.2))
	# ggsave("imag.jpg", plot = mean_sd) 
	
	# dev.print(png, filename = paste(output_dir, "mean_sd2.png", sep=""),  height = 400)
	
	# TODO: for each row in the data_filtered rasterbrick (ie: each TS of a pixel) => subtract the mean of the control points from that vector
	data_diff <<- data_filtered
	data_diff[] <<- NA
	for( rl in 1:nrow(data_filtered) ) { 
		v <- getValues(data_filtered, rl, 1)
		data_diff[rl,] <<- v - final$mean
	}
	
	obs <- names(final[1])
	first_point <- names(final[2])
	plot(final[,obs], final[,first_point], type="l", main=c("cleared", names(final[2]), " signal"), xaxt='n', xlab='Month', ylab='NDVI')
	axis(side=1, at = seq(0,276,12))
	# dev.print(png, file = c(output_dir, "first_contol_point.png"), width = 1024, height = 768)
	
	# ggplot(data=final, mapping=aes(x=obs, y=diff)) + geom_line() + geom_point()
	
	final2 <- reshape::melt(final, id.var="months") # reshape to long format
	values <- NULL
	if (farm == 5) {
		values <- c('red', 'green', "blue", "black", "yellow", "magenta","pink")
	} else if (farm == 19) {
		values <- c('red', 'green', "blue", "black", "yellow", "orange", "grey", "pink")
	} else if (farm == 11) {
		values <- c('red', 'green', "blue", "black", "yellow", "orange", "grey", "pink", "turquoise", "purple")
	} else if (farm == 14) {
		values <- c('red', 'green', "blue", "black", "yellow", "orange", "grey", "pink", "turquoise", "purple", "plum")
	} else if (farm == 1) {
		values <- c('red', 'green', "blue", "black", "yellow", "orange", "grey", "pink", "turquoise", "purple", "plum")
	}
	print(ggplot(data=final2, mapping=aes(x=months, y=value, color=variable)) + 
		geom_line() +
		scale_color_manual(values=values) +
			scale_x_continuous(name="Year (12 month intervals from 2000-2022)", limits=c(0, 276), breaks=seq(0,276,12)) +
			ylab("NDVI Value (Rescaled 0-2") + 
			ggtitle(label=sprintf("Time Series of Control Points in Farm %s with Mean and Standard Deviation", control))) #+
	# geom_point()
	# dev.print(png, file = c(output_dir, "all_control_points"), width = 1024, height = 768)
	
}

####################################################################################
###########     Plot the different stages    ###########
####################################################################################

visualise_stages <- function(output_dir, farm) {
	
	# plot original raster layers
	# par(oma = c(2, 2, 2, 2))
	# plot(data_stack)
	# dev.print(png, file = c(output_dir, "raw_data.png"), width = 1024, height = 768)
	
	# plot gap filled raster layers
	# par(oma = c(2, 2, 2, 2))
	# plot(data_gapfill)
	# dev.print(png, file = c(output_dir, "gapfill_data.png"), width = 1024, height = 768)
	
	# plot filtered raster layers
	# par(oma = c(2, 2, 2, 2))
	# plot(data_filtered)
	# dev.print(png, file = c(output_dir, "filtered_data.png"), width = 1024, height = 768)
	
	# plot raster layers minus the mean of control points
	# par(oma = c(2, 2, 2, 2))
	# plot(data_diff)
	# dev.print(png, file = c(output_dir, "diff_data.png"), width = 1024, height = 768)
	
	
	# plot some points from the study area
	df <- NULL
	points <- NULL
	if (farm == 5) {
		# farm 5
		points <- c("2", "3", "5", "6", "8",
					"9", "10", '11', "12", "13",
					"14", "15")
		df <- data.frame(point = points,
						 lat = c(-30.30210591,-30.30210591,-30.30237541,-30.3026449,-30.30237541,
						 		-30.3026449,-30.30210591,-30.30237541,-30.3026449,-30.30210591,
						 		-30.30237541,-30.3026449),
						 lon = c(21.76613442,21.76640391,21.76613442,21.76613442,21.76640391,
						 		21.76640391,21.76667341,21.76667341,21.76667341,
						 		21.7669429,21.7669429,21.7669429))
	} else if (farm == 19) {
		# farm 19
		points <- c("1", "2", "3", "4", "5",
					"6", "7", "8", "9", "10",
					"11", "12")
		df <- data.frame(point = points,
						 lat = c(-30.85079689,-30.85079689,-30.85079689,-30.85079689,-30.85106638,
						 		-30.85106638,-30.85133588,-30.85133588,-30.85133588,-30.85106638,
						 		-30.85079689,-30.85079689),
						 lon = c(22.09922973,22.09949922,22.09976871,22.10003821,22.09922973,
						 		22.09896023,22.09922973,22.09896023,22.09869074,22.09869074,
						 		22.09869074,22.09896023))
	} else if (farm == 11) {
		points <- c("11", "12", "13", "14", "15", "16", "17", "18", "19")
		df <- data.frame(point = points,
						 lat = c(-30.68775266,-30.68775266,-30.68775266,-30.68802216,-30.68802216,
						 		-30.68802216,-30.68829165,-30.68829165,-30.68829165),
						 lon = c(22.02673568,22.02646619,22.02619669,22.02619669,22.02646619,
						 		22.02673568,22.02673568,22.02646619,22.02619669))
	} else if (farm == 14) {
		points <- c("1", "2", "3", "4", "5", "6")
		df <- data.frame(point = points,
						 lat = c(-30.69880194,-30.69853245,-30.69853245,-30.69826295,-30.69826295,-30.69799346),
						 lon = c(22.22265825,22.22292774,22.22265825,22.22265825,22.22292774,22.22265825))
	} else if (farm == 1) {
		points <- c("1", "3", "4", "5", "6", "7", "11", "12")
		df <- data.frame(point = points,
						 lat = c(-30.60178389,-30.601245,-30.601514,-30.602053,-30.60286187,
						 		-30.60205338,-30.60259237,-30.60232288),
						 lon = c(22.47571366,22.475714,22.475714,22.475983,22.47544417,
						 		22.47598316,22.47571366,22.47598316))
	}
	
	df_sf <- st_as_sf(x = df,
					  coords = c("lon", "lat"),
					  crs = 4326)
	
	# extract the point values
	values <- raster::extract(data_diff, df_sf)
	
	# TODO: combine with coordinates: not sure if this is neccessary?
	# complete <- cbind(df_sf, values)
	
	final <- as.data.frame(t(values)) # transpose to have the observations for each point as rows instead of cols
	row_av <- apply(final, 1, mean)
	# row_sd <- apply(final, 1, sd)
	final <- cbind(final, mean=row_av)#, sd=row_sd)
	row.names(final) <- NULL # remove the filenames as rows from the df
	final <- cbind(rownames(final), data.frame(final, row.names=NULL)) # add the row numbers as it's own col
	colnames(final) <- c("months", points, "SD", "mean")
	final <- data.frame(lapply(final,as.numeric)) # convert all cols to numeric format and return as a dataframe
	# final$diff <- (final$cleared512 - final$control521)
	var <- cbind(obs=final$obs, mean=final$mean, sd = final$sd) # create list with variables for plotting
	
	# TODO: plot mean and SD for averaged control sites
	# print(ggplot(final, aes(x=months, y=mean, colour="red")) + 
		  	# geom_line()+
		  	# geom_linerange(aes(ymin=mean-sd, ymax=mean+sd, colour="grey"), width=.2))

	# obs <- names(final[1])
	# first_point <- names(final[2])
	# plot(final$months, final$X5.1.8, type="l", main=c("Difference with control points mean"), xaxt='n', xlab='Month', ylab='NDVI')
	# axis(side=1, at = seq(0,276,12))
	# dev.print(png, file = c(output_dir, "first_contol_point.png"), width = 1024, height = 768)
	
	# ggplot(data=final, mapping=aes(x=obs, y=diff)) + geom_line() + geom_point()
	
	final2 <- reshape::melt(final, id.var="months") # reshape to long format
	values <- NULL
	if (farm == 5) {
		values <- c('red', 'green', "blue", "chartreuse", "yellow",
					"orange", "grey","pink", "turquoise", "purple",
					"magenta", "deepskyblue", "forestgreen", "black")
	} else if (farm == 19) {
		values <- c('red', 'green', "blue", "chartreuse", "yellow",
					"orange", "grey","pink", "turquoise", "purple",
					"magenta", "deepskyblue", "forestgreen", "black")
	} else if (farm == 11) {
		values <- c('red', 'green', "blue", "chartreuse", "yellow",
					"orange", "grey","pink", "turquoise", "purple",
					"magenta", "forestgreen", "black")
	} else if (farm == 14) {
		values <- c('red', 'green', "blue", "chartreuse", "yellow",
					"orange", "forestgreen", "black")
	} else if (farm == 1) {
		values <- c('red', 'green', "blue", "chartreuse", "yellow",
					"orange", "grey","pink", "forestgreen", "black")
	}
	print(ggplot(data=final2, mapping=aes(x=months, y=value, color=variable)) + 
		geom_line() +
		scale_color_manual(values=values) +
		scale_x_continuous(name="Year (12 month intervals from 2000-2022)", limits=c(0, 276), breaks=seq(0,276,12)) +
		ylab("NDVI (Rescaled 0-2) Difference with control") + 
		ggtitle(label=sprintf("Time Series of Data Points in Farm %s with Mean", farm))) 
		# geom_line(data=final$mean, size = 3))
		# scale_size_manual(values=))
	#+
	# geom_point()
	# dev.print(png, file = c(output_dir, "all_control_points"), width = 1024, height = 768)
	
}

####################################################################################
###########     Attempt Greenbrown     ###########
####################################################################################

gb <- function(ouptut_dir, farm) {
	BP <- NULL
	if (farm == 5) {
		BP <- 3
	} else if (farm == 19) {
		BP <- 4
	} else if (farm == 11) {
		BP <- 4
	} else if (farm == 14) {
		BP <- 2
	} else if (farm == 1) {
		BP <- 3
	}
	# data <- data_filtered
	
	# # unfiltered and no gaps filled
	# # also try method SeasonalAdjusted  - seems to give more breakpoints (this method attempts to remove seasonality from TS which is needed if using monthly obs)
	trendmap_raw <<- TrendRaster(data_stack, start=c(2000, 1), freq=12, method="SeasonalAdjusted", breaks=BP, mosum.pval = 0.05, funSeasonalCycle = MeanSeasonalCycle)
	par(oma = c(2, 2, 1, 2)+0.1)
	plot(trendmap_raw, col=brgr.colors(20), legend.width=4)
	mtext(sprintf("Farm %s - SeasonalAdjusted (MeanSeasonalCycle) on Raw Data", farm), side = 3, adj = 0,
		  outer = TRUE)	# 
	# # gap filled
	# trendmap_gapfill <- TrendRaster(data_gapfill, start=c(2000, 1), freq=12, method="SeasonalAdjusted", breaks=BP, mosum.pval = 0.05, funSeasonalCycle = MeanSeasonalCycle)
	# par(oma = c(2, 2, 2, 2))
	# plot(trendmap_gapfill, col=brgr.colors(20), legend.width=4)
	# # dev.print(png, file = c(output_dir, "gb_gapfill.png"), width = 1024, height = 768)
	# 
	# # filtered and gap filled
	trendmap_filtered <- TrendRaster(data_filtered, start=c(2000, 1), freq=12, method="SeasonalAdjusted", breaks=BP, mosum.pval = 0.05, funSeasonalCycle = MeanSeasonalCycle)
	par(oma = c(2, 2, 1, 2)+0.1)
	plot(trendmap_filtered, col=brgr.colors(20), legend.width=4)
	mtext(sprintf("Farm %s - SeasonalAdjusted (MeanSeasonalCycle) on Filtered Data", farm), side = 3, adj = 0,
		  outer = TRUE)	
	
	# run each farm with different methods
	# for (i in 5:6) {
		
		trendmap_diff_sams <<- TrendRaster(data_diff, start=c(2000, 1), freq=12, method="SeasonalAdjusted", breaks=BP, mosum.pval = 0.05, funSeasonalCycle = MeanSeasonalCycle)
		par(oma = c(2, 2, 1, 2)+0.1)
		plot(trendmap_diff_sams, col=brgr.colors(20), legend.width=4)
		mtext(sprintf("Farm %s - SeasonalAdjusted (MeanSeasonalCycle) on Differenced Data", farm), side = 3, adj = 0,
			  outer = TRUE)


		print(paste("Extracting BP (SAMS) from final data_diff trendmap for farm: ", farm))
		n <- names(trendmap_diff_sams)
		breaks <- n[grepl("BP*", n)] # get names of BP* layers
		print(breaks)
		for (i in 1:length(breaks)) {
			b <- breaks[i]
			print(paste("Printing layer of: ", b))
			layer <- subset(trendmap_diff_sams, b)
			for (j in 1:nrow(layer)) {
				print(layer[j,])
			}
		}
		print(cat("...\n..."))


		# trendmap_diff_salm <<- TrendRaster(data_diff, start=c(2000, 1), freq=12, method="SeasonalAdjusted", breaks=BP, mosum.pval = 0.05, funSeasonalCycle = LmSeasonalCycle)
		# par(oma = c(2, 2, 1, 2)+0.1)
		# plot(trendmap_diff_salm, col=brgr.colors(20), legend.width=4)
		# mtext(sprintf("Farm %s - SeasonalAdjusted (LmSeasonalCycle) on Differenced Data", farm), side = 3, adj = 0,
		# 	  outer = TRUE)
		# 
		# print(paste("Extracting BP (SALM) from final data_diff trendmap for farm: ", farm))
		# n <- names(trendmap_diff_salm)
		# breaks <- n[grepl("BP*", n)] # get names of BP* layers
		# print(breaks)
		# for (i in 1:length(breaks)) {
		# 	b <- breaks[i]
		# 	print(paste("Printing layer of: ", b))
		# 	layer <- subset(trendmap_diff_sams, b)
		# 	for (j in 1:nrow(layer)) {
		# 		print(layer[j,])
		# 	}
		# }
		# print(cat("...\n..."))
		# 
		# 
		# trendmap_diff_stm <<- TrendRaster(data_diff, start=c(2000, 1), freq=12, method="STM", breaks=BP, mosum.pval=0.05)
		# par(oma = c(2, 2, 1, 2)+0.1)
		# plot(trendmap_diff_stm, col=brgr.colors(20), legend.width=4)
		# mtext(sprintf("Farm %s - STM on Differenced Data", farm), side = 3, adj = 0,
		# 	  outer = TRUE)
		# 
		# 	print(paste("Extracting BP (STM) from final data_diff trendmap for farm: ", farm))
		# n <- names(trendmap_diff_stm)
		# breaks <- n[grepl("BP*", n)] # get names of BP* layers
		# print(breaks)
		# for (i in 1:length(breaks)) {
		# 	b <- breaks[i]
		# 	print(paste("Printing layer of: ", b))
		# 	layer <- subset(trendmap_diff_stm, b)
		# 	for (j in 1:nrow(layer)) {
		# 		print(layer[j,])
		# 	}
		# }
		# print(cat("...\n..."))
# 
	# }

	
	

	
	
	# TODO: try greenbrown with different methods/params and then do a compareclassification 
	# to compare the greening/browning classifications 
	
	# trendmap_diff_sams <<- TrendRaster(data_diff, start=c(2000, 1), freq=12, method="SeasonalAdjusted", breaks=BP, mosum.pval = 0.05, funSeasonalCycle = MeanSeasonalCycle)
	# par(oma = c(2, 2, 1, 2)+0.1)
	# plot(trendmap_diff_sams, col=brgr.colors(20), legend.width=4)
	# mtext(sprintf("Farm %s - SeasonalAdjusted (MeanSeasonalCycle) on Differenced Data", farm), side = 3, adj = 0,
	# 	  outer = TRUE)
	# 
	# trendmap_diff_salm <<- TrendRaster(data, start=c(2000, 1), freq=12, method="SeasonalAdjusted", breaks=BP, mosum.pval = 0.05, funSeasonalCycle = LmSeasonalCycle)
	# par(oma = c(2, 2, 1, 2)+0.1)
	# plot(trendmap_diff_salm, col=brgr.colors(20), legend.width=4)
	# mtext("SeasonalAdjusted (LmSeasonalCycle) on Differenced Data", side = 3, adj = 0,
	# 	  outer = TRUE)
	# 
	# trendmap_diff_stm <<- TrendRaster(data, start=c(2000, 1), freq=12, method="STM", breaks=BP, mosum.pval=0.05)
	# par(oma = c(2, 2, 1, 2)+0.1)
	# plot(trendmap_diff_stm, col=brgr.colors(20), legend.width=4)
	# mtext("STM on Differenced Data", side = 3, adj = 0,
	# 	  outer = TRUE)
	# 
	# trendmap_diff_aat <<- TrendRaster(data, start=c(2000, 1), freq=12, method="AAT", breaks=BP, mosum.pval=1)
	# par(oma = c(2, 2, 1, 2)+0.1)
	# plot(trendmap_diff_aat, col=brgr.colors(20), legend.width=4)
	# mtext("AAT on Differenced Data", side = 3, adj = 0,
	# 	  outer = TRUE)
	
	# trendclassmap_sams <<- TrendClassification(trendmap_diff_sams, min.length=2)
	# plot(trendclassmap_sams, col=brgr.colors(3), legend.width=4)
	# 
	# trendclassmap_stm <<- TrendClassification(trendmap_diff_stm, min.length=2)
	# plot(trendclassmap_stm, col=brgr.colors(3), legend.width=4)
	# 
	# compare_sams_stm <<- CompareClassification(x=trendclassmap_sams, y=trendclassmap_stm,
											   # names=list('SAMS'=c("Br", "No", "Gr"), 'STM'=c("Br", "No", "Gr")))
	
	
}


####################################################################################
###########     Reading and writing raster bricks     ###########
####################################################################################

save_bricks <- function(farm) {
	# write data bricks
	outfile1 <- writeRaster(data_stack, filename=sprintf('bricks/%s/data_stack_%s.tif', farm, farm), format="GTiff", overwrite=TRUE,options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
	outfile2 <- writeRaster(data_gapfill, filename=sprintf('bricks/%s/data_gapfill_%s.tif', farm, farm), format="GTiff", overwrite=TRUE,options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
	outfile3 <- writeRaster(data_filtered, filename=sprintf('bricks/%s/data_filtered_%s.tif', farm, farm), format="GTiff", overwrite=TRUE,options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
	outfile4 <- writeRaster(data_diff, filename=sprintf('bricks/%s/data_diff_%s.tif', farm, farm), format="GTiff", overwrite=TRUE,options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
	
	# write control bricks	
	outfile5 <- writeRaster(control_stack, filename=sprintf('bricks/%s/control_stack_%s.tif', farm, farm), format="GTiff", overwrite=TRUE,options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
	outfile6 <- writeRaster(control_gapfill, filename=sprintf('bricks/%s/control_gapfill_%s.tif', farm, farm), format="GTiff", overwrite=TRUE,options=c("INTERLEAVE=BAND","COMPRESS=LZW"))	
	outfile7 <- writeRaster(control_filtered, filename=sprintf('bricks/%s/control_filtered_%s.tif', farm, farm), format="GTiff", overwrite=TRUE,options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
}


read_bricks <- function(farm) {
	data_diff_name <- NULL
	data_filtered_name <- NULL
	data_stack_name <- NULL
	control_filtered_name <- sprintf('bricks/%s/control_filtered_%s.tif', farm, farm)
	
	if (farm == 5) {
		data_diff_name <- 'bricks/5/data_diff_5.tif'
		data_filtered_name <- 'bricks/5/data_filtered_5.tif'
		data_stack_name <- 'bricks/5/data_stack_5.tif'
	} else if (farm == 19) {
		data_diff_name <- 'bricks/19/data_diff_19.tif'
		data_filtered_name <- 'bricks/19/data_filtered_19.tif'
		data_stack_name <- 'bricks/19/data_stack_19.tif'
	} else if (farm == 11) {
		data_diff_name <- 'bricks/11/data_diff_11.tif'
		data_filtered_name <- 'bricks/11/data_filtered_11.tif'
		data_stack_name <- 'bricks/11/data_stack_11.tif'
	} else if (farm == 14) {
		data_diff_name <- 'bricks/14/data_diff_14.tif'
		data_filtered_name <- 'bricks/14/data_filtered_14.tif'
		data_stack_name <- 'bricks/14/data_stack_14.tif'
	} else if (farm == 1) {
		data_diff_name <- 'bricks/1/data_diff_1.tif'
		data_filtered_name <- 'bricks/1/data_filtered_1.tif'
		data_stack_name <- 'bricks/1/data_stack_1.tif'
	}
	
	
	
	data_diff <<- brick(x=data_diff_name)
	data_filtered <<- brick(x=data_filtered_name)
	data_stack <<- brick(x=data_stack_name)
	# control_filtered <<- brick(x=control_filtered_name)
}

##################################################
#######      RUN ALL AREAS       #################
##################################################
# create a list of raster .tif files from the working directory
# data_path <- 'farm5_actual_ndvi_2005-2022'
# data_path <- 'farm5_area1_2000-2022'

# control <- 19
# data_path <- 'farm19/data'
# control_path <- 'farm19/control'

# control <- 11
# data_path <- 'farm11_self/data'
# control_path <- 'farm11_self/control'

data_list <- c("farm1_self/data", "farm5_self/data", "farm11_self/data", "farm14_WFW/data", "farm19_WFW/data")
control_list <- c("farm1_self/control", "farm5_self/control", "farm11_self/control", "farm14_WFW/control", "farm19_WFW/control")
control_farm <- c(1, 5, 11, 14, 19)
outputdir_list <- c("farm1_self", "farm5_self", "farm11_self", "farm14_WFW", "farm19_WFW")
# dev.off()
# run through all the areas of interest and output results
for (i in 2:2) {
	print(paste("iteration:     ", i))
	data_path <- data_list[i]
	control_path <- control_list[i]
	farm <- control_farm[i]
	output_dir <- outputdir_list[i]
	
	print(sprintf("data path:    %s",  data_path))
	print(sprintf("control path: %s", control_path))
	print(sprintf("farm:      %d", farm))
	print(sprintf("output dir:   %s", output_dir))
	
	# print("Loading Data")
	# load_data(data_path, control_path)
	
	# read in bricks from disk instead of computing them each time
	print("Reading brick")
	read_bricks(farm)
	
	# check for seasonality - need a univariate time series
	# print("Checking for seasonality")
	# Seasonality()
	
	# print("Filling gaps and smoothing with other thingy from GB")
	# original_filter()
	
	# print("Filling gaps")
	# fill_gaps()
	# print("Applying SG filter")
	# filter()
	# print("Differencing with control points")
	# difference(output_dir, farm)
	# print("Visualising stages")
	# visualise_stages(output_dir, farm)
	# print("Running greenbrown on stages")
	# gb(output_dir, farm)
	
	
	
	# save all proccessed raster bricks to disk so quicker to load
	# print("Saving bricks to disk")
	# save_bricks(farm)
	
	# nxt <- raster(ncols=5, nrows=10)
	# values(nxt) <- 1
	# plot(nxt, main=data_path)
	# print(cat(".\n"))
}

