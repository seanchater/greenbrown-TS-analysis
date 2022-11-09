library(raster)
library(reshape2)
library(rgdal)
library(scales)
library(dplyr)
library(greenbrown)
library(ggplot2)
library(signal) # for sgolayfilt function for gap filling
library(sf)






data_list <- c("farm5_self/data", "farm11_self/data", "farm14_WFW/data", "farm19_WFW/data")
control_list <- c("farm5_self/control", "farm11_self/control", "farm14_WFW/control", "farm19_WFW/control")
control_farm <- c(5, 11, 14, 19)
outputdir_list <- c("farm5_self", "farm11_self", "farm14_WFW", "farm19_WFW")

	print(paste("iteration:     ", 3))
	data_path <- data_list[3]
	control_path <- control_list[3]
	control <- control_farm[3]
	output_dir <- outputdir_list[3]
	
	print(sprintf("data path:    %s",  data_path))
	print(sprintf("control path: %s", control_path))
	print(sprintf("control:      %d", control))
	print(sprintf("output dir:   %s", output_dir))
	
	# print("Loading Data")
	# load_data(data_path, control_path)
	# print("Filling gaps")
	# fill_gaps()
	# print("Applying SG filter")
	# filter()
	# print("Differencing with control points")
	# difference(output_dir)
	# print("Visualising stages")
	# visualise_stages(output_dir)
	# print("Running greenbrown on stages")
	# gb(output_dir)



####################################################################################
###########      Load in Data      ###########
####################################################################################

# load_data <- function(data_path, control_path) {
	data_files = list.files(path=data_path, pattern='.tif$', full.names=TRUE)
	control_files = list.files(path=control_path, pattern='.tif$', full.names=TRUE)
	data_files
	
	# create a stack from the file list
	data_stack <- stack(data_files)
	control_stack <- stack(control_files)
	crs(data_stack)
	
	# Add 1 to all values in each raster stack
	# This is to make all NDVI values positive to simplify calculations further on such as with differencing section
	for (i in 1:nlayers(data_stack)) {
		data_stack[[i]] <- data_stack[[i]] + 1
	}
	for (i in 1:nlayers(control_stack)) {
		control_stack[[i]] <- control_stack[[i]] + 1
	}
# }

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

# fill_gaps <- function() {
	data_gapfill <- data_stack
	data_gapfill[] <- NA
	for( rl in 1:nrow(data_stack) ) { 
		v <- getValues(data_stack, rl, 1)
		data_gapfill[rl,] <- as.matrix(t(apply(v, MARGIN=1, FUN=impute.loess)))
	}
	
	control_gapfill <- control_stack
	control_gapfill[] <- NA
	for( rl in 1:nrow(control_stack) ) { 
		v <- getValues(control_stack, rl, 1)
		control_gapfill[rl,] <- as.matrix(t(apply(v, MARGIN=1, FUN=impute.loess)))
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
		data_gapfill[rl,] <- v
	}
	
	for( rl in 1:nrow(control_gapfill) ) { 
		v <- getValues(control_gapfill, rl, 1)
		v[is.na(v)] <- rowMeans(v, na.rm=TRUE)
		control_gapfill[rl,] <- v
	}
# }

####################################################################################
###########     Filtering      ###########
####################################################################################

sg_filter <- function(x) {
	v=as.vector(x)
	ndvi_ts = ts(v, start=c(2005,1), end=c(2022,12), frequency=12)
	x=sgolayfilt(ndvi_ts, p=1, n=3, ts=30)
}

# filter <- function() {
	# https://matinbrandt.wordpress.com/2014/12/02/smoothingfiltering-a-ndvi-time-series-using-a-savitzky-golay-filter-and-r/
	# check for gaps and fill - SG filtering https://www.sciencedirect.com/science/article/pii/S0034425716304692#s0120
	# data_stack_sg <- sgolayfilt(data_gapfill[[2:nlayers(data_gapfill)]])
	
	
	data_filtered <- data_gapfill
	data_filtered[] <- NA
	for( rl in 1:nrow(data_gapfill) ) { 
		v <- getValues(data_gapfill, rl, 1)
		data_filtered[rl,] <- as.matrix(t(apply(v, MARGIN=1, FUN=sg_filter)))
	}
	
	control_filtered <- control_gapfill
	control_filtered[] <- NA
	for( rl in 1:nrow(control_gapfill) ) { 
		v <- getValues(control_gapfill, rl, 1)
		control_filtered[rl,] <- as.matrix(t(apply(v, MARGIN=1, FUN=sg_filter)))
	}
# }

# TODO: rename layers in rasterbricks to appropriate names

####################################################################################
###########     Control difference     ###########
####################################################################################

# difference <- function(output_dir) {
	# set the points where you want to extract values from the raster stack
	df <- NULL
	points <- NULL
	if (control == 5) {
		# farm 5
		points <- c("5.2.1", "5.2.2", "5.2.3", "5.2.4", "5.2.5")
		df <- data.frame(point = points,
						 lat = c( -30.3039923733089, -30.3037228787237, -30.3037228787237, -30.3034533841385, -30.3034533841385),
						 lon = c(21.7666734076222, 21.7666734076222, 21.7669429022075, 21.7666734076222, 21.7669429022075))
	} else if (control == 19) {
		# farm 19
		points <- c("17.1.1", "17.1.2", "17.1.3", "17.1.4", "17.1.5", "17.1.6")
		df <- data.frame(point = points,
						 lat = c(-30.89768894, -30.89768894, -30.89795844, -30.89795844, -30.89795844, -30.89795844),
						 lon = c(22.08683297, 22.08656348, 22.08656348, 22.08683297, 22.08710247, 22.08737196))
	} else if (control == 11) {
		points <- c("11.4.1", "11.4.2", "11.4.3", "11.4.4", "11.4.5", "11.4.6", "11.4.7", "11.4.8")
		df <- data.frame(point = points,
						 lat = c(-30.6909866, -30.6907171, -30.6907171, -30.69044761, -30.69044761, -30.69017811, -30.69017811, -30.69017811),
						 lon = c(22.0555716, 22.0555716, 22.0558411, 22.0558411, 22.0555716, 22.0555716, 22.0558411, 22.05611059))
	} else if (control == 14) {
		points <- c("18.1.1", "18.1.2", "18.1.3", "18.1.4", "18.1.9", "18.1.10", "18.1.11", "18.1.12", "18.1.13")
		df <- data.frame(point = points,
						 lat = c(-30.6610727, -30.6610727, -30.6608032, -30.66134219, -30.66053371, -30.66026421, -30.65999472, -30.65999472, -30.66026421),
						 lon = c(22.22643117, 22.22616168, 22.22616168, 22.22643117, 22.22643117, 22.22643117, 22.22643117, 22.22616168, 22.22616168))
	}
	
	df_sf <- st_as_sf(x = df,
					  coords = c("lon", "lat"),
					  crs = 4326)
	
	# extract the point values
	values <- raster::extract(control_filtered, df_sf)
	
	# TODO: combine with coordinates: not sure if this is neccessary?
	# complete <- cbind(df_sf, values)
	
	final <- as.data.frame(t(values)) # transpose to have the observations for each point as rows instead of cols
	# need to average the pixel values of (minimum) 5 different control points per layer/month 
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
	data_diff <- data_filtered
	data_diff[] <- NA
	for( rl in 1:nrow(data_filtered) ) { 
		v <- getValues(data_filtered, rl, 1)
		data_diff[rl,] <- v - final$mean
	}
	
	obs <- names(final[1])
	first_point <- names(final[2])
	plot(final[,obs], final[,first_point], type="l", main=c("cleared", names(final[2]), " signal"), xaxt='n', xlab='Month', ylab='NDVI')
	axis(side=1, at = seq(0,276,12))
	# dev.print(png, file = c(output_dir, "first_contol_point.png"), width = 1024, height = 768)
	
	# ggplot(data=final, mapping=aes(x=obs, y=diff)) + geom_line() + geom_point()
	
	final2 <- reshape::melt(final, id.var="months") # reshape to long format
	values <- NULL
	if (control == 5) {
		values <- c('red', 'green', "blue", "black", "yellow")
	} else if (control == 19) {
		values <- c('red', 'green', "blue", "black", "yellow", "orange", "grey", "pink")
	} else if (control == 11) {
		values <- c('red', 'green', "blue", "black", "yellow", "orange", "grey", "pink", "turquoise", "purple")
	} else if (control == 14) {
		values <- c('red', 'green', "blue", "black", "yellow", "orange", "grey", "pink", "turquoise", "purple", "plum", "wheat", "magenta", "slategray2", "maroon1")
	}
	print(ggplot(data=final2, mapping=aes(x=months, y=value, color=variable)) + 
		  	geom_line() +
		  	scale_color_manual(values=values)) #+
	# geom_point()
	# dev.print(png, file = c(output_dir, "all_control_points"), width = 1024, height = 768)
	
# }

####################################################################################
###########     Plot the different stages    ###########
####################################################################################

# visualise_stages <- function(output_dir) {
	# plot original raster layers
	plot(data_stack)
	# dev.print(png, file = c(output_dir, "raw_data.png"), width = 1024, height = 768)
	
	# plot gap filled raster layers
	plot(data_gapfill)
	# dev.print(png, file = c(output_dir, "gapfill_data.png"), width = 1024, height = 768)
	
	# plot filtered raster layers
	plot(data_filtered)
	# dev.print(png, file = c(output_dir, "filtered_data.png"), width = 1024, height = 768)
	
	# plot raster layers minus the mean of control points
	plot(data_diff)
	# dev.print(png, file = c(output_dir, "diff_data.png"), width = 1024, height = 768)
# }

####################################################################################
###########     Attempt Greenbrown     ###########
####################################################################################

# gb <- function(ouptut_dir) {
	BP <- NULL
	if (control == 5) {
		BP <- 2
	} else if (control == 19) {
		BP <- 3
	} else if (control == 11) {
		BP <- 4
	} else if (control == 14) {
		BP <- 1
	}
	
	# unfiltered and no gaps filled
	# also try method SeasonalAdjusted  - seems to give more breakpoints (this method attempts to remove seasonality from TS which is needed if using monthly obs)
	trendmap2_non <- TrendRaster(data_stack, start=c(2000, 1), freq=12, method="SeasonalAdjusted", breaks=BP, mosum.pval = 0.05, funSeasonalCycle = MeanSeasonalCycle)
	plot(trendmap2_non, col=brgr.colors(20), legend.width=2)
	# dev.print(png, file = c(output_dir, "gb_raw.png"), width = 1024, height = 768)
	
	# gap filled
	trendmap2_gapfill <- TrendRaster(data_gapfill, start=c(2000, 1), freq=12, method="SeasonalAdjusted", breaks=BP, mosum.pval = 0.05, funSeasonalCycle = MeanSeasonalCycle)
	plot(trendmap2_gapfill, col=brgr.colors(20), legend.width=2)
	# dev.print(png, file = c(output_dir, "gb_gapfill.png"), width = 1024, height = 768)
	
	# filtered and gap filled
	trendmap2_all <- TrendRaster(data_filtered, start=c(2000, 1), freq=12, method="SeasonalAdjusted", breaks=BP, mosum.pval = 0.05, funSeasonalCycle = MeanSeasonalCycle)
	plot(trendmap2_all, col=brgr.colors(20), legend.width=2)
	# dev.print(png, file = c(output_dir, "gb_filtered.png"), width = 1024, height = 768)
	
	# filtered and gap filled and differenced with control mean
	trendmap4 <- TrendRaster(data_diff, start=c(2000, 1), freq=12, method="SeasonalAdjusted", breaks=BP, mosum.pval = 0.05, funSeasonalCycle = MeanSeasonalCycle)
	plot(trendmap4, col=brgr.colors(20), legend.width=2)
	# dev.print(png, file = c(output_dir, "gb_diff.png"), width = 1024, height = 768)
	
	# trendclassmap1 <- TrendClassification(data_filtered, min.length=8, max.pval=0.05)
	# plot(trendclassmap1, col=brgr.colors(3), legend.width=2, main="TrendClass1")
# }

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

# data_list <- c("farm5_self/data", "farm11_self/data", "farm14_WFW/data", "farm19_WFW/data")
# control_list <- c("farm5_self/control", "farm11_self/control", "farm14_WFW/control", "farm19_WFW/control")
# control_farm <- c(5, 11, 14, 19)
# outputdir_list <- c("farm5_self", "farm11_self", "farm14_WFW", "farm19_WFW")
# # dev.off()
# # run through all the areas of interest and output results
# for (i in 3:3) {
# 	print(paste("iteration:     ", i))
# 	data_path <- data_list[i]
# 	control_path <- control_list[i]
# 	control <- control_farm[i]
# 	output_dir <- outputdir_list[i]
# 	
# 	print(sprintf("data path:    %s",  data_path))
# 	print(sprintf("control path: %s", control_path))
# 	print(sprintf("control:      %d", control))
# 	print(sprintf("output dir:   %s", output_dir))
# 	
# 	print("Loading Data")
# 	load_data(data_path, control_path)
# 	print("Filling gaps")
# 	fill_gaps()
# 	print("Applying SG filter")
# 	filter()
# 	print("Differencing with control points")
# 	difference(output_dir)
# 	print("Visualising stages")
# 	visualise_stages(output_dir)
# 	print("Running greenbrown on stages")
# 	gb(output_dir)
	
	
	
# 	nxt <- raster(ncols=5, nrows=10)
# 	values(nxt) <- 1
# 	plot(nxt, main=data_path)
# }

