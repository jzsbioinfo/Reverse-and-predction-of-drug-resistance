
# 2021 01 11
# 31th_20X
# calculate growth rate for segmentation
# analyze cell cluster fluorescence intensity


setwd("E:/Project/2019__MitochondriaHeterogeneityHeritability/Software/YeastColonyTracking__Matlab_pipeline/A_20201110_new_analysis/31th_20X/")

# time interval = 60  20201110

# plot for one well, mCherry and PDR5 --------------------------------------
data_dir <- "Z:/MD/Analysis_data_Zhisheng/31th_20X/analysis/"


######
"D06 824"
"D09 476"

######

one_well <- "B02"

growth_rate_one_well <- read.csv(paste0(data_dir,"growth_rate_data/",one_well,"_all_valid_area.csv"))
mCherry_one_well <- read.csv(paste0(data_dir,"RFP_data/",one_well,"_all_valid_area_RFP.csv"))
YFP_one_well <- read.csv(paste0(data_dir,"GFP_data/",one_well,"_all_valid_area_GFP.csv"))


# 1. function to process one well----------------------------------

# what about NA
process_one_well <- function(all_valid_area) {
  
  
  
  time_point = 1:nrow(all_valid_area)
  
  # 1. filter out unvalid data for time points
  all_valid_area_filter <- all_valid_area[,colSums(all_valid_area>0)>=10] # filter 1 : more than 10 point 
  
  # avoid data.frame to vector
  if ( is.vector(all_valid_area_filter) ) {
    all_valid_area_filter <- data.frame(all_valid_area_filter)
    colnames(all_valid_area_filter) <- colnames(all_valid_area)[colSums(all_valid_area>0)>=10]
  }
  ncol(all_valid_area_filter)
  
  # get max and no_zero_min of columns
  
  
  col_no_zero_mins <-  function(x) {
    
    no_zero_min <-  function(x) {
      x <- x[x>0]
      x <- min(x)
      return(x)
    }
    
    apply(x, 2, no_zero_min)
  }
  
  col_last_not_zero <-  function(x) {
    
    last_not_zero <- function(x) {
      x <- x[x>0]
      x <- x[length(x)]
    }
    
    apply(x, 2, last_not_zero)
  }
  
  
  colmaxs <-  function(x) {
    apply(x, 2, max)
  }
  
  all_valid_area_filter <- all_valid_area_filter[,col_last_not_zero(all_valid_area_filter)/col_no_zero_mins(all_valid_area_filter) >3]  #filter 2 : 3*double
  
  # avoid data.frame to vector
  if ( is.vector(all_valid_area_filter) ) {
    all_valid_area_filter <- data.frame(all_valid_area_filter)
    colnames(all_valid_area_filter) <- colnames(all_valid_area)[colSums(all_valid_area>0)>=10]
  }
  # 2. five point sliding window to get max grwoth rate with R square > 0.9
  
  three_point_growth_rate <- function(data) {
    # function to process one cluster get max grwoth rate with R square > 0.9
    
    
    
    max_growth_rate <- function(data) {
      
      temp_max_growth_rate <- data.frame(matrix(ncol = 2, nrow = 0))
      colnames(temp_max_growth_rate) <- c("sw_Growth_rate","sw_R_square")
      
      index <- data>0
      data <- data[index]
      
      # start from time point 6 -----
      for (i in 6:(length(data)-4)) {
        
        
        gr <- summary(lm(log2(data[c(i:(i+4))]) ~ time_point[index][c(i:(i+4))]))
        temp_max_growth_rate[i-5,1] <- gr$coefficients[2]
        temp_max_growth_rate[i-5,2] <- ifelse(is.nan(gr$adj.r.squared),0,gr$adj.r.squared)
      }
      
      # have > 0.9 return all > 0.9, else return the max R square one
      if (max(temp_max_growth_rate[,2])>0.9) {
        temp_max_growth_rate <- temp_max_growth_rate[temp_max_growth_rate[,2]>0.9,] 
      } else {
        temp_max_growth_rate <- temp_max_growth_rate[temp_max_growth_rate[,2]==max(temp_max_growth_rate[,2]),]
      }
      
      
      
      temp_max_growth_rate <- temp_max_growth_rate[temp_max_growth_rate[,1]==max(temp_max_growth_rate[,1]),]
      
      return(temp_max_growth_rate)
    }
    
    
    
    
    
    if ( is.data.frame(data)) {
      three_point_growth_rate <-  sapply(data, max_growth_rate,simplify = TRUE, USE.NAMES = F)
      three_point_growth_rate <- t(three_point_growth_rate)
    } else {
      three_point_growth_rate <-  data.frame("sw_Growth_rate"=c(0,0), "sw_R_square"=c(0,0))
    }
    
    
    
    return(three_point_growth_rate)
  }
  
  three_point_growth_rate_result <- three_point_growth_rate(all_valid_area_filter)
  
  
  # 3. lm fit for all point >0 and after 6 time points to get the growth rate and R2  
  
  growth_rate <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(growth_rate) <- c("Growth_rate","R_square")
  
  if (is.data.frame(all_valid_area_filter)) {
    for (i in 1:ncol(all_valid_area_filter)) {
      index <- all_valid_area_filter[,i] > 0
      gr <- lm(log2(all_valid_area_filter[,i][index])[6:length(index)] ~ time_point[index][6:length(index)])
      gr <- summary(gr)
      growth_rate[i,1] <- gr$coefficients[2]
      growth_rate[i,2] <- gr$adj.r.squared
    }
  } else {
    growth_rate <-  data.frame("Growth_rate"=c(0,0), "R_square"=c(0,0))
  }
  
  
  # tidy the result table
  colname_temp <- colnames(three_point_growth_rate_result)
  rowname_temp <- rownames(three_point_growth_rate_result)
  three_point_growth_rate_result <- as.data.frame(matrix(unlist(three_point_growth_rate_result),
                                                         nrow = nrow(three_point_growth_rate_result),
                                                         byrow = F))
  colnames(three_point_growth_rate_result) <- colname_temp
  rownames(three_point_growth_rate_result) <- rowname_temp
  
  growth_rate_table <- as.data.frame(cbind(three_point_growth_rate_result,growth_rate))
  
  growth_rate_table$max_growth_rate <- 0 # not valid
  
  for (i in 1:nrow(growth_rate_table)) {
    if (growth_rate_table[i,2] >=0.9 & growth_rate_table[i,4] >=0.9) {
      growth_rate_table$max_growth_rate[i] = max(growth_rate_table[i,1], growth_rate_table[i,3])
    }
    
    if (growth_rate_table[i,2] >=0.9 & growth_rate_table[i,4] <0.9) {
      growth_rate_table$max_growth_rate[i] = growth_rate_table[i,1]
    }
    
    if (growth_rate_table[i,2] <0.9 & growth_rate_table[i,4] >=0.9) {
      growth_rate_table$max_growth_rate[i] = growth_rate_table[i,3]
    }
  }
  
  
  # filter all point fit R>0
  all_valid_area_filter <- all_valid_area_filter[,growth_rate_table$Growth_rate > 0]
  growth_rate_table <- growth_rate_table[growth_rate_table$Growth_rate > 0,]
  
  # filter R>0.9 and growth rate>0
  all_valid_area_filter <- all_valid_area_filter[,growth_rate_table$max_growth_rate > 0]
  growth_rate_table <- growth_rate_table[growth_rate_table$max_growth_rate > 0,]
  
  all_result <- list(growth_rate_table,all_valid_area_filter)
  return(all_result)
  
}

valid_one_well <- process_one_well(growth_rate_one_well)

plot(density(valid_one_well[[1]]$max_growth_rate))  # max_growth_rate
hist(valid_one_well[[1]]$max_growth_rate)
plot(density(valid_one_well[[1]]$Growth_rate))  # Growth_rate
hist(valid_one_well[[1]]$Growth_rate)

cell_clusters <-  rownames(valid_one_well[[1]])


growth_rate_one_well <- growth_rate_one_well[,colnames(growth_rate_one_well) %in% cell_clusters]

mCherry_one_well <- mCherry_one_well[,colnames(mCherry_one_well) %in% cell_clusters]

YFP_one_well <- YFP_one_well[,colnames(YFP_one_well) %in% cell_clusters]


mCherry_one_well_normbyarea <- mCherry_one_well/growth_rate_one_well
YFP_one_well_normbyarea <- YFP_one_well/growth_rate_one_well

# library("plot3D")
# 
# scatter3D(mCherry_one_well$cell_cluster8, 
#           YFP_one_well$cell_cluster8,
#           1:30)
# 
# 
# colnames(growth_rate_one_well)



#-------- mCherry -----------------------
mCherry_col1 <- filtered_mCherry_data[[1]]
for ( i in 2:6) {
  temp <- filtered_mCherry_data[[i]]
  mCherry_col1 <- cbind(mCherry_col1,temp)
}

plot(sort(as.numeric(mCherry_col1[1,])))
plot(sort(colSums(mCherry_col1[c(1:3),])))


plot(sort(colSums(mCherry_one_well[c(1:5),])))
plot(sort(colSums(mCherry_one_well[c(1:3),])))

plot(sort(as.numeric(mCherry_one_well[1,])))





YFP_col1 <- filtered_YFP_data[[1]]
for ( i in 2:6) {
  temp <- filtered_YFP_data[[i]]
  YFP_col1 <- cbind(YFP_col1,temp)
}


plot(sort(as.numeric(mCherry_col1[1,])/as.numeric(YFP_col1[1,])))

plot(sort(colSums(mCherry_col1[c(1:3),])/colSums(YFP_col1[c(1:3),])))

plot(sort(colSums(mCherry_col1[c(1:5),])/colSums(YFP_col1[c(1:5),])))

mCherry_col1
YFP_col1


#----------------------------------------


png(paste0(data_dir,"/","20201114 mix_colony_analysis/",one_well,"_mCherry.png"), width = 500, height = 300, units="mm",res=300)



oldpar <- par()
par(mfcol=c(6,6))
par(mar=c(1,1,1,1))
for ( i in 1:ncol(mCherry_one_well)) {
  temp_data <- mCherry_one_well[,i]
  
  temp_data <- temp_data[temp_data>0]
  
  
  plot(1:length(temp_data),temp_data,xlim=c(0,30),axes=FALSE, frame=TRUE, pch=22, col="pink",bg="purple")
  axis(side=1,  tck= -0.02)
  axis(side=2,  tck= -0.02)
  text(0,(max(temp_data)),cell_clusters[i])
}
par(oldpar)


dev.off()

#------- YFP-------------------------------
png(paste0(data_dir,"/","20201114 mix_colony_analysis/",one_well,"_YFP.png"), width = 500, height = 300, units="mm",res=300)

oldpar <- par()
par(mfcol=c(6,11))
par(mar=c(1,1,1,1))
for ( i in 1:ncol(YFP_one_well_normbyarea)) {
  temp_data <- YFP_one_well[,i]
  
  temp_data <- temp_data[temp_data>0]
  
  
  plot(1:length(temp_data),temp_data,xlim=c(0,30),axes=FALSE, frame=TRUE, pch=22, col="black",bg="green")
  axis(side=1,  tck= -0.02)
  axis(side=2,  tck= -0.02)
  text(0,(max(temp_data)),cell_clusters[i])
}
par(oldpar)

dev.off()

#----------- mCherry + YFP --------------------------------

png(paste0(data_dir,"/","20201114 mix_colony_analysis/",one_well,"_mCherry_YFP.png"), width = 500, height = 300, units="mm",res=300)

oldpar <- par()
par(mfcol=c(6,11))
par(mar=c(1,1,1,1))
for ( i in 1:ncol(mCherry_one_well)) {
  
  i= which(cell_clusters=="cell_cluster476")
  temp_data_mCherry <- mCherry_one_well[,i]
  
  temp_data_mCherry <- temp_data_mCherry[temp_data_mCherry>0]
  
  
  temp_data_YFP <- YFP_one_well[,i]
  
  temp_data_YFP <- temp_data_YFP[temp_data_YFP>0]
  
  
  plot(1:length(temp_data_mCherry),temp_data_mCherry,xlim=c(0,30),axes=FALSE, frame=TRUE, pch=22, col="pink",bg="purple",ylim = c(min(temp_data_mCherry,temp_data_YFP),max(temp_data_mCherry,temp_data_YFP)))
  points(1:length(temp_data_YFP),temp_data_YFP,xlim=c(0,30),axes=FALSE, frame=TRUE, pch=22, col="black",bg="green")
  axis(side=1,  tck= -0.02)
  axis(side=2,  tck= -0.02)
  text(0,(max(temp_data_mCherry,temp_data_YFP)),cell_clusters[i])
}
par(oldpar)

dev.off()


#----------mCherry + YFP  same Y axis --------------------------

png(paste0(data_dir,"/","20201114 mix_colony_analysis/",one_well,"_mCherry_YFP_sameYaxis.png"), width = 500, height = 300, units="mm",res=300)


oldpar <- par()
par(mfcol=c(6,11))
par(mar=c(1,1,1,1))
for ( i in 1:ncol(mCherry_one_well)) {
  temp_data_mCherry <- mCherry_one_well[,i]
  
  temp_data_mCherry <- 5*temp_data_mCherry[temp_data_mCherry>0]
  
  
  temp_data_YFP <- YFP_one_well[,i]
  
  temp_data_YFP <- temp_data_YFP[temp_data_YFP>0]
  
  
  plot(1:length(temp_data_mCherry),temp_data_mCherry,xlim=c(0,30),axes=FALSE, frame=TRUE, pch=22, col="pink",bg="purple",ylim = c(min(mCherry_one_well,YFP_one_well),max(mCherry_one_well,YFP_one_well)))
  points(1:length(temp_data_YFP),temp_data_YFP,xlim=c(0,30),axes=FALSE, frame=TRUE, pch=22, col="black",bg="green")
  axis(side=1,  tck= -0.02)
  axis(side=2,  tck= -0.02)
  text(0,(max(mCherry_one_well,YFP_one_well)),cell_clusters[i])
}
par(oldpar)

dev.off()



# grwoth curve --------------------------------------

# same Y axis
png(paste0(data_dir,"/","20201114 mix_colony_analysis/",one_well,"_growth_curve_sameYaxis.png"), width = 500, height = 300, units="mm",res=300)


oldpar <- par()
par(mfcol=c(6,6))
par(mar=c(1,1,1,1))
for ( i in 1:ncol(growth_rate_one_well)) {
  temp_data <- growth_rate_one_well[,i]
  
  temp_data <- temp_data[temp_data>0]
  
  
  plot(1:length(temp_data),temp_data,xlim=c(0,30),ylim = c(0,max(growth_rate_one_well)),
       axes=FALSE, frame=TRUE, pch=22, col="#fdbb84",bg="#3182bd")
  axis(side=1,  tck= -0.02)
  axis(side=2,  tck= -0.02)
  text(0,(max(growth_rate_one_well)),cell_clusters[i])
}
par(oldpar)

dev.off()



# different Y axis
png(paste0(data_dir,"/","20201114 mix_colony_analysis/",one_well,"_growth_curve.png"), width = 500, height = 300, units="mm",res=300)

oldpar <- par()
par(mfcol=c(6,11))
par(mar=c(1,1,1,1))
for ( i in 1:ncol(growth_rate_one_well)) {
  temp_data <- growth_rate_one_well[,i]
  
  temp_data <- temp_data[temp_data>0]
  
  
  plot(1:length(temp_data),temp_data,xlim=c(0,30),#ylim = c(0,max(growth_rate_one_well)),
       axes=FALSE, frame=TRUE, pch=22, col="#fdbb84",bg="#3182bd")
  axis(side=1,  tck= -0.02)
  axis(side=2,  tck= -0.02)
  text(0,(max(temp_data)),cell_clusters[i])
}
par(oldpar)

dev.off()


# same Y axis, log2
png(paste0(data_dir,"/","20201114 mix_colony_analysis/",one_well,"_growth_curve_sameYaxislog2.png"), width = 500, height = 300, units="mm",res=300)


oldpar <- par()
par(mfcol=c(6,11))
par(mar=c(1,1,1,1))
for ( i in 1:ncol(growth_rate_one_well)) {
  temp_data <- growth_rate_one_well[,i]
  
  temp_data <- temp_data[temp_data>0]
  
  
  plot(1:length(temp_data),log2(temp_data),xlim=c(0,30),ylim = c(0,max(log2(growth_rate_one_well))),
       axes=FALSE, frame=TRUE, pch=22, col="#fdbb84",bg="#3182bd")
  axis(side=1,  tck= -0.02)
  axis(side=2,  tck= -0.02)
  text(0,(max(log2(growth_rate_one_well))),cell_clusters[i])
}
par(oldpar)

dev.off()



# classify class1 and class2 --------------------

# 1.separate survival and non_survival
onewell_grwoth_rate <- valid_one_well[[1]]
onewell_grwoth_rate <- onewell_grwoth_rate[order(onewell_grwoth_rate$max_growth_rate),]
onewell_grwoth_rate <- onewell_grwoth_rate[order(onewell_grwoth_rate$Growth_rate),]
plot(onewell_grwoth_rate$Growth_rate)

big_colony <- onewell_grwoth_rate[onewell_grwoth_rate$Growth_rate>0.15,]


big_colony_name <- rownames(big_colony)
big_colony_name <- big_colony_name[order(as.numeric(gsub("cell_cluster","",big_colony_name)))]


# mCherry+ YFP ----
oldpar <- par()
par(mfcol=c(4,7))
par(mar=c(1,1,1,1))
for ( i in 1:length(big_colony_name)) {
  
  temp_data_mCherry <- mCherry_one_well[,colnames(mCherry_one_well)==big_colony_name[i]]
  
  temp_data_mCherry <- temp_data_mCherry[temp_data_mCherry>0]
  
  
  temp_data_YFP <- YFP_one_well[,colnames(YFP_one_well)==big_colony_name[i]]
  
  temp_data_YFP <- temp_data_YFP[temp_data_YFP>0]
  
  
  plot(1:length(temp_data_mCherry),temp_data_mCherry,xlim=c(0,30),axes=FALSE, frame=TRUE, pch=22, col="pink",bg="purple",ylim = c(min(temp_data_mCherry,temp_data_YFP),max(temp_data_mCherry,temp_data_YFP)))
  points(1:length(temp_data_YFP),temp_data_YFP,xlim=c(0,30),axes=FALSE, frame=TRUE, pch=22, col="black",bg="green")
  axis(side=1,  tck= -0.02)
  axis(side=2,  tck= -0.02)
  text(0,(max(temp_data_mCherry,temp_data_YFP)),big_colony_name[i])
}
par(oldpar)


# mCherry + YFP sameYaxis---------


oldpar <- par()
par(mfcol=c(4,7))
par(mar=c(1,1,1,1))
for ( i in 1:length(big_colony_name)) {
  temp_data_mCherry <- mCherry_one_well[,colnames(mCherry_one_well)==big_colony_name[i]]
  
  temp_data_mCherry <- temp_data_mCherry[temp_data_mCherry>0]
  
  
  temp_data_YFP <- YFP_one_well[,colnames(YFP_one_well)==big_colony_name[i]]
  
  temp_data_YFP <- temp_data_YFP[temp_data_YFP>0]
  
  
  plot(1:length(temp_data_mCherry),temp_data_mCherry,xlim=c(0,30),axes=FALSE, frame=TRUE, pch=22, col="pink",bg="purple",ylim = c(min(mCherry_one_well,YFP_one_well),max(mCherry_one_well,YFP_one_well)))
  points(1:length(temp_data_YFP),temp_data_YFP,xlim=c(0,30),axes=FALSE, frame=TRUE, pch=22, col="black",bg="green")
  axis(side=1,  tck= -0.02)
  axis(side=2,  tck= -0.02)
  text(0,(max(mCherry_one_well,YFP_one_well)),big_colony_name[i])
}
par(oldpar)

# grwoth curve --------------------------------------

# same Y axis

oldpar <- par()
par(mfcol=c(6,6))
par(mar=c(1,1,1,1))
for ( i in 1:length(big_colony_name)) {
  temp_data <- growth_rate_one_well[,colnames(growth_rate_one_well)==big_colony_name[i]]
  
  temp_data <- temp_data[temp_data>0]
  
  
  plot(1:length(temp_data),temp_data,xlim=c(0,30),ylim = c(0,max(growth_rate_one_well)),
       axes=FALSE, frame=TRUE, pch=22, col="#fdbb84",bg="#3182bd")
  axis(side=1,  tck= -0.02)
  axis(side=2,  tck= -0.02)
  text(0,(max(growth_rate_one_well)),cell_clusters[i])
}
par(oldpar)








#####################################################################


#-------- mCherry avarage  -----------------------

oldpar <- par()
par(mfcol=c(6,10))
par(mar=c(1,1,1,1))
for ( i in 1:ncol(mCherry_one_well_normbyarea)) {
  temp_data <- mCherry_one_well_normbyarea[,i]
  
  temp_data <- temp_data[temp_data>0]
  
  
  plot(1:length(temp_data),temp_data,xlim=c(0,30),axes=FALSE, frame=TRUE, pch=22, col="pink",bg="purple")
  axis(side=1,  tck= -0.05)
  axis(side=2,  tck= -0.05)
  text(0,(max(temp_data)),cell_clusters[i])
}
par(oldpar)



#------- YFP  avarage-------------------------------

oldpar <- par()
par(mfcol=c(6,10))
par(mar=c(1,1,1,1))
for ( i in 1:ncol(YFP_one_well_normbyarea)) {
  temp_data <- YFP_one_well_normbyarea[,i]
  
  temp_data <- temp_data[temp_data>0]
  
  
  plot(1:length(temp_data),temp_data,xlim=c(0,30),axes=FALSE, frame=TRUE, pch=22, col="black",bg="green")
  axis(side=1,  tck= -0.05)
  axis(side=2,  tck= -0.05)
  text(0,(max(temp_data)),cell_clusters[i])
}
par(oldpar)


# mCherry + YFP avarage --------------------------------------


oldpar <- par()
par(mfcol=c(6,10))
par(mar=c(1,1,1,1))
for ( i in 1:ncol(mCherry_one_well_normbyarea)) {
  temp_data_mCherry <- mCherry_one_well_normbyarea[,i]
  
  temp_data_mCherry <- temp_data_mCherry[temp_data_mCherry>0]
  
  
  temp_data_YFP <- YFP_one_well_normbyarea[,i]
  
  temp_data_YFP <- temp_data_YFP[temp_data_YFP>0]
  
  ylim_low <- min(na.omit(temp_data_mCherry),na.omit(temp_data_YFP))
  ylim_up <- max(na.omit(temp_data_mCherry),na.omit(temp_data_YFP))
  
  
  plot(1:length(temp_data_mCherry),temp_data_mCherry,xlim=c(0,30),axes=FALSE, frame=TRUE, pch=22, col="pink",bg="purple",ylim = c(ylim_low,ylim_up))
  points(1:length(temp_data_YFP),temp_data_YFP,xlim=c(0,30),axes=FALSE, frame=TRUE, pch=22, col="black",bg="green")
  axis(side=1,  tck= -0.05)
  axis(side=2,  tck= -0.05)
  text(0,(max(temp_data_mCherry)),cell_clusters[i])
}
par(oldpar)











# Process bright field data#######################################################################
# read in segmentation area data-----------------
data_dir <- "Z:/MD/Analysis_data_Zhisheng/31th_20X/analysis/"
growth_rate_filenames <- list.files(paste0(data_dir,"growth_rate_data"), pattern=".*.csv", full.names=TRUE)

growth_rate_filenames_well_label <- gsub("_all_valid_area.*","",gsub(".*growth_rate_data/","",growth_rate_filenames))

plot_order <- order(substr(growth_rate_filenames_well_label,2,3))

growth_rate_file_list <- lapply(growth_rate_filenames,read.csv)

growth_rate_file_list[[1]][1:30,1:5]

#############################

library(purrr)

# 1. function to process one well----------------------------------

# what about NA
process_one_well <- function(all_valid_area) {
  
  
  
  time_point = 1:nrow(all_valid_area)
  
  # 1. filter out unvalid data for time points
  all_valid_area_filter <- all_valid_area[,colSums(all_valid_area>0)>=10] # filter 1 : more than 10 point 
  
  # avoid data.frame to vector
  if ( is.vector(all_valid_area_filter) ) {
    all_valid_area_filter <- data.frame(all_valid_area_filter)
    colnames(all_valid_area_filter) <- colnames(all_valid_area)[colSums(all_valid_area>0)>=10]
  }
  ncol(all_valid_area_filter)
  
  # get max and no_zero_min of columns
  
  
  col_no_zero_mins <-  function(x) {
    
    no_zero_min <-  function(x) {
      x <- x[x>0]
      x <- min(x)
      return(x)
    }
    
    apply(x, 2, no_zero_min)
  }
  
  col_last_not_zero <-  function(x) {
    
    last_not_zero <- function(x) {
      x <- x[x>0]
      x <- x[length(x)]
    }
    
    apply(x, 2, last_not_zero)
  }
  
  
  colmaxs <-  function(x) {
    apply(x, 2, max)
  }
  
  all_valid_area_filter <- all_valid_area_filter[,col_last_not_zero(all_valid_area_filter)/col_no_zero_mins(all_valid_area_filter) >3]  #filter 2 : 3*double
  
  # avoid data.frame to vector
  if ( is.vector(all_valid_area_filter) ) {
    all_valid_area_filter <- data.frame(all_valid_area_filter)
    colnames(all_valid_area_filter) <- colnames(all_valid_area)[colSums(all_valid_area>0)>=10]
  }
  # 2. five point sliding window to get max grwoth rate with R square > 0.9
  
  three_point_growth_rate <- function(data) {
    # function to process one cluster get max grwoth rate with R square > 0.9
    
    
    
    max_growth_rate <- function(data) {
      
      temp_max_growth_rate <- data.frame(matrix(ncol = 2, nrow = 0))
      colnames(temp_max_growth_rate) <- c("sw_Growth_rate","sw_R_square")
      
      index <- data>0
      data <- data[index]
      
      # start from time point 6 -----
      for (i in 6:(length(data)-4)) {
        
        
        gr <- summary(lm(log2(data[c(i:(i+4))]) ~ time_point[index][c(i:(i+4))]))
        temp_max_growth_rate[i-5,1] <- gr$coefficients[2]
        temp_max_growth_rate[i-5,2] <- ifelse(is.nan(gr$adj.r.squared),0,gr$adj.r.squared)
      }
      
      # have > 0.9 return all > 0.9, else return the max R square one
      if (max(temp_max_growth_rate[,2])>0.9) {
        temp_max_growth_rate <- temp_max_growth_rate[temp_max_growth_rate[,2]>0.9,] 
      } else {
        temp_max_growth_rate <- temp_max_growth_rate[temp_max_growth_rate[,2]==max(temp_max_growth_rate[,2]),]
      }
      
      
      
      temp_max_growth_rate <- temp_max_growth_rate[temp_max_growth_rate[,1]==max(temp_max_growth_rate[,1]),]
      
      return(temp_max_growth_rate)
    }
    
    
    
    
    
    if ( is.data.frame(data)) {
      three_point_growth_rate <-  sapply(data, max_growth_rate,simplify = TRUE, USE.NAMES = F)
      three_point_growth_rate <- t(three_point_growth_rate)
    } else {
      three_point_growth_rate <-  data.frame("sw_Growth_rate"=c(0,0), "sw_R_square"=c(0,0))
    }
    
    
    
    return(three_point_growth_rate)
  }
  
  three_point_growth_rate_result <- three_point_growth_rate(all_valid_area_filter)
  
  
  # 3. lm fit for all point >0 and after 6 time points to get the growth rate and R2  
  
  growth_rate <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(growth_rate) <- c("Growth_rate","R_square")
  
  if (is.data.frame(all_valid_area_filter)) {
    for (i in 1:ncol(all_valid_area_filter)) {
      index <- all_valid_area_filter[,i] > 0
      gr <- lm(log2(all_valid_area_filter[,i][index])[6:length(index)] ~ time_point[index][6:length(index)])
      gr <- summary(gr)
      growth_rate[i,1] <- gr$coefficients[2]
      growth_rate[i,2] <- gr$adj.r.squared
    }
  } else {
    growth_rate <-  data.frame("Growth_rate"=c(0,0), "R_square"=c(0,0))
  }
  
  
  # tidy the result table
  colname_temp <- colnames(three_point_growth_rate_result)
  rowname_temp <- rownames(three_point_growth_rate_result)
  three_point_growth_rate_result <- as.data.frame(matrix(unlist(three_point_growth_rate_result),
                                                         nrow = nrow(three_point_growth_rate_result),
                                                         byrow = F))
  colnames(three_point_growth_rate_result) <- colname_temp
  rownames(three_point_growth_rate_result) <- rowname_temp
  
  growth_rate_table <- as.data.frame(cbind(three_point_growth_rate_result,growth_rate))
  
  growth_rate_table$max_growth_rate <- 0 # not valid
  
  for (i in 1:nrow(growth_rate_table)) {
    if (growth_rate_table[i,2] >=0.9 & growth_rate_table[i,4] >=0.9) {
      growth_rate_table$max_growth_rate[i] = max(growth_rate_table[i,1], growth_rate_table[i,3])
    }
    
    if (growth_rate_table[i,2] >=0.9 & growth_rate_table[i,4] <0.9) {
      growth_rate_table$max_growth_rate[i] = growth_rate_table[i,1]
    }
    
    if (growth_rate_table[i,2] <0.9 & growth_rate_table[i,4] >=0.9) {
      growth_rate_table$max_growth_rate[i] = growth_rate_table[i,3]
    }
  }
  

  # filter all point fit R>0
  all_valid_area_filter <- all_valid_area_filter[,growth_rate_table$Growth_rate > 0]
  growth_rate_table <- growth_rate_table[growth_rate_table$Growth_rate > 0,]
  
  # filter R>0.9 and growth rate>0
  all_valid_area_filter <- all_valid_area_filter[,growth_rate_table$max_growth_rate > 0]
  growth_rate_table <- growth_rate_table[growth_rate_table$max_growth_rate > 0,]
  
  all_result <- list(growth_rate_table,all_valid_area_filter)
  return(all_result)
  
}


#--this is for one well--------------------------

all_valid_area <- growth_rate_file_list[[1]]  # D02

one_well <- process_one_well(all_valid_area)

one_well[[1]]$max_growth_rate <- tt[[1]]$max_growth_rate/time_interval
one_well[[1]]$sw_Growth_rate <- tt[[1]]$sw_Growth_rate/time_interval
one_well[[1]]$Growth_rate <- tt[[1]]$Growth_rate/time_interval


one_well[[1]]
one_well[[2]]

# library(openxlsx)
# write.xlsx(one_well[[1]],"./20191011_20191003_TL25_YFP_H01_2/H01_calculated_grwoth_rate.xlsx",col.names=T,row.names=T)
# write.xlsx(one_well[[2]],"./20191011_20191003_TL25_YFP_H01_2/H01_raw_area_size.xlsx",col.names=T)



# 2. apply the function to all wells to get result------------------------------

#debugonce(process_one_well) # this is used to debug function
growth_rate_file_list_result <- map(growth_rate_file_list,possibly(process_one_well, NA_real_))

# growth_rate_file_list_result <- lapply(growth_rate_file_list, process_one_well)

growth_rate_file_list_result <- growth_rate_file_list_result[c(plot_order)]


#############################



# 1. function to process one well----------------------------------


process_one_well <- function(all_valid_area) {


  
  time_point = 1:nrow(all_valid_area)
  
  # 1. filter out unvalid data for time points
  all_valid_area_filter <- all_valid_area[,colSums(all_valid_area>0)>=8] # filter 1 : more than 8 point 
  
  ncol(all_valid_area_filter)
  
  # get max and no_zero_min of columns
 
  
  col_no_zero_mins <-  function(x) {
    
    no_zero_min <-  function(x) {
      x <- x[x>0]
      x <- min(x)
      return(x)
    }
    
    apply(x, 2, no_zero_min)
  }
  
  col_last_not_zero <-  function(x) {
    
    last_not_zero <- function(x) {
      x <- x[x>0]
      x <- x[length(x)]
    }
    
    apply(x, 2, last_not_zero)
  }
  
  
  colmaxs <-  function(x) {
    apply(x, 2, max)
  }
  
  all_valid_area_filter <- all_valid_area_filter[,col_last_not_zero(all_valid_area_filter)/col_no_zero_mins(all_valid_area_filter) >3]  #filter 2 : 3*double
  

  # 2. five point sliding window to get max grwoth rate with R square > 0.9

  three_point_growth_rate <- function(data) {
    # function to process one cluster get max grwoth rate with R square > 0.9
   
    
    
    max_growth_rate <- function(data) {
      
      temp_max_growth_rate <- data.frame(matrix(ncol = 2, nrow = 0))
      colnames(temp_max_growth_rate) <- c("sw_Growth_rate","sw_R_square")
      
      index <- data>0
      data <- data[index]
      
      # time point 4 -----
      for (i in 4:(length(data)-4)) {
        

        gr <- summary(lm(log2(data[c(i:(i+4))]) ~ time_point[index][c(i:(i+4))]))
        temp_max_growth_rate[i-3,1] <- gr$coefficients[2]
        temp_max_growth_rate[i-3,2] <- ifelse(is.nan(gr$adj.r.squared),0,gr$adj.r.squared)
      }
      
      # have > 0.9 return all > 0.9, else return the max R square one
      if (max(temp_max_growth_rate[,2])>0.9) {
        temp_max_growth_rate <- temp_max_growth_rate[temp_max_growth_rate[,2]>0.9,] 
        } else {
          temp_max_growth_rate <- temp_max_growth_rate[temp_max_growth_rate[,2]==max(temp_max_growth_rate[,2]),]
        }
          
    
      
      temp_max_growth_rate <- temp_max_growth_rate[temp_max_growth_rate[,1]==max(temp_max_growth_rate[,1]),]
      
      return(temp_max_growth_rate)
    }
    
   
    
    
    
    if ( is.data.frame(data)) {
      three_point_growth_rate <-  sapply(data, max_growth_rate,simplify = TRUE, USE.NAMES = F)
      three_point_growth_rate <- t(three_point_growth_rate)
    } else {
      three_point_growth_rate <-  data.frame("sw_Growth_rate"=c(0,0), "sw_R_square"=c(0,0))
    }
    
    
    
    return(three_point_growth_rate)
  }

  three_point_growth_rate_result <- three_point_growth_rate(all_valid_area_filter)
  
  
  # 3. lm fit for all point >0 to get the growth rate and R2  

  growth_rate <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(growth_rate) <- c("Growth_rate","R_square")
  
  if (is.data.frame(all_valid_area_filter)) {
    for (i in 1:ncol(all_valid_area_filter)) {
      index <- all_valid_area_filter[,i] > 0
      gr <- lm(log2(all_valid_area_filter[,i][index]) ~ time_point[index])
      gr <- summary(gr)
      growth_rate[i,1] <- gr$coefficients[2]
      growth_rate[i,2] <- gr$adj.r.squared
    }
  } else {
    growth_rate <-  data.frame("Growth_rate"=c(0,0), "R_square"=c(0,0))
  }
  
  
  # tidy the result table
  colname_temp <- colnames(three_point_growth_rate_result)
  rowname_temp <- rownames(three_point_growth_rate_result)
  three_point_growth_rate_result <- as.data.frame(matrix(unlist(three_point_growth_rate_result),
                                                         nrow = nrow(three_point_growth_rate_result),
                                                         byrow = F))
  colnames(three_point_growth_rate_result) <- colname_temp
  rownames(three_point_growth_rate_result) <- rowname_temp
  
  growth_rate_table <- as.data.frame(cbind(three_point_growth_rate_result,growth_rate))
  
  growth_rate_table$max_growth_rate <- 0 # not valid
  
  for (i in 1:nrow(growth_rate_table)) {
    if (growth_rate_table[i,2] >=0.9 & growth_rate_table[i,4] >=0.9) {
      growth_rate_table$max_growth_rate[i] = max(growth_rate_table[i,1], growth_rate_table[i,3])
    }
    
    if (growth_rate_table[i,2] >=0.9 & growth_rate_table[i,4] <0.9) {
      growth_rate_table$max_growth_rate[i] = growth_rate_table[i,1]
    }
    
    if (growth_rate_table[i,2] <0.9 & growth_rate_table[i,4] >=0.9) {
      growth_rate_table$max_growth_rate[i] = growth_rate_table[i,3]
    }
  }
  
  # filter all point fit R>0
  all_valid_area_filter <- all_valid_area_filter[,growth_rate_table$Growth_rate > 0]
  growth_rate_table <- growth_rate_table[growth_rate_table$Growth_rate > 0,]
  
  # filter R>0.9 and growth rate>0
  all_valid_area_filter <- all_valid_area_filter[,growth_rate_table$max_growth_rate > 0]
  growth_rate_table <- growth_rate_table[growth_rate_table$max_growth_rate > 0,]
  
  all_result <- list(growth_rate_table,all_valid_area_filter)
  return(all_result)
  
}


#--this is for one well--------------------------
all_valid_area <- growth_rate_file_list[[85]]

one_well <- process_one_well(all_valid_area)

one_well[[1]]$max_growth_rate <- tt[[1]]$max_growth_rate/time_interval
one_well[[1]]$sw_Growth_rate <- tt[[1]]$sw_Growth_rate/time_interval
one_well[[1]]$Growth_rate <- tt[[1]]$Growth_rate/time_interval


one_well[[1]]
one_well[[2]]

# library(openxlsx)
# write.xlsx(one_well[[1]],"./20191011_20191003_TL25_YFP_H01_2/H01_calculated_grwoth_rate.xlsx",col.names=T,row.names=T)
# write.xlsx(one_well[[2]],"./20191011_20191003_TL25_YFP_H01_2/H01_raw_area_size.xlsx",col.names=T)



# 2. apply the function to all wells to get result------------------------------

#debugonce(process_one_well) # this is used to debug function

growth_rate_file_list_result <- lapply(growth_rate_file_list, process_one_well)

growth_rate_file_list_result <- growth_rate_file_list_result[c(plot_order)]




# 3. plot all wells growth rate distribution----------------------------------------

time_interval <- 100/60  # 100 mins
  
oldpar <- par()
par(mfcol=c(6,6))
par(mar=c(1,1,1,1))
for (i in 1:length(growth_rate_file_list_result)) {
  
  tryCatch(
    { p <- growth_rate_file_list_result[[i]][[1]]$max_growth_rate/time_interval
    plot(density(p),
         main = "",col="#3182bd",lwd=2)
    abline(v=median(p),lty=2,lwd=2,col="#FF007F")
    text(x=median(p),y=max(density(p)$y)/2,
         paste0("med_GR=","\n",round(median(p),2),"\n","cluster_num=","\n",length(p)),cex=2,face="bold") },
    #warning = function(w) { warning-handler-code },
    error = function(e) { message('Error @ ',x) ; plot(0,0) }
    
    #finally = { cleanup-code }
  )
}

par(oldpar)


# 4. plot all wells growth curve---------------------------------------

library(reshape2)
library(RColorBrewer)
display.brewer.all()
color_use <- brewer.pal(12,"Paired")
color_use[11] <- "#dd1c77"


time_interval <- 100/60

oldpar <- par()
par(mfcol=c(6,6))
par(mar=c(2,2,2,2))

for (i in 1:length(growth_rate_file_list_result)) {
  
  
  
  
  if (is.na(growth_rate_file_list_result[[i]])) { 
    plot(1,1,col="white",main = i)
    next
  } else {
    p <- growth_rate_file_list_result[[i]][[2]]
  }
  
  if (is.vector(p)) {
    plot(x=c(1:length(p)),y=log2(p)-log2(p[1]),
         xlim = c(0,36),ylim = c(0,10),
         #xlab = "Time(hour)", ylab="log2(area increase)",
         type = "o",pch=21,cex=1.5,lwd=2,col="black",main=i
         #,bg="black"
    )
  } else if ( is.data.frame(p) & ncol(p)>1 ) {
    #p <- sweep(p,2,FUN="-",STATS = colmaxs(p))
    
    p$Time <- (1:nrow(p))*time_interval
    
    
    # need to smooth?
    y_value <- log2(p[,1][p[,1]>0])-log2(p[,1][1])
    
    
    plot(x=p$Time[p[,1]>0],y=y_value,
         xlim = c(0,36),ylim = c(0,10),
         #xlab = "Time(hour)", ylab="log2(area increase)",
         type = "o",pch=21,cex=1.5,lwd=2,col="black",main=i
         #,bg="black"
    )
    
    for (j in 2:(ncol(p)-1)) {
      lines(x=p$Time[p[,j]>0],y=log2(p[,j][p[,j]>0])-log2(p[,j][1]),type = "o",pch=21,cex=1.5,lwd=2,
            #bg=color_use[j%%12+1],
            col=color_use[j%%12+1]
      )
    }
    
    text(x=20,y=9.5,
         paste0("cluster num=",ncol(p)-1),cex=1.5)
    
    #legend("topleft",paste0("cluster number=",ncol(p)-1,cex=1))
    
  } else {
    plot(x=1,y=1,col="white", main=i)
  }

}
  


par(oldpar)




# 5. plot all wells growth curve to pdf----------------------------------------

grDevices::cairo_pdf("./20191009_20191003_TL25_YFP_all/growthcurver.pdf", width = 12*2, height = 8*2)

old_par <- graphics::par(mfcol = c(8, 12), mar = c(0.5,0.5,0.5,0.5))

time_interval <- 100/60

for (i in 1:length(growth_rate_file_list_result)) {
  
  p <- growth_rate_file_list_result[[i]][[2]]
  
  
  #p <- sweep(p,2,FUN="-",STATS = colmaxs(p))
  
  p$Time <- (1:nrow(p))*time_interval
  
  
  # need to smooth?
  y_value <- log2(p[,1][p[,1]>0])-log2(p[,1][1])
  
  
  plot(x=p$Time[p[,1]>0],y=y_value,
       xlim = c(0,36),ylim = c(0,10),
       #xlab = "Time(hour)", ylab="log2(area increase)",
       type = "o",pch=21,cex=1.5,lwd=2,col="black"
       #,bg="black"
  )
  
  for (j in 2:(ncol(p)-1)) {
    lines(x=p$Time[p[,j]>0],y=log2(p[,j][p[,j]>0])-log2(p[,j][1]),type = "o",pch=21,cex=1.5,lwd=2,
          #bg=color_use[j%%12+1],
          col=color_use[j%%12+1]
    )
  }
  
  text(x=20,y=9.5,
       paste0("cluster num=",ncol(p)-1),cex=1.5,face="bold")
  
  #legend("topleft",paste0("cluster number=",ncol(p)-1,cex=1))
}

grDevices::dev.off()




#######################################
# Process mCherry data###############################################################################

# read in mCherry data------------------------------------------------------------

mCherry_filenames <- list.files("Z:/MD/Analysis_data_Zhisheng/31th_20X//analysis/RFP_data/", pattern=".*.csv", full.names=TRUE)

mCherry_filenames_well_label <- gsub("_all_valid_area_RFP.*","",gsub(".*RFP_data/","",mCherry_filenames))

plot_order <- order(substr(mCherry_filenames_well_label,2,3))

mCherry_file_list <- lapply(mCherry_filenames,read.csv)

mCherry_file_list <- mCherry_file_list[c(plot_order)]

mCherry_file_list[[1]][1:5,1:10]


# 1. get the filtered cell cluster according to the growth rate data------------

filter_FI_data_with_growth_rate_data <- function(FI_data,growth_rate_data) {
  
  for (i in 1:length(growth_rate_data)) {
    FI_data[[i]] <- FI_data[[i]][,colnames(FI_data[[i]]) %in% rownames(growth_rate_data[[i]][[1]])]
    
    #FI_data[[i]] <- FI_data[[i]][order(c(1,10:19,2,20,3:9)) ,]  # this is because of the order of RFP is not correct
  }
  return(FI_data)
}


filtered_mCherry_data <- filter_FI_data_with_growth_rate_data(mCherry_file_list,growth_rate_file_list_result)


# 2. plot all wells FI vs. time---------------------------------------

library(reshape2)
library(RColorBrewer)
display.brewer.all()
color_use <- brewer.pal(12,"Paired")
color_use[11] <- "#dd1c77"


time_interval <- 100/60

oldpar <- par()
par(mfcol=c(8,12))
par(mar=c(1,1,1,1))

for (i in 1:length(filtered_mCherry_data)) {
  
  p <- filtered_mCherry_data[[i]]
  
  
  #p <- sweep(p,2,FUN="-",STATS = colmaxs(p))
  
  p$Time <- (1:nrow(p))*time_interval
  
  
  # need to smooth?
  
  
  
  plot(x=p$Time[p[,1]>0],y=log2(p[,1][p[,1]>0]+1),
       xlim = c(0,36),ylim = c(6,10),
       #xlab = "Time(hour)", ylab="log2(area increase)",
       type = "o",pch=21,cex=1,lwd=2,col="black"
       #,bg="black"
  )
  
  for (j in 2:(ncol(p)-1)) {
    lines(x=p$Time[p[,j]>0],y=log2(p[,j][p[,j]>0]+1),type = "o",pch=21,cex=1,lwd=2,
          #bg=color_use[j%%12+1],
          col=color_use[j%%12+1]
    )
  }
  
  text(x=20,y=4.5,
       paste0("cluster num=",ncol(p)-1),cex=1.5
       #,font=list(family="sans", face=2)
  )      
  # font: The font "face" (1=plain, 2=bold, 3=italic, 4=bold-italic)
  # font=list(family=string, face=numberorstring)
  
  #legend("topleft",paste0("cluster number=",ncol(p)-1,cex=1))
}

par(oldpar)


# 3. plot all wells area vs. FI---------------------------------------

library(reshape2)
library(RColorBrewer)
display.brewer.all()
color_use <- brewer.pal(12,"Paired")
color_use[11] <- "#dd1c77"


time_interval <- 100/60

oldpar <- par()
par(mfcol=c(8,12))
par(mar=c(1,1,1,1))

for (i in 1:length(filtered_mCherry_data)) {

  p <- filtered_mCherry_data[[i]]
  
  area <- growth_rate_file_list_result[[i]][[2]]
  #p <- sweep(p,2,FUN="-",STATS = colmaxs(p))
  
  #p$Time <- (1:nrow(p))*time_interval
  
  
  # need to smooth?
  # n=2
  # plot(x=log2(area[,n][area[,n]>0]+1),y=log2(p[,n][p[,n]>0]+1),
  #      xlim = c(7,13),ylim = c(7,9),
  #      main =n,xlab = "log2(area size)", ylab="log2(mCherry signal)",
  #      type = "o",pch=21,cex=1,lwd=2,col="black",cex.lab=1.5
  #      #,bg="black"
  # )
  # lines(x=log2(area[,n][area[,n]>0]+1),y=log2(p[,n][p[,n]>0]+1),type = "o",pch=21,cex=1,lwd=2,col="red")
  # 
  
  plot(x=log2(area[,1][area[,1]>0]+1),y=log2(p[,1][p[,1]>0]+1),
       xlim = c(6,16),ylim = c(6,10),
       #xlab = "Time(hour)", ylab="log2(area increase)",
       type = "o",pch=21,cex=1,lwd=2,col="black"
       #,bg="black"
  )
  
  for (j in 1:(ncol(p)-1)) {
    lines(x=log2(area[,j][area[,j]>0]+1),y=log2(p[,j][p[,j]>0]+1),type = "o",pch=21,cex=1,lwd=2,
          #bg=color_use[j%%12+1],
          col=color_use[j%%12+1]
          #col="grey"
          #col=alpha("red",0.4)
    )
    #if (j>=12) {break}
  }
  
  text(x=15,y=9,
       paste0("cluster num=",ncol(p)-1),cex=1.5
       #,font=list(family="sans", face=2)
  )      
  # font: The font "face" (1=plain, 2=bold, 3=italic, 4=bold-italic)
  # font=list(family=string, face=numberorstring)
  
  #legend("topleft",paste0("cluster number=",ncol(p)-1,cex=1))
}

par(oldpar)


# 4. plot single wells comparison of area vs. FI---------------------------------------

display.brewer.all()
color_use <- brewer.pal(12,"Paired")
color_use[11] <- "#dd1c77"

time_interval <- 100/60


# single cell vs. single cell----------------------

n=50 # cluster number

i=8  # well


p <- filtered_mCherry_data[[i]]

area <- growth_rate_file_list_result[[i]][[2]]


plot(x=log2(area[,n][area[,n]>0]+1),y=log2(p[,n][p[,n]>0]+1),
     xlim = c(7,13),ylim = c(7,9),
     main =n,xlab = "log2(area size)", ylab="log2(mCherry signal)",
     type = "o",pch=21,cex=1,lwd=2,col="black",cex.lab=1.5
     #,bg="black"
)


i=1 # well

p <- filtered_mCherry_data[[i]]

area <- growth_rate_file_list_result[[i]][[2]]


lines(x=log2(area[,n][area[,n]>0]+1),y=log2(p[,n][p[,n]>0]+1),type = "o",pch=21,cex=1,lwd=2,col="red")
legend("topleft",col = c("black","red"),legend = c("H01","A01"),lty = 1,lwd=2)



# single well vs. single well---------------

i=4

p <- filtered_mCherry_data[[i]]

area <- growth_rate_file_list_result[[i]][[2]]

plot(x=log2(area[,1][area[,1]>0]+1),y=log2(p[,1][p[,1]>0]+1),
     xlim = c(6,16),ylim = c(6,10),
     #xlab = "Time(hour)", ylab="log2(area increase)",
     type = "o",pch=21,cex=1,lwd=2,col="grey"
     #,bg="black"
)

for (j in 1:(ncol(p)-1)) {
  lines(x=log2(area[,j][area[,j]>0]+1),y=log2(p[,j][p[,j]>0]+1),type = "o",pch=21,cex=1,lwd=2,
        #bg=color_use[j%%12+1],
        #col=color_use[j%%12+1]
        #col="grey"
        col=alpha("grey",0.4)
  )
  #if (j>=12) {break}
}


i=8

p <- filtered_mCherry_data[[i]]

area <- growth_rate_file_list_result[[i]][[2]]

for (j in 1:(ncol(p)-1)) {
  lines(x=log2(area[,j][area[,j]>0]+1),y=log2(p[,j][p[,j]>0]+1),type = "o",pch=21,cex=1,lwd=2,
        #bg=color_use[j%%12+1],
        #col=color_use[j%%12+1]
        #col="grey"
        col=alpha("red",0.4)
  )
  #if (j>=12) {break}
}


# single well to 96 line-----------
oldpar <- par()
par(mfcol=c(8,12))
par(mar=c(1,1,1,1))

i=1

p <- filtered_mCherry_data[[i]]

area <- growth_rate_file_list_result[[i]][[2]]


for (j in 1:max(96,(ncol(p)-1))) {
#for (j in 50:70) {  
  plot(x=log2(area[,j][area[,j]>0]+1),y=log2(p[,j][p[,j]>0]+1),type = "o",pch=21,cex=1,lwd=2,
       xlim = c(6,16),ylim = c(6,10),
        #bg=color_use[j%%12+1],
        #col=color_use[j%%12+1]
        #col="grey"
        col=alpha("red",0.4)
  )
  abline(h=8,lty=2)
  #if (j>=12) {break}
}

par(oldpar)


# 5. get the average FI of mCherry across all well----------------------

average_FI_across_all_well <- function(filtered_FI_data,growth_rate_data) {
  
  for (i in 1:length(filtered_FI_data)) {
    filtered_FI_data[[i]] <- filtered_FI_data[[i]]/growth_rate_data[[i]][[2]]
  }
  return(filtered_FI_data)
}


average_mCherry_data <- average_FI_across_all_well(filtered_mCherry_data[well_num],growth_rate_file_list_result[well_num])

growth_rate_file_list_result[[1]][[2]][,1:5]
filtered_mCherry_data[[1]][,1:5]
average_mCherry_data[[1]][,1:5]

# 6. plot all wells average FI vs. time---------------------------------------

library(reshape2)
library(RColorBrewer)
display.brewer.all()
color_use <- brewer.pal(12,"Paired")
color_use[11] <- "#dd1c77"


time_interval <- 100/60

oldpar <- par()
par(mfcol=c(8,12))
par(mar=c(1,1,1,1))

for (i in 1:length(average_mCherry_data)) {

  p <- average_mCherry_data[[i]]
  
  
  #p <- sweep(p,2,FUN="-",STATS = colmaxs(p))
  
  p$Time <- (1:nrow(p))*time_interval
  
  
  # need to smooth?

  
  
  plot(x=p$Time[p[,1]>0],y=log2(p[,1][p[,1]>0]+1),
       xlim = c(0,36),ylim = c(0,5),
       #xlab = "Time(hour)", ylab="log2(area increase)",
       type = "o",pch=21,cex=1,lwd=2,col="black"
       #,bg="black"
  )
  
  for (j in 2:(ncol(p)-1)) {
    lines(x=p$Time[p[,j]>0],y=log2(p[,j][p[,j]>0]+1),type = "o",pch=21,cex=1,lwd=2,
          #bg=color_use[j%%12+1],
          col=color_use[j%%12+1]
    )
  }
  
  text(x=20,y=4.5,
       paste0("cluster num=",ncol(p)-1),cex=1.5
       #,font=list(family="sans", face=2)
       )      
  # font: The font "face" (1=plain, 2=bold, 3=italic, 4=bold-italic)
  # font=list(family=string, face=numberorstring)
  
  #legend("topleft",paste0("cluster number=",ncol(p)-1,cex=1))
}

par(oldpar)


#to pdf---------------------------
grDevices::cairo_pdf("./20191007_20190929_TL25_mCherry_all/average_mCherry_curve.pdf", width = 12*2, height = 8*2)

old_par <- graphics::par(mfcol = c(8, 12), mar = c(0.5,0.5,0.5,0.5))

time_interval <- 100/60

for (i in 1:length(average_mCherry_data)) {
  
  p <- average_mCherry_data[[i]]
  
  
  #p <- sweep(p,2,FUN="-",STATS = colmaxs(p))
  
  p$Time <- (1:nrow(p))*time_interval
  
  
  # need to smooth?
  
  
  
  plot(x=p$Time[p[,1]>0],y=log2(p[,1][p[,1]>0]+1),
       xlim = c(0,36),ylim = c(0,5),
       #xlab = "Time(hour)", ylab="log2(area increase)",
       type = "o",pch=21,cex=1,lwd=2,col="black"
       #,bg="black"
  )
  
  for (j in 2:(ncol(p)-1)) {
    lines(x=p$Time[p[,j]>0],y=log2(p[,j][p[,j]>0]+1),type = "o",pch=21,cex=1,lwd=2,
          #bg=color_use[j%%12+1],
          col=color_use[j%%12+1]
    )
  }
  
  text(x=20,y=4.5,
       paste0("cluster num=",ncol(p)-1),cex=1.5
       #,font=list(family="sans", face=2)
  )      
  # font: The font "face" (1=plain, 2=bold, 3=italic, 4=bold-italic)
  # font=list(family=string, face=numberorstring)
  
  #legend("topleft",paste0("cluster number=",ncol(p)-1,cex=1))
}

grDevices::dev.off()




# 7. plot all wells average(average FI) vs. growth rate---------------------------------------

# filter first two time points, as the cell cluster area is too small, it can lead to bias. And use the --mean-- of next three time points to get the mean(mCherry/area) vs. time 




library(reshape2)
library(RColorBrewer)
display.brewer.all()
color_use <- brewer.pal(12,"Paired")
color_use[11] <- "#dd1c77"


# is.nan doesn't actually have a method for data frames, unlike is.na. So, let's fix that!
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

# get the xlim and ylim---------------------
time_interval <- 100/60

max_gr <- function(x) {
  x <- max(x[[1]]$Growth_rate)
  return(x)
}

min_gr <- function(x) {
  x <- min(x[[1]]$Growth_rate)
  return(x)
}

x_max <- max(unlist(lapply(growth_rate_file_list_result, max_gr)))/time_interval

x_min <- min(unlist(lapply(growth_rate_file_list_result, min_gr)))/time_interval

max_FI <- function(y) {
  y <- max(colMeans(y[4:8,]))
  return(y)
}

y_max <- max(unlist(lapply(average_mCherry_data, max_FI)))

y_min = 0


# plot 1   with xlim and ylim----------------------------------
oldpar <- par()
par(mfcol=c(8,12))
par(mar=c(2,2,2,2))

for (i in 1:length(average_mCherry_data)) {
  

  p <- average_mCherry_data[[i]][4:8,]
  
  
  #p <- sweep(p,2,FUN="-",STATS = colmaxs(p))
  
  p <- colMeans(p)
  
  
  # need to smooth?
  
  
  
  plot(x=growth_rate_file_list_result[[i]][[1]]$Growth_rate/time_interval,y=p,
       xlim = c(x_min,x_max),ylim = c(y_min,y_max),
       #xlab = "Time(hour)", ylab="log2(area increase)",
       type = "p",pch=21,cex=1,lwd=2,col="black",cex.axis=1.5
       #,bg="black"
  )
  
  # for (j in 2:(ncol(p)-1)) {
  #   lines(x=p$Time[p[,j]>0],y=log2(p[,j][p[,j]>0]),type = "o",pch=21,cex=1,lwd=2,
  #         #bg=color_use[j%%12+1],
  #         col=color_use[j%%12+1]
  #   )
  # }
  
  # text(x=20,y=4.5,
  #      paste0("cluster num=",ncol(p)-1),cex=1.5
  #      #,font=list(family="sans", face=2)
  # )      
  # font: The font "face" (1=plain, 2=bold, 3=italic, 4=bold-italic)
  # font=list(family=string, face=numberorstring)
  
  #legend("topleft",paste0("cluster number=",ncol(p)-1,cex=1))
}

par(oldpar)

# plot 2   without xlim and ylim--------------------------------
oldpar <- par()
par(mfcol=c(8,12))
par(mar=c(2,2,2,2))

for (i in 1:length(average_mCherry_data)) {
  
  
  p <- average_mCherry_data[[i]][4:8,]
  
  
  #p <- sweep(p,2,FUN="-",STATS = colmaxs(p))
  
  p <- colMeans(p)
  
  
  # need to smooth?
  
  
  x_value = growth_rate_file_list_result[[i]][[1]]$Growth_rate/time_interval
  plot(x=x_value,y=p,
       #xlim = c(x_min,x_max),ylim = c(y_min,y_max),
       #xlab = "Time(hour)", ylab="log2(area increase)",
       type = "p",pch=21,cex=1,lwd=2,col="black",cex.axis=1.5
       #,bg="black"
  )
  abline(lm(p~x_value),col="red")
  text(x=max(x_value)*0.8,y=max(p)*0.95, paste0("R=",round(cor(x=x_value,y=p,method = "pearson"),2)),cex = 1.5)
  #text(x=,y=,)
  # for (j in 2:(ncol(p)-1)) {
  #   lines(x=p$Time[p[,j]>0],y=log2(p[,j][p[,j]>0]),type = "o",pch=21,cex=1,lwd=2,
  #         #bg=color_use[j%%12+1],
  #         col=color_use[j%%12+1]
  #   )
  # }
  
  # text(x=20,y=4.5,
  #      paste0("cluster num=",ncol(p)-1),cex=1.5
  #      #,font=list(family="sans", face=2)
  # )      
  # font: The font "face" (1=plain, 2=bold, 3=italic, 4=bold-italic)
  # font=list(family=string, face=numberorstring)
  
  #legend("topleft",paste0("cluster number=",ncol(p)-1,cex=1))
}

par(oldpar)








# 8. this is to check a single well growth rate--------------------------
library(ggplot2)
library(reshape2)

data <- growth_rate_file_list_result[[12]][[2]]
data$time <- c(1:20)*time_interval
melt_data <- melt(data,id.vars = "time")

melt_data$value <- log2(melt_data$value)

melt_data$value[melt_data$value<=0] <- NA

melt_data$color <- "grey"



i_data <- c("cell_cluster278","cell_cluster309","cell_cluster317")

melt_data$color[melt_data$variable %in% i_data] <- c(rep("cell_cluster278",20),
                                                     rep("cell_cluster309",20),
                                                     rep("cell_cluster317",20)
                                                     )

ggplot(data=melt_data,aes(x=time,y=value,group=variable)) +
  geom_point(color="grey") +
  geom_line(color="grey") +
  
  geom_line(data=subset(melt_data,variable %in% i_data),aes(color=color),size=2) +
  labs(y="log2(area)") +
  theme_bw()



# 9.Classification of survive one and non-survive one--------------------

# 9.1 hole palte, use colmean of 4 to 8, D 32ug/ml--------------
library(matrixStats)
oldpar <- par()
par(mfcol=c(8,12))
par(mar=c(1,1,1,1))



for (i in 1:96) {
  

  p <- filtered_mCherry_data[[i]]
  
  if (ncol(p)<1) {
    plot(-1,-1,main=i)
  } else {
    #area <- growth_rate_file_list_result[[i]][[2]]
    
    p <- log2(p[4:8,]+1)
    #area <- log2(area[1:7,]+1)
    
    # slope
    # slope <- c()
    # 
    # for (k in 1:ncol(p)) {
    #   slope[k] <- lm(p[,k]~area[,k])$coefficients[2]
    # }
    
    # colmax
    colmaxs <- function(x) {
      x <- apply(x, 2, max)
    }
    
    max_p <- colmaxs(p)
    
    
    mu <- mean(max_p)
    sigma <- sqrt(sum((max_p-mu)^2)/length(max_p))
    
    down_bound <- qnorm(0.995, mean = mu, sd = sigma, lower.tail = F, log.p = FALSE)
    # range
    
    # colmins <- function(x) {
    #   x <- apply(x, 2, min)
    # }
    # 
    # range_p <- colmaxs(p) - colmins(p)
    # 
    # slope <- as.numeric(slope)
    # 
    # 
    # final_value <- slope/mean(slope)+max_p/mean(max_p)+range_p/mean(range_p)
    # 
    plot(sort(max_p), main = i)
    abline(h=down_bound, lty=2)
  }
 
}


par(oldpar)

# 9.2 use colmean of 4 to 8, all well-------------

oldpar <- par()
par(mfcol=c(2,9))
par(mar=c(1,1,1,1))

well_num <- sort(c(seq(1,41,5),seq(2,42,5)))

for (i in well_num) {
  
  
  p <- filtered_mCherry_data[[i]]
  
  if (ncol(p)<1) {
    plot(-1,-1,main=i)
  } else {
    #area <- growth_rate_file_list_result[[i]][[2]]
    
    # p <- log2(p[4:8,]+1)
    #area <- log2(area[1:7,]+1)
    
    # slope
    # slope <- c()
    # 
    # for (k in 1:ncol(p)) {
    #   slope[k] <- lm(p[,k]~area[,k])$coefficients[2]
    # }
    
    # colmax
    # colmaxs <- function(x) {
    #   x <- apply(x, 2, max)
    # }
    # 

    max_p <- colMeans(p[4:8,])

    
    mu <- mean(max_p)
    sigma <- sqrt(sum((max_p-mu)^2)/length(max_p))
    
    down_bound <- qnorm(0.995, mean = mu, sd = sigma, lower.tail = F, log.p = FALSE)
    # range
    
    # colmins <- function(x) {
    #   x <- apply(x, 2, min)
    # }
    # 
    # range_p <- colmaxs(p) - colmins(p)
    # 
    # slope <- as.numeric(slope)
    # 
    # 
    # final_value <- slope/mean(slope)+max_p/mean(max_p)+range_p/mean(range_p)
    # 
    plot(sort(max_p), main = i)
    abline(h=220, lty=2)
  }
  
}


par(oldpar)


# 9.3 use colmean of 4 to 8, D 32ug/ml--------------

oldpar <- par()
par(mfcol=c(2,9))
par(mar=c(2,2,2,2))

well_num <- sort(c(seq(1,41,5),seq(2,42,5)))



n=0
for (i in well_num) {
  
  n=n+1
  p <- filtered_mCherry_data[[i]]
  
  if (ncol(p)<1) {
    plot(-1,-1,main=i)
  } else {
    #area <- growth_rate_file_list_result[[i]][[2]]
    
    # p <- log2(p[4:8,]+1)
    #area <- log2(area[1:7,]+1)
    
    # slope
    # slope <- c()
    # 
    # for (k in 1:ncol(p)) {
    #   slope[k] <- lm(p[,k]~area[,k])$coefficients[2]
    # }
    
    # colmax
    # colmaxs <- function(x) {
    #   x <- apply(x, 2, max)
    # }
    # 
    
    max_p <- colMeans(p[4:6,])
    
    print(max(max_p))
    
    mu <- mean(max_p)
    sigma <- sqrt(sum((max_p-mu)^2)/length(max_p))
    
    down_bound <- qnorm(0.995, mean = mu, sd = sigma, lower.tail = F, log.p = FALSE)
    # range
    
    # colmins <- function(x) {
    #   x <- apply(x, 2, min)
    # }
    # 
    # range_p <- colmaxs(p) - colmins(p)
    # 
    # slope <- as.numeric(slope)
    # 
    # 
    # final_value <- slope/mean(slope)+max_p/mean(max_p)+range_p/mean(range_p)
    # 
    plot(sort(max_p), main = mCherry_filenames_well_label[plot_order][well_num]
[n])
    # abline(h=down_bound, lty=2)
    abline(h=220, lty=2)
  }
  
}


par(oldpar)


# 10. Grwoth rate vs. mCherry signal----------
  


# detect outliers----------

# library(extremevalues)
# getOutliersII(max_p, alpha=c(0.01, 0.01), FLim=c(0, 1),  distribution="normal", returnResiduals=TRUE)

colmaxs <- function(x) {
  x <- apply(x, 2, max)
}


# 10. After decide the time interval used to predict, we need to calculate the down_bound--------

# 10.1 time points and survive_in_data----
SlidingWindow <- function(data,step_length) {
  length_list = length(data)-step_length+1
  result <- vector(mode = "list", length=length_list)
  for (i in 1:length_list) {
    result[[i]] <- data[i:(i+step_length-1)]
  }
  return(result)
}


data <- 1:6

SlidingWindow(data,6)

all_try_time_points <- c(SlidingWindow(data,1),SlidingWindow(data,2),
                            SlidingWindow(data,3),SlidingWindow(data,4),
                            SlidingWindow(data,5),SlidingWindow(data,6)
                            )



# read in survive cell

library(openxlsx)
survive <- read.xlsx("Z:/MD/Analysis_data_Zhisheng/15th_20191117_for_prediction_10X/analysis/20191121mannully chosed resistant colony.XLSX")

survive <- survive[,order(substr(colnames(survive),2,3))]

survive_in_data <- vector(mode = "list", length = ncol(survive))

well_num <- sort(c(seq(1,41,5),seq(2,42,5)))

for (i in 1:ncol(survive)) {
  survive_in_data[[i]] <- survive[,i][survive[,i] %in% gsub("cell_cluster","",rownames(growth_rate_file_list_result[[well_num[i]]][[1]]))]
  
}


# 10.2 check all time points------------------------------------

F1_score_table <- data.frame(matrix(ncol = 4,nrow = 0))
colnames(F1_score_table) <- c("Precision","Recall","F1_score","time_point_used")

dir_to_save <- "Z:/MD/Analysis_data_Zhisheng/15th_20191117_for_prediction_10X/analysis/"
dir_to_save <- paste0(dir_to_save,"absolute_mCherry2/")


# average_mCherry_data or filtered_mCherry_data[well_num], # average_mCherry_data does not work----
FI_data_used_here <- filtered_mCherry_data[well_num]    

for (k in 1:length(all_try_time_points)) {
  
  
  
  
  time_period <- unlist(all_try_time_points[k])
  
  #data <- filtered_YFP_data[well_num]
  Get_colmean_of_choosen_time_period <- function(data,time_period) {
    
    Colmean_of_choosen_time_period <- function(data) {
      data <- colMeans(data[time_period,])
    }
    
    # Colmean_of_choosen_time_period function has only one input
    data <- lapply(data, Colmean_of_choosen_time_period)
    
    data <- unlist(data)
    
    return(data)
  }
  
  
  
  filtered_mCherry_data_all <- Get_colmean_of_choosen_time_period(FI_data_used_here,time_period)
  
  length(filtered_mCherry_data_all)
  
  
  
  title_plot <- paste("Time point",min(time_period),"to",max(time_period))
  
  png(paste0(dir_to_save,"/",title_plot," density.png"), width = 100, height = 100, units="mm",res=300)
  plot(density(filtered_mCherry_data_all),main = title_plot)
  
  dev.off()
  
  png(paste0(dir_to_save,"/",title_plot," hist.png"), width = 100, height = 100, units="mm",res=300)
  hist(filtered_mCherry_data_all,breaks =1000,main = title_plot)
  dev.off()
  
  #density_data <- density(filtered_mCherry_data_all,bw =7)
  
  
  # fit a binomal distribution model to find two main peak--------------
  
  library(flexmix)
  
  kde <- density(filtered_mCherry_data_all)
  m1 <- FLXMRglm(family = "gaussian")
  m2 <- FLXMRglm(family = "gaussian")
  fit <- flexmix(filtered_mCherry_data_all ~ 1, data = as.data.frame(filtered_mCherry_data_all), k = 2, model = list(m1, m2))
  c1 <- parameters(fit, component=1)[[1]]
  c2 <- parameters(fit, component=2)[[1]]
  
  
  png(paste0(dir_to_save,"/",title_plot," binormal.png"), width = 100, height = 100, units="mm",res=300)
  plot(kde)
  abline(v=c1[[1]], lty=2, col='blue')
  abline(v=c2[[1]], lty=2, col='red')
  dev.off()
  
  c1_peak <- min(c1[[1]],c2[[1]])
  c2_peak <- max(c1[[1]],c2[[1]])
  # density_data$x[which(diff(sign(diff(density_data$y)))==-2)+1]
  # density_data$y[which(diff(sign(diff(density_data$y)))==-2)+1]
  # 
  # minimum_point <- optimize(approxfun(density_data$x,density_data$y),interval=c(300,500))$minimum
  # 
  # density_data$x[which.min(abs(diff(density_data$y)))]
  
  
  # find a point in the range and give smallest difference (L1)----------------
  filtered_mCherry_data_all <- sort(filtered_mCherry_data_all)
  # how many points in this range
  length_range <- length(filtered_mCherry_data_all[filtered_mCherry_data_all>c1_peak & filtered_mCherry_data_all<c2_peak])
  
  # whole range, 1/5 of length_range to calculate distance
  
  c1_point <- which((filtered_mCherry_data_all>c1_peak)==1)[1]
  c2_point <- which((filtered_mCherry_data_all>c2_peak)==1)[1]-1
  
  # xrange
  c1_point_x <- c1_point+trunc(length_range/10)  # trunc ceiling and round etc.
  c2_point_x <- c2_point-trunc(length_range/10)
  
  #index <- which.max(as.numeric(filtered_mCherry_data_all > minimum_point))
  
  filtered_mCherry_data_all_limit <- filtered_mCherry_data_all[c1_point:c2_point]
  loss_function <- function(x,data=filtered_mCherry_data_all_limit) {
    x <- which.min(abs(x-data))
    data <- filtered_mCherry_data_all_limit[(x-trunc(length_range/10)-1):(x+trunc(length_range/10)-1)]
    loss <- sum((x-data)^2)
    return(loss)
  }
  
  minimum_point <- optim(par=filtered_mCherry_data_all[c1_point_x], fn=loss_function, gr=NULL, 
                         data=filtered_mCherry_data_all_limit,method =  "Brent",
                         lower = filtered_mCherry_data_all[c1_point_x], 
                         upper = filtered_mCherry_data_all[c2_point_x])$par
  
  
  
  # optimise(loss_function, interval=c(min(filtered_mCherry_data_all),max(filtered_mCherry_data_all)),
  #          maximum = FALSE,
  #          tol = .Machine$double.eps^0.25)
  ###################################################
  
  
  get_outlier_cell_mCherry <- function(p) {
    
    
    max_p <- colMeans(p[time_period,])
    
    down_bound <- minimum_point
    
    result <- max_p[max_p<down_bound]
    
    return(result)
    # plot(density(max_p))
    # 
    # qnorm(0.995, mean = mu, sd = sigma, lower.tail = F, log.p = FALSE)
  }
  
  all_outlier_cell <- lapply(FI_data_used_here, get_outlier_cell_mCherry)
  
  
  
  
  # 11. check precision and recall----------
  
  summary_prediction <- data.frame(matrix(ncol = 9,nrow = 0))
  colnames(summary_prediction) <- c("A:predicted_survive","B:Real_survive","Cells_in_A_not_B","Cells_in_B_not_A","#Overlap_of_A_and_B","#Cells_in_A_not_B","#Cells_in_B_not_A","Precision","Recall")
  
  for (i in 1:length(survive_in_data)) {
    index1 <- gsub("cell_cluster","",names(all_outlier_cell[[i]])) %in% survive_in_data[[i]]
    index2 <- survive_in_data[[i]] %in% gsub("cell_cluster","",names(all_outlier_cell[[i]]))
    
    
    summary_prediction[i,1] <- length(all_outlier_cell[[i]])
    summary_prediction[i,2] <- length(sort(unique(survive_in_data[[i]])))
    summary_prediction[i,3] <- paste(gsub("cell_cluster","",names(all_outlier_cell[[i]]))[!index1],sep = ",",collapse = ",")
    summary_prediction[i,4] <- paste(unique(survive_in_data[[i]][!index2]), sep=",", collapse = ",")
    
    summary_prediction[i,5] <- sum(index1)
    summary_prediction[i,6] <- length(gsub("cell_cluster","",names(all_outlier_cell[[i]]))[!index1])
    summary_prediction[i,7] <- length(survive_in_data[[i]][!index2])
    
  }
  summary_prediction$Precision <- summary_prediction$`#Overlap_of_A_and_B`/(summary_prediction$`A:predicted_survive`)
  
  summary_prediction$Recall <- summary_prediction$`#Overlap_of_A_and_B`/summary_prediction$`B:Real_survive`
  
  
  rownames(summary_prediction) <- mCherry_filenames_well_label[plot_order][well_num]
  
  summary_prediction[19,] <- c(sum(summary_prediction$`A:predicted_survive`),sum(summary_prediction$`B:Real_survive`),"","",sum(summary_prediction$`#Overlap_of_A_and_B`),sum(summary_prediction$`#Cells_in_A_not_B`),sum(summary_prediction$`#Cells_in_B_not_A`),sum(summary_prediction$`#Overlap_of_A_and_B`)/sum(summary_prediction$`A:predicted_survive`),sum(summary_prediction$`#Overlap_of_A_and_B`)/sum(summary_prediction$`B:Real_survive`))
  
  
  
  rownames(summary_prediction)[19] <- "Total"
  
  #summary_prediction
  
  write.xlsx(summary_prediction,paste(dir_to_save,"20191123 prediction", paste("Time point",min(time_period),"to",max(time_period)),"h.xlsx"),row.names=T)
  
  Precision <- as.numeric(summary_prediction$Precision[nrow(summary_prediction)])
  Recall <- as.numeric(summary_prediction$Recall[nrow(summary_prediction)])
  F1_score <- (2*Precision*Recall)/(Precision+Recall)
  time_point_used <- ifelse(max(time_period)==min(time_period),
                            as.character(max(time_period)),
                            paste0(min(time_period),":",max(time_period)))
  F1_score_table[k,] <- c(Precision,Recall,F1_score,time_point_used)
  
}


# 10.3 plot F1_score_table figure-----------
library(ggplot2)
library(reshape2)
library(cowplot)
F1_score_table$time_point_used <- factor(F1_score_table$time_point_used, level=unique(F1_score_table$time_point_used))
F1_score_table_melt <- melt(F1_score_table,id.vars = "time_point_used")
F1_score_table_melt$value <- as.numeric(F1_score_table_melt$value)


ggplot(data = F1_score_table_melt, aes(x=time_point_used,y=value,group=variable)) +
  geom_point(aes(col=variable),size=5) +
  geom_line(aes(col=variable),size=3) +
  scale_y_continuous(breaks = seq(0.75,1, 0.05),limits = c(0.75,1)) +
  geom_hline(yintercept=0.95,linetype="dashed",color = "red", size=2) +
  labs(x="Time point used (mCherry data)")






growth_rate_file_list_result[[well_num[1]]][[2]]



filtered_mCherry_data[well_num[1]]

#########################################
# Process YFP  data###############################################################################

dir_to_save <- "Z:/MD/Analysis_data_Zhisheng/31th_20X//analysis/"

YFP_filenames <- list.files("Z:/MD/Analysis_data_Zhisheng/31th_20X//analysis/GFP_data/", pattern=".*.csv", full.names=TRUE)

YFP_filenames_well_label <- gsub("_all_valid_area_GFP.*","",gsub(".*GFP_data/","",YFP_filenames))

plot_order <- order(substr(YFP_filenames_well_label,2,3))

YFP_file_list <- lapply(YFP_filenames,read.csv)

YFP_file_list <- YFP_file_list[c(plot_order)]

YFP_file_list[[1]][1:5,1:10]


# 1. get the filtered cell cluster according to the growth rate data------------

filter_FI_data_with_growth_rate_data <- function(FI_data,growth_rate_data) {
  
  for (i in 1:length(growth_rate_data)) {
    FI_data[[i]] <- FI_data[[i]][,colnames(FI_data[[i]]) %in% rownames(growth_rate_data[[i]][[1]])]
    
    #FI_data[[i]] <- FI_data[[i]][order(c(1,10:19,2,20,3:9)) ,]  # this is because of the order of RFP is not correct
  }
  return(FI_data)
}


filtered_YFP_data <- filter_FI_data_with_growth_rate_data(YFP_file_list,growth_rate_file_list_result)

# 5. get the average FI of YFP across all well----------------------

average_FI_across_all_well <- function(filtered_FI_data,growth_rate_data) {
  
  for (i in 1:length(filtered_FI_data)) {
    filtered_FI_data[[i]] <- filtered_FI_data[[i]]/growth_rate_data[[i]][[2]]
  }
  return(filtered_FI_data)
}


average_YFP_data <- average_FI_across_all_well(filtered_YFP_data[well_num],growth_rate_file_list_result[well_num])

growth_rate_file_list_result[[1]][[2]][,1:5]
filtered_YFP_data[[1]][,1:5]
average_YFP_data[[1]][,1:5]


# 10. After decide the time interval used to predict, we need to calculate the down_bound--------

# 10.1 time points and survive_in_data----
SlidingWindow <- function(data,step_length) {
  length_list = length(data)-step_length+1
  result <- vector(mode = "list", length=length_list)
  for (i in 1:length_list) {
    result[[i]] <- data[i:(i+step_length-1)]
  }
  return(result)
}


data <- 1:6

SlidingWindow(data,6)

all_try_time_points <- c(SlidingWindow(data,1),SlidingWindow(data,2),
                         SlidingWindow(data,3),SlidingWindow(data,4),
                         SlidingWindow(data,5),SlidingWindow(data,6)
)



# read in survive cell

library(openxlsx)
survive <- read.xlsx("Z:/MD/Analysis_data_Zhisheng/15th_20191117_for_prediction_10X/analysis/20191121mannully chosed resistant colony.XLSX")

survive <- survive[,order(substr(colnames(survive),2,3))]

survive_in_data <- vector(mode = "list", length = ncol(survive))

well_num <- sort(c(seq(1,41,5),seq(2,42,5)))

for (i in 1:ncol(survive)) {
  survive_in_data[[i]] <- survive[,i][survive[,i] %in% gsub("cell_cluster","",rownames(growth_rate_file_list_result[[well_num[i]]][[1]]))]
  
}


# 10.2 check all time points------------------------------------

F1_score_table_YFP <- data.frame(matrix(ncol = 4,nrow = 0))
colnames(F1_score_table_YFP) <- c("Precision","Recall","F1_score","time_point_used")

dir_to_save <- "Z:/MD/Analysis_data_Zhisheng/15th_20191117_for_prediction_10X/analysis/"
dir_to_save <- paste0(dir_to_save,"absolute_YFP/")


# average_YFP_data or filtered_YFP_data[well_num], # average_YFP_data does not work----
FI_data_used_here <- filtered_YFP_data[well_num]    

for (k in 1:length(all_try_time_points)) {



  
  time_period <- unlist(all_try_time_points[k])
  
  #data <- filtered_YFP_data[well_num]
  Get_colmean_of_choosen_time_period <- function(data,time_period) {
    
    Colmean_of_choosen_time_period <- function(data) {
      data <- colMeans(data[time_period,])
    }
    
    # Colmean_of_choosen_time_period function has only one input
    data <- lapply(data, Colmean_of_choosen_time_period)
    
    data <- unlist(data)
    
    return(data)
  }
  
  
  
  filtered_YFP_data_all <- Get_colmean_of_choosen_time_period(FI_data_used_here,time_period)
  
  length(filtered_YFP_data_all)
  
  
  
  title_plot <- paste("Time point",min(time_period),"to",max(time_period))
  
  png(paste0(dir_to_save,"/",title_plot," density.png"), width = 100, height = 100, units="mm",res=300)
  plot(density(filtered_YFP_data_all),main = title_plot)
  
  dev.off()
  
  png(paste0(dir_to_save,"/",title_plot," hist.png"), width = 100, height = 100, units="mm",res=300)
  hist(filtered_YFP_data_all,breaks =1000,main = title_plot)
  dev.off()
  
  #density_data <- density(filtered_YFP_data_all,bw =7)
  
  
  # fit a binomal distribution model to find two main peak--------------
  
  library(flexmix)

  kde <- density(filtered_YFP_data_all)
  m1 <- FLXMRglm(family = "gaussian")
  m2 <- FLXMRglm(family = "gaussian")
  fit <- flexmix(filtered_YFP_data_all ~ 1, data = as.data.frame(filtered_YFP_data_all), k = 2, model = list(m1, m2))
  c1 <- parameters(fit, component=1)[[1]]
  c2 <- parameters(fit, component=2)[[1]]
  

  png(paste0(dir_to_save,"/",title_plot," binormal.png"), width = 100, height = 100, units="mm",res=300)
  plot(kde)
  abline(v=c1[[1]], lty=2, col='blue')
  abline(v=c2[[1]], lty=2, col='red')
  dev.off()
  
  c1_peak <- min(c1[[1]],c2[[1]])
  c2_peak <- max(c1[[1]],c2[[1]])
  # density_data$x[which(diff(sign(diff(density_data$y)))==-2)+1]
  # density_data$y[which(diff(sign(diff(density_data$y)))==-2)+1]
  # 
  # minimum_point <- optimize(approxfun(density_data$x,density_data$y),interval=c(300,500))$minimum
  # 
  # density_data$x[which.min(abs(diff(density_data$y)))]
  
  
  # find a point in the range and give smallest difference (L1)----------------
  filtered_YFP_data_all <- sort(filtered_YFP_data_all)
  # how many points in this range
  length_range <- length(filtered_YFP_data_all[filtered_YFP_data_all>c1_peak & filtered_YFP_data_all<c2_peak])
  
  # whole range, 1/5 of length_range to calculate distance
  
  c1_point <- which((filtered_YFP_data_all>c1_peak)==1)[1]
  c2_point <- which((filtered_YFP_data_all>c2_peak)==1)[1]-1
  
  # xrange
  c1_point_x <- c1_point+trunc(length_range/10)  # trunc ceiling and round etc.
  c2_point_x <- c2_point-trunc(length_range/10)
  
  #index <- which.max(as.numeric(filtered_YFP_data_all > minimum_point))
  
  filtered_YFP_data_all_limit <- filtered_YFP_data_all[c1_point:c2_point]
  loss_function <- function(x,data=filtered_YFP_data_all_limit) {
    x <- which.min(abs(x-data))
    data <- filtered_YFP_data_all_limit[(x-trunc(length_range/10)-1):(x+trunc(length_range/10)-1)]
    loss <- sum((x-data)^2)
    return(loss)
  }

  minimum_point <- optim(par=filtered_YFP_data_all[c1_point_x], fn=loss_function, gr=NULL, 
                         data=filtered_YFP_data_all_limit,method =  "Brent",
                         lower = filtered_YFP_data_all[c1_point_x], 
                         upper = filtered_YFP_data_all[c2_point_x])$par
  
  

  # optimise(loss_function, interval=c(min(filtered_YFP_data_all),max(filtered_YFP_data_all)),
  #          maximum = FALSE,
  #          tol = .Machine$double.eps^0.25)
  ###################################################
  
  
  get_outlier_cell_YFP <- function(p) {
    
    
    max_p <- colMeans(p[time_period,])
    
    up_bound <- minimum_point
    
    result <- max_p[max_p>up_bound]
    
    return(result)
    # plot(density(max_p))
    # 
    # qnorm(0.995, mean = mu, sd = sigma, lower.tail = F, log.p = FALSE)
  }
  
  all_outlier_cell <- lapply(FI_data_used_here, get_outlier_cell_YFP)
  
  
  
  
  # 11. check precision and recall----------
  
  summary_prediction <- data.frame(matrix(ncol = 9,nrow = 0))
  colnames(summary_prediction) <- c("A:predicted_survive","B:Real_survive","Cells_in_A_not_B","Cells_in_B_not_A","#Overlap_of_A_and_B","#Cells_in_A_not_B","#Cells_in_B_not_A","Precision","Recall")
  
  for (i in 1:length(survive_in_data)) {
    index1 <- gsub("cell_cluster","",names(all_outlier_cell[[i]])) %in% survive_in_data[[i]]
    index2 <- survive_in_data[[i]] %in% gsub("cell_cluster","",names(all_outlier_cell[[i]]))
    
    
    summary_prediction[i,1] <- length(all_outlier_cell[[i]])
    summary_prediction[i,2] <- length(sort(unique(survive_in_data[[i]])))
    summary_prediction[i,3] <- paste(gsub("cell_cluster","",names(all_outlier_cell[[i]]))[!index1],sep = ",",collapse = ",")
    summary_prediction[i,4] <- paste(unique(survive_in_data[[i]][!index2]), sep=",", collapse = ",")
    
    summary_prediction[i,5] <- sum(index1)
    summary_prediction[i,6] <- length(gsub("cell_cluster","",names(all_outlier_cell[[i]]))[!index1])
    summary_prediction[i,7] <- length(survive_in_data[[i]][!index2])
    
  }
  summary_prediction$Precision <- summary_prediction$`#Overlap_of_A_and_B`/(summary_prediction$`A:predicted_survive`)
  
  summary_prediction$Recall <- summary_prediction$`#Overlap_of_A_and_B`/summary_prediction$`B:Real_survive`
  
  
  rownames(summary_prediction) <- YFP_filenames_well_label[plot_order][well_num]
  
  summary_prediction[19,] <- c(sum(summary_prediction$`A:predicted_survive`),sum(summary_prediction$`B:Real_survive`),"","",sum(summary_prediction$`#Overlap_of_A_and_B`),sum(summary_prediction$`#Cells_in_A_not_B`),sum(summary_prediction$`#Cells_in_B_not_A`),sum(summary_prediction$`#Overlap_of_A_and_B`)/sum(summary_prediction$`A:predicted_survive`),sum(summary_prediction$`#Overlap_of_A_and_B`)/sum(summary_prediction$`B:Real_survive`))
  
  
  
  rownames(summary_prediction)[19] <- "Total"
  
  #summary_prediction
  
  write.xlsx(summary_prediction,paste(dir_to_save,"20191123 prediction", paste("Time point",min(time_period),"to",max(time_period)),"h.xlsx"),row.names=T)
  
  Precision <- as.numeric(summary_prediction$Precision[nrow(summary_prediction)])
  Recall <- as.numeric(summary_prediction$Recall[nrow(summary_prediction)])
  F1_score <- (2*Precision*Recall)/(Precision+Recall)
  time_point_used <- ifelse(max(time_period)==min(time_period),
                            as.character(max(time_period)),
                            paste0(min(time_period),":",max(time_period)))
  F1_score_table_YFP[k,] <- c(Precision,Recall,F1_score,time_point_used)
  
}


# 10.3 plot F1_score_table figure-----------
library(ggplot2)
library(reshape2)
library(cowplot)
F1_score_table_YFP$time_point_used <- factor(F1_score_table_YFP$time_point_used, level=unique(F1_score_table_YFP$time_point_used))
F1_score_table_YFP_melt <- melt(F1_score_table_YFP,id.vars = "time_point_used")
F1_score_table_YFP_melt$value <- as.numeric(F1_score_table_YFP_melt$value)


ggplot(data = F1_score_table_YFP_melt, aes(x=time_point_used,y=value,group=variable)) +
  geom_point(aes(col=variable),size=5) +
  geom_line(aes(col=variable),size=3) +
  scale_y_continuous(breaks = seq(0.5,1, 0.05),limits = c(0.5,1)) +
  geom_hline(yintercept=0.95,linetype="dashed",color = "red", size=2) +
  labs(x="Time point used (YFP data)")






#####################################
# Combine mCherry and YFP, to Predict survive#################################




# mCherry and YFP data
FI_data_used_here_1 <- filtered_mCherry_data[well_num]    
FI_data_used_here_2 <- filtered_YFP_data[well_num]


# combine mCherry and PDR5--------
F1_score_table_combine <- data.frame(matrix(ncol = 7,nrow = 0))
colnames(F1_score_table_combine) <- c("Precision_test","Recall_test","F1_score_test",
                                      "Precision_train","Recall_train","F1_score_train",
                                      "time_point_used")

dir_to_save <- "Z:/MD/Analysis_data_Zhisheng/15th_20191117_for_prediction_10X/analysis/"
dir_to_save <- paste0(dir_to_save,"combine_mCherry_and_YFP2/")


for (k in 1:length(all_try_time_points)) {
  
  k=2
  
  time_period <- unlist(all_try_time_points[k])
  
  #data <- filtered_YFP_data[well_num]
  Get_colmean_of_choosen_time_period <- function(data,time_period) {
    
    Colmean_of_choosen_time_period <- function(data) {
      data <- colMeans(data[time_period,])
    }
    
    # Colmean_of_choosen_time_period function has only one input
    data <- lapply(data, Colmean_of_choosen_time_period)
    
    data <- unlist(data)
    
    return(data)
  }
  
  
  filtered_mCherry_data_all <- Get_colmean_of_choosen_time_period(FI_data_used_here_1,time_period)
  filtered_YFP_data_all <- Get_colmean_of_choosen_time_period(FI_data_used_here_2,time_period)
  
  


  # create input data for logistic regression------------
  survive_or_not <- function(data,survive_in_data) {
    data <- gsub("cell_cluster","",names(data)) %in% survive_in_data
  }
  
  survive_or_not <- unlist(map2(FI_data_used_here_1, survive_in_data, survive_or_not))
  


  logistic_input <- data.frame("mCherry"=filtered_mCherry_data_all, 
                               "PDR5"=filtered_YFP_data_all,
                               "Survive"= survive_or_not)
  logistic_input$Survive <- factor(as.numeric(logistic_input$Survive), level=c("1","0"))
  # plot mCherry vs. PDR5 and correlation----------
  title_plot <- paste("Time point",min(time_period),"to",max(time_period)," (Raw data)")
  
  logistic_input_survive <- logistic_input[logistic_input$Survive=="1",]
  
  png(paste0(dir_to_save,"/",title_plot," combine.png"), width = 200, height = 200, units="mm",res=300)
  plot(filtered_mCherry_data_all,filtered_YFP_data_all,
       main = title_plot, xlab="mCherry",ylab="PDR5",cex.lab=1.5,cex.main=1.5,cex.axis=1.5,
       pch=21,bg="red",col="purple")
  points(logistic_input_survive$mCherry,logistic_input_survive$PDR5,pch=21,bg="blue",col="purple")
  text(x=0.8*max(filtered_mCherry_data_all),y=0.8*max(filtered_YFP_data_all),
       paste("r = ",round(cor(filtered_mCherry_data_all,
                              filtered_YFP_data_all,method = "pearson"),2)),
       cex = 2)
  legend("topright",legend = c("Reistance","Sensitive"),pch=21,col="purple",pt.bg=c("blue","red"),cex=1.5)

  dev.off()
  
  

  # fit a logistic model to predict survive--------------
  
  # 1. Pre-processing data using Caret

  # 2. Split the data into training set:test set = 8:2----------
  set.seed(123)
  index <- createDataPartition(logistic_input$Survive, p=0.8, list=FALSE)
  trainSet <- logistic_input[ index,]
  testSet <- logistic_input[-index,]
  
  str(trainSet)
  
  #grep("logistic",names(getModelInfo()))
  
  # 3. feature selection --------------
  #predictors <- c("mCherry")
  
  predictors <- c("mCherry","PDR5")
  
  outcomeName <- c("Survive")
  
  # 4. Training models using Caret-----------------
  
  #model_glm_logistic <- train(trainSet[,predictors],trainSet[,outcomeName],method='glm',family = binomial)
  
  #model_glm_logistic <- train(Survive ~ mCherry + PDR5, data=trainSet, method='glm',family = binomial)
  
  
  # use k-fold cross-validation to train the model
  fitControl <- trainControl(
    method = "repeatedcv",
    number = 5,
    repeats = 10)
  
  # 5. Parameter tuning using Caret------------
  #modelLookup(model='glm')
  
  # training the model
  model_glm_logistic <- train(trainSet[,predictors],trainSet[,outcomeName],
                            method='glm',family = binomial,
                            trControl=fitControl)
                            # tuneGrid=grid 
  # model_glm_logistic2 <- train(Survive ~ mCherry + PDR5, data=trainSet,
  #                             method='glm',family = binomial,
  #                             trControl=fitControl)
  # summarizing the model
  print(model_glm_logistic)
  
  summary(model_glm_logistic)
  
  #plot(model_glm_logistic)
  # 6. Variable importance estimation using caret----------------
  #Variable Importance
  # varImp(object=model_glm_logistic)
  # plot(varImp(object=model_glm_logistic),main="logistic - Variable Importance")
  # 
  # 7. Predictions using Caret------------------------
  
  predictions_train <- predict.train(object=model_glm_logistic,trainSet[,predictors],type="raw")
  
  F_train_table <- table(predictions_train,trainSet[,outcomeName])
  
  Precision_train <- precision(F_train_table, relevant = rownames(F_train_table)[1])
  Recall_train <- recall(F_train_table, relevant = rownames(F_train_table)[1])
  F1_train <- F_meas(F_train_table, relevant = rownames(F_train_table)[1], beta = 1)
  
  
  
  predictions<-predict.train(object=model_glm_logistic,testSet[,predictors],type="raw")
  
  
  
  F_table <- table(predictions,testSet[,outcomeName])
  
  # Accuracy is percent of correctly classified samples
  
  confusionMatrix(predictions,testSet[,outcomeName])
  
  
  Precision <- precision(F_table, relevant = rownames(F_table)[1])
  Recall <- recall(F_table, relevant = rownames(F_table)[1])
  F1 <- F_meas(F_table, relevant = rownames(F_table)[1], beta = 1)
  
  time_point_used <- ifelse(max(time_period)==min(time_period),
                            as.character(max(time_period)),
                            paste0(min(time_period),":",max(time_period)))
  
  F1_score_table_combine[k,] <- c(Precision,Recall,F1,
                                  Precision_train,Recall_train,F1_train,
                                  time_point_used)
  
  
  
  
  predictions_all<-predict.train(object=model_glm_logistic,logistic_input,type="raw")
  
  logistic_input$Survive <- predictions_all
  
  
  
  logistic_input_survive_predicted <- logistic_input[logistic_input$Survive=="1",]
  
  title_plot <- paste("Time point",min(time_period),"to",max(time_period)," (Predicted data)")
  
  png(paste0(dir_to_save,"/",title_plot," combine_predicted.png"), width = 200, height = 200, units="mm",res=300)
  plot(filtered_mCherry_data_all,filtered_YFP_data_all,
       main = title_plot, xlab="mCherry",ylab="PDR5",cex.lab=1.5,cex.main=1.5,cex.axis=1.5,
       pch=21,bg="red",col="purple")
  points(logistic_input_survive_predicted$mCherry,logistic_input_survive_predicted$PDR5,pch=21,bg="blue",col="purple")
  text(x=0.8*max(filtered_mCherry_data_all),y=0.8*max(filtered_YFP_data_all),
       paste("r = ",round(cor(filtered_mCherry_data_all,
                              filtered_YFP_data_all,method = "pearson"),2)),
       cex = 2)
  legend("topright",legend = c("Reistance","Sensitive"),pch=21,col="purple",pt.bg=c("blue","red"),cex=1.5)
  
  dev.off()
  
}


write.xlsx(F1_score_table_combine,paste0(dir_to_save,"F1_score_table_combine.xlsx"))

# only mCherry--------
F1_score_table_mCherry <- data.frame(matrix(ncol = 7,nrow = 0))
colnames(F1_score_table_mCherry) <- c("Precision_test","Recall_test","F1_score_test",
                                      "Precision_train","Recall_train","F1_score_train",
                                      "time_point_used")
dir_to_save <- "Z:/MD/Analysis_data_Zhisheng/15th_20191117_for_prediction_10X/analysis/"
dir_to_save <- paste0(dir_to_save,"mCherry_predict/")


for (k in 1:length(all_try_time_points)) {
  
  

  
  
  time_period <- unlist(all_try_time_points[k])
  
  #data <- filtered_YFP_data[well_num]
  Get_colmean_of_choosen_time_period <- function(data,time_period) {
    
    Colmean_of_choosen_time_period <- function(data) {
      data <- colMeans(data[time_period,])
    }
    
    # Colmean_of_choosen_time_period function has only one input
    data <- lapply(data, Colmean_of_choosen_time_period)
    
    data <- unlist(data)
    
    return(data)
  }
  
  
  filtered_mCherry_data_all <- Get_colmean_of_choosen_time_period(FI_data_used_here_1,time_period)
  filtered_YFP_data_all <- Get_colmean_of_choosen_time_period(FI_data_used_here_2,time_period)
  
  
  
  
  # create input data for logistic regression------------
  survive_or_not <- function(data,survive_in_data) {
    data <- gsub("cell_cluster","",names(data)) %in% survive_in_data
  }
  
  survive_or_not <- unlist(map2(FI_data_used_here_1, survive_in_data, survive_or_not))
  
  
  
  logistic_input <- data.frame("mCherry"=filtered_mCherry_data_all, 
                               "PDR5"=filtered_YFP_data_all,
                               "Survive"= survive_or_not)
  logistic_input$Survive <- factor(as.numeric(logistic_input$Survive), level=c("1","0"))
  # plot mCherry vs. PDR5 and correlation----------
  title_plot <- paste("Time point",min(time_period),"to",max(time_period)," (Raw data)")
  
  logistic_input_survive <- logistic_input[logistic_input$Survive=="1",]
  
  png(paste0(dir_to_save,"/",title_plot," mCherry.png"), width = 200, height = 200, units="mm",res=300)
  plot(filtered_mCherry_data_all,filtered_YFP_data_all,
       main = title_plot, xlab="mCherry",ylab="PDR5",cex.lab=1.5,cex.main=1.5,cex.axis=1.5,
       pch=21,bg="red",col="purple")
  points(logistic_input_survive$mCherry,logistic_input_survive$PDR5,pch=21,bg="blue",col="purple")
  text(x=0.8*max(filtered_mCherry_data_all),y=0.8*max(filtered_YFP_data_all),
       paste("r = ",round(cor(filtered_mCherry_data_all,
                              filtered_YFP_data_all,method = "pearson"),2)),
       cex = 2)
  legend("topright",legend = c("Reistance","Sensitive"),pch=21,col="purple",pt.bg=c("blue","red"),cex=1.5)
  
  dev.off()
  
  
  logistic_input <- logistic_input[,c(1,3)]
  class(logistic_input$Survive)
  # fit a logistic model to predict survive--------------
  
  # 1. Pre-processing data using Caret
  
  # 2. Split the data into training set:test set = 8:2----------
  set.seed(123)
  index <- createDataPartition(logistic_input$Survive, p=0.8, list=FALSE)
  trainSet <- logistic_input[ index,]
  testSet <- logistic_input[-index,]
  
  str(trainSet)
  
  #grep("logistic",names(getModelInfo()))
  
  # 3. feature selection --------------
  predictors <- c("mCherry")
  
  #predictors <- c("mCherry","PDR5")
  
  outcomeName <- c("Survive")
  
  # 4. Training models using Caret-----------------
  
  #model_glm_logistic <- train(trainSet[,predictors],trainSet[,outcomeName],method='glm',family = binomial)
  
  #model_glm_logistic <- train(Survive ~ mCherry + PDR5, data=trainSet, method='glm',family = binomial)
  
  
  # use k-fold cross-validation to train the model
  fitControl <- trainControl(
    method = "repeatedcv",
    number = 5,
    repeats = 10)
  
  # 5. Parameter tuning using Caret------------
  #modelLookup(model='glm')
  
  # training the model
  model_glm_logistic <- train(trainSet[,predictors,drop=FALSE],trainSet[,outcomeName],
                              method='glm',family = binomial,
                              trControl=fitControl)
  # tuneGrid=grid 
  # model_glm_logistic2 <- train(Survive ~ mCherry + PDR5, data=trainSet,
  #                             method='glm',family = binomial,
  #                             trControl=fitControl)
  # summarizing the model
  print(model_glm_logistic)
  
  summary(model_glm_logistic)
  
  #plot(model_glm_logistic)
  # 6. Variable importance estimation using caret----------------
  #Variable Importance
  # varImp(object=model_glm_logistic)
  # plot(varImp(object=model_glm_logistic),main="logistic - Variable Importance")
  # 
  # 7. Predictions using Caret------------------------
  
  predictions_train <- predict.train(object=model_glm_logistic,trainSet[,predictors,drop=FALSE],type="raw")
  
  F_train_table <- table(predictions_train,trainSet[,outcomeName])
  
  Precision_train <- precision(F_train_table, relevant = rownames(F_train_table)[1])
  Recall_train <- recall(F_train_table, relevant = rownames(F_train_table)[1])
  F1_train <- F_meas(F_train_table, relevant = rownames(F_train_table)[1], beta = 1)
  
  
  
  predictions<-predict.train(object=model_glm_logistic,testSet[,predictors,drop=FALSE],type="raw")
  
  
  
  F_table <- table(predictions,testSet[,outcomeName])
  
  # Accuracy is percent of correctly classified samples
  
  confusionMatrix(predictions,testSet[,outcomeName])
  
  
  Precision <- precision(F_table, relevant = rownames(F_table)[1])
  Recall <- recall(F_table, relevant = rownames(F_table)[1])
  F1 <- F_meas(F_table, relevant = rownames(F_table)[1], beta = 1)
  
  time_point_used <- ifelse(max(time_period)==min(time_period),
                            as.character(max(time_period)),
                            paste0(min(time_period),":",max(time_period)))
  
  F1_score_table_mCherry[k,] <- c(Precision,Recall,F1,
                                  Precision_train,Recall_train,F1_train,
                                  time_point_used)
  
  
  
  
  predictions_all<-predict.train(object=model_glm_logistic,logistic_input,type="raw")
  
  logistic_input$Survive <- predictions_all
  
  logistic_input$PDR5 <- filtered_YFP_data_all
  
  
  
  
  logistic_input_survive_predicted <- logistic_input[logistic_input$Survive=="1",]
  
  title_plot <- paste("Time point",min(time_period),"to",max(time_period)," (Predicted data)")
  
  png(paste0(dir_to_save,"/",title_plot," mCherry_predicted.png"), width = 200, height = 200, units="mm",res=300)
  plot(filtered_mCherry_data_all,filtered_YFP_data_all,
       main = title_plot, xlab="mCherry",ylab="PDR5",cex.lab=1.5,cex.main=1.5,cex.axis=1.5,
       pch=21,bg="red",col="purple")
  points(logistic_input_survive_predicted$mCherry,logistic_input_survive_predicted$PDR5,pch=21,bg="blue",col="purple")
  text(x=0.8*max(filtered_mCherry_data_all),y=0.8*max(filtered_YFP_data_all),
       paste("r = ",round(cor(filtered_mCherry_data_all,
                              filtered_YFP_data_all,method = "pearson"),2)),
       cex = 2)
  legend("topright",legend = c("Reistance","Sensitive"),pch=21,col="purple",pt.bg=c("blue","red"),cex=1.5)
  
  dev.off()
  
}

write.xlsx(F1_score_table_mCherry,paste0(dir_to_save,"F1_score_table_mCherry.xlsx"))

# only PDR5-------------
F1_score_table_PDR5 <- data.frame(matrix(ncol = 7,nrow = 0))
colnames(F1_score_table_PDR5) <- c("Precision_test","Recall_test","F1_score_test",
                                      "Precision_train","Recall_train","F1_score_train",
                                      "time_point_used")
dir_to_save <- "Z:/MD/Analysis_data_Zhisheng/15th_20191117_for_prediction_10X/analysis/"
dir_to_save <- paste0(dir_to_save,"PDR5_predict/")


for (k in 1:length(all_try_time_points)) {
  
  

  
  time_period <- unlist(all_try_time_points[k])
  
  #data <- filtered_YFP_data[well_num]
  Get_colmean_of_choosen_time_period <- function(data,time_period) {
    
    Colmean_of_choosen_time_period <- function(data) {
      data <- colMeans(data[time_period,])
    }
    
    # Colmean_of_choosen_time_period function has only one input
    data <- lapply(data, Colmean_of_choosen_time_period)
    
    data <- unlist(data)
    
    return(data)
  }
  
  
  filtered_mCherry_data_all <- Get_colmean_of_choosen_time_period(FI_data_used_here_1,time_period)
  filtered_YFP_data_all <- Get_colmean_of_choosen_time_period(FI_data_used_here_2,time_period)
  
  
  
  
  # create input data for logistic regression------------
  survive_or_not <- function(data,survive_in_data) {
    data <- gsub("cell_cluster","",names(data)) %in% survive_in_data
  }
  
  survive_or_not <- unlist(map2(FI_data_used_here_1, survive_in_data, survive_or_not))
  
  
  
  logistic_input <- data.frame("mCherry"=filtered_mCherry_data_all, 
                               "PDR5"=filtered_YFP_data_all,
                               "Survive"= survive_or_not)
  logistic_input$Survive <- factor(as.numeric(logistic_input$Survive), level=c("1","0"))
  # plot mCherry vs. PDR5 and correlation----------
  title_plot <- paste("Time point",min(time_period),"to",max(time_period)," (Raw data)")
  
  logistic_input_survive <- logistic_input[logistic_input$Survive=="1",]
  
  png(paste0(dir_to_save,"/",title_plot," PDR5.png"), width = 200, height = 200, units="mm",res=300)
  plot(filtered_mCherry_data_all,filtered_YFP_data_all,
       main = title_plot, xlab="mCherry",ylab="PDR5",cex.lab=1.5,cex.main=1.5,cex.axis=1.5,
       pch=21,bg="red",col="purple")
  points(logistic_input_survive$mCherry,logistic_input_survive$PDR5,pch=21,bg="blue",col="purple")
  text(x=0.8*max(filtered_mCherry_data_all),y=0.8*max(filtered_YFP_data_all),
       paste("r = ",round(cor(filtered_mCherry_data_all,
                              filtered_YFP_data_all,method = "pearson"),2)),
       cex = 2)
  legend("topright",legend = c("Reistance","Sensitive"),pch=21,col="purple",pt.bg=c("blue","red"),cex=1.5)
  
  dev.off()
  
  
  logistic_input <- logistic_input[,c(2,3)]
  class(logistic_input$Survive)
  # fit a logistic model to predict survive--------------
  
  # 1. Pre-processing data using Caret
  
  # 2. Split the data into training set:test set = 8:2----------
  set.seed(123)
  index <- createDataPartition(logistic_input$Survive, p=0.8, list=FALSE)
  trainSet <- logistic_input[ index,]
  testSet <- logistic_input[-index,]
  
  str(trainSet)
  
  #grep("logistic",names(getModelInfo()))
  
  # 3. feature selection --------------
  predictors <- c("PDR5")
  
  #predictors <- c("mCherry","PDR5")
  
  outcomeName <- c("Survive")
  
  # 4. Training models using Caret-----------------
  
  #model_glm_logistic <- train(trainSet[,predictors],trainSet[,outcomeName],method='glm',family = binomial)
  
  #model_glm_logistic <- train(Survive ~ mCherry + PDR5, data=trainSet, method='glm',family = binomial)
  
  
  # use k-fold cross-validation to train the model
  fitControl <- trainControl(
    method = "repeatedcv",
    number = 5,
    repeats = 10)
  
  # 5. Parameter tuning using Caret------------
  #modelLookup(model='glm')
  
  # training the model
  model_glm_logistic <- train(trainSet[,predictors,drop=FALSE],trainSet[,outcomeName],
                              method='glm',family = binomial,
                              trControl=fitControl)
  # tuneGrid=grid 
  # model_glm_logistic2 <- train(Survive ~ mCherry + PDR5, data=trainSet,
  #                             method='glm',family = binomial,
  #                             trControl=fitControl)
  # summarizing the model
  print(model_glm_logistic)
  
  summary(model_glm_logistic)
  
  #plot(model_glm_logistic)
  # 6. Variable importance estimation using caret----------------
  #Variable Importance
  # varImp(object=model_glm_logistic)
  # plot(varImp(object=model_glm_logistic),main="logistic - Variable Importance")
  # 
  # 7. Predictions using Caret------------------------
  
  predictions_train <- predict.train(object=model_glm_logistic,trainSet[,predictors,drop=FALSE],type="raw")
  
  F_train_table <- table(predictions_train,trainSet[,outcomeName])
  
  Precision_train <- precision(F_train_table, relevant = rownames(F_train_table)[1])
  Recall_train <- recall(F_train_table, relevant = rownames(F_train_table)[1])
  F1_train <- F_meas(F_train_table, relevant = rownames(F_train_table)[1], beta = 1)
  
  
  
  predictions<-predict.train(object=model_glm_logistic,testSet[,predictors,drop=FALSE],type="raw")
  
  
  
  F_table <- table(predictions,testSet[,outcomeName])
  
  # Accuracy is percent of correctly classified samples
  
  confusionMatrix(predictions,testSet[,outcomeName])
  
  
  Precision <- precision(F_table, relevant = rownames(F_table)[1])
  Recall <- recall(F_table, relevant = rownames(F_table)[1])
  F1 <- F_meas(F_table, relevant = rownames(F_table)[1], beta = 1)
  
  time_point_used <- ifelse(max(time_period)==min(time_period),
                            as.character(max(time_period)),
                            paste0(min(time_period),":",max(time_period)))
  
  F1_score_table_PDR5[k,] <- c(Precision,Recall,F1,
                                  Precision_train,Recall_train,F1_train,
                                  time_point_used)
  
  
  
  
  predictions_all<-predict.train(object=model_glm_logistic,logistic_input,type="raw")
  
  logistic_input$Survive <- predictions_all
  
  logistic_input$mCherry <- filtered_mCherry_data_all
  
  
  
  
  logistic_input_survive_predicted <- logistic_input[logistic_input$Survive=="1",]
  
  title_plot <- paste("Time point",min(time_period),"to",max(time_period)," (Predicted data)")
  
  png(paste0(dir_to_save,"/",title_plot," PDR5_predicted.png"), width = 200, height = 200, units="mm",res=300)
  plot(filtered_mCherry_data_all,filtered_YFP_data_all,
       main = title_plot, xlab="mCherry",ylab="PDR5",cex.lab=1.5,cex.main=1.5,cex.axis=1.5,
       pch=21,bg="red",col="purple")
  points(logistic_input_survive_predicted$mCherry,logistic_input_survive_predicted$PDR5,pch=21,bg="blue",col="purple")
  text(x=0.8*max(filtered_mCherry_data_all),y=0.8*max(filtered_YFP_data_all),
       paste("r = ",round(cor(filtered_mCherry_data_all,
                              filtered_YFP_data_all,method = "pearson"),2)),
       cex = 2)
  legend("topright",legend = c("Reistance","Sensitive"),pch=21,col="purple",pt.bg=c("blue","red"),cex=1.5)
  
  dev.off()
  
}

write.xlsx(F1_score_table_PDR5,paste0(dir_to_save,"F1_score_table_PDR5.xlsx"))





#----------------------------plot-----------------------------------------------

library(ggplot2)
library(reshape2)
library(cowplot)
F1_score_table_combine$time_point_used <- factor(F1_score_table_combine$time_point_used, level=unique(F1_score_table_combine$time_point_used))
F1_score_table_combine_melt <- melt(F1_score_table_combine,id.vars = "time_point_used")
F1_score_table_combine_melt$value <- as.numeric(F1_score_table_combine_melt$value)


ggplot(data = F1_score_table_combine_melt, aes(x=time_point_used,y=value,group=variable)) +
  geom_point(aes(col=variable),size=5) +
  geom_line(aes(col=variable),size=3) +
  scale_y_continuous(breaks = seq(0.8,1, 0.025),limits = c(0.8,1)) +
  geom_hline(yintercept=0.95,linetype="dashed",color = "red", size=2) +
  labs(x="Time point used (mCherry and YFP data)")

#------compare diffeert predictions------------------------
F1_score_all <- data.frame(matrix(ncol = 4,nrow = 21))
colnames(F1_score_all) <- c("mcherry","PDR5","combine","time_point_used")

F1_score_all$mcherry <- F1_score_table_mCherry$F1_score_test
F1_score_all$PDR5 <- F1_score_table_PDR5$F1_score_test
F1_score_all$combine <- F1_score_table_combine$F1_score_test
F1_score_all$time_point_used <- factor(F1_score_table_combine$time_point_used, level=unique(F1_score_table_combine$time_point_used))


F1_score_all_melt <- melt(F1_score_all,id.vars = "time_point_used")
F1_score_all_melt$value <- as.numeric(F1_score_all_melt$value)

ggplot(data = F1_score_all_melt, aes(x=time_point_used,y=value,group=variable)) +
  geom_point(aes(col=variable),size=5) +
  geom_line(aes(col=variable),size=3) +
  scale_y_continuous(breaks = seq(0.8,1, 0.025),limits = c(0.8,1)) +
  geom_hline(yintercept=0.95,linetype="dashed",color = "red", size=2) +
  labs(x="Time point used",y="F1 score")



#------find abnormal cells----------------------------------

num_cells_in_well <- unlist(lapply(filtered_mCherry_data[well_num], length))

names(num_cells_in_well) <- mCherry_filenames_well_label[plot_order][well_num]

create_names <- function(num_cells_in_well) {
  result <- c()
  for (i in 1:length(num_cells_in_well)) {
    result <- c(result,rep(names(num_cells_in_well)[i],num_cells_in_well[i]))
  }
  return(result)
}

num_cells_in_well_all <- create_names(num_cells_in_well)

logistic_input <- data.frame("mCherry"=filtered_mCherry_data_all, 
                             "PDR5"=filtered_YFP_data_all,
                             "Survive"= survive_or_not)


rownames(logistic_input) <- paste(names(filtered_mCherry_data_all),num_cells_in_well_all)


logistic_input[logistic_input$Survive==TRUE & logistic_input$mCherry > 230,]












