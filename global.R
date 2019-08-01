#############################################################################

# PHYSICAL HABITAT CALCULATOR
# DESCRIPTION: This code calculates Physical Habitat parameters from Virginia DEQ Probablistic 
#    Monitoring Data read in from a Microsoft Access (.mdb) file and outputs results to a comma 
#    separated value (.csv) file and a table within the Access file names FinalSummary
#    This version outputs three Relative Bed Stability (LRBS) results, utilizing the most recent 
#    Kaufmann et.al. 2008 publication for residual pool and wood volume corrections (LRBS_fin)
# AUTHORS: Emma Jones, Jason Hill, Larry Willis: 
#    Virginia DEQ Water Monitoring Division, Blue Ridge Regional Office
# See README.txt for additional help and troubleshooting information
# LAST UPDATED: 11/30/2018 by Emma Jones (emma.jones@deq.virginia.gov)
#               07/29/2019 by David Winters (David_Winters@comcast.net)

#############################################################################

start.time <- Sys.time()
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken



rm(list=ls())


# Install packages (only if first time using tool)
#install.packages("reshape2")  
#install.packages("plyr")
#install.packages("dplyr")
#install.packages("RODBC")
#install.packages("magrittr")

# Load packages
library(reshape2) #restructures data with melt(),dcast()
library(plyr) # data management functions
library(dplyr) # even better dataframe manipulation functions
#replace plyr and dplyr with tibbles and tidyr
library(RODBC) # connects to databases
library(magrittr) # piping package
#use purrr

library(tibble)
library(tidyr)
library(purrr)
library(stringr)
library(svSocket)

# Restrict number of decimal places reported in output
options(digits=3)

# Set your working directory, calculator will not run without correct file path
setwd("C:/Users/kkn73995/Downloads/PhabTools/")

####################################
#David's Functions

#This provides a simpler method of initializing an empty matrix
initMat <- function(row, col, rownames, colnames) {
  return(data.frame(matrix(vector(), row, col, dimnames=list(rownames, colnames)),
                    stringsAsFactors=F))
}

#Returns the average of the logarithms of the two inputs
logAvg <- function(input1, input2) {
  return((log10(input1)+log10(input2))/2)
}

#This is a common equation found below, self explanatory
stacksEq <- function(slope, interval) {
  return((0.12 + 0.25*slope)*interval)
}

#Another commonly found equation in the code
areaEq <- function(depth, interval) {
  return(depth*interval/100)
}

accessData <- function(filePath) {
  #Names of the input tables
  list1 <- c("Reach", "Thalweg", "Wet", "BankFull", "Wood", "Incision", "Embeddedness", "Substrate")
  
  #New names assigned to the tables as used in the code
  list2 <- c("reach", "thalweg", "wet", "bankfull", "wood", "incision", "embed", "substrate")
  
  #For each pair of names, apply a function to input the data from Microsoft Access
  #This is used to fix an issue between moving between 32 and 64 bit computers
  for(i in 1:length(list1)) {
    access_query_32_EVJ(filePath, list1[[i]], list2[[i]])
  }
}

create1 <- function(thalwegX) {
  #This creates a new table by cutting out the center X-7 lines of data,
  #just leaving the identificating data and adding 4 new columns of information
  #calculated from the cut out data
  len <- ncol(thalwegX)
  long <- nrow(thalwegX)
  output <- initMat(long, 11, c(), 
            c(names(thalwegX)[c(1:3,(len-3):len)], "Xdepth","Vdepth","Xdep_Vdep","Stacks_Equation"))
  output[1:7] <- thalwegX[,-(4:(len-4))]

  
  #Remove the first 3 and last 4 columns of the input table
  Ttest<-thalwegX[,-(c(1:3,(len-3):len))]
  
  #The final columns are row: means, standard deviation, difference between mean and 
  #standard deviation, and an application of the Stacks equation to columns 6 and 7 
  #(see above)
  output[,8] <- apply(Ttest, 1, mean)
  output[,9] <- apply(Ttest, 1, sd)
  output[,10] <- output[,8] - output[,9]
  output[,11] <- stacksEq(output[,7], output[,6])
  
  return(output)
}

pickmin <- function(input) {
  #A very simple function, but it allows us to use apply to return the 
  #minimum values from two columns in conjunction with apply.
  #Ex. min(input[,1:2]) would return the minimum value out of ALL of the table values
  #This is used to return a vector of values instead of a single value
  return(min(input[1],input[2]))
}

expandThalweg <- function(input) {
  #The input to this function is a single value of the input title followed by 3 sets of data
  #stacked on top of one another. The first value is used to create the output titles.
  #Ex. 'AB0_Residual', 'AB0_Depth', 'AB0_Area'
  
  #The length of 1 of the 3 data inputs
  correctcolnum <- (length(input)-1)/3
  
  #We initialize the matrices to make it faster by allocating space all at once
  output <- initMat(correctcolnum, 3, c(), c(paste(c('Residual','Depth','Area'),input[1],sep='_')))
  quickmat <- initMat(correctcolnum, 3, c(), c())
  
  #This takes the data and breaks it into a matrix with 3 columns  
  quickmat[,1] <- as.numeric(input[2:(correctcolnum+1)]) #Current
  quickmat[,2] <- as.numeric(input[(correctcolnum+2):(2*correctcolnum+1)]) #Residual
  quickmat[,3] <- as.numeric(input[(2*correctcolnum+2):(3*correctcolnum+1)]) #Interval
  
  #This reformats the data and prepares it for output
  output[,1] <- quickmat[,1]
  output[,2] <- quickmat[,2] - output[,1]
  output[,3] <- as.numeric(output[,2])*as.numeric(quickmat[,3])/100
  
  return(output)
}

loopfunction <- function(input) {
  
  #The output of this function begins with the rows of its input, but then also includes
  #Residual, Depth, and Area columns for each substrate so it will have
  # 4*(number of substrates) + 11 identifying columns
  #Note that a lot of the weird definitions for having values in terms of the numbers of
  #columns of the input allow the function to be applied to multiple tables
  
  output <- initMat(nrow(input), 4*(ncol(input)-11)+11,c(), c())
  output[,1:ncol(input)] <- input
  names(output)[1:ncol(input)] <- names(input)
  output <- arrange(output, SampleID)
  
  #A temporary matrix used in creating the residuals for each substrate
  temp3 <- initMat(nrow(input), 2,c(), c())
  
  #This is a matrix of residuals. I had to use a for loop for this and could not use an 
  #apply function because each residual uses the previous residual as an input
  residualmat <- initMat(nrow(input), (ncol(input)-11),c(), c())
  residualmat[,1] <- apply(output[,c('Xdep_Vdep','AB0')], 1, pickmin)
  
  for(i in 2:ncol(residualmat)) {
    temp <- output['Stacks_Equation']
    temp2 <- output[7+i]
    temp3[,1] <- residualmat[,i-1]+temp[,1]
    temp3[,2] <- temp2[,1]
    residualmat[,i] <- apply( temp3[,1:2], 1, pickmin)
  }
  
  #This right here is the product of trying to over optimize. Basically I stack 1 row of names and
  #3 sets of data on top of one another so that I can use an apply function and remove one
  #of the 'for' loops.
  crazyinputmatrix <- initMat((3*length(output[[1]])+1), ((ncol(output)-11)/4+1), c(), c())
  NC <- ncol(crazyinputmatrix)
  
  #Row of substrate names
  crazyinputmatrix[1,] <- names(output[,8:(.25*(ncol(output)-11)+8)])
  
  #Rows of Residuals
  crazyinputmatrix[2:(length(output[[1]])+1),
                   1:(NC-1)] <- residualmat
  
  #Rows of Substrate Values
  crazyinputmatrix[(length(output[[1]])+2):(2*length(output[[1]])+1),
                   1:(NC-1)] <- output[8:(.25*(ncol(output)-11)+7)]
  
  #Rows of Interval
  crazyinputmatrix[(2*length(output[[1]])+2):(3*length(output[[1]])+1),
                   1:(NC-1)] <- output[,6]
  
  #This separator is here because R doesn't like it if you directly put
  #the data into the columns for some reason
  new <-  t(apply(crazyinputmatrix[,1:(ncol(crazyinputmatrix)-1)], 2, expandThalweg))
  
  output[,(ncol(input)+1):ncol(output)] <- new
  output[,(ncol(input)+1):ncol(output)] <-
    apply(output[,(ncol(input)+1):ncol(output)], c(1,2), as.numeric)
  
  #A for loop to put in the right names. It is short so I didn't bother making it an apply function
  for(i in 8:(ncol(input)-4)) {
    j <- ncol(input) -23 + 3*i
    names(output)[j:(j+2)] <- c(paste('Residual',names(output)[i],sep='_'),
                                        paste('Depth',names(output)[i],sep='_'),
                                        paste('Area',names(output)[i],sep='_'))
  }
  return(output)
}

sArea1 <- function(thalwegloop) {
  #An adjustment that allows for the same method to be used on multiple tables
  c <- (ncol(thalwegloop)-11)/4
  
  #Takes 7 rows from the input plus one random column to hold space for data
  sumArea <- thalwegloop[,c(1:3,(c+8):(c+11),8)]
  
  #Column 8 will be where we put the sum of the areas for each sample
  names(sumArea)[8] <- "AreaSum"
  
  #From the input, we extract the area columns starting at c+12 and then every third row
  tTest <- thalwegloop[,-(1:(c+11))]
  
  #Sneaky method of removing the 2nd and 3rd columns out of sets of 3 for 100
  for(i in 1:c) {
    tTest <- tTest[,-i]
    tTest <- tTest[,-i]
  }
  #Sum the areas for each sample and place in column 8
  sumArea[,8] <- apply(tTest, 1, sum)
  
  return(sumArea)
}

removeBadValues <- function(input, badValues) {
  #if a value is in the list badValues, then return NA
  if(input %in% badValues) {
    input <- NA
  }
  return(input)
}

removeNA <- function(input) {
  #Sets NA values to 0
  if(is.na(input)) {
    input <- 0
  }  
  return(input)
}

populateTempVec <- function(substrate, colnames, end) {
  #This method is used in create Wide to count the number of each substrate
  #in one row/ sample. It returns a row vector
  tempvector <- initMat(1, length(colnames), c(), colnames)

  for(j in 4:end) {
    #If a substrate exists, add 1 to its column. Substrates not found in the
    #sample will retrun as NA and will be changed later to 0
    if(!is.na(substrate[j])) {
      #The values in substrate are the same as the names of the columns of tempvector
      tempvector[[1,substrate[j]]] <- 
        ifelse(is.na(tempvector[1,substrate[j]]), 1, tempvector[1,substrate[j]] + 1)
    }
  }
  
  return(t(tempvector[1,1:(length(colnames)-4)]))
}

createWide <- function(substrate, remove) {
  #This suppresses a useless error from the 'gather' function
  oldw <- getOption("warn")
  options(warn = -1) 
  
  #This function skips over a few steps that Emma used and is used for multiple tables
  
  #Depending on the desired output indicated by 'remove', certain substrates will not be counted
  if(remove == 1) {
    badValues <- c('MM','OT','WD','','MISSING','NA','RC')
  } else {
    badValues <- c('MM','','MISSING','NA')
  }
  
  end <- ncol(substrate)
  range <- 4:end
  rem <- removeBadValues
  
  #This chunk expands the data frame (the same as melt from reshape2) and
  #records the unique substrates found in alphabetical order
  intermed <- as.data.frame(apply(substrate[,range], c(1,2), rem, badValues<-badValues)) 
  long <- gather(intermed, names(intermed), key='Useless', value='Values')
  uniquesubs <- unique(long[,2], na.rm=TRUE) %>% sort()
  
  #The substrate names are used as the column names for the matrices
  col1 <- c(uniquesubs, c('totalcount', 'SampleID', 'StationID', 'Date'))
  col2 <- c(c('SampleID', 'StationID', 'Date'), paste(uniquesubs, 'PCT', sep='_'))
  columns <- length(col1)
  
  #Ininitalize an empty vector to count the number of each substrate
  tempvector <- initMat(nrow(substrate), length(col1), c(), col1)
  
  #Initialize the output matrix in the format we want
  output <- initMat(nrow(substrate), length(col2), c(), col2)
  
  #This will mark any "Bad" value (Missing, blank, etc.) as 'xd'
  substrate[,4:ncol(substrate)] <- intermed
  
  #Of the acceptable values (see above) count the number of each substrate per SampleID
  tempvector[,1:(columns-4)] <- t(apply(substrate, 1, populateTempVec, colnames<-col1, end<-end))
  
  
  #Takes the SampleID, StationID, and Date for each sample from the input table
  tempvector[,ncol(tempvector)-2] <- substrate[,1]
  tempvector[,ncol(tempvector)-1] <- as.character(substrate[,2])
  tempvector[,ncol(tempvector)] <- as.Date(substrate[,3],origin = "1970-01-01")
  
  #Sets all NA values to 0
  tempvector[,1:(ncol(tempvector)-4)] <- apply(tempvector[,1:(ncol(tempvector)-4)], c(1,2), removeNA)
  
  #creates a new column in tempvector of the total number of substrates per row
  tempvector[,'totalcount'] <- apply(tempvector[,1:(ncol(tempvector)-4)], 1, sum)
  
  #Places the SampleID, StationID, and Date for each sample in the correct position 
  #in the output table
  output[,'SampleID'] <- tempvector[,'SampleID']
  output[,'StationID'] <- tempvector[,'StationID']
  output[,'Date'] <- tempvector[,'Date']
  
  #Calculates percentages of each substrate
  #This works because of the order we named the columns in tempvector
  for(j in 1:(ncol(tempvector)-4)) {
    output[,j+3] <- tempvector[,j]/tempvector[,ncol(tempvector)-3]
  }
  
  #Changes data types to the format we want
  output[,2] <- as.factor(output[,2])
  output[,3] <- as.Date(output[,3], origin = "1970-01-01")
  output <- arrange(output, SampleID)
  output[,3] <- as.character(output[,3])
  output[,3] <- as.POSIXct(output[,3], origin = "1970-01-01")
  
  #Returns warnings to normal
  options(warn = oldw)
  return(output)
}

createSubSum <- function(sw) {
  
  #This summarizes the substrates by creating a new data frame and combining them
  #in a sum and rescaled by some particular values I assume are specific to each substrate
  substrate_summary <- data.frame(SampleID = sw['SampleID'], 
                                  StationID = sw['StationID'], Date = sw['Date'],
                                  cb_PCT= sw['cb_PCT'], CB_PCT= sw['CB_PCT'], fn_PCT= sw['fn_PCT'],
                                  FN_PCT= sw['FN_PCT'], GC_PCT= sw['GC_PCT'], GF_PCT= sw['GF_PCT'],
                                  hp_PCT= sw['hp_PCT'], HP_PCT= sw['HP_PCT'], rr_PCT= sw['rr_PCT'],
                                  RR_PCT= sw['RR_PCT'], RS_PCT= sw['RS_PCT'], sa_PCT= sw['sa_PCT'],
                                  SA_PCT= sw['SA_PCT'], SB_PCT= sw['SB_PCT'], XB_PCT= sw['XB_PCT'],
                                  BL_PCT= 1, RRprop=1, BLprop=1, CBprop=1, GCprop=1,
                                  GFprop=1, SAprop=1, FNprop=1, LSUB_DMM=1, SUB_DMM=1)
  
  substrate_summary['RR_PCT'] <- (sw['RR_PCT']+sw['RS_PCT']+sw['HP_PCT'])
  substrate_summary['BL_PCT'] <- (sw['XB_PCT']+sw['SB_PCT'])
  substrate_summary['RRprop'] <- ((sw['RR_PCT']+sw['RS_PCT']+sw['HP_PCT'])*logAvg(4000, 8000))
  substrate_summary['BLprop'] <- ((sw['XB_PCT']+sw['SB_PCT'])*logAvg(250, 4000))
  substrate_summary['CBprop'] <- (sw['CB_PCT']*logAvg(64, 250))
  substrate_summary['GCprop'] <- (sw['GC_PCT']*logAvg(16, 64))
  substrate_summary['GFprop'] <- (sw['GF_PCT']*logAvg(2, 16))
  substrate_summary['SAprop'] <- (sw['SA_PCT']*logAvg(0.06, 2))
  substrate_summary['FNprop'] <- (sw['FN_PCT']*logAvg(0.001, 0.06))
  substrate_summary['LSUB_DMM'] <- (substrate_summary['RRprop'] + substrate_summary['BLprop']
                                    + substrate_summary['CBprop'] + substrate_summary['GCprop']
                                    + substrate_summary['GFprop'] + substrate_summary['SAprop']
                                    + substrate_summary['FNprop']) 
  substrate_summary['SUB_DMM'] <- (10^substrate_summary['LSUB_DMM'])
  
  return(substrate_summary)
}

createSubPCT <- function(substrate.wide2.2) {
  #This converts the percentages in substrate.wide2.2 from decimals to integers
  #and rearranges the columns. Also it sums some percentages together.
  
  substratePCT <- data.frame(SampleID = substrate.wide2.2['SampleID'],
                             StationID = substrate.wide2.2['StationID'],
                             Date = substrate.wide2.2['Date'],
                             BR_PCT=1,
                             HP_PCT=substrate.wide2.2['HP_PCT']*100,
                             RC_PCT=substrate.wide2.2['RC_PCT']*100,
                             BL_PCT=1,
                             CB_PCT=substrate.wide2.2['CB_PCT']*100,
                             GC_PCT=substrate.wide2.2['GC_PCT']*100,
                             GF_PCT=substrate.wide2.2['GF_PCT']*100,
                             SA_PCT=substrate.wide2.2['SA_PCT']*100,
                             FN_PCT=substrate.wide2.2['FN_PCT']*100,
                             WD_PCT=substrate.wide2.2['WD_PCT']*100,
                             OT_PCT=substrate.wide2.2['OT_PCT']*100,
                             BL_CB_GR_PCT=1,
                             SA_FN_PCT=1,
                             TotSubstrate_PCT=1)
  substratePCT['BR_PCT'] <- (substrate.wide2.2['RR_PCT']+substrate.wide2.2['RS_PCT'])*100
  substratePCT['BL_PCT'] <- (substrate.wide2.2['XB_PCT']+substrate.wide2.2['SB_PCT'])*100
  substratePCT['BL_CB_GR_PCT'] <- (substrate.wide2.2['XB_PCT']*100+substrate.wide2.2['SB_PCT']*100
                                   +substrate.wide2.2['CB_PCT']*100
                                   +substrate.wide2.2['GC_PCT']*100
                                   +substrate.wide2.2['GF_PCT']*100)
  substratePCT['SA_FN_PCT'] <-(substrate.wide2.2['SA_PCT']*100+substrate.wide2.2['FN_PCT']*100)
  substratePCT['TotSubstrate_PCT'] <- ((substrate.wide2.2['RR_PCT']*100+substrate.wide2.2['RS_PCT']*100)
                                       +substrate.wide2.2['HP_PCT']*100
                                       +substrate.wide2.2['RC_PCT']*100+(substrate.wide2.2['XB_PCT']
                                                                         +substrate.wide2.2['SB_PCT'])*100+substrate.wide2.2['CB_PCT']*100
                                       +substrate.wide2.2['GC_PCT']*100+substrate.wide2.2['GF_PCT']*100
                                       +substrate.wide2.2['SA_PCT']*100+substrate.wide2.2['FN_PCT']*100
                                       +substrate.wide2.2['WD_PCT']*100+substrate.wide2.2['OT_PCT']*100)
  return(substratePCT)
}

createssalt <- function(sPCT, ss) {
  #This merges particular columns together from the two inputs
  intersum <- ss[c('SampleID', 'StationID', 'Date', 'LSUB_DMM', 'SUB_DMM')]
  intersum <- arrange(intersum, SampleID)
  sPCT <- arrange(sPCT, SampleID)
  #This helps keep the dataframes in the correct format
  output <- initMat(nrow(sPCT),(ncol(sPCT)+2),c(),c()) 
  output <- merge(sPCT, intersum)
  return(output)
}

createINssalt <- function(sPCT, ss) {
  #This merges particular columns together from the two inputs
  intersum <- ss[c('SampleID', 'StationID', 'Date', 'IN_LSUB_DMM', 'IN_SUB_DMM')]
  intersum <- arrange(intersum, SampleID)
  sPCT <- arrange(sPCT, SampleID)
  #This helps keep the dataframes in the correct format
  output <- initMat(nrow(sPCT),(ncol(sPCT)+2),c(),c()) 
  output <- merge(sPCT, intersum)
  return(output)
}

createOUTssalt <- function(sPCT, ss) {
  #This merges particular columns together from the two inputs
  intersum <- ss[c('SampleID', 'StationID', 'Date', 'OUT_LSUB_DMM', 'OUT_SUB_DMM')]
  intersum <- arrange(intersum, SampleID)
  sPCT <- arrange(sPCT, SampleID)
  #This helps keep the dataframes in the correct format
  output <- initMat(nrow(sPCT),(ncol(sPCT)+2),c(),c()) 
  output <- merge(sPCT, intersum)
  return(output)
}

createISS <- function(sw) {
  #Takes the in substrate columns and rearranges, sums, and multiplies them by different constants
  ISS <- data.frame(SampleID = sw['SampleID'],
                    StationID = sw['StationID'],
                    Date = sw['Date'],
                    cb_PCT= sw['cb_PCT'], CB_PCT= sw['CB_PCT'],
                    FN_PCT= sw['FN_PCT'], GC_PCT= sw['GC_PCT'], 
                    GF_PCT= sw['GF_PCT'], HP_PCT= sw['HP_PCT'], 
                    rr_PCT= sw['rr_PCT'], RR_PCT= sw['RR_PCT'], 
                    RS_PCT= sw['RS_PCT'], sa_PCT= sw['sa_PCT'],
                    SA_PCT= sw['SA_PCT'], SB_PCT= sw['SB_PCT'],
                    XB_PCT= sw['XB_PCT'], BL_PCT= 1, IN_RRprop=1,
                    IN_BLprop=1, IN_CBprop=1, IN_GCprop=1, IN_GFprop=1,
                    IN_SAprop=1, IN_FNprop=1, IN_LSUB_DMM=1, IN_SUB_DMM=1)
  
  ISS['RR_PCT'] <- (sw['RR_PCT']+sw['RS_PCT']+sw['HP_PCT'])
  ISS['BL_PCT'] <- (sw['XB_PCT']+sw['SB_PCT'])
  ISS['IN_RRprop'] <- ((sw['RR_PCT']+sw['RS_PCT']+sw['HP_PCT'])*logAvg(4000, 8000))
  ISS['IN_BLprop'] <- ((sw['XB_PCT']+sw['SB_PCT'])*logAvg(250, 4000))
  ISS['IN_CBprop'] <- (sw['CB_PCT']*logAvg(64, 250))
  ISS['IN_GCprop'] <- (sw['GC_PCT']*logAvg(16, 64))
  ISS['IN_GFprop'] <- (sw['GF_PCT']*logAvg(2, 16))
  ISS['IN_SAprop'] <- (sw['SA_PCT']*logAvg(0.06, 2))
  ISS['IN_FNprop'] <- (sw['FN_PCT']*logAvg(0.001, 0.06))
  ISS['IN_LSUB_DMM'] <- (ISS['IN_RRprop'] + ISS['IN_BLprop'] + ISS['IN_CBprop'] + ISS['IN_GCprop']
                         + ISS['IN_GFprop'] + ISS['IN_SAprop'] + ISS['IN_FNprop']) 
  ISS['IN_SUB_DMM'] <- (10^ISS['IN_LSUB_DMM'])
  
  return(ISS)
}

createOSS <- function(sw) {
  #Takes the out substrate columns and rearranges, sums, and multiplies them by different constants
  OSS <- data.frame(SampleID = sw['SampleID'],
                    StationID = sw['StationID'],
                    Date = sw['Date'],
                    cb_PCT= sw['cb_PCT'], CB_PCT= sw['CB_PCT'],
                    fn_PCT= sw['fn_PCT'], FN_PCT= sw['FN_PCT'],
                    GC_PCT= sw['GC_PCT'], GF_PCT= sw['GF_PCT'],
                    hp_PCT= sw['hp_PCT'], HP_PCT= sw['HP_PCT'],
                    RR_PCT= sw['RR_PCT'], RS_PCT= sw['RS_PCT'],
                    SA_PCT= sw['SA_PCT'], SB_PCT= sw['SB_PCT'],
                    XB_PCT= sw['XB_PCT'], BL_PCT= 1, OUT_RRprop=1,
                    OUT_BLprop=1, OUT_CBprop=1, OUT_GCprop=1, OUT_GFprop=1,
                    OUT_SAprop=1, OUT_FNprop=1, OUT_LSUB_DMM=1, OUT_SUB_DMM=1)
  
  OSS['RR_PCT'] <- (sw['RR_PCT']+sw['RS_PCT']+sw['HP_PCT'])
  OSS['BL_PCT'] <- (sw['XB_PCT']+sw['SB_PCT'])
  OSS['OUT_RRprop'] <- ((sw['RR_PCT']+sw['RS_PCT']+sw['HP_PCT'])*logAvg(4000, 8000))
  OSS['OUT_BLprop'] <- ((sw['XB_PCT']+sw['SB_PCT'])*logAvg(250, 4000))
  OSS['OUT_CBprop'] <- (sw['CB_PCT']*logAvg(64, 250))
  OSS['OUT_GCprop'] <- (sw['GC_PCT']*logAvg(16, 64))
  OSS['OUT_GFprop'] <- (sw['GF_PCT']*logAvg(2, 16))
  OSS['OUT_SAprop'] <- (sw['SA_PCT']*logAvg(0.06, 2))
  OSS['OUT_FNprop'] <- (sw['FN_PCT']*logAvg(0.001, 0.06))
  OSS['OUT_LSUB_DMM'] <- (OSS['OUT_RRprop'] + OSS['OUT_BLprop'] + OSS['OUT_CBprop'] + OSS['OUT_GCprop']
                          + OSS['OUT_GFprop'] + OSS['OUT_SAprop'] + OSS['OUT_FNprop']) 
  OSS['OUT_SUB_DMM'] <- (10^OSS['OUT_LSUB_DMM'])
  
  return(OSS)
}

createINPCT <- function(sw) {
  #This combines substrate percentages and multiplies by 100 to get integer forms of the percents
  IPCT <- data.frame(SampleID = sw['SampleID'], StationID = sw['StationID'], Date = sw['Date'],
                     IN_BR_PCT=1, IN_HP_PCT=1, IN_RC_PCT=1, IN_BL_PCT=1,
                     IN_CB_PCT=1, IN_GC_PCT=1, IN_GF_PCT=1, IN_SA_PCT=1,
                     IN_FN_PCT=1, IN_WD_PCT=1, IN_OT_PCT=1, IN_BL_CB_GR_PCT=1,
                     IN_SA_FN_PCT=1, IN_TotSubstrate_PCT=1)
  
  IPCT['IN_BR_PCT'] <- (sw['RR_PCT']+sw['RS_PCT'])*100
  IPCT['IN_HP_PCT'] <- (sw['HP_PCT']*100)
  IPCT['IN_RC_PCT'] <- (sw['RC_PCT']*100)
  IPCT['IN_BL_PCT'] <- (sw['XB_PCT']+sw['SB_PCT'])*100
  IPCT['IN_CB_PCT'] <- (sw['CB_PCT']*100)
  IPCT['IN_GC_PCT'] <- (sw['GC_PCT']*100)
  IPCT['IN_GF_PCT'] <- (sw['GF_PCT']*100)
  IPCT['IN_SA_PCT'] <- (sw['SA_PCT']*100)
  IPCT['IN_FN_PCT'] <- (sw['FN_PCT']*100)
  IPCT['IN_WD_PCT'] <- (sw['WD_PCT']*100)
  IPCT['IN_OT_PCT'] <- (sw['OT_PCT']*100)
  IPCT['IN_BL_CB_GR_PCT'] <- (IPCT['IN_BL_PCT']+IPCT['IN_CB_PCT']+IPCT['IN_GC_PCT']+IPCT['IN_GF_PCT'])
  IPCT['IN_SA_FN_PCT'] <- (IPCT['IN_SA_PCT'] + IPCT['IN_FN_PCT']) 
  IPCT['IN_TotSubstrate_PCT'] <-(IPCT['IN_BR_PCT']+IPCT['IN_HP_PCT']+IPCT['IN_RC_PCT']+IPCT['IN_BL_PCT']+
                                   IPCT['IN_CB_PCT']+IPCT['IN_GC_PCT']+IPCT['IN_GF_PCT']+IPCT['IN_SA_PCT']+
                                   IPCT['IN_FN_PCT']+IPCT['IN_WD_PCT']+IPCT['IN_OT_PCT'])
  
  return(IPCT)
  
}

createOUTPCT <- function(sw) {
  #This combines substrate percentages and multiplies by 100 to get integer forms of the percents
  OPCT <- data.frame(SampleID = sw['SampleID'], StationID = sw['StationID'], Date = sw['Date'],
                     OUT_BR_PCT=1, OUT_HP_PCT=1, OUT_RC_PCT=1, OUT_BL_PCT=1,
                     OUT_CB_PCT=1, OUT_GC_PCT=1, OUT_GF_PCT=1, OUT_SA_PCT=1,
                     OUT_FN_PCT=1, OUT_WD_PCT=1, OUT_OT_PCT=1, OUT_BL_CB_GR_PCT=1,
                     OUT_SA_FN_PCT=1, OUT_TotSubstrate_PCT=1)
  
  OPCT['OUT_BR_PCT'] <- (sw['RR_PCT']+sw['RS_PCT'])*100
  OPCT['OUT_HP_PCT'] <- (sw['HP_PCT']*100)
  OPCT['OUT_RC_PCT'] <- (sw['RC_PCT']*100)
  OPCT['OUT_BL_PCT'] <- (sw['XB_PCT']+sw['SB_PCT'])*100
  OPCT['OUT_CB_PCT'] <- (sw['CB_PCT']*100)
  OPCT['OUT_GC_PCT'] <- (sw['GC_PCT']*100)
  OPCT['OUT_GF_PCT'] <- (sw['GF_PCT']*100)
  OPCT['OUT_SA_PCT'] <- (sw['SA_PCT']*100)
  OPCT['OUT_FN_PCT'] <- (sw['FN_PCT']*100)
  OPCT['OUT_WD_PCT'] <- (sw['WD_PCT']*100)
  OPCT['OUT_OT_PCT'] <- (sw['OT_PCT']*100)
  OPCT['OUT_BL_CB_GR_PCT'] <- (OPCT['OUT_BL_PCT']+OPCT['OUT_CB_PCT']+OPCT['OUT_GC_PCT']+OPCT['OUT_GF_PCT'])
  OPCT['OUT_SA_FN_PCT'] <- (OPCT['OUT_SA_PCT'] + OPCT['OUT_FN_PCT']) 
  OPCT['OUT_TotSubstrate_PCT'] <-(OPCT['OUT_BR_PCT']+OPCT['OUT_HP_PCT']+OPCT['OUT_RC_PCT']+OPCT['OUT_BL_PCT']+
                                    OPCT['OUT_CB_PCT']+OPCT['OUT_GC_PCT']+OPCT['OUT_GF_PCT']+OPCT['OUT_SA_PCT']+
                                    OPCT['OUT_FN_PCT']+OPCT['OUT_WD_PCT']+OPCT['OUT_OT_PCT'])
  
  return(OPCT)
}

check1 <- function(input) {
  #Take an input of a long row and break it into a column of substrates and a column of percentages 
  rlen <- length(input)/2
  tempvector <- initMat(rlen, 2, c(),c())
  tempvector[,1] <- input[1:rlen]
  tempvector[,2] <- input[(rlen+1):(2*rlen)]
  
  #Then return only the percentages we desire by filtering using the check2 function.
  output <- initMat(rlen, 1,c(),c())
  output <- apply(tempvector, 1, check2)
  return(output)
}

check2 <- function(input) {
  #The input to this function is a vector of a substrate and a percentage
  #If the substrate is in the list belowand the percentage does not exist,
  #return pi (3.14). If the input does exist, return the input. If the data selected is for 
  #a substrate not found in the list of gravel, boulders, etc. then return NA.
  if(isTRUE(input[1] %in% c('XB','SB','CB','GC','GF'))) {
    if(is.na(input[2])){
      input[2] <- 3.14
    }
    output <- input[2]
  } else {output <- NA}
  return(output)
}

finalcheck <- function(input) {
  #For a sample as input, if there is a single percentage present, return the sample
  #Else this means the entire row is filled with NA so it returns a row of 'remove'
  for(i in 1:length(input)){
    if(!is.na(input[i])){
      return(input)
    }
  }
  input[] <- 'remove'
  return(input)
}

removepi <- function(input) {
  #If the input is pi, return NA, else return the input.
  if(isTRUE(input==3.14)) {
    input<-NA
  }
  return(input)
}

simplify <- function(input) {
  #Used to count the elements that are neither 'MISSING' nor 'NA'
  if(isTRUE(input %in% c('MISSING', 'NA'))) {
    return(0)
  } else { return(1)}
}

simplify2 <- function(input) {
  #Used to count the number of elements in the list 'XB','SB','CB','GC','GF'
  if(isTRUE(input %in% c('XB','SB','CB','GC','GF'))) {
    return(1)
  } else { return(0)}
}

createes1 <- function(embed2, sub.wide2) {
  
  #This counts the total number of elements of embed in each row that are not 'MISSING' or NA
  embed.sub1.1 <- initMat(nrow(embed2), 4, c(), c('SampleID', 'StationID', 'Date', 'totalcount'))
  embed.sub1.1[,1:3] <- embed2[1:3]
  sub.wide3 <- apply(sub.wide2, c(1,2), simplify)
  embed.sub1.1[,4] <- apply(sub.wide3[,4:58],1,sum)
  return(embed.sub1.1)
}

createes1.5 <- function(embed2, sub.wide2) {
  
  #This counts the number of elements per row in embed that are 'XB', 'SB', 'CB', 'GC', or 'GF'
  #It also only includes the rows that have nonzero counts
  embed.sub1.1 <- initMat(nrow(embed2), 4, c(), c('SampleID', 'StationID', 'Date', 'BL_CB_GRCount'))
  embed.sub1.1[,1:3] <- embed2[1:3]
  sub.wide3 <- apply(sub.wide2, c(1,2), simplify2)
  embed.sub1.1[,4] <- apply(sub.wide3[,4:58],1,sum)
  
  embed.sub1.2 <- initMat(0, 4, c(), c('SampleID', 'StationID', 'Date', 'BL_CB_GRCount'))
  counter <- 1
  #This only adds a row to embed.sub1.2 if there is a nonzero count
  for(i in 1:nrow(embed.sub1.1)) {
    if(!isTRUE(embed.sub1.1[i,4]==0)) {
      embed.sub1.2[counter,] <- apply(embed.sub1.1[i,], c(1,2), as.character)
      counter <- counter + 1
    }
  }
  
  #General data management
  embed.sub1.2[,1] <- as.numeric(embed.sub1.2[,1])
  embed.sub1.2[,2] <- as.factor(embed.sub1.2[,2])
  embed.sub1.2[,3] <- as.POSIXct(embed.sub1.2[,3], origin='1970-01-01')
  embed.sub1.2[,4] <- as.numeric(embed.sub1.2[,4])
  
  return(embed.sub1.2)
}

createes2 <- function(embed2, sub.wide2) {
  
  #The goal of this function is to filter the inputs of embed2 using sub.wide2 as the criteria 
  #embed2 contains numbers and sub.wide2 contains substrate names (characters)
  # entry (i,j) in embed 2 corresponds with its name (i,j) of sub.wide2
  #If the substrate in sub.wide2 is 'XB','SB','CB','GC', or 'GF', then we keep the value
  # even if the value is NA. If the substrate is not in the above list, then the entry will be
  #NA. If the entire row is nothing but NA, then we remove the row, but not if one of the NA's
  # was kept from the previous table. That is the reason for this function being so messy.
  
  #This stacks the values of embed2 and sub.wide2 on top of each other so we can use an apply function
  sw3 <- initMat(2*nrow(embed2), 55, c(), names(sub.wide2)[4:length(sub.wide2)])
  sw3[1:nrow(embed2),] <- sub.wide2[,4:ncol(embed2)]
  for(i in 1:ncol(sub.wide2)){
    sw3[(nrow(sub.wide2)+1):(2*nrow(sub.wide2)),names(sub.wide2)[i]] <- embed2[names(sub.wide)[i]]
  }
  
  #This returns only the values we want, but it returns 3.14 instead of NA to differentiate
  #the two 'types' of NA
  sw4 <- initMat(nrow(embed2),ncol(embed2),c(),names(sub.wide2))
  sw4[,1:3]<-embed2[,1:3]
  sw4[,4:ncol(sw4)] <- apply(sw3[,1:55], 2, check1)
  
  #Marks a row for removal if only filled with NA
  sw4[,4:ncol(sw4)] <- t(apply(sw4[,4:ncol(sw4)], 1, finalcheck))
  
  #This is the same dataframe as above, but with all of the empty rows removed
  sw5 <- data.frame(matrix(nrow = nrow(sw4), ncol = ncol(sw4)))
  names(sw5) <- names(sw4)
  
  #Only adds the proper rows to sw5
  counter <- 1
  for(i in 1:nrow(sw4)) {
    if(!isTRUE(sw4[i,4]=='remove')){
      sw5[counter,c(1,4:ncol(sw4))] <- sw4[i,c(1,4:ncol(sw4))]
      sw5[counter,2] <- as.character(sw4[i,2])
      sw5[counter,3] <- as.character(sw4[i,3])
      counter <- counter + 1
    }
  }
  
  #General data management including removing the pis
  sw5 <- sw5[-(counter:nrow(sw5)),]
  sw5[,4:ncol(sw5)] <- apply(sw5[,4:ncol(sw5)], c(1,2), removepi)
  sw5[,1] <- as.numeric(sw5[,1])
  sw5[,2] <- as.factor(sw5[,2])
  sw5[,3] <- as.POSIXct(sw5[,3], origin='1970')
  sw5[,4:ncol(sw5)] <- apply(sw5[,4:ncol(sw5)], c(1,2), as.numeric)
  
  return(sw5)
}

createBCG <- function(embed2, sub.wide2, embed.sub1, embed.sub1.5, embed.sub2) {
  
  #This function is a convoluted way of merging 3 data frames: embed.sub1, embed.sub1.5, and embed.sub2
  #There are missing rows in two of the data frames so a simple merge does not produce the correct result
  #This method transposes the matrices and makes the sampleIDs into column names to match the 
  #data frames together
  
  #A vector of means
  BL_CB_GRmeanEmbed <- apply(embed.sub2[,4:58],1,mean,na.rm=T)
  
  #e1 is embed.sub.1 transposed
  e1 <- initMat(5, nrow(embed.sub1), c('SampleID','StationID','Date','BL_CB_GRmeanEmbed','BL_CB_GR_transectPCT'),
                c(1:848))
  e1[1:3,] <- t(embed.sub1)[1:3,]
  names(e1) <- 1:ncol(e1)
  
  #e1.5 is embed.sub.1.5 transposed
  e1.5 <- initMat(4, nrow(embed.sub1.5), c('SampleID','StationID','Date','BL_CB_GRCount'),
                  embed.sub1.5[,1])
  e1.5 <- t(embed.sub1.5)
  
  #e2.1 is embed.sub.2 and the vector of means transposed
  e2.0 <- data.frame(cbind(embed.sub2[,1:3],BL_CB_GRmeanEmbed))
  e2.1 <- initMat(4, 1:nrow(e2.0), c('SampleID','StationID','Date','BL_CB_GRmeanEmbed'),
                  c())
  e2.1 <- t(e2.0)
  
  #Set the column names of the vectors to SampleID
  names(e1) <- e1[1,]
  names(e1.5) <- e1.5[1,]
  colnames(e2.1) <- as.character(as.numeric(e2.1[1,]))
  
  #Append the means to e1
  for(i in colnames(e2.1)[1:ncol(e2.1)]){
    e1[4, as.numeric(i)] <- e2.1[4,as.character(as.numeric(i))]
  }
  
  #Calculate and add percentages of embed to e1
  counter <- 1
  for(i in names(e1.5)[1:ncol(e1.5)]){
    e1[5,as.numeric(i)] <- 100*as.numeric(e1.5[4,as.character(counter)])/as.numeric(embed.sub1[as.numeric(i),4])
    counter <- counter+1
  }
  
  #General management of data types - revert to nontransposed
  ef <- t(e1)
  ef[,5] <- round(as.numeric(ef[,5]), digits=3)
  ef[,3] <- as.character(ef[,3])
  ef <- as.data.frame(ef)
  ef[,1] <- as.numeric(ef[,1])
  ef[,3] <- as.POSIXct(ef[,3], origin='1970-01-01')
  ef[,4] <- as.numeric(as.character(ef[,4]))
  ef[,'BL_CB_GR_transectPCT'] <- as.numeric(as.character(ef[,'BL_CB_GR_transectPCT']))
  
  return(ef)
}

wsALT <- function(wood, BL_CB_GRmeanEmbed.df) {
  
  #This rescales some of the columns of wood, selects columns, appends vectors from the
  #environment, and merges with a second data frame
  wood_summary <- mutate(wood, XBKF_W=XBKF_W, ReachLength=ReachLength, wetsdsl=wetsdsl*0.058
                         ,wetsdml=wetsdml*0.182, wetsdll=wetsdll*0.438, wetmdsl=wetmdsl*0.333
                         ,wetmdml=wetmdml*1.042, wetmdll=wetmdll*2.501, wetldsl=wetldsl*0.932
                         ,wetldml=wetldml*2.911, wetldll=wetldll*6.988, wetxdsl=wetxdsl*3.016
                         ,wetxdml=wetxdml*9.421, wetxdll=wetxdll*22.620
                         ,VLW=wetsdsl+wetsdml+wetsdll+wetmdsl+wetmdml+wetmdll+wetldsl+wetldml+wetldll
                         +wetxdsl+wetxdml+wetxdll, VLW_msq=VLW/(ReachLength*XBKF_W))%>%
    select(SampleID,StationID,Date,ReachLength,VLW,VLW_msq)
  wood_summary.alt <- mutate(wood_summary,XBKF_H=XBKF_H,Xwid=Xwid,INC_H=INC_H,Xembed=Xembed,Vembed=Vembed) %>%
    merge(BL_CB_GRmeanEmbed.df,by=c('SampleID','StationID','Date')) %>%
    arrange(SampleID)
  return(wood_summary.alt)
}

wsIO <- function(wood_summary.alt) {
  #Appends vectors in the environment to the end of wood_summary.alt
  output <- mutate(wood_summary.alt,IN_Xembed=IN_Xembed,IN_Vembed=IN_Vembed,OUT_Xembed=OUT_Xembed
                   ,OUT_Vembed=OUT_Vembed)
  return(output)
}

create10QA <- function(thalweg, keepNames) {
  
  #This counts and calculates percentage of NA in thalweg where the thalweg number is 10
  thalweg_10QA.5 <- subset(thalweg, ThalwegN=='10', select=keepNames)
  thalweg_10QA.5[,4:103] <- ifelse(is.na(thalweg_10QA.5[,4:103]),1,0)
  thalweg_10QA.5['thalweg_NAcount'] <- sum(thalweg_10QA.5[,4:103])
  thalweg_10QA.5 <-thalweg_10QA.5[,-(4:107)]
  thalweg_10QA.5['thalwegPCT_NA'] <- (thalweg_10QA.5['thalweg_NAcount']/100)*100
  thalweg_10QA.5 <- arrange(thalweg_10QA.5,SampleID)
  
  return(thalweg_10QA.5)
}

createQA_thal <- function(thalweg,thalweg_10QA) {
  
  #This counts and calculates percentage of NA in thalweg where the thalweg number is
  # 15 and then merges with thalweg_10QA
  thalweg_10QA.5 <- subset(thalweg, ThalwegN=='15')
  thalweg_10QA.5[,4:153] <- ifelse(is.na(thalweg_10QA.5[,4:153]),1,0)
  thalweg_10QA.5['thalweg_NAcount'] <- sum(thalweg_10QA.5[,4:153])
  thalweg_10QA.5 <- thalweg_10QA.5[,-(4:157)]
  thalweg_10QA.5['thalwegPCT_NA'] <- (thalweg_10QA.5['thalweg_NAcount']/150)*100
  thalweg_10QA.5 <- arrange(thalweg_10QA.5,SampleID)
  QA_thalweg2 <- rbind(thalweg_10QA.5,thalweg_10QA) %>% arrange(SampleID)
  
  return(QA_thalweg2)
}

createQAgen <- function(input, name) {
  
  #This function counts and calculates the percentage of NA elements of the input
  #Used for 'wet', 'incision', 'bankfull', 'embed
  thalweg_10QA.5 <- input
  thalweg_10QA.5[,4:ncol(input)] <- ifelse(is.na(thalweg_10QA.5[,4:ncol(input)]),1,0)
  
  thalweg_10QA.5[paste(name, 'NAcount', sep='_')] <- apply(thalweg_10QA.5[,4:ncol(input)], 1, sum)
  thalweg_10QA.5 <-thalweg_10QA.5[,-(4:ncol(input))]
  thalweg_10QA.5[paste(name, 'PCT_NA', sep='')] <- (thalweg_10QA.5[paste(name, 'NAcount', sep='_')]/(ncol(input)-3))*100
  QA_embed2 <- arrange(thalweg_10QA.5,SampleID)
  
  return(QA_embed2)
}

createQA_sub <- function(substrate) {
  
  #This function counts and calculates the percentage of NA, 'MISSING' and '' elements of substrate
  thalweg_10QA.5 <- substrate
  thalweg_10QA.5[,4:108] <- ifelse(is.na(thalweg_10QA.5[,4:108]),1,
                                   ifelse(thalweg_10QA.5[,4:108]=='MISSING',1,
                                          ifelse(thalweg_10QA.5[,4:108]=='',1,0)))
  
  thalweg_10QA.5['substrateDataIssues'] <- apply(thalweg_10QA.5[,4:108], 1, sum)
  thalweg_10QA.5 <-thalweg_10QA.5[,-(4:108)]
  thalweg_10QA.5['substratePCT_NA'] <- (thalweg_10QA.5['substrateDataIssues']/105)*100
  QA_substrate2 <- arrange(thalweg_10QA.5,SampleID)
  
  return(QA_substrate2)
}

createQAsummary <- function(QA_thalweg,QA_bankfull,QA_wet,QA_incision,QA_embed,QA_substrate) {
  
  #This function merges the previously created summaries together and
  # if any parameters (except wood) have >50% missing data flag site
  output <- merge(QA_thalweg,QA_bankfull,by=c('SampleID','StationID','Date')) %>% 
    merge(QA_wet,by=c('SampleID','StationID','Date')) %>%
    merge(QA_incision,by=c('SampleID','StationID','Date')) %>%
    merge(QA_embed,by=c('SampleID','StationID','Date')) %>%
    merge(QA_substrate,by=c('SampleID','StationID','Date')) %>%
    mutate(SiteFlag=ifelse(thalwegPCT_NA>50|bankfullPCT_NA>50|wetPCT_NA>50|incisionPCT_NA>50
                           |embedPCT_NA>50|substratePCT_NA>50,"Flag","Site Accepted")) %>%
    arrange(SampleID)
  
  return(output)
}

createTMDLSummary <- function(summary) {
  rho=998; rhosed=2650; g=9.807
  
  #This function takes a data frame as input, applies functions to certain columns,
  # then selects certain columns and rearranges them into a final output.
  output <- mutate(summary,radiusBKF=((Xdepth/100)+XBKF_H)/2,Dcbf=13.7*radiusBKF*Slope,RBS1=SUB_DMM/Dcbf
                   ,LRBS1=log10(SUB_DMM/Dcbf),RW=VLW_msq*1000, RP=(RP100*0.5)*10
                   ,RBFmm=((Xdepth*10)+(XBKF_H*1000))/2,R_bf=RBFmm-RW-RP,D_cbf=(13.7*R_bf)*(Slope/100)
                   ,LDMB=log10(D_cbf),LRBS2=LSUB_DMM-LDMB,XBKF_W=XBKF_W,BKF_depth_in_meters=(Xdepth/100)+XBKF_H
                   ,BKFW_BKFD=XBKF_W/BKF_depth_in_meters,incised_depth=(Xdepth/100)+INC_H) %>%
    select(-c(Vdepth,Xdep_Vdep,Stacks_Equation,AreaSum,SUB_DMM,VLW,XBKF_H,INC_H,Vembed,radiusBKF,Dcbf
              ,RBS1,LRBS1,RW,RP,RBFmm,R_bf,D_cbf,LDMB)) %>% # Remove select variables
    select(SampleID,StationID,Date,SiteFlag,ThalwegN,ReachLength,Interval,Slope,RP100,BR_PCT,HP_PCT,RC_PCT
           ,BL_PCT,CB_PCT,GC_PCT,GF_PCT,SA_PCT,FN_PCT,WD_PCT,OT_PCT,BL_CB_GR_PCT,SA_FN_PCT,TotSubstrate_PCT
           ,LSUB_DMM,VLW_msq,Xdepth,Xwid,XBKF_W,BKF_depth_in_meters,BKFW_BKFD,incised_depth,Xembed
           ,BL_CB_GR_transectPCT,BL_CB_GRmeanEmbed,LRBS2)
  return(output)
}

createMasterSummary <- function(summary) {
  #This code applies a lot of calculations to the columns of the summary to create a master summary
  
  rho=998; rhosed=2650; g=9.807
  output <- mutate(summary,radiusBKF=((Xdepth/100)+XBKF_H)/2,Dcbf=13.7*radiusBKF*Slope
                   ,RBS1=SUB_DMM/Dcbf,LRBS1=log10(SUB_DMM/Dcbf) # End LRBS1 calculations 
                   ,RW=VLW_msq*1000,RP=(RP100*0.5)*10,RBFmm=((Xdepth*10)+(XBKF_H*1000))/2
                   ,R_bf=RBFmm-RW-RP,D_cbf=(13.7*R_bf)*(Slope/100),LDMB=log10(D_cbf) 
                   ,LRBS2=LSUB_DMM-LDMB,XBKF_W=XBKF_W,BKF_depth_in_meters=(Xdepth/100)+XBKF_H
                   ,BKFW_BKFD=XBKF_W/BKF_depth_in_meters,incised_depth=(Xdepth/100)+INC_H
                   ,incised_height_BKFD=INC_H/BKF_depth_in_meters
                   ,BKFW_incised_depth=XBKF_W/incised_depth
                   ,incised_D_BKFD=incised_depth/BKF_depth_in_meters # End LRBS2 calculations
                   ,LS=log10(0.01+Slope),LRP=log10(0.1+RP100),LW=log10(0.0001+VLW_msq)
                   ,Dgm=10^LSUB_DMM,Dgm_m=(10^LSUB_DMM)/1000,Dbf_th=XBKF_H+(Xdepth/100)
                   ,Rb3=0.65*Dbf_th
                   ,Ct_rpwd=(1.21*(((RP100/100)^1.08)*(((RP100/100)+VLW_msq)^0.638)))/(Dbf_th^3.32)
                   ,LCt_rpwd=log10(Ct_rpwd),Cp3_mill=(1/8)*((2.03*log10(12.2*(Rb3/((10^LSUB_DMM)/1000))))^-2)
                   ,Cp3_mill_2=ifelse(Cp3_mill<0.002,0.002,Cp3_mill),LCp3_mill=log10(Cp3_mill)
                   ,Ct_rpwd_cl=Ct_rpwd,Ct_rpwd_cl_2=ifelse(Ct_rpwd<Cp3_mill_2,Cp3_mill_2,Ct_rpwd)
                   ,Ct_rpwd_cl_3= ifelse(Ct_rpwd_cl_2=='NA',Cp3_mill,Ct_rpwd_cl_2)
                   ,Cp3Ctrpwd_rat=Cp3_mill_2/Ct_rpwd
                   ,Cp3Ctrpwd_rat=ifelse(Cp3Ctrpwd_rat>1.0,1.0,Cp3Ctrpwd_rat)
                   ,Rrpw3=Rb3*(Cp3Ctrpwd_rat^0.3333),Rrpw3_RtRat=Rrpw3/Rb3
                   ,Reyp3=(((g*Rb3*(Slope/100))^0.5)*((10^LSUB_DMM)/1000))/0.00000102
                   ,LReyp3=log10(Reyp3)
                   ,Shld_Px3= 0.5*((0.22*(Reyp3^(-0.6)))+(0.06*(10^(-7.7*(Reyp3^(-0.6))))))
                   ,Shld_Px3=ifelse(is.na(LReyp3),NA,Shld_Px3),LShld_Px3=log10(Shld_Px3)
                   ,LShld_Px3=ifelse(LReyp3<=1.4166,(-1.3991-(0.24419*LReyp3)),LShld_Px3)
                   ,Shld_Px3=10^(LShld_Px3)
                   ,Dcbf_fin=1000*((rho*g*Rrpw3*(Slope/100))/(Shld_Px3*(rhosed-rho)*g))
                   ,LDcbf_fin=log10(Dcbf_fin),LRBS_fin=LSUB_DMM-LDcbf_fin) 
  return(output)
}

updateQASummary <- function(QAsummary, MasterSummary) {
  #This merges the two inputs, arranges the data by SampleID, and selects certain columns
  output <- merge(QAsummary,MasterSummary,by=c('SampleID','StationID','Date','SiteFlag')) %>% 
    arrange(SampleID) %>%  
    select(SampleID,StationID,Date,thalweg_NAcount,thalwegPCT_NA,bankfull_NAcount,bankfullPCT_NA
           ,wet_NAcount,wetPCT_NA,incision_NAcount,incisionPCT_NA,embed_NAcount,embedPCT_NA
           ,substrateDataIssues,substratePCT_NA,SiteFlag,RBFmm,LS,LRP,LW,Dgm,Dgm_m,Dbf_th,Rb3
           ,Ct_rpwd,LCt_rpwd,Cp3_mill,Cp3_mill_2,LCp3_mill,Ct_rpwd_cl,Ct_rpwd_cl_2,Ct_rpwd_cl_3
           ,Cp3Ctrpwd_rat,Rrpw3,Rrpw3_RtRat,Reyp3,LReyp3,Shld_Px3,LShld_Px3,Dcbf_fin,LDcbf_fin)
  return(output)
}

updateMasterSummary <- function(MasterSummary) {
  #This selects certain columns from the previous version
  output <- select(MasterSummary,SampleID,StationID,Date,SiteFlag,ThalwegN,ReachLength,Interval
                   ,Slope,Xdepth,Vdepth,Xdep_Vdep,Stacks_Equation,AreaSum,RP100,BR_PCT,HP_PCT
                   ,RC_PCT,BL_PCT,CB_PCT,GC_PCT,GF_PCT,SA_PCT,FN_PCT,WD_PCT,OT_PCT,BL_CB_GR_PCT
                   ,SA_FN_PCT,TotSubstrate_PCT,LSUB_DMM,SUB_DMM,VLW,VLW_msq,XBKF_H,Xwid,INC_H
                   ,Xembed,Vembed,BL_CB_GRmeanEmbed,BL_CB_GR_transectPCT,radiusBKF,Dcbf,RBS1
                   ,LRBS1,RW,RP,R_bf,D_cbf,LDMB,LRBS2,XBKF_W,BKF_depth_in_meters,BKFW_BKFD,incised_depth
                   ,incised_height_BKFD,BKFW_incised_depth,incised_D_BKFD,LRBS_fin)
  return(output)
}

createINOUTMasterSummary <- function(INOUTsummary) {
  rho=998; rhosed=2650; g=9.807
  
  #I did not write this code so I do not know exactly what it does, but basically it just selects
  #certain columns and applies a lot of math operations to them
  output <- INOUTsummary %>%
    select(c(SampleID,StationID,Date,SiteFlag,ThalwegN,ReachLength,Interval,Slope,Xdepth
             ,RP100,VLW_msq,IN_BR_PCT:OUT_Vembed)) %>%
    mutate(IN_radiusBKF=((Xdepth/100)+XBKF_H)/2,IN_Dcbf=13.7*IN_radiusBKF*Slope
           ,IN_RBS1=IN_SUB_DMM/IN_Dcbf, IN_LRBS1=log10(IN_SUB_DMM/IN_Dcbf) # End IN_LRBS1 calculations 
           ,RW=VLW_msq*1000, RP=(RP100*0.5)*10, RBFmm=((Xdepth*10)+(XBKF_H*1000))/2
           ,R_bf=RBFmm-RW-RP, D_cbf=(13.7*R_bf)*(Slope/100), LDMB=log10(D_cbf) 
           ,IN_LRBS2=IN_LSUB_DMM-LDMB # End IN_LRBS2 calculations
           ,LS=log10(0.01+Slope), LRP=log10(0.1+RP100), LW=log10(0.0001+VLW_msq)
           ,IN_Dgm= 10^IN_LSUB_DMM, IN_Dgm_m=(10^IN_LSUB_DMM)/1000
           ,Dbf_th=XBKF_H+(Xdepth/100),Rb3=0.65*Dbf_th
           ,Ct_rpwd=(1.21*(((RP100/100)^1.08)*(((RP100/100)+VLW_msq)^0.638)))/(Dbf_th^3.32)
           ,LCt_rpwd=log10(Ct_rpwd) 
           ,IN_Cp3_mill=(1/8)*((2.03*log10(12.2*(Rb3/((10^IN_LSUB_DMM)/1000))))^-2)
           ,IN_Cp3_mill_2=ifelse(IN_Cp3_mill<0.002,0.002,IN_Cp3_mill), IN_LCp3_mill=log10(IN_Cp3_mill)
           ,IN_Ct_rpwd_cl=Ct_rpwd, IN_Ct_rpwd_cl_2= ifelse(Ct_rpwd<IN_Cp3_mill_2,IN_Cp3_mill_2,Ct_rpwd)
           ,IN_Ct_rpwd_cl_3= ifelse(IN_Ct_rpwd_cl_2=='NA',IN_Cp3_mill,IN_Ct_rpwd_cl_2)
           ,IN_Cp3Ctrpwd_rat=IN_Cp3_mill_2/Ct_rpwd
           ,IN_Cp3Ctrpwd_rat=ifelse(IN_Cp3Ctrpwd_rat>1.0,1.0,IN_Cp3Ctrpwd_rat)
           ,IN_Rrpw3=Rb3*(IN_Cp3Ctrpwd_rat^0.3333), IN_Rrpw3_RtRat=IN_Rrpw3/Rb3
           ,IN_Reyp3=(((g*Rb3*(Slope/100))^0.5)*((10^IN_LSUB_DMM)/1000))/0.00000102
           ,IN_LReyp3=log10(IN_Reyp3)
           ,IN_Shld_Px3= 0.5*((0.22*(IN_Reyp3^(-0.6)))+(0.06*(10^(-7.7*(IN_Reyp3^(-0.6))))))
           ,IN_Shld_Px3=ifelse(is.na(IN_LReyp3),NA,IN_Shld_Px3), IN_LShld_Px3=log10(IN_Shld_Px3)
           ,IN_LShld_Px3=ifelse(IN_LReyp3<=1.4166,(-1.3991-(0.24419*IN_LReyp3)),IN_LShld_Px3)
           ,IN_Shld_Px3=10^(IN_LShld_Px3)
           ,IN_Dcbf_fin=1000*((rho*g*IN_Rrpw3*(Slope/100))/(IN_Shld_Px3*(rhosed-rho)*g))
           ,IN_LDcbf_fin=log10(IN_Dcbf_fin), IN_LRBS_fin=IN_LSUB_DMM-IN_LDcbf_fin 
           ,OUT_radiusBKF=((Xdepth/100)+XBKF_H)/2, OUT_Dcbf=13.7*OUT_radiusBKF*Slope
           ,OUT_RBS1=OUT_SUB_DMM/OUT_Dcbf, OUT_LRBS1=log10(OUT_SUB_DMM/OUT_Dcbf) # End OUT_LRBS1 calculations 
           ,OUT_LRBS2=OUT_LSUB_DMM-LDMB # End OUT_LRBS2 calculations
           ,OUT_Dgm= 10^OUT_LSUB_DMM, OUT_Dgm_m=(10^OUT_LSUB_DMM)/1000
           ,OUT_Cp3_mill=(1/8)*((2.03*log10(12.2*(Rb3/((10^OUT_LSUB_DMM)/1000))))^-2)
           ,OUT_Cp3_mill_2=ifelse(OUT_Cp3_mill<0.002,0.002,OUT_Cp3_mill)
           ,OUT_LCp3_mill=log10(OUT_Cp3_mill), OUT_Ct_rpwd_cl=Ct_rpwd
           ,OUT_Ct_rpwd_cl_2= ifelse(Ct_rpwd<OUT_Cp3_mill_2,OUT_Cp3_mill_2,Ct_rpwd)
           ,OUT_Ct_rpwd_cl_3= ifelse(OUT_Ct_rpwd_cl_2=='NA',OUT_Cp3_mill,OUT_Ct_rpwd_cl_2)
           ,OUT_Cp3Ctrpwd_rat=OUT_Cp3_mill_2/Ct_rpwd
           ,OUT_Cp3Ctrpwd_rat=ifelse(OUT_Cp3Ctrpwd_rat>1.0,1.0,OUT_Cp3Ctrpwd_rat)
           ,OUT_Rrpw3=Rb3*(OUT_Cp3Ctrpwd_rat^0.3333), OUT_Rrpw3_RtRat=OUT_Rrpw3/Rb3
           ,OUT_Reyp3=(((g*Rb3*(Slope/100))^0.5)*((10^OUT_LSUB_DMM)/1000))/0.00000102
           ,OUT_LReyp3=log10(OUT_Reyp3)
           ,OUT_Shld_Px3= 0.5*((0.22*(OUT_Reyp3^(-0.6)))+(0.06*(10^(-7.7*(OUT_Reyp3^(-0.6))))))
           ,OUT_Shld_Px3=ifelse(is.na(OUT_LReyp3),NA,OUT_Shld_Px3)
           ,OUT_LShld_Px3=log10(OUT_Shld_Px3)
           ,OUT_LShld_Px3=ifelse(OUT_LReyp3<=1.4166,(-1.3991-(0.24419*OUT_LReyp3)),OUT_LShld_Px3)
           ,OUT_Shld_Px3=10^(OUT_LShld_Px3)
           ,OUT_Dcbf_fin=1000*((rho*g*OUT_Rrpw3*(Slope/100))/(OUT_Shld_Px3*(rhosed-rho)*g))
           ,OUT_LDcbf_fin=log10(OUT_Dcbf_fin), OUT_LRBS_fin=OUT_LSUB_DMM-OUT_LDcbf_fin)%>%
    select(SampleID,StationID,Date,IN_BR_PCT,IN_HP_PCT,IN_RC_PCT,IN_BL_PCT,IN_CB_PCT,IN_GC_PCT,IN_GF_PCT
           ,IN_SA_PCT,IN_FN_PCT,IN_WD_PCT,IN_OT_PCT,IN_BL_CB_GR_PCT,IN_SA_FN_PCT,IN_TotSubstrate_PCT
           ,IN_LSUB_DMM,IN_SUB_DMM,OUT_BR_PCT,OUT_HP_PCT,OUT_RC_PCT,OUT_BL_PCT,OUT_CB_PCT,OUT_GC_PCT
           ,OUT_GF_PCT,OUT_SA_PCT,OUT_FN_PCT,OUT_WD_PCT,OUT_OT_PCT,OUT_BL_CB_GR_PCT,OUT_SA_FN_PCT
           ,OUT_TotSubstrate_PCT,OUT_LSUB_DMM,OUT_SUB_DMM,IN_Xembed,IN_Vembed,OUT_Xembed,OUT_Vembed
           ,IN_LRBS1,IN_LRBS2,IN_LRBS_fin,OUT_LRBS1,OUT_LRBS2,OUT_LRBS_fin)
  return(output)
}

createsummary <- function(tab1, tab2, tab3, tab4) {
  #Arrange all the inputs by the SampleID
  tab1 <- arrange(tab1, SampleID)
  tab2 <- arrange(tab2, SampleID)
  tab3 <- arrange(tab3, SampleID)
  tab4 <- arrange(tab4, SampleID)
  
  #Then merge them together in order
  summary <- merge(tab1, tab2, by=c('SampleID', 'StationID', 'Date'))
  summary <- merge(summary, tab3, by=c('SampleID', 'StationID', 'Date'))
  summary <- merge(summary, tab4, by=c('SampleID', 'StationID', 'Date'))
  
  #One final rearrangement because it is out of order and remove a column
  summary <- arrange(summary, SampleID)
  summary <- summary[,-(31)]
  names(summary)[6] <- 'ReachLength'
  return(summary)
}

createINOUTsummary <- function(tab1, tab2, tab3) {
  #Arrange the inputs
  tab1 <- arrange(tab1, SampleID)
  tab2 <- arrange(tab2, SampleID)
  tab3 <- arrange(tab3, SampleID)
  
  #Merge the tables together based on their common columns then arrange by SampleID
  summary <- merge(tab1, tab2, by=c('SampleID', 'StationID', 'Date', 'BR_PCT', 'HP_PCT', 'BL_PCT', 'CB_PCT',
                                    'GC_PCT', 'GF_PCT', 'SA_PCT', 'FN_PCT', 'WD_PCT', 'OT_PCT', 'BL_CB_GR_PCT',
                                    'SA_FN_PCT', 'TotSubstrate_PCT', 'LSUB_DMM', 'SUB_DMM', 'RC_PCT'))
  summary <- merge(summary, tab3, by=c('SampleID', 'StationID', 'Date', 'VLW', 'VLW_msq', 'XBKF_H', 
                                       'Xwid', 'INC_H', 'Xembed', 'Vembed', 'BL_CB_GRmeanEmbed', 'BL_CB_GR_transectPCT'))
  s <- arrange(summary, SampleID)
  
  #The only thing this does is it reorders the columns
  summary2 <- cbind(s[,c("SampleID", "StationID", "Date", "SiteFlag", "ThalwegN",            
                         "ReachLength", "Interval", "Slope", "Xdepth", "Vdepth",              
                         "Xdep_Vdep", "Stacks_Equation", "AreaSum", "RP100", "BR_PCT",              
                         "HP_PCT", "RC_PCT", "BL_PCT", "CB_PCT", "GC_PCT",              
                         "GF_PCT", "SA_PCT", "FN_PCT", "WD_PCT", "OT_PCT",              
                         "BL_CB_GR_PCT", "SA_FN_PCT", "TotSubstrate_PCT", "LSUB_DMM", "SUB_DMM",             
                         "VLW", "VLW_msq", "XBKF_H", "Xwid", "INC_H",               
                         "Xembed", "Vembed", "BL_CB_GRmeanEmbed", "BL_CB_GR_transectPCT", "IN_BR_PCT",           
                         "IN_HP_PCT", "IN_RC_PCT", "IN_BL_PCT", "IN_CB_PCT", "IN_GC_PCT",
                         "IN_GF_PCT", "IN_SA_PCT", "IN_FN_PCT", "IN_WD_PCT", "IN_OT_PCT",           
                         "IN_BL_CB_GR_PCT", "IN_SA_FN_PCT", "IN_TotSubstrate_PCT", "IN_LSUB_DMM", "IN_SUB_DMM",          
                         "OUT_BR_PCT", "OUT_HP_PCT", "OUT_RC_PCT", "OUT_BL_PCT", "OUT_CB_PCT",          
                         "OUT_GC_PCT", "OUT_GF_PCT", "OUT_SA_PCT", "OUT_FN_PCT", "OUT_WD_PCT",          
                         "OUT_OT_PCT", "OUT_BL_CB_GR_PCT", "OUT_SA_FN_PCT", "OUT_TotSubstrate_PCT", "OUT_LSUB_DMM",        
                         "OUT_SUB_DMM", "IN_Xembed", "IN_Vembed", "OUT_Xembed", "OUT_Vembed")])
  return(summary2)
}

finalCharandOrder <- function(input) {
  #Change the Date to be a character and arrange by the SampleID
  input$Date <- as.character(input$Date)
  output <- input[order(input$SampleID),]
  return(output)
}

outputSummariesCSV <- function(TMDLSummary, MasterSummary, INOUT_MasterSummary, QAsummary) {
  #This outputs each of the 4 summaries to a csv file in the same folder as the script
  #The titles are the name of the summary followed by the date
  write.csv(TMDLSummary, file=paste("TMDLsummary_",Sys.Date(),".csv",sep=""), row.names=FALSE)
  write.csv(MasterSummary, file=paste("MasterSummary_",Sys.Date(),".csv",sep=""), row.names=FALSE)
  write.csv(INOUT_MasterSummary, file=paste("INOUT_MasterSummary_",Sys.Date(),".csv",sep=""), row.names=FALSE)
  write.csv(QAsummary, file=paste("QAsummary_",Sys.Date(),".csv",sep=""),row.names=FALSE)
}

sqlOutput <- function(channel, table, tableName) {
  #This function currently does not work as the data is outputted as csv files not into a SQL database
  sqlSave(channel, table, tablename=paste(tableName,Sys.Date(),sep=""), append=FALSE, rownames=FALSE, colnames=FALSE)
}


####################################
#Emma's Data Access Function

access_query_32_EVJ <- function(db_path, db_table, table_out ) {
  
  # Emma Jones adaptation (5/6/2019) of manotheshark Stack Overflow solution
  # https://stackoverflow.com/questions/13070706/how-to-connect-r-with-access-database-in-64-bit-window
  
  # variables to make values uniform
  sock_port <- 8642L
  sock_con <- "sv_con"
  ODBC_con <- "a32_con"
  #db_path <- "data/EDASxp_Family_090117.accdb"
  
  if (file.exists(db_path)) {
    
    # build ODBC string
    ODBC_str <- local({
      s <- list()
      s$path <- paste0("DBQ=", gsub("(/|\\\\)+", "/", path.expand(db_path)))
      s$driver <- "Driver={Microsoft Access Driver (*.mdb, *.accdb)}"
      s$threads <- "Threads=4"
      s$buffer <- "MaxBufferSize=4096"
      s$timeout <- "PageTimeout=5"
      paste(s, collapse=";")
    })
    
    # start socket server to transfer data to 32 bit session
    startSocketServer(port=sock_port, server.name="access_query_32", local=TRUE)
    
    # build expression to pass to 32 bit R session
    expr <- "library(svSocket)"
    expr <- c(expr, "library(RODBC)")
    expr <- c(expr, sprintf("%s <- odbcDriverConnect('%s')", ODBC_con, ODBC_str))
    expr <- c(expr, sprintf("if('%1$s' %%in%% sqlTables(%2$s)$TABLE_NAME) {%1$s <- sqlFetch(%2$s, '%1$s')} else {%1$s <- 'table %1$s not found'}", db_table, ODBC_con))
    expr <- c(expr, sprintf("%s <- socketConnection(port=%i)", sock_con, sock_port))
    expr <- c(expr, sprintf("evalServer(%s, %s, %s)", sock_con, table_out, db_table))
    expr <- c(expr, "odbcCloseAll()")
    expr <- c(expr, sprintf("close(%s)", sock_con))
    expr <- paste(expr, collapse=";")
    
    # launch 32 bit R session and run expressions
    prog <- file.path(R.home(), "bin", "i386", "Rscript.exe")
    system2(prog, args=c("-e", shQuote(expr)), stdout=NULL, wait=TRUE, invisible=TRUE)
    
    # stop socket server
    stopSocketServer(port=sock_port)
    
    # display table fields
    message("retrieved: ", table_out, " - ", paste(colnames(get(table_out)), collapse=", "))
  } else {
    warning("database not found: ", db_path)
  }
}

####################################