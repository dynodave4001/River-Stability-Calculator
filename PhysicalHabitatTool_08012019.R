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
#Takes the data from the file denoted by its file path and name
accessData("C:/Users/kkn73995/Downloads/PhabTools/PHab_12-3-2018.mdb")

## Residual Pool Calculations
# Combine columns from the Reach and Thalweg tables 
thalweg <- merge(thalweg, reach[,c('SampleID','StationID','Date','ThalwegN','ReachLength','Interval'
                ,'Slope')], by=c('SampleID','StationID','Date')) %>% arrange(SampleID)

# Identify columns to keep in 10 thalweg dataframe by determining which to drop and taking the inverse
keepNames <- grep('10|11|12|13|14',names(thalweg),invert=TRUE,value=TRUE)
# Divide thalweg dataframe into 10 and 15 thalweg dataframes
thalweg_10 <- subset(thalweg, ThalwegN=='10', select=keepNames) 
thalweg_15 <- subset(thalweg, ThalwegN=='15')

# Thalweg Dataframe Formatting Prior to Calculations
thalweg_10.1 <- create1(thalweg_10)
thalweg_15.1 <- create1(thalweg_15)

# Append new calculations to thalweg_10 and create new dataframe to alter in Residual Pool Calculation loop
thalweg_10loop <- merge(thalweg_10,thalweg_10.1,by=c('SampleID','StationID','Date','ThalwegN'
                          ,'ReachLength','Interval','Slope')) %>% loopfunction()
thalweg_15loop <- merge(thalweg_15,thalweg_15.1,by=c('SampleID','StationID','Date','ThalwegN'
                          ,'ReachLength','Interval','Slope')) %>% loopfunction()

# delete unnecessary dataframes
rm(thalweg_15.1, thalweg_10.1) 

# Calculate summary Area metrics
sumArea10 <- sArea1(thalweg_10loop)
sumArea15 <- sArea1(thalweg_15loop)
sumArea <- rbind(sumArea10, sumArea15) %>% arrange(SampleID)

# Summarize the data of Thalweg
thalweg_summary <- merge(thalweg,sumArea, by=c('SampleID', 'StationID', 'Date')) %>% arrange(SampleID) %>%
                           mutate(RP100=AreaSum/ReachLength*100)
thalweg_summary.alt <- thalweg_summary[,c(1:3,154:163)]

# Clean up workspace before going to substrate calculations
rm(sumArea, sumArea10, sumArea15, thalweg_10, thalweg_10loop, thalweg_15, thalweg_15loop)

# Substrate Calculations
substrate[,4:108] <- lapply(substrate[,4:108],as.character) # Convert all station data to character format
substrate.wide1 <- createWide(substrate, 1)
substrate.wide2 <- createWide(substrate, 2)

# Calculate median particle size from substrate percentages
substrate_summary <- createSubSum(substrate.wide1)

# Calculate substrate percentages including OT, RC, & WD
substratePCT <- createSubPCT(substrate.wide2)

substrate_summary.alt <- arrange(createssalt(substratePCT, substrate_summary),SampleID)

## IN/OUT channel substrate calculations
# Subset substrate dataframe to separate inside and outermost channel observations
IN_substrate <- subset(substrate, select=grep('SampleID|StationID|Date|Left Middle|Middle|Right Middle'
                                              ,names(substrate),value=TRUE))
OUT_substrate <- subset(substrate, select=grep('Left Middle|Middle|Right Middle',names(substrate)
                                               ,invert=TRUE,value=TRUE))

# Transform to wide format for QA step
IN_substrate.wide1  <- createWide(IN_substrate,  1)
OUT_substrate.wide1 <- createWide(OUT_substrate, 1)

IN_substrate.wide2  <- createWide(IN_substrate,  2)
OUT_substrate.wide2 <- createWide(OUT_substrate, 2)

# Calculate median particle size from substrate percentages
IN_substrate_summary <- createISS(IN_substrate.wide1)
OUT_substrate_summary <- createOSS(OUT_substrate.wide1)

# Final Substrate Percentage report
IN_substratePCT <- createINPCT(IN_substrate.wide2)
OUT_substratePCT <- createOUTPCT(OUT_substrate.wide2)

IN_substrate_summary.alt  <- arrange( createINssalt( IN_substratePCT, IN_substrate_summary), SampleID)
OUT_substrate_summary.alt <- arrange(createOUTssalt(OUT_substratePCT,OUT_substrate_summary), SampleID)

# Combine all substrate summaries
INOUTsubstrate_summary.alt <- merge(substrate_summary.alt,IN_substrate_summary.alt)
INOUTsubstrate_summary.alt <- merge(INOUTsubstrate_summary.alt,OUT_substrate_summary.alt)
INOUTsubstrate_summary.alt <- arrange(INOUTsubstrate_summary.alt, SampleID)


#Clean up workspace before proceeding
rm(substrate_summary, substrate.wide1, substrate.wide2, substratePCT, IN_substrate_summary,
   IN_substrate.wide1, IN_substrate.wide2, IN_substrate, IN_substratePCT,
   IN_substrate_summary.alt, OUT_substrate_summary, OUT_substrate.wide1,
   OUT_substrate.wide2, OUT_substrate, OUT_substratePCT, OUT_substrate_summary.alt)


## Bankfull
XBKF_H <- apply(bankfull[, 4:14], 1, mean, na.rm=T)  # Average Bankfull Height
# na.rm=T adjusts denominator for missing data
XBKF_W <- apply(bankfull[,15:25], 1, mean, na.rm=T)  # Average Bankfull Width

## Wetted Width
Xwid <- apply(wet[, 4:24], 1, mean, na.rm=T) # Mean Wetted Width

## Incision
INC_H <- apply(incision[, 4:14], 1, mean, na.rm=T) # Average Incision Height

## Embeddedness
Xembed <- apply(embed[, 4:58], 1, mean, na.rm=T)
Vembed <- apply(embed[, 4:58], 1, sd, na.rm=T)

# IN/OUT Embeddedness- subset embed dataframe to separate inside and outermost channel observations
IN_embed <- subset(embed, select=grep('SampleID|StationID|Date|LM|M|RM',names(embed),value=TRUE))
OUT_embed <- subset(embed, select=grep('LM|M|RM',names(embed),invert=TRUE,value=TRUE))
# Embeddedness, inside and outside of channel
IN_Xembed <- apply(IN_embed[, 4:36], 1, mean, na.rm=T)
IN_Vembed <- apply(IN_embed[, 4:36], 1, sd, na.rm=T)
OUT_Xembed <- apply(OUT_embed[, 4:25], 1, mean, na.rm=T)
OUT_Vembed <- apply(OUT_embed[, 4:25], 1, sd, na.rm=T)


# Mean embeddedness of BL(XB+SB),CB,GR (GC+GF)
embed.long <- melt(embed,id.vars=c('SampleID','StationID','Date'),variable.name='Location'
                   ,value.name='EmbedPCT')
# sub <- unique(substrate.long1[,c(1:6)])# Retrieve particle count and totalcount from substrate
# sub2 <- unique(substrate.long2[,c(1:6)])

sub.wide <- select(substrate,-matches("2")) # Remove mid thalweg substrate measures
names(sub.wide) <- c('SampleID','StationID','Date'
                     ,'AL','BL','CL','DL','EL','FL','GL','HL','IL','JL','KL'
                     ,'ALM','BLM','CLM','DLM','ELM','FLM','GLM','HLM','ILM','JLM','KLM'
                     ,'AM','BM','CM','DM','EM','FM','GM','HM','IM','JM','KM'
                     ,'ARM','BRM','CRM','DRM','ERM','FRM','GRM','HRM','IRM','JRM','KRM'
                     ,'AR','BR','CR','DR','ER','FR','GR','HR','IR','JR','KR')

embed2 <- arrange(embed, SampleID)
sub.wide2 <- arrange(sub.wide,SampleID)

embed.sub1 <- createes1(embed2, sub.wide2)

embed.sub1.5 <- createes1.5(embed2, sub.wide2)

embed.sub2 <- createes2(embed2, sub.wide2)

BL_CB_GRmeanEmbed.df <- createBCG(embed2, sub.wide2, embed.sub1, embed.sub1.5, embed.sub2)

## Wood
attach(reach)
wood_summary.alt <- wsALT(wood, BL_CB_GRmeanEmbed.df)
wood_summaryINOUT <- wsIO(wood_summary.alt)

# Clean up workspace
rm(embed.long, embed.sub1, embed.sub1.5, embed.sub2, sub.wide, sub.wide2, embed2, BL_CB_GRmeanEmbed.df)


## QA STEP: Search previous parameters for missing data
thalweg_10QA <- create10QA(thalweg,keepNames)

QA_thalweg <- createQA_thal(thalweg,thalweg_10QA)

QA_bankfull <- createQAgen(bankfull, 'bankfull')

QA_wet <- createQAgen(wet, 'wet')

QA_incision <- createQAgen(incision, 'incision')

QA_embed <- createQAgen(embed, 'embed')

QA_substrate <- createQA_sub(substrate)

QAsummary <- createQAsummary(QA_thalweg, QA_bankfull, QA_wet, QA_incision, QA_embed, QA_substrate)


# Clean Workspace
rm(QA_bankfull, QA_embed, QA_incision, QA_substrate, QA_thalweg, QA_wet, thalweg_10QA, thalweg_summary)


# Final Data Summary
summary <- createsummary(QAsummary[,c(1:3,16)], thalweg_summary.alt,substrate_summary.alt, wood_summary.alt)
INOUTsummary <- createINOUTsummary(summary,INOUTsubstrate_summary.alt,wood_summaryINOUT[,-(4)])

# Establish constants for additional summary parameters
rho=998; rhosed=2650; g=9.807

# Change Date format and sort dataframe by SampleID for final output
TMDLSummary <- createTMDLSummary(summary) %>% finalCharandOrder()

MasterSummary <- createMasterSummary(summary)

QAsummary <- updateQASummary(QAsummary, MasterSummary) %>% finalCharandOrder()

MasterSummary <- updateMasterSummary(MasterSummary) %>% finalCharandOrder()

INOUT_MasterSummary <- createINOUTMasterSummary(INOUTsummary) %>% finalCharandOrder()


#Final Cleanup
rm(bankfull, embed, IN_embed, incision, INOUTsubstrate_summary.alt, INOUTsummary, OUT_embed, reach,
   substrate, substrate_summary.alt, summary, thalweg, thalweg_summary.alt, wet, wood, wood_summary.alt,
   wood_summaryINOUT)

# Output summary data frame to CSV
outputSummariesCSV(TMDLSummary, MasterSummary, INOUT_MasterSummary, QAsummary)

# channel <- "C:/Users/kkn73995/Downloads/PhabTools/"
# #Output summary data frame to Access
# #First ensure all previous table versions (or tables with identical names) are deleted of Access database
# sqlOutput(channel, TMDLSummary, "TMDLsummary_")
# sqlOutput(channel, MasterSummary, "MasterSummary_")
# sqlOutput(channel, INOUT_MasterSummary, "INOUT_MasterSummary_")
# sqlOutput(channel, QAsummary, "QAsummary_")
# odbcClose(channel)