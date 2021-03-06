# LocalAdaptation_RangeShifts
Source code for the AoB paper "Simulating the effects of local adaptation and life history on the ability of plants to track climate shifts"

Emily V. Moran, 2020
University of California, Merced

Similar to range_shift_model, but with improved dispersal calculations and allowing for multiple life histories besides trees.

##Locally adapted tree model

TreeMod_LA_Setup.r - Intializes the basic locally adapted tree model, setting up landscape, adding in first individuals, calculating dispersal probabilities for all seed and pollen dispersal kernels, creating a seed dispersal image (if desired) and saving an output file "TreeMod_LA_Setup.RData"

TreeMod_LA_40_Stabil.r - Inputs "TreeMod_LA_Setup.RData" and runs 40 m average dispersal model for 200 years in original stable climate. Outputs a map of allele distribution (if desires), calculated occupancy, and saves output file "TreeMod_LA_40_Stabil.RData".  Can be changed to 80 m dispersal version by altering line commenting.

TreeMod_LA_40_Shift.r - Inputs "TreeMod_LA_40_Stabil.RData" and runs 40 m average dispersal model for 20 years of climate change. Outputs a map of allele distribution (if desires), calculated occupancy, and saves output file "TreeMod_LA_40_Shift.RData".  Can be changed to 80 m dispersal version by altering line commenting.

TreeMod_LA_40_PostShift.r - Inputs "TreeMod_LA_40_Shift.RData" and runs 40 m average dispersal model for 200 years of new stable climate. Outputs a map of allele distribution (if desires), calculates occupancy and saves output file "TreeMod_LA_40_PostShift.RData".  Can be changed to 80 m dispersal version by altering line commenting.

##Plastic tree model

TreeMod_P_Setup.r - Same as TreeMod_LA_Setup.r, except for plastic species. Saves output file "TreeMod_P_Setup.RData"

TreeMod_P_40_Stabil.r - Same as TreeMod_LA_40_Stabil.r except for plastic species. Inputs "TreeMod_P_Setup.RData". Saves output file "TreeMod_P_40_Stabil.RData".  

TreeMod_P_40_Shift.r - Same as TreeMod_LA_40_Shift.r except for plastic species. Inputs "TreeMod_P_40_Stabil.RData". Saves output file "TreeMod_P_40_Shift.RData".  

TreeMod_P_40_PostShift.r - Same as TreeMod_LA_40_PostShift.r except for plastic species. Inputs "TreeMod_P_40_Shift.RData". Saves output file "TreeMod_P_40_PostShift.RData".  

##Locally adapted perennial model

PerenMod_LA_Setup.r - Same as TreeMod_LA_Setup.r, except with perennial's demographic parameters. Saves output file "PerenMod_LA_Setup.RData"

PerenMod_LA_20_Stabil.r - Same as TreeMod_LA_20_Stabil.r except with perennial's demographic parameters, and with a 100 year stabilization period. Inputs "PerenMod_LA_Setup.RData". Saves output file "PerenMod_LA_20_Stabil.RData".  

PerenMod_LA_20_Shift.r - Same as TreeMod_LA_20_Shift.r except with perennial's demographic parameters. Inputs "PerenMod_LA_20_Stabil.RData". Saves output file "PerenMod_LA_20_Shift.RData".  

PerenMod_LA_20_PostShift.r - Same as TreeMod_LA_20_PostShift.r except with perennial's demographic parameters, and with a 100 year stabilization period. Inputs "PerenMod_LA_20_Shift.RData". Saves output file "PerenMod_LA_20_PostShift.RData".  

##Plastic perennial model

PerenMod_P_Setup.r - Same as PerenMod_LA_Setup.r, except for plastic species. Saves output file "PerenMod_P_Setup.RData"

PerenMod_P_20_Stabil.r - Same as PerenMod_LA_20_Stabil.r except for plastic species. Inputs "PerenMod_P_Setup.RData". Saves output file " PerenMod_P_20_Stabil.RData".  

PerenMod_P_20_Shift.r - Same as PerenMod_LA_20_Shift.r except for plastic species. Inputs "PerenMod_P_20_Stabil.RData". Saves output file "PerenMod_P_20_Shift.RData".  

PerenMod_P_20_PostShift.r - Same as PerenMod_LA_20_PostShift.r except for plastic species. Inputs "PerenMod_P_20_Shift.RData". Saves output file " PerenMod_P_20_PostShift.RData".  

##Locally adapted annual model

AnnMod_LA_Setup.r - Same as TreeMod_LA_Setup.r, except with annual's demographic parameters. Saves output file "AnnMod_LA_Setup.RData"

AnnMod_LA_20_Stabil.r - Same as TreeMod_LA_20_Stabil.r except with annual's demographic parameters, and with a 20 year stabilization period. Inputs "AnnMod_LA_Setup.RData". Saves output file "AnnMod_LA_20_Stabil.RData".  

AnnMod_LA_20_Shift.r - Same as TreeMod_LA_20_Shift.r except with annual's demographic parameters. Inputs "AnnMod_LA_20_Stabil.RData". Saves output file "AnnMod_LA_20_Shift.RData".  

AnnMod_LA_20_PostShift.r - Same as TreeMod_LA_20_PostShift.r except with annual's demographic parameters, and with a 20 year stabilization period. Inputs "AnnMod_LA_20_Shift.RData". Saves output file "AnnMod_LA_20_PostShift.RData".  

##Plastic annual model

AnnMod_P_Setup.r - Same as AnnMod_LA_Setup.r, except for plastic species. Saves output file " AnnMod_P_Setup.RData"

AnnMod_P_20_Stabil.r - Same as P AnnMod_LA_20_Stabil.r except for plastic species. Inputs " AnnMod_P_Setup.RData". Saves output file "AnnMod_P_20_Stabil.RData".  

AnnMod_P_20_Shift.r - Same as AnnMod_LA_20_Shift.r except for plastic species. Inputs "AnnMod_P_20_Stabil.RData". Saves output file " AnnMod_P_20_Shift.RData".  

AnnMod_P_20_PostShift.r - Same as AnnMod_LA_20_PostShift.r except for plastic species. Inputs "AnnMod_P_20_Shift.RData". Saves output file "AnnMod_P_20_PostShift.RData".  

##Model variants

TreeMod_LA_40_Fst2_Shift.r - Same as TreeMod_LA_40_Shift.r, except that climate gradient shifts at twice the rate (40 m/yr)

TreeMod_LA_Lg_Setup.r - Same as TreeMod_LA_Setup.r, except that climate gradient is 4x as long and climate bands are wider.

TreeMod_LA_Lg_40_Shift.r - The climate shift period file corresponding to TreeMod_LA_Lg_Setup.r

PerenMod_PW_Setup.r - Setup file for the "broadly plastic perennial" scenarios.

PerenMod_LA_Fit_Setup.r - Setup file for the locally adapted perennial scenario in which the effects of sub-optimal climate on annual demographic rates is halved.

