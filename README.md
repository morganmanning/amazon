Timeframe:
- Zabalo: 2018-07-03 - 2018-10-09
- Siekopai: 2022-07-14 - 2022-11-24
- Siona: 2022-07-27 - 2022-11-25
- Sinangoe: 2022-08-02 - 2022-12-01



Data:
- 'CollaredPeccary.csv': detection matrix for Pecari tajacu by Michael Esbach pulled from the Zabalo Google Drive
- 'Deer.csv': detection matrix for Mazama americana by Michael Esbach pulled from the Zabalo Google Drive
- 'Paca.csv': detection matrix for Cuniculus paca by Michael Esbach pulled from the Zabalo Google Drive
- 'HuntingData2018.csv': records of what was hunted, where it was hunted, and how many were hunted surrounding Zabalo
- Any .csv with 'Stations' in the name: data on camera trap location and operation
- 'RecordTable2018.csv': raw camera trap data (non-independent data not removed)
- Any .csv with 'IndependentRecords': camera trap data with only independent data (filtered with 'assessTemporalIndependence.R'
- Any .csv with 'StationsFormatted' or 'RecordsFormatted': data passed through 'formattingRawDates.R' to standardize dates
- all data within Zabalo/FPDistData folder: from https://github.com/lbeaudrot/Elevational-Shifts/tree/master 
- within Zabalo/FPDistData folder: Pantheria data from https://esapubs.org/archive/ecol/E090/184/metadata.htm 



Scripts:
- 'assessTemporalIndependence.R': function to assess temporal independence between camera trap detections (written by Michael Esbach)
- 'formattingRawDates.R': pulls in raw camera trap records (e.g. ZABIndependentRecords.csv) and camera station (e.g. ZABStations.csv) information to standardize date format
- 'rawToDetectionDF.R': imports camera trap data after date standardization and outputs detection history matrices per species using 'camtrapR' (imported data must include camera operability)
- 'DataExploration.R': script looking at 'Deer.csv', 'Paca.csv', and 'CollaredPeccary.csv' and running preliminary models with site covariates (siteCov2018.csv)
- 'allSpeciesOccupancy.R': mostly self-driving script where you input the names of detection matrices (in .csv format) and covariate data and it outputs detection plots, the best models for each species, clumps based on optimal grouping function (optimalClumping.R), and plots predictions with the selected occupancy/detection covariates (only really useful for Zabalo data without adaptation; meant to replace separate scripts for each animal)
- 'allSpeciesPlotting.R': uses 'allSpeciesOccupancy.R' .RData objects to plot predictions across range of selected covariates for each species established in 'allSpeciesOccupancy.R'
- 'basicOccupancyExample.R': exactly what it sounds like... the blueprint for occupancy modeling and plotting
- 'brocketOccupancy.R'/'pacaOccupancy.R'/'peccaryOccupancy.R': runs occupancy models and plots predictions for stated species; can now use 'allSpeciesOccupancy.R' instead
- 'convergence_test.R': quick function to test convergence of null model when you input the detection matrix and grouping factor
- 'dataFormatting.R': scales all continuous covariates, plots all hunting events versus camera locations, and creates the hunting intensity variable
- 'FunctionalPhylogeneticDistance.R': uses Lydia Beaudrot methods (https://github.com/lbeaudrot/Elevational-Shifts/tree/master) to calculate functional diversity between collared peccary, lowland paca, and red brocket
- 'hunting.R': plotting hunting occurrences versus camera locations, looks at initial method of assigning each hunt to a camera, and eventually weights hunts for every camera (redundant)
- 'optimalClumping.R': creates function to determine the best clumping/grouping factor to reduce the margin of error
- 'occupancySGE.R': input names of detection .csv matrices and outputs occupancy model predictions
- 'SinangoeVsZabalo.Rmd': report looking at the differences/similarities in diversity indices and occupancy between Sinangoe and Zabalo
  



***Chapter 3 Schedule***
| Week | Date | Task Due | Meeting Topic |
| --- | --- | --- | --- |
| 1 | February 8 | Pull questions from literature | Yay or nay questions/question topics |
| 2 | February 15 | First very rough draft to MA | Discuss draft |
| 3 | February 22 | Second less rough draft to MA | Discuss draft |
| 4 | February 29 | Possible third least rough draft to MA | Discuss draft |
| 5 | March 7 | First round of edits to GR | |
| 6 | March 14 | Revise per GR edits and pass to MA | Discuss GR edits |
| 7 | March 21 | Revise per MA edits | |
| 8 | March 28 | Second round of edits to GR and MA | |
| 9 | April 4 | Revise | Discuss final GR edits |
| 10 | April 11 | Final touches and formatting | Final sign off on format |
| 11 | April 18 | Give the survey | |
| 12 | April 25 | Remind about the survey | | 







