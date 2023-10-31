# from the R package camtrapR

assessTemporalIndependence <- function(intable,
                                       deltaTimeComparedTo,
                                       columnOfInterest,     # species/individual column
                                       cameraCol,
                                       camerasIndependent,
                                       stationCol,
                                       minDeltaTime,
                                       removeNonIndependentRecords = TRUE)
{
  # check if all Exif DateTimeOriginal tags were read correctly
  if(any(is.na(intable$DateTimeOriginal))){
    which.tmp <- which(is.na(intable$DateTimeOriginal))
    if(length(which.tmp) == nrow(intable)) stop("Could not read any Exif DateTimeOriginal tag at station: ", paste(unique(intable[which.tmp, stationCol])), " Consider checking for corrupted Exif metadata.")
    warning(paste("Could not read Exif DateTimeOriginal tag of", length(which.tmp),"image(s) at station", paste(unique(intable[which.tmp, stationCol]), collapse = ", "), ". Will omit them. Consider checking for corrupted Exif metadata. \n",
                  paste(file.path(intable[which.tmp, "Directory"],
                                  intable[which.tmp, "FileName"]), collapse = "\n")), call. = FALSE, immediate. = TRUE)
    intable <- intable[-which.tmp ,]
    rm(which.tmp)
  }
  
  # prepare to add time difference between observations columns
  intable <- data.frame(intable,
                        delta.time.secs  = NA,
                        delta.time.mins  = NA,
                        delta.time.hours = NA,
                        delta.time.days  = NA)
  
  # introduce column specifying independence of records
  if(minDeltaTime == 0) {
    intable$independent <- TRUE    # all independent if no temporal filtering
  } else {
    intable$independent <- NA
  }
  
  
  for(xy in 1:nrow(intable)){     # for every record
    
    # set independent = TRUE if it is the 1st/only  record of a species / individual
    
    if(camerasIndependent == TRUE){
      if(intable$DateTimeOriginal[xy]  == min(intable$DateTimeOriginal[which(intable[, columnOfInterest] == intable[xy, columnOfInterest] &
                                                                             intable[, stationCol]       == intable[xy, stationCol] &
                                                                             intable[, cameraCol]        == intable[xy, cameraCol]) ])){    # cameras at same station assessed independently
        intable$independent[xy]       <- TRUE
        intable$delta.time.secs[xy]   <- 0
      }
    } else {
      if(intable$DateTimeOriginal[xy]  == min(intable$DateTimeOriginal[which(intable[, columnOfInterest] == intable[xy, columnOfInterest] &
                                                                             intable[, stationCol]       == intable[xy, stationCol]) ])){
        intable$independent[xy]       <- TRUE
        intable$delta.time.secs[xy]   <- 0
      }
    }
    
    if(is.na(intable$delta.time.secs[xy])) {   # if not the 1st/only record, calculate time difference to previous records of same species at this station
      
      if(deltaTimeComparedTo == "lastIndependentRecord"){
        
        if(camerasIndependent == TRUE){
          which_time2 <- which(intable[, columnOfInterest]       == intable[xy, columnOfInterest] &    # same species/individual
                                 intable[, stationCol]              == intable[xy, stationCol] &          # at same station
                                 intable[, cameraCol]               == intable[xy, cameraCol] &           # at same camera
                                 intable$independent                == TRUE &                             # independent (first or only record of a species at a station)
                                 intable$DateTimeOriginal           <  intable$DateTimeOriginal[xy])      # earlier than record xy
        } else {
          which_time2 <- which(intable[, columnOfInterest]       == intable[xy, columnOfInterest] &
                                 intable[, stationCol]             == intable[xy, stationCol] &
                                 intable$independent               == TRUE &
                                 intable$DateTimeOriginal          <  intable$DateTimeOriginal[xy])
        }
      }  else {
        if(camerasIndependent  == TRUE){
          which_time2 <- which(intable[, columnOfInterest]       == intable[xy, columnOfInterest] &
                                 intable[, stationCol]             == intable[xy, stationCol] &
                                 intable[, cameraCol]              == intable[xy, cameraCol] &
                                 intable$DateTimeOriginal          <  intable$DateTimeOriginal[xy])
        } else {
          which_time2 <- which(intable[, columnOfInterest]       == intable[xy, columnOfInterest] &
                                 intable[, stationCol]             == intable[xy, stationCol] &
                                 intable$DateTimeOriginal          <  intable$DateTimeOriginal[xy])
        }
      }
      
      # time difference to last (independent) record
      diff_tmp <- min(na.omit(difftime(time1 = intable$DateTimeOriginal[xy],            # delta time to last independent record
                                       time2 = intable$DateTimeOriginal[which_time2],
                                       units = "secs")))
      
      # save delta time in seconds
      intable$delta.time.secs[xy] <-  diff_tmp
      if(intable$delta.time.secs[xy] >= (minDeltaTime * 60) | intable$delta.time.secs[xy] == 0){
        intable$independent[xy] <- TRUE
      } else {
        intable$independent[xy] <- FALSE
      }
      
    }   # end   if(intable$DateTimeOriginal[xy] == min(...)} else {...}
  }     # end for(xy in 1:nrow(intable))
  
  
  if(removeNonIndependentRecords){
  # keep only independent records
  outtable <- intable[intable$delta.time.secs >= (minDeltaTime * 60) |
                        intable$delta.time.secs == 0,]
  } else {
    outtable <- intable
  }
  
  return(outtable)
}
