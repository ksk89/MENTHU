#' calculateMENTHU2
#'
#' @param inData A DNA sequence/window/context
#' @param cutSite The expected index of the nucleotide to the left of the cut; -1 if cut site is in middle of wildtype sequence
#' @param weight A weight factor for caculating pattern score; default is 20, as in Bae paper
#'
#' @return


#For each sequence in the list, calculate the slope competition
calculateMENTHU2 <- function(inData, cutSite = -1, weight = 20){
  
  # If the cut site isn't specified, assume it is in the middle of the input
  if(cutSite == -1){
    cutSite <- floor(nchar(inData)/2)
  }
  
  # Get ordered pattern scores for every deletion pattern in the wildtype sequence
  patternScoreDF <- calculateBae(inData, cutSite, weight)
  patternScoreDF <- patternScoreDF[order(-patternScoreDF$patternScore),]
  
  # Filter out MHs shorter than 3 from patternScoreDF
  threePlus <- patternScoreDF[which(nchar(patternScoreDF$microhomology) >= 3),]
  
  # Filter out entries with more than 5bp separation between MHs
  critOne <- threePlus[which(threePlus$delLength - nchar(patternScoreDF$microhomology) <= 5),]
  
  # Generate output values
  if (nrow(critOne) == 0) { # If no entries comply with criterion
    ratio <- 0
    fShift <- "NA"
  } else if (nrow(critOne) == 1) { # If only one entry complies with criterion
    ratio <- "Max"
    fShift <- (if(critOne[1, 3] %% 3 == 0){"No"} else {"Yes"})
  } else { # Calculation of MENTHU score 2.0
    ratio <- critOne[1,4]/critOne[2,4]
    fShift <- (if(critOne[1, 3] %% 3 == 0){"No"} else {"Yes"})
  }
  
  # Organize output values as data frama
  outFrame <- data.frame(seq = inData,
                          MENTHUscore = ratio,
                          frameShift = fShift,
                          stringsAsFactors = FALSE)
  
  return(outFrame)
}