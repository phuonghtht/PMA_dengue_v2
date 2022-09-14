#'Description of the method
#' @param a First number
#' @param b second number
#' @return sum of a and b
#' @export
tong <- function(a,b){
  return(a+b)
}


#' Titer computtation
#' @description This is a function that performs titer computation for an antigen, edited from "robustComputeTiter{aatk"}
#' @return EC50 value; if luminescence at dilution 20 < top*0.5, EC50 is assigned into 2 groups:10 or 20.
#' @export
robustComputeTiter.new <-function (dilutionvalues, log2.transform = FALSE, top = 65535,
                                   bottom = 3000, monotonicity.cutoff = 0.5, log2.output = TRUE,
                                   generate.plot = FALSE)
{
  dilutionvalues <- dilutionvalues[order(dilutionvalues$DilutionSteps), ]
  if (log2.transform) {
    dilutionvalues$Value = log2(dilutionvalues$Value)
    top = log2(top)
    bottom = log2(bottom)
  }
  high = sum(dilutionvalues$Value > 0.5 * top, na.rm = T)# is there any luminescence signal > EC50 ( 32767.5)
  if (high == nrow(dilutionvalues)) {
    EC50 = max(dilutionvalues$DilutionSteps)# assign upper value
    names(EC50) = "EC50"
    if (generate.plot)
      plot(dilutionvalues$DilutionSteps, dilutionvalues$Value,
           log = "y")
  }


  if (high == 0){
    if (sum(!is.na(dilutionvalues$Value))==0){#
      EC50=0
      names(EC50) = "EC50"
    }else{
      if (mean(dilutionvalues$Value[1:2], na.rm = T)<=20000){
        EC50 = min(dilutionvalues$DilutionSteps)/2 # assign lower value : luminescence<20000
        names(EC50) = "EC50"
        if (generate.plot)
          plot(dilutionvalues$DilutionSteps, dilutionvalues$Value,
               log = "y")
      }else{
        EC50 = min(dilutionvalues$DilutionSteps) # assign lower value: 20000<luminescence<32767.5
        names(EC50) = "EC50"
        if (generate.plot)
          plot(dilutionvalues$DilutionSteps, dilutionvalues$Value,
               log = "y")
      }
    }
  }



  if (high > 0 & high < nrow(dilutionvalues)) {
    mean.values = aggregate(dilutionvalues$Value, list(dilutionvalues$DilutionSteps),
                            mean, na.rm = TRUE)
    number.pairs.required = sum((nrow(mean.values) - 1):1) *
      monotonicity.cutoff
    number.monotonic = sum(as.dist(outer(mean.values[, 2],
                                         mean.values[, 2], "<")), na.rm = TRUE)
    if (number.monotonic >= number.pairs.required) {
      EC50 = tryCatch(computeTiter(dilutionvalues = dilutionvalues,
                                   top = top, bottom = bottom, generate.plot = generate.plot),
                      error = function(e) {
                        message("Fit did not converge!")
                        NA
                      })
    }
    else {
      EC50 = NA
    }
  }
  if (log2.output)
    return(log2(EC50))
  else return(EC50)
}



#'Creating an agArraySet object
#' @description  This function is edited from the function "createAgArraySet", package "aatk".agArrayList is kept as original regardless of Flag label.
#' @return An agArraySet object
#' @export
createAgArraySet.new<- function (metadata, agArrayList, Value = "Median", Func = robustComputeTiter.new,
                                 parallel = TRUE, pos.contr.label = NULL, exclude.trays = NULL,
                                 exclude.antigens = NULL)
{
  if (parallel & requireNamespace("parallel")) {
    no_cores = max(1, parallel::detectCores() - 1)
    mylapply = function(...) parallel::mclapply(mc.cores = no_cores,
                                                ...)
    message("Using parallel computation...")
  }
  else {
    mylapply = lapply
    message("Not using parallel computation...")
  }
  message("Reading in arrays...")
  agArrayListF = agArrayList
  valcols = c("Median", "Mean", "SD", "Median.minus.B", "Mean.minus.B",
              "SignalNoiseRatio")
  # agArrayListF = lapply(agArrayListF, function(x) {
  #   x$data[!x$data$Flags %in% c(1, 3), valcols] = NA
  #   x
  # }) # this step is not nescessary!Be cautious when reading the slides!
  message("Computing titers...")
  agArrayTitersList = mylapply(agArrayListF, function(singleArray) {
    sapply(spltSample(singleArray, Value = Value), function(x) fitAntigens(x,
                                                                           Func = Func))
  })
  agArrayTitersR = t(Reduce(cbind, agArrayTitersList))
  if (sum(!rownames(metadata) %in% rownames(agArrayTitersR)) >
      0) {
    stop("Worklistfile does not match csv.files! Please check if the names of your scanfiles match the `filenames` column in your worklistfile.")
  }
  agArrayTitersRaw = agArrayTitersR[rownames(metadata), ]
  agTiters = apply(agArrayTitersRaw[rownames(metadata), ],
                   2, function(x) x + log2(metadata$startdilution))
  colnames(agTiters) = gsub("EC50", "", colnames(agTiters))
  output = list(meta = metadata, rawdata = lapply(agArrayList,
                                                  function(x) x[["data"]]), rawtiters = agArrayTitersRaw,
                titers = agTiters)
  class(output) = "agArraySet"
  if (!is.null(exclude.antigens)) {
    message("Removing antigen(s) ", paste(exclude.antigens,
                                          collapse = " ,"))
    selAntigens = antigens(output)[!antigens(output) %in%
                                     exclude.antigens]
    output = selectAntigens(output, selAntigens)
  }
  if (!is.null(exclude.trays)) {
    output = subset(output, which(!output$meta$tray %in%
                                    exclude.trays))
    message("Removing tray(s) ", paste(exclude.trays, collapse = " ,"))
  }
  if (!is.null(pos.contr.label)) {
    message("Performing normalization...")
    output = normalizeAgArrayLoess(output, pos.contr.label = pos.contr.label)
  }
  return(output)
}


#'Test SD
#' @description  This function is to test SD of possitive control.
#' @return within 2SD
#' @export
isInRange <- function(numArray) {
  result <- NULL
  arraySize <- length(numArray)
  for (i in 1 : arraySize) {
    mean <- mean(numArray[-i])
    sd <- sd(numArray[-i])
    if (numArray[i] > mean + 2 * sd || numArray[i] < mean - 2 * sd) {
      result <- c(result,FALSE)
    } else {
      result <- c(result,TRUE)
    }
  }
  return(result)}


#'Adjusting Titers
#' @description  This function is to adjust titers to the limit of detection from log2(10) to log2(1280)
#' @return Titers which have no value below log2(10) or above log2(1280)
#' @param col_vec Columns that indicate titers
#' @param min_thres minimum threshold
#' @param max_thres maximum threshold
#' @param ...
#' @export
AdjustTiter <- function(df, col_vec, min_thres = NA, val_min = log2(10), max_thres = NA, val_max = log2(1280)){
  if (!is.na(min_thres)){
    for (idx in col_vec){
      df[[idx]][df[[idx]] <= min_thres] <- val_min
      df[[idx]][is.na(df[[idx]])] <- NA
    }
  }

  if (!is.na(max_thres)){
    for (idx in col_vec){
      df[[idx]][df[[idx]] >= max_thres] <- val_max
    }
  }

  return (df)
}


#'Transforming ELISA_PanbioUnit into ELISA Status c("Neg","Pos","Uncertain")
#' @return df with an extra col called " E_Status"
#' @param df data.frame
#' @param var a vector of ELISA result in PanbioUnit format
#' @export
E_Status<- function(df,var){
  df$E_Status<-NA
  for (i in 1: length(df[[var]])){
    if(df[[var]][i] < 9){
      df$E_Status[i] <- "Neg"
    } else if (df[[var]][i] > 11){
      df$E_Status[i] <- "Pos"
    }else{
      df$E_Status[i] <- "Uncertain"
    }
  }
  return(df)
}




#'Check out whether a sample is homotypic infection or not?
#' @return df with an extra col called "NoSerotyeInfected"
#' @param df data.frame
#' @param col_vec vector of Ags variables
#' @param cutoff default is 5
#' @export
IsHomotypic <- function(df,col_vec, cutoff=5){
  df$NoSerotyeInfected <-NA
  if(sum(is.na(df[,col_vec]))>=1){ # count NAs value
    stop("Check if there is NAs value!Then crossing out all of them!")
  }else
    for( i in (1: length(df$sampleID))){
      if(sum(df[i,col_vec]<=cutoff)==(length(col_vec)-1)){
        df$NoSerotyeInfected[i]<-1
      }
      else if (sum(df[i,col_vec]<=cutoff)==(length(col_vec)-2)){
        df$NoSerotyeInfected[i]<-2
      }

      else if (sum(df[i,col_vec]<=cutoff)==(length(col_vec)-3)){
        df$NoSerotyeInfected[i]<-3
      }
      else{
        df$NoSerotyeInfected[i]<-"more than 3"
      }
    }

  return(df)
}

# > IsHomotypic
# function(df,col_vec, cutoff=5){
#   df$homotypic <-NA
#   for( i in 1: length(df$sampleID)){
#     if(sum(df[i,col_vec]<=cutoff)==7){
#       df$homotypic[i]<-"Yes"
#     }else{
#       df$homotypic[i]<-"No"
#     }
#   }
#   return(df)
# }

#'Assigning PMA status
#' @return df with an extra col called "cutoff_5"
#' @param df data.frame
#' @param col_vec vector of Ags variables
#' @param cutoff default is 5
#' @export
PMAStatus<- function( df, col_vec, cutoff=5){
  df$cutoff_5<- NA
  if(sum(is.na(df[,col_vec]))>=1){
    stop("Make sure having no NA value!")
  }else
    for( i in (1:length(df$sampleID))){
      if(sum(df[i,col_vec]<=cutoff)==length(col_vec)){
        df$cutoff_5[i] <- "Neg"
      }else{
        df$cutoff_5[i] <- "Pos"
      }
    }
  return(df)
}

#'Assigning serotype of homotypic infection samples.
#' @description Which antigen has titer that is greater than cutoff ( default = 5)
#' is assigned such antigen as an infecting serotype
#' @return df with an extra col called pred.serotype
#' @param df data.frame
#' @param col_vec vector of Ags variables
#' @param cutoff default is 5
#' @export
AssignSerotypeHomo <- function(df,col_vec,cutoff = 5){
  df$pred.serotype <- NA
  if(sum(is.na(df[,col_vec]))>=1){ # count NAs value
    stop("Check if there is NAs value!Then crossing out all of them!")
  }else
    for(i in 1: length(df$sampleID)){
      idx <- which(df[i,col_vec] >=cutoff) # idx of col_vec
      idx <- col_vec[1] - 1 + idx # idx in a df
      df$pred.serotype[i] <- colnames(df[idx])
    }
  return(df)
}
