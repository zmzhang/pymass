PIT <- function(mzXMLFile){
  library(PITracer)
  PIC = findPureIonChromatogram(
    mzXMLFile,
    maxGap = 5,
    RangeEstimatedPPMTolerance=c(5, 50),
    MaxMassCalibrationTolerance=100,        #numeric or NA (NA: without calibration)
    SaturatedIntensity=Inf,
    SaturatedPPM=20)

  Chromatogram = PIC$Chromatogram 

  ChromInfo = t(sapply(Chromatogram, function(x) cbind(mz = x[which.max(x[,"intensity"]),"mz"],
                                                     mzmin = min(x[,"mz"]),
                                                     mzmax = max(x[,"mz"]),
                                                     rt = x[which.max(x[,"intensity"]),"rt"],
                                                     rtmin = min(x[,"rt"]),
                                                     rtmax = max(x[,"rt"]),
                                                     intensity = max(x[,"intensity"]),
                                                     length = diff(range(x[,"scan"]))+1
                                                     )))
  colnames(ChromInfo) = c("mz", "mzmin", "mzmax", "sc", "scmin", "scmax", "intensity", "length")


  ChromInfo_refine = ChromInfo[ChromInfo[,"length"] >= 10 & ChromInfo[,"intensity"]>300,]

  return(ChromInfo_refine)
}


XC <- function(mzXMLFile){
  library(xcms)
  xset <- xcmsSet(mzXMLFile)
  return(xset@peaks[,c(1,2,3,4,5,6,9)])
}
