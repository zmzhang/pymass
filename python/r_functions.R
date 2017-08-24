PIT <- function(mzXMLFile, maxGap,CalibrationTolerance, plen, intensity){
  library(PITracer)
  PIC = findPureIonChromatogram(
    mzXMLFile,
    maxGap = maxGap,
    RangeEstimatedPPMTolerance=c(5, 50),
    MaxMassCalibrationTolerance=CalibrationTolerance,        #numeric or NA (NA: without calibration)
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


  ChromInfo_refine = ChromInfo[ChromInfo[,"length"] >= plen & ChromInfo[,"intensity"]>intensity,]

  return(ChromInfo_refine)
}


XC <- function(mzXMLFile, w1, w2, snr, intensity){
  library(xcms)
  xraw <- xcmsRaw(mzXMLFile)
  peaks = findPeaks.centWave(xraw, peakwidth=c(w1,w2), snthresh=snr, prefilter=c(1,intensity))
  return(peaks[,c(1,2,3,4,5,6,9)])
}
