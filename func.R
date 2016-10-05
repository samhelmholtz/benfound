# a function to get the first digit of an non-zero integer
getFirstDigit <- function(k){
    stringSplit <- (head(strsplit(as.character(abs(k)),'')))[[1]]
    if(stringSplit[1]=="0" && length(stringSplit)>2){
        return(as.numeric(stringSplit[3]))
    }
    return(as.numeric(stringSplit[1]))
}

## Assumes k is an integer > 10
getTwoSignificantDigits <- function(k){
    return(as.integer(substring(gsub("\\.","",as.character(abs(k))),1,2)))
}

## Scales the input vector to the range [0,1]
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

getNumberOfDigitsInPositiveInteger <- function(positiveInteger, integerBase){
    return(floor(log(positiveInteger,base=integerBase)+1))
}

## Returns the vector of first digits from values
getFirstDigits <- function(values){    
    result <- c()
    for(v in values){
        result <- c(result, getFirstDigit(v))
    }
    return(result)
}

## Filters zeros from the given vector
stripZeros <- function(values){
    return(values[values != 0])
}

## Basic diagnostic plot of a given first digits vector
plotDigitDistribution <- function(firstDigitsVector){
    hist(firstDigitsVector, breaks=seq(-0.5,9.5,1), xlim=c(0.5,9.5))
}

getRemainderAfterRemovingFirstDigit <- function(integerValue, integerBase){
    logValue <- log(integerValue, base=integerBase)
    floorOfLogValue <- floor(logValue)
    return(integerValue-integerBase^floorOfLogValue)
}

getTwoDigitBenfordDensities <- function(){
    return(pbenf(2))
}

calculateMAD <- function(firstTwoDigitActualProportions, firstTwoDigitExpectedProportions){
    return(sum(abs(firstTwoDigitActualProportions-firstTwoDigitExpectedProportions))/length(firstTwoDigitActualProportions))
}

calculateKS <- function(vectorData){
    ## From the vector data we calculate the mantissae of the logarithms
    mantissae <- log(vectorData, 10)-floor(log(vectorData,10))

    ## Test for conformity to a uniform distribution
    ## One-sample KS test
    ## Note that this does the full test, so this is not the optimized version.
    ksr <- ks.test(mantissae,"punif")

    return(ksr)
}

## Given a vector of raw values, a window size and an increment (increment = 1 to shift the window by one value each time), returns the vector of window starts. For example, given x=1:10, windowSize=3, increment=2, it will return 1,3,5,7 (note that the window starting at 7 also includes the values 8 and 9)
calculateWindowStarts <- function(x,windowSize,increment){
    return(seq(1,length(x)-windowSize,increment))
}

## Given a vector of raw values (e.g. tweeting user follower counts), a window size and an increment, this function first calculates the locations for the windows. If increment=1, then the window locations are just the positions of the vector 1,2,3,4,...,(length(x)-windowSize).
## For each window, it then calculates our Benfordness conformity measure p-value and returns it in a vector of p-values. The length of the returned vector is |windowStarts|
extractKSSignal <- function(x, windowSize, increment){

    windowStarts <- calculateWindowStarts(x,windowSize,increment)

    ksValues <- c();
    for(s in windowStarts){
        window <- x[s:(s+windowSize-1)]
        ksr <- calculateKS(window)
        ksValues <- c(ksValues,ksr$p.value)
    }
    
    return(ksValues)
}

## A utility function to generate the data for the digit histograms for each window in a vector. Each row represents a window. Each row has 9 entries, representing the first-digit histogram frequencies for that window. The result is a matrix with dimensions |windowStarts| x 9
extractDigitHistograms <- function(x, windowSize, increment){
    windowStarts <- calculateWindowStarts(x,windowSize,increment)

    digitHistograms <- matrix(,0,9)
    for(s in windowStarts){
        window <- x[s:(s+windowSize-1)]
        firstDigits <- as.numeric(table(c(1:9,getFirstDigits(window))))-1
        digitHistograms <- rbind(digitHistograms,matrix(firstDigits,1,9))
    }
    
    return(digitHistograms)
}

## Computes the mean for each window of size windowSize and increment
extractMeansSignal <- function(x,windowSize,increment){
    windowStarts <- calculateWindowStarts(x,windowSize,increment)

    meanValues <- c();
    for(s in windowStarts){
        window <- x[s:(s+windowSize-1)]
        meanValue <- mean(window)
        meanValues <- c(meanValues,meanValue)
    }
    
    return(meanValues)
}

## Computes the Sd for each window of size windowSize and increment
extractSdSignal <- function(x, windowSize, increment){
    windowStarts <- calculateWindowStarts(x,windowSize,increment)

    sdValues <- c();
    for(s in windowStarts){
        window <- x[s:(s+windowSize-1)]
        sdValue <- sd(window)
        sdValues <- c(sdValues,sdValue)
    }
    
    return(sdValues)
}

## Computes the skewness for each window of size windowSize and increment
extractSkewnessSignal <- function(x, windowSize, increment){
    windowStarts <- calculateWindowStarts(x,windowSize,increment)

    skewnessValues <- c();
    for(s in windowStarts){
        window <- x[s:(s+windowSize-1)]
        skewnessValue <- skewness(window)
        skewnessValues <- c(skewnessValues,skewnessValue)
    }
    
    return(skewnessValues)
}

## Computes the kurtosis for each window of size windowSize and increment
extractKurtosisSignal <- function(x, windowSize, increment){
    windowStarts <- calculateWindowStarts(x,windowSize,increment)

    kurtosisValues <- c();
    for(s in windowStarts){
        window <- x[s:(s+windowSize-1)]
        kurtosisValue <- kurtosis(window)
        kurtosisValues <- c(kurtosisValues,kurtosisValue)
    }
    
    return(kurtosisValues)
}

## Computes the MAD (see Nigrini) for each window of size windowSize and increment
extractMADSignal <- function(x, windowSize, increment){
    ## Note that we have to filter out entries with values less than 10, because we're doing a two-digit analysis
    xgt10 <- x[x>=10]

    benfDensities <- getTwoDigitBenfordDensities()
    
    ## Now getting the first two digits for the entire vector
    firstDigits <- getTwoSignificantDigits(xgt10)
    
    ## Now calculate the initial frequency table
    firstPortion <- firstDigits[1:windowSize]
    freq <- hist(firstPortion,seq(9.5,99.5,1),plot=FALSE)
    runningCounts <- freq$counts

    ## And calculate the first MAD value
    actualProportions <- runningCounts/windowSize
    madValues <- c(calculateMAD(actualProportions, benfDensities))

    totalLength <- length(firstDigits)
    nextIndex <- 1+increment
    while(nextIndex<=totalLength){
        ## Find the flux
        fallenOut <- firstDigits[(nextIndex-increment):(nextIndex-1)]
        comeIn <- firstDigits[(nextIndex+windowSize-increment):(nextIndex+windowSize-1)]
        comeIn <- comeIn[!is.na(comeIn)]

        fallenOutFreq <- hist(fallenOut,seq(9.5,99.5,1),plot=FALSE)
        comeInFreq <- hist(comeIn,seq(9.5,99.5,1),plot=FALSE)
        fluxCounts <- comeInFreq$counts-fallenOutFreq$counts
        runningCounts <- runningCounts+fluxCounts
        actualProportions <- runningCounts/windowSize
        madValues <- c(madValues,calculateMAD(actualProportions, benfDensities))

        ## Advance
        nextIndex <- nextIndex+increment
    }

    return(madValues)
}

## plots the first two digit histogram for the given raw vector data. Ignores values < 10
firstTwoDigitHistogram <- function(rawVectorData){
    gt10 <- rawVectorData[rawVectorData>=10]

    firstTwoDigits <- getTwoSignificantDigits(gt10)

    benfDensities <- getTwoDigitBenfordDensities()

    hist(firstTwoDigits,seq(9.5,100.5,1))
    points(10:99,benfDensities*length(gt10))    
}

## This is a function used for the optimized version of our algorithm, for online scenarios. See the code in the section "Experiments on Synthetic Data" of the paper.
shiftSortedWindow <- function(oldWindowUnsorted, oldWindowSorted, newValue){
    w <- length(oldWindowUnsorted)
    allIndices <- 1:w
    valueToEject <- oldWindowUnsorted[1]
    indexOfValueToEjectInSortedWindow <- which(oldWindowSorted==valueToEject)[1]
    newWindowUnsorted <- c(oldWindowUnsorted[2:w],newValue)
    newWindowSorted <- oldWindowSorted[setdiff(allIndices,indexOfValueToEjectInSortedWindow)]
    ## Now inject the new value into the sorted window
    ## First we have to find where it belongs
    if(newValue>=newWindowSorted[w-1]){
        newWindowSorted <- c(newWindowSorted,newValue)        
    } else{
        insertIndex <- which(newValue<newWindowSorted)[1]
        if(insertIndex==1){
            newWindowSorted <- c(newValue,newWindowSorted[insertIndex:(w-1)])            
        } else{
            newWindowSorted <- c(newWindowSorted[seq(1,(insertIndex-1),1)],newValue,newWindowSorted[insertIndex:(w-1)])
        }
    }
    return(list(newWindowUnsorted=newWindowUnsorted,newWindowSorted=newWindowSorted))
}

## This is used in the built-in ks.test function to calculate the p-value. In our case, if we want an optimized algorithm, we can work directly with the KS statistic (not the p-value) if we have a constant window size and alpha threshold. This function just calculates that statistic. Assumes the input is sorted. It finds the maximum deviation from to the uniform ECDF. See the source of the built-in function ks.test for more (it motivated this function).
calculateKsOneSidedStatisticUniform <- function(sortedX){
    n <- length(sortedX)
    x <- punif(sortedX) - (0:(n - 1))/n
    return(max(c(x, 1/n - x)))
}
