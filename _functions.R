#function for plotting the variograms estimates as fuzzy variogram
plotVariograms <- function(data.vgm, vgm.low, vgm.mid, vgm.high){
  #### plotting all variograms together ####
  # creating data frame for the variances
  Fitted <- data.frame(dist = seq(0.01, max(data.vgm$dist), length = 101))
  Fitted$low <- variogramLine(vgm.low, dist_vector = Fitted$dist)$gamma
  Fitted$mid <- variogramLine(vgm.mid, dist_vector = Fitted$dist)$gamma
  Fitted$high <- variogramLine(vgm.high, dist_vector = Fitted$dist)$gamma
  
  #convert the dataframes to a long format
  Empirical <- data.vgm
  Modeled <- melt(Fitted, id.vars = "dist", measure.vars = c("low", "mid", "high"))
  colnames(Modeled)[2] <- "variogram"
  
  ggplot(Empirical, aes(x = dist, y = gamma)) +  geom_point(size = 3.5) + 
    geom_line(data = Modeled, aes(x = dist, y=value, group = variogram, linetype = variogram, size = variogram)) +
    scale_size_manual(values=c(1,1,1)) +
    scale_linetype_manual(values = c("dotted", "solid", "dashed")) +
    labs(#title = "Fuzzy variogram of PM10 in the Czech Republic \n",
         x = "distance (°)", y = "semivariance \n") + 
    theme(legend.position="none", plot.title = element_text(lineheight=0.8, face="bold", size= 20, family = "sans"), axis.title = element_text(size= 15), axis.text = element_text(size= 13), panel.background = element_rect(fill = 'white', colour = 'black'), panel.grid.major = element_line(colour = "gray"), panel.grid.minor  = element_line(colour = "gray"))
}
  

# returns estimates of minimal and maximal values based on limits of sill, range, nugget a models
# vgm_start - is modal variogram selected for the data
# krigingModel - model of the kriging used in the calculations
# dataSet - the data for modelling, gridPred - grid of locations, at which the values should be predicted
# logResults - was the dataset logaritmized prior to the calculation?

# data format of result is list with two vectors - result$min_data, result$max_data

calculateMinMaxPoints <- function(psills, ranges, nuggets, models, vgm_start, 
                                  krigingModel, dataSet, gridPred, logResults=FALSE){
  
  numberOfCalculations = length(psills) * length(ranges) * length(nuggets) * length(models)
  
  #initial prediction is done based on modal model
  pred = krige(krigingModel, dataSet, gridPred, model = vgm_start)
  
  #were the data logaritmized prior to the calculation?
  if(logResults){
    min_data = exp(pred$var1.pred)
    max_data = exp(pred$var1.pred)
  }
  else{
    min_data = pred$var1.pred
    max_data = pred$var1.pred
  }
  
  calculationNumber = 1
  
  # for all combinations of the sills, ranges, nuggers and also models(if there is more than one)
  for (psill in psills){
    for (range in ranges){
      for(nugget in nuggets){
        for(model in models){
          
          #prepare the varigogram
          vgm <- vgm(psill = psill, model= model,range= range, nugget=nugget)
          
          #calculated the krigging
          pred <- krige(krigingModel, dataSet, gridPred, model = vgm)
          
          #if the dataset was logaritmized we need to exponentiate the outcome
          if(logResults){
            temp_data = exp(pred$var1.pred)
          }
          else{
            temp_data = pred$var1.pred
          }
          
          #compare the obtained preditions to minimal and maximal values that we allready have
          #if the value is outside of the range, than adjust the range
          for(i in 1:length(temp_data)){
            if(temp_data[i] < min_data[i]){
              min_data[i] = temp_data[i]
            }
            
            if(max_data[i] < temp_data[i]){
              max_data[i] = temp_data[i]
            }
          }
          
          print(paste("Calculation",calculationNumber,"of",numberOfCalculations,"done.",(calculationNumber/numberOfCalculations)*100,"%.", sep = " "))
          calculationNumber = calculationNumber + 1
        }
      } 
    }
  }
  
  #prepare the result and return it from the function
  result <- list(min.data=min_data,max.data=max_data)
  return(result)
}

#determine possibility of exceedance of threshold based on modal and maximal value
possibilityExceedance <- function(modal, max, threshold){
  data = cbind(modal, max)
  result = ifelse(data[,2] <= threshold, 0, ifelse( data[,1] >= threshold,  1, 1-((threshold-data[,1])/(data[,2]-data[,1])) ))
  return(result)
}

#determine necessity of exceedance of threshold based on modal and maximal value
necessityExceedance <- function(min, modal, threshold){
  data = cbind(min, modal)
  result = ifelse(data[,2] <= threshold, 0, ifelse( data[,1] >= threshold,  1, 1-((threshold-data[,1])/(data[,2]-data[,1])) ))
  return(result)
}
