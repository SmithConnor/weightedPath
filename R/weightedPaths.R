

weighted_paths = function(x, y, family, B = 10, min, max, c = 0.5){

  n = x %>%
    nrow()

  # Generate weights
  weights = stats::rexp(n = n*B,
                 rate = 1) %>%
    matrix(data = .,
           nrow = n,
           ncol = B)

  # Fit LASSO
  modelFit = weights %>%
    base::apply(X = .,
          MARGIN = 2,
          FUN = glmnet::glmnet,
          x = x,
          y = y,
          family = family
    )

  restrictFit = lapply(X = modelFit,
                       FUN = L1_range,
                       min = min,
                       max = max)

  lambdaRange = matrix(data = NA_integer_, nrow = B, ncol = 2)
  for (b in 1:B){
    lambdaRange[b,] = restrictFit[[b]]$lambdaRange
  }

  lambdaVector = lapply(X = restrictFit,
                    FUN = pick_lambda) %>%
    unlist(x = .) %>%
    unique(x = .) %>%
    sort(x = .)

  L1Vector = lapply(X = restrictFit,
                    FUN = pick_L1) %>%
    unlist(x = .) %>%
    unique(x = .) %>%
    sort(x = .)

  modelCoef = base::lapply(X = modelFit,
                     FUN = calc_coef,
                     lambdaVector = lambdaVector,
                     min = min,
                     max = max)

  bootCoef = list()

  for(i in 1:length(lambdaVector)){
    coefMatrix = matrix(data = NA_integer_, nrow = (ncol(x)+1), ncol = B)
    for (b in 1:B){
      for(j in 1:(ncol(x)+1)){
        coefMatrix[j, b] = modelCoef[[b]][j,i]
      }
    }
    bootCoef[[i]] = coefMatrix
  }

  coefInt = lapply(bootCoef, conf_int, c = c)

  coefIntMat = list()
  coefIntMat[[1]] = matrix(NA_integer_, (ncol(x)+1), length(lambdaVector))
  coefIntMat[[2]] = matrix(NA_integer_, (ncol(x)+1), length(lambdaVector))

  for(i in 1:2){
    for (k in 1:length(lambdaVector)){
      for(j in 1:(ncol(x)+1)){
        coefIntMat[[i]][j, k] = coefInt[[k]][i,j]
      }
    }
  }


  averageCoef = base::Reduce(`+`, modelCoef)/B

  avgCoef = averageCoef %>%
    as.matrix(x = .) %>%
    t(x = .) %>%
    data.frame(.)


  avgCoef$id = lambdaVector


  for(i in 1:2){
    coefIntMat[[i]] = coefIntMat[[i]] %>%
      t(.) %>%
      data.frame(.)
    coefIntMat[[i]]$id = lambdaVector
    colnames(coefIntMat[[i]]) = colnames(avgCoef)
  }


  plotData = reshape2::melt(data = avgCoef,
                            id.var = "id",
                            factorsAsStrings = TRUE)

  plotDataLower = reshape2::melt(data = coefIntMat[[1]],
                                 id.var = "id",
                                 factorsAsStrings = TRUE)

  plotDataUpper = reshape2::melt(data = coefIntMat[[2]],
                            id.var = "id",
                            factorsAsStrings = TRUE)



  colnames(plotData) = c("Lambda", "Variable", "Coef")
  colnames(plotDataLower) = c("Lambda", "Variable", "Coef")
  colnames(plotDataUpper) = c("Lambda", "Variable", "Coef")

  plotData$Lower = plotDataLower$Coef
  plotData$Upper = plotDataUpper$Coef

  return(list(plotData = plotData,
              modelCoef = modelCoef))

}

##########################################################


wPath = weighted_paths(x = x,
                          y = y,
                          family = "binomial",
                          B = 10,
                          min = 1,
                          max = 5,
                       c = 0.6)

plot = ggplot(data = wPath$plotData, mapping = aes(x = Lambda, y = Coef, colour = Variable)) +
  geom_line() + geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.3, linetype = "blank")

ggplotly(plot)


##########################################################

model_range = function(vector, min, max){

  # adjust for no model of size

  while(min %in% vector == FALSE){
    min = min + 1
  }

  while(max %in% vector == FALSE){
    max = max - 1
  }

  minPos = vector %>%
    purrr::detect_index(.x = .,
                 .f = is_value,
                 value = min) - 1

  maxPos = vector %>%
    base::rev(x = .) %>%
    purrr::detect_index(.x = .,
                 .f = is_value,
                 value = max)

  # Adjust Value
  maxPos = vector %>%
    base::length(x = .) - maxPos + 1

  return(list(lowerBound = minPos,
              upperBound = maxPos,
              min = min,
              max = max))
}

##########################################################

is_value = function(x, value){
  x == value
}

##########################################################

L1_range = function(model, min, max){

  rangePos = model$df %>%
    model_range(vector = .,
               min = min,
               max = max)

  lambdaRange = c(model$lambda[rangePos$lowerBound],
                  model$lambda[rangePos$upperBound])

  restrictedLambda = (model$lambda[rangePos$lowerBound:rangePos$upperBound] - min(lambdaRange))/(max(lambdaRange) - min(lambdaRange))

  coef = rbind(model$a0, model$beta)

  rownames(coef)[1] = "int"

  restrictedCoef = coef[,rangePos$lowerBound:rangePos$upperBound]

  L1 = apply(X = coef,
             MARGIN = 2,
             FUN = L1_norm)

  L1Range = c(L1[rangePos$lowerBound],
              L1[rangePos$upperBound])

  restrictedL1 = (L1[rangePos$lowerBound:rangePos$upperBound] - min(L1Range))/(max(L1Range) - min(L1Range))

  return(list(coef = restrictedCoef,
              mappedLambda = restrictedLambda,
              mappedL1 = restrictedL1,
              min = rangePos$min,
              max = rangePos$max,
              lambdaRange = lambdaRange ))
}

##########################################################

L1_norm = function(vector){
  vector %>%
    abs(x = .) %>%
    sum(x = .)
}

##########################################################

pick_L1 = function(list){
  list$mappedL1
}

##########################################################

pick_lambda = function(list){
  list$mappedLambda
}

##########################################################

calc_coef = function(model, lambdaVector, min, max){

  nLambda = lambdaVector %>%
    length(x = .)

  rangePos = model$df %>%
    model_range(vector = .,
                min = min,
                max = max)

  lambdaRange = c(model$lambda[rangePos$lowerBound],
                  model$lambda[rangePos$upperBound])

  lambdaMin = min(lambdaRange)
  lambdaMax = max(lambdaRange)

  lambdaNew = lambdaVector * (lambdaMax - lambdaMin) + lambdaMin

  newCoef = coef(object = model,
                 s = lambdaNew)

  return(newCoef)

}

##########################################################

conf_int = function(matrix, c){
  output = apply(matrix, 1, conf_int_vector, c = c)
  return(output)
}

##########################################################

conf_int_vector = function(vector, c){
  length = vector %>%
    length(.)
  vector = vector %>%
    sort(.)
  which = ceiling(length * c)
  if (which %% 2 == 1){
    which = which + 1
  }
  first = (length - which)/2 + 1
  int = c(first, length - first + 1)
  return(vector[int])
}

