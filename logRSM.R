# based on regRSM package (http://cran.r-project.org/web/packages/regRSM/)

options(width = 140)

compute_initial_weights = function(y,x){
# The function returns initial weights.

  initial_weights = as.numeric(cor(y,x))^2
  initial_weights = initial_weights/(sum(initial_weights))
  return(initial_weights)
}


is.wholenumber = function(x, tol = .Machine$double.eps^0.5){
# The function checks if the given number is whole number.

  abs(x - round(x)) < tol
}

compute_scores = function(y,x,m,B,initial_weights=NULL){
# The function returns RSM scores.

  p = ncol(x)
  scores = numeric(p)
  ns = numeric(p)

  for(k in 1:B){
    submodel = sample(1:p,size=m,replace=FALSE,prob=initial_weights)
    lm1 = glm(y~x[,submodel], family = binomial)
    weights = as.numeric((summary(lm1)$coef[-1,3])^2)
    scores[submodel] =  scores[submodel] + weights
    ns[submodel] = ns[submodel] + 1
  }
  ns = ifelse(ns!=0,ns,1)
  scores = scores/ns

  return(scores)
}

select_finalmodel_bic = function(y,x,order1,thrs,penalty){
# The function returns final model using Bayesian Information Criterion.

  n = length(y)
  lm0 = lm(y~1)
  beta0 = as.numeric(coef(lm0))
  rss0 = sum(lm0$residuals^2)
  bic0 = n*log(rss0/n)+1*penalty

  xo = x[,order1[1:thrs]]
  xo = cbind(1,xo)
  xo = as.matrix(xo)


  qr1 = qr(xo)

  Q = qr.Q(qr1)
  R = qr.R(qr1)
  R_inv = solve(R)

  tQy = t(Q) %*% y

  betabj = vector("list",thrs)
  rssj =  numeric(thrs)
  bic = numeric(thrs)


  RSthrs = y-Q %*% tQy
  rssj[thrs] = t(RSthrs) %*% RSthrs
  for(j in thrs:2){
    rssj[j-1] = rssj[j]+ tQy[(j+1)]^2
    bic[j] = n*log(rssj[j]/n)+(j+1)*penalty
  }

  bic[1] =  n*log(rssj[1]/n)+2*penalty

  sel = which.min(bic)
  betab = as.numeric(R_inv[1:(sel+1),1:(sel+1)] %*% tQy[1:(sel+1)])
  model_sel = order1[1:sel]

  if(bic0<bic[sel]){
    betab = beta0
    model_sel = 0
  }

  Result = list(model=model_sel,informationCriterion=bic,coefficients=betab)
  return(Result)
}

new.logRSM <- function()
{

  logRSM=list(scores=NULL,model=NULL,time=list(user=0,system=0,elapsed=0),
      data_transfer=list(user=0,system=0,elapsed=0),
      coefficients=NULL, predError=NULL,input_data=list(x=NULL,y=NULL),
      control=list(useGIC=NULL,selval=NULL,screening=NULL,init_weights=FALSE,parallel=NULL,m=NULL,B=NULL))

  attr(logRSM,"class")="logRSM"
  return(logRSM)
}

logRSM = function(y,x,yval=NULL,xval=NULL,m=NULL,B=NULL, parallel="NO",nslaves=c(4),
    store_data=FALSE,screening=NULL,init_weights=FALSE,useGIC=TRUE,thrs=NULL,penalty=NULL,initial_weights=NULL)
{
  if (init_weights) {
    if (!is.null(initial_weights)) {
      stop('init_weights cannot be TRUE if initial_weigths are provided')
    }
  }
  data_x = x;
  x = as.matrix(x)
  y = as.numeric(y)
  n = length(y)
  p = ncol(x)
  scores = NULL

  startTime <- proc.time()

  # Set default values of m and B
  if(is.null(m)){
    m = floor(min(n-1,p)/2)
  }else{
    if(m>(n-2)) stop("Parameter m cannot be larger than the number of observations minus two!")
    if(m<=0) stop("Parameter m must be a positive number!")
    if(!is.wholenumber(m)) stop("Parameter m must be a whole number!")
  }
  if(is.null(B)){
    B = 1000
  }else{
    if(B<=0) stop("Parameter B must be a positive number!")
    if(!is.wholenumber(B)) stop("Parameter B must be a whole number!")
  }

  #Check for screeneing
  if(!is.null(screening))
  {
    if((screening>=1)||(screening<=0)) stop("screening must be in (0,1)")

    iw =  compute_initial_weights(y,x)
    sel = which(iw>=quantile(iw,screening))
    if(m>length(sel)) stop('Parameter m cannot be larger than the number of attributes remaining after screening procedure!')
    x = x[,sel]
  }
  #Check for initial_weights
  if(init_weights){
    initial_weights = compute_initial_weights(y,x)
  }

  #RSM method esence
  d1=d2=proc.time()
  scores = compute_scores(y,x,m,B,initial_weights)

  #Set score 0, when variable is not selected by screeneing
  if(!is.null(screening)){
    scores1 = numeric(ncol(data_x))
    scores1[sel] = scores
    scores = scores1
  }


  selval = ifelse(!is.null(yval) && !is.null(xval),TRUE,FALSE)
  if(selval) useGIC = FALSE

  if(useGIC){
    if(is.null(penalty)){
      penalty = log(length(y))
    }else{
      if(penalty<0) stop("Penalty must be positive!")
    }
    if(is.null(thrs)){
      thrs = ifelse(p<=floor(n/2),p,floor(n/2))
    }else{
      if(thrs>min(p,(n-2))) stop("Parameter thrs cannot be larger than min(p,(n-2))!")
      if(thrs<=1) stop("Parameter thrs must be greater than one!")
      if(!is.wholenumber(thrs)) stop("Parameter thrs must be a whole number!")
    }
    order1 = sort(scores,decreasing=TRUE,index.return=TRUE)$ix
    selected_model = select_finalmodel_bic(y,data_x,order1,thrs,penalty)
    model = selected_model$model
    coefficients =  as.numeric(selected_model$coefficients)
    predError = NULL
    informationCriterion = selected_model$informationCriterion
  }else{
    if(selval==TRUE){
      order1 = sort(scores,decreasing=TRUE,index.return=TRUE)$ix
      selected_model = select_finalmodel_qr(y,data_x,yval,xval,order1)
      model = selected_model$model
      coefficients =  as.numeric(selected_model$coefficients)
      predError = selected_model$predError
      informationCriterion = NULL
    }else{
      model = NULL
      coefficients =  NULL
      predError = NULL
      informationCriterion = NULL
    }
  }
  stopTime <- proc.time()

  logRSM = new.logRSM()
  logRSM$scores = scores
  logRSM$model = model
  logRSM$time = stopTime-startTime
  logRSM$coefficients = coefficients
  logRSM$predError = predError
  logRSM$informationCriterion = informationCriterion
  logRSM$data_transfer = d2-d1
  if(store_data) { logRSM$input_data$x=data_x; logRSM$input_data$y=y }

  logRSM$control$useGIC = useGIC
  logRSM$control$selval = selval
  logRSM$control$screening = screening
  logRSM$control$init_weights =  init_weights
  logRSM$control$parallel = parallel
  logRSM$control$m = m
  logRSM$control$B = B

  return(logRSM)
}

generateTestData = function(columns = 10000, observations = 100, shifted = 10, shift = 1.0) {
  if (shifted > columns) {
    stop ("shifted cannot be more then columns")
  }

  data = matrix(rnorm(columns * observations), nrow = observations, ncol = columns)
  colnames(data) <-  c(paste('V_N01_' , 1:(ncol(data) - shifted), sep=''), paste('V_S_' , 1:shifted, sep='') )
  row.names(data) <- paste('', 1:nrow(data), sep='')

  shiftCol = numeric(observations)
  class = numeric(observations)
  for (i in 1:observations ) {
    class[i] <- i %% 2
    if (i %% 2 == 0) {
      shiftCol[i] <- shift
    } else {
      shiftCol[i] <- -shift
    }
  }

  for (i in (columns - shifted + 1):columns) {
    data[, i] = data[, i] + shiftCol;
  }

  result = list(data = data, class = class)
  return(result)
}

calculateLogRSMWithRandomData = function(columns = 10000) {
  for (i in 1:1) {
    cat ('calculating ', i, '\n')
    data = generateTestData(columns = columns, observations = 100, shifted = 100)
    reg = logRSM(data$class, data$data, m = 10)
    result = rev(colnames(data$data)[order(reg$scores)])
    print (result[1:10])
    write.table(t(result), file = "c:/workspaces/dmlab/regresja/result2.csv", col.names=FALSE, row.names=FALSE, sep = ",", quote = FALSE, append = TRUE)
  }
}

calculate_initial_weigths_with_t = function(class, data) {
  result = numeric(ncol(data));
  for (i in 1:ncol(data)) {
    x = numeric();
    y = numeric();
    for (j in 1:length(class)) {
      if (class[j] == 0) {
        x[length(x) + 1] <- data[j, i];
      } else {
        y[length(y) + 1] <- data[j, i];
      }
    }
    tResult = t.test(x, y, paired = FALSE, alternative = 'two.sided', var.equal = TRUE)$statistic;
    result[i] = abs(tResult);
  }
  result = result / sum(result)
  return (result)
}

test_calculate_initial_weights_with_t = function() {

  dane = read.csv('testfile.csv', header = FALSE)
  class = dane[, 1]
  data = dane[, 2:ncol(dane)]

  result = calculate_initial_weigths_with_t(class, data);
  print (result)
}

# test_calculate_initial_weights_with_t()

subexperiment = function(subsets, i, b, shifted, calculateWeigths = TRUE) {
  print(paste('licze', i, 'dla', subsets, 'i', b))

  inputpath = paste('dane_', i, '.csv', sep = '')
  dane = read.csv(inputpath, header = FALSE)
  class = dane[, 1]
  data = dane[, 2:ncol(dane)]

  weigthsAA = NULL
  fileName = 'log_'
  if (calculateWeigths) {
    calculate_initial_weigths_with_t(class, data);
    fileName = 'log_wg_'
  }

  colnames(data) <- c(paste('V_S_' , 1:shifted, sep=''), paste('V_N_' , (shifted+1):ncol(data), sep='') )

  reg = logRSM(class, data, m = subsets, B = b, initial_weights = weigthsAA)
  result = rev(colnames(data)[order(reg$scores)])
  print (result[1:10])

  resultpath = paste(fileName, b, '_', subsets, '.csv', sep = '')
  write.table(t(result), file = resultpath, col.names=FALSE, row.names=FALSE, sep = ",", quote = FALSE, append = TRUE)
}

experimentWithWeights = function(subsets, b, shifted = 20) {
  for (i in 1:100) {
    print(system.time(subexperiment(subsets, i, b, shifted, calculateWeigths = TRUE)))
  }
}

experimentWithoutWeights = function(subsets, b, shifted = 20) {
  for (i in 1:100) {
    print(system.time(subexperiment(subsets, i, b, shifted, calculateWeigths = FALSE)))
  }
}

