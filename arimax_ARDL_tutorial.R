install.packages('forecast')
install.packages('tseries')
install.packages('ggplot2')
install.packages('TSA')
install.packages('Ecdat')
install.packages('Hmisc')
install.packages('astsa')

library(forecast)
library(tseries)
library(ggplot2)
library(TSA)
library(Ecdat)
library(Hmisc)
library(astsa)

## learning to model one series as its own ARIMA process

daily_data = read.csv('/home/rstudio/dynamic_regression/day.csv', header=TRUE, stringsAsFactors=FALSE)
daily_data$Date = as.Date(daily_data$dteday)

#quick plot
ggplot(daily_data, 
       aes(Date, cnt)) + geom_line() + scale_x_date('month') + ylab("Daily Bike Checkouts") + xlab("")

#create ts object
count_ts = ts(daily_data[, c('cnt')])
daily_data$clean_cnt = tsclean(count_ts)

#plot again
ggplot() + 
  geom_line(data = daily_data, aes(x = Date, y = clean_cnt)) + ylab('Cleaned Bicycle Count') 

#smoothing it manually with moving average of order 7 (week based) order 30 (montsh based)
daily_data$cnt_ma = ma(daily_data$clean_cnt, order=7) # using the clean count with no outliers
daily_data$cnt_ma30 = ma(daily_data$clean_cnt, order=30)

#plot
ggplot() + 
  geom_line(data = daily_data, aes(x = Date, y = clean_cnt, colour = "Counts")) +
  geom_line(data = daily_data, aes(x = Date, y = cnt_ma,   colour = "Weekly Moving Average"))  +
  ylab('Bicycle Count') 

#frequency 30 means a unit of measure is a month, 30 obs each month
count_ma = ts(na.omit(daily_data$cnt_ma), frequency=30) 
decomp = stl(count_ma, s.window="periodic") #seasonal decomposition of ts by loess
deseasonal_cnt <- seasadj(decomp)
plot(decomp) 

count_d1 = diff(deseasonal_cnt, differences = 1)
plot(count_d1)
adf.test(count_d1, alternative = "stationary")
Acf(count_d1, main='ACF for Differenced Series')
Pacf(count_d1, main='PACF for Differenced Series') 

#feed data into auto.arima. Note only deseasonal is enough. DO not manually difference it.

fit<-auto.arima(deseasonal_cnt, seasonal=FALSE)
tsdisplay(residuals(fit), lag.max=45, main='(2,1,0) Model Residuals') 

fit2 = auto.arima(deseasonal_cnt)
fit2
tsdisplay(residuals(fit2), lag.max=15, main='Seasonal Model Residuals')

fcast <- forecast(fit2, h=30)
plot(fcast)
fcast

#testing GOF with holdout sample
hold <- window(ts(deseasonal_cnt), start=701)
fit_no_holdout = auto.arima(ts(deseasonal_cnt[-c(701:725)]), seasonal = TRUE, max.p = 7, max.q = 7)
fit_no_holdout

fcast_no_holdout <- forecast(fit_no_holdout,h=25)
plot(fcast_no_holdout, main=" Forecast with holdout sample")
lines(ts(deseasonal_cnt))

fit_w_seasonality = auto.arima(deseasonal_cnt, seasonal=TRUE, max.p = 7, max.q = 7)
fit_w_seasonality
ggtsdisplay(residuals(fit_w_seasonality), lag.max=15, main='Seasonal Model Residuals')
#use ggtsdisplay to get a quick plot of a ts
checkresiduals(fit)
#this function is a Q-statistic test to see if the residual ts is white noise + tsdisplay


#Next, transfer function models in tsa package (the added advantage is that tsa allows
#for lagged effects of the covariates, inaddition to lagged effects of response series)

data(airmiles)
ggtsdisplay(log(airmiles),ylab='Log(airmiles)',xlab='Year',main='Airmiles Logged')

acf(diff(diff(window(log(airmiles),end=c(2001,8)),12)),lag.max=48,main='First-order differenced and Yearly Differences')

#I will need to come back to this later ... cannot really understand it
air.m1=arimax(log(airmiles),
              order=c(0,1,1),
              seasonal=list(order=c(0,1,1),period=12),
              xtransf=data.frame(I911=1*(seq(airmiles)==69),I911=1*(seq(airmiles)==69)),
              transfer=list(c(0,0),c(1,0)),
              xreg=data.frame(Dec96=1*(seq(airmiles)==12),
                              Jan97=1*(seq(airmiles)==13),
                              Dec02=1*(seq(airmiles)==84)),
              method='ML')

######################################
##from another tutorial using tsa package
data(airquality)
#there are missing values. so only subset rows with complete records
ozone <- subset(na.omit(airquality))

######################################
#build trainig set and test set
set.seed(123)
#calculate 70% of rows = 78 rows
#ceiling means if the value is 77.8, return 78
#floor is the opposite operator
N.train <- ceiling(0.7 * nrow(ozone)) 
#total rows - training rows = Num test rows
N.test <- nrow(ozone) - N.train
#get sequence from 1:N train rows and the rest as test set
trainset <- seq(1:N.train)
testset <- setdiff(seq(1:nrow(ozone)), trainset)

######################################
# assign each date a decimal day in a year
# (not sure why this step is necessary?)

# only m and d is available, so we manually add a year
year <- 1973 # known from data documentation
dates <- as.Date(paste0(year, "-", ozone$Month, "-", ozone$Day))
dates

#%j	Decimal day of the year	
min.date <- as.numeric(format(min(dates), "%j"))
max.date <- as.numeric(format(max(dates), "%j"))

#start and end are the times of the first and last observation
#frequency is the number of observations per unit time (1=annual, 4=quartly, 12=monthly if the unit is a year, etc.).
ozone.ts <- ts(ozone$Ozone, start = min.date, end = max.date, frequency = 1)
ozone.ts <- window(ozone.ts, 121, 231) # deal with repetition  due to missing time values
# assumes that measurements are consecutive although they are not，due to missing value
ozone$t <- seq(start(ozone.ts)[1], end(ozone.ts)[1]) 

######################################
#now we have a ts objecdt ozone
ggtsdisplay(ozone$Ozone) # we could get this without creating demical days, right??

ozone.ts <- ts(ozone$Ozone, frequency = 30 ) #daily data, unit = month
decomp = stl(ozone.ts, s.window="periodic") #seasonal decomposition of ts by loess
plot(decomp) 
#visual inspection shows no obvious seasonal trend

#manually remove the seasonal trend or not
deseasonal_cnt <- seasadj(decomp)
ggtsdisplay(deseasonal_cnt)

######################################
#creating ARIMAX regression

auto.arima(ozone$Ozone[trainset], seasonal = TRUE)
#the best model is ARIMA(0,1,2)

features <- c("Solar.R", "Wind", "Temp") # exogenous features
A <- arimax(x = ozone$Ozone[trainset], 
            xreg = ozone[trainset,features], #a datafrme containing covariates
            order = c(0,1,2))
A

#predict(model to predict with, dataframe with testset covariates)
preds.temporal <- predict(A, newxreg = ozone[testset, features])
preds.temporal
#hand calculate R squared
Rsquared.temporal <- cor(preds.temporal$pred, ozone[testset, "Ozone"])^2


preds2 <- predict(A, newxreg = ozone[trainset, features])
preds2
#hand calculate R squared
Rsquared2<- cor(preds2$pred, ozone[trainset, "Ozone"])^2


##############################
#ARIMAX on the airquality data set
######

data(Icecream)

# create a time-series object
library(lubridate) # for week function
wk <- week(c(as.Date("1951-03-18"), as.Date("1953-07-11"))) #decimal week of a year
# manually put into 3 years of months 
# from March 1951 - July 1953
months <- c(seq(3,12), seq(1,12), seq(1,7)) 
#get decimal week of each four-weekly observation
wks <- c(seq(wk[1], 52, 4), seq(1, 52, 4), seq(1, 52, 4))
# these two outputs were not used in creating ts!

#start week, end week of that year, freq = each year how many four weekly obs?
ice.ts <- ts(Icecream$cons, start = c(1951, 3), end = c(1953, 7), frequency = 52/4)
plot(decompose(ice.ts))
#visual inspection shows obvious seasonal trend
ggtsdisplay(ice.ts)


train <- 1:20
test <- 21:30
features <- c("income", "temp")
A <- arimax(x = window(ice.ts, c(1951,3), c(1951,3+20-1)),
           xreg = Icecream[train, features],
           order = c(1,1,0))
A
fitted(A)
plot(fitted(A),lty =3)
lines(window(ice.ts, c(1951,3), c(1951,22), lty=3))


preds <- predict(A, newxreg = Icecream[test, features]) 
preds
plot(preds$pred)
lines(window(ice.ts, c(1951, 10), c(1953,7)),lty=3) # actual values in dashed black


##############################
#ARIMAX on sales dataset
######

data(sales)
data(lead)
#use lead to predict sales
sales_140<-window(sales,end=140)
lead_140<-window(lead,end=140)
mock140 <- 

lead_140_Z <- lead_140 - mean(lead_140)
sales_140_D <- diff(sales_140)
lead_140_D <- diff(lead_140_Z)

mod <- arimax(sales_140,
              order=c(0,0,1),
              include.mean = TRUE,
              xtransf=lead_140_Z, #note it is not xreg
              #fixed=c(NA,NA,0,0,0,0,0,NA,NA),
              transfer=list(c(0,1)), #if there are multiple covariate, list them all as ARMA order p,q
              method="ML")
mod


                  
