# This script contains 4 functions used to model ET0 and water balance:
# 1. 'snowmod' estimates snowfall and snowpack and net moisture input as a function of temperature, precip, and existing
#   snowpack.  It also outputs a vector of albedo values, generally 0.2 if there is no snow, or 0.8 if there is snow.
# 2. 'monthlyETO' for calculating monthly reference evapotranspiration
# 3. 'dailyET0' for calculating daily reference evapotranspiration
# 4. 'aetmod' estimates actual et, deficit, soil moisture and runoff as a function of moisture input, existing
#   soil moisture, and soil water capacity. 
#
# Author: Alan Swanson 2012
###############################################################################

snowmod <- function(tmean,ppt,radiation=NULL,snowpack_prev=NULL,albedo=0.23,albedo_snow=0.8){
  # This function computes monthly estimated snowfall and snowmelt. Output includes end-of-month snowpack,
  # water "input" (snowmelt plus rain), and albedo.
  # Arguments:
  #  tmean - vector of mean monthly temperatures
  #  radiation - vector of shortwave solar radiation in MJ/m^2/day.
  #  snowpack_prev - vector of snowpack at the beginning of the month.  If NULL this is 
  #    taken to be zero.
  #  albedo - a single value for albedo in the absence of snow cover.  
  #  albedo_snow - single value of albedo given snow cover
  #
  # Value:  dataframe with three columns for end-of-month snowpack, H2O input (rain plus snowmelt),
  #  and albedo.
  
  
  N <- length(tmean)
  if(is.null(snowpack_prev)) snowpack_prev <- rep(0,N)
  snowpack <- rep(NA,N)
  input <- rep(NA,N)
  
  # this is for radiation in MJ/m^2/day
  mf <- function(t,t0,t1) pmin(pmax((t-t0)/(t1-t0),0),1)
  linrmelt <- function(temp,radiation,b0,b1,b2) pmax((b0+temp*b1+radiation*b2),0)
  parvec <- c(-4.604,6.329,-398.4,81.75,25.05)
  mfsnow <- mf(tmean,parvec[1],parvec[2])
  mfmelt <- linrmelt(tmean,radiation,parvec[3],parvec[4],parvec[5])
  
  # calculate values
  snow <- (1-mfsnow)*ppt
  rain <- mfsnow*ppt	
  melt <- pmin(mfmelt,snow+snowpack_prev) 
  snowpack <- snowpack_prev+snow-melt 
  input <-rain+melt
  
  # make vector of albedo values
  albedo <- rep(albedo,N)
  albedo[snowpack>0 | (snowpack_prev>0)] <- albedo_snow
  
  return(data.frame(snowpack=snowpack,input=input,albedo=albedo))
}


monthlyET0 <- function(radiation,tmax,tmin,wind,lat,elev,dpt,tmean_prev,albedo=0.23,month){
  # This function runs Reference ET estimates for monthly timesteps using methods based on
  # the Penman-Montieth equation as presented in Allen et al (1998).  It incorporates a 
  # modification which adjusts stomatal conductance downwards at temperatures below 5 C.
  #
  # Arguments:
  # radiation: vector of monthly average shortwave radiation in MJ/m^2/day
  # tmax, tmin: vectors of monthly average maximum and minimum temperatures in C, 
  # wind: vector of monthly average wind speed in m/s at 10m above ground, 
  # tmean_prev: vector of mean temp for the previous month, 
  # lat: vector of latitude in degrees 
  # elev: vector of elevation in meters, 
  # dpt: vector of dewpoint temperature in C.
  # tmean_prev: vector of mean temp of previous month in C
  # albedo: vector or scalar of albedo values, 
  # month: scalar 1-12.
  
  #
  # Value: 
  # Returns a vector of ET0 values.
  
  t0<-unclass(Sys.time())	
  daysinmonth=c(31,28,31,30,31,30,31,31,30,31,30,31)
  d2=c(31,59,90,120,151,181,212,243,273,304,334,365)
  d1=c(1,32,60,91,121,152,182,213,244,274,305,335)
  DoY <- (d1[month]+d2[month])/2 # use middle day of month to represent monthly average. 
  n_days <- daysinmonth[month]
  
  # calculate soil heat flux (total for the month) using change in temperature from previous month
  tmean <- (tmax+tmin)/2 
  G <- 0.14*(tmean-tmean_prev) # fixed from previous version
  
  # convert to wind height at 2m
  hw=10 # height of wind measurements 
  wind <- wind*(4.87/log(67*hw-5.42))  # convert to wind height at 2m
  
  # stomatal conductance adjustment for low temperatures
  sr=100 # stomatal resistance sec/m
  ks_min=.01 # minimum value for temps below T1
  Tl=-10       # minimum temp (sc goes to ks_min below this temp)
  T0=5		# optimal temp
  Th=100     # maximum temp (sc goes to zero above this)
  thresh=5   # temperature threshold below which to apply Jarvis equation (ks=1 above this temp)
  b4 <- (Th-T0)/(Th-Tl)
  b3 <- 1/((T0-Tl)*(Th-T0)^b4)
  ks <- pmax(pmin(b3*(tmean-Tl)*(Th-tmean)^b4,1),ks_min)
  ks[is.na(ks)] <- ks_min
  ks[tmean>=thresh] <- 1
  
  # convert to stomatal resistance.
  sr  <- sr/ks
  
  # ra is aerodynamic resistance, rs is bulk surface resistance
  ra  <- 208/wind #(log((2-2/3*0.12)/(0.123*0.12))*log((2-2/3*0.12)/(0.1*0.123*0.12)))/(0.41^2*wind) # equal to 208/wind for hh=hw=2.
  #ra <- 208/wind
  rs <- sr/(0.5*24*0.12) # value of 70 when sr=100
  
  # Saturation vapor pressure , 
  es <- 0.6108*exp(tmin*17.27/(tmin+237.3))/2+0.6108*exp(tmax*17.27/(tmax+237.3))/2     
  ea <- 0.6108*exp((dpt)*17.27/((dpt)+237.3))
  vpd <- es - ea
  vpd[vpd<0] <- 0    # added because this can be negative if dewpoint temperature is greater than mean temp (implying vapor pressure greater than saturation).
  
  # delta - Slope of the saturation vapor pressure vs. air temperature curve at the average hourly air temperature 
  delta  <- (4098 * es)/(tmean + 237.3)^2  
  
  P <- 101.3*((293-0.0065*elev)/293)^5.26  # Barometric pressure in kPa
  lambda <- 2.501-2.361e-3*tmean # latent heat of vaporization    
  cp  <- 1.013*10^-3 # specific heat of air
  gamma <- cp*P/(0.622*lambda) # Psychrometer constant (kPa C-1)
  pa <- P/(1.01*(tmean+273)*0.287) # mean air density at constant pressure
  
  # Calculate potential max solar radiation or clear sky radiation	
  GSC=0.082      # MJ m -2 min-1 (solar constant)
  phi <- pi*lat/180 
  dr <- 1+0.033*cos(2*pi/365*DoY)      
  delt <- 0.409*sin(2*pi/365*DoY-1.39)     
  omegas <- acos(-tan(phi)*tan(delt)) 
  Ra <- 24*60/pi*GSC*dr*(omegas*sin(phi)*sin(delt) +cos(phi)*cos(delt)*sin(omegas))    # Daily extraterrestrial radiation
  Rso <- Ra*(0.75+2e-5*elev)     #For a cloudless day, Rs is roughly 75% of extraterrestrial radiation (Ra)
  
  
  # radfraction is a measure of relative shortwave radiation, or of
  # possible radiation (cloudy vs. clear-sky)
  radfraction <- radiation/Rso
  radfraction[radfraction>1] <- 1
  
  # longwave  and net radiation
  longw <- 4.903e-9*n_days*((tmax+273.15)^4+(tmin+273.15)^4)/2*(.34-.14*sqrt(ea))*(1.35*radfraction-.35)     
  netrad <- radiation*n_days*(1-albedo)-longw     
  
  # ET0
  et0 <- .408*((delta*(netrad-G))+(pa*cp*vpd/ra*3600*24*n_days))/(delta+gamma*(1+rs/ra))
  return(et0)
} 



dailyET0 <- function(radiation,tmax,tmin,wind,lat,elev,albedo=0.23,dpt,doy){
  # This is a ET0 function designed for daily inputs.  
  
  # Arguments:
  # radiation: vector of monthly average shortwave radiation in MJ/m^2/day
  # tmax, tmin: vectors of monthly average maximum and minimum temperatures in C, 
  # wind: vector of monthly average wind speed in m/s, 
  # tmean_prev: vector of mean temp for the previous month, 
  # lat: vector of latitude in degrees 
  # elev: vector of elevation in meters, 
  # albedo: scaler or vector of albedo values, 
  # doy: scalar day of year 1-365,
  
  #
  # Value: 
  # Returns a vector of ET0 values.
  tmean <- (tmin+tmax)/2
  n_days <- 1 
  G <- 0 # assume soil heat flux to be zero
  
  # wind adjustment to 2m from 10m output
  hw=10 # height at which wind is measured
  wind <- wind*(4.87/log(67*hw-5.42))  
  
  # stomatal conductance adjustment for low temperatures
  sr=100 # stomatal resistance sec/m
  ks_min=.01 # minimum value for temps below T1
  Tl=-10       # minimum temp (sc goes to ks_min below this temp)
  T0=5		# optimal temp
  Th=100     # maximum temp (sc goes to zero above this)
  thresh=5   # temperature threshold below which to apply Jarvis equation (ks=1 above this temp)
  b4=(Th-T0)/(Th-Tl)  # from Jarvis 1978
  b3=1/((T0-Tl)*(Th-T0)^b4)
  ks=pmax(pmin(b3*(tmean-Tl)*(Th-tmean)^b4,1),ks_min)
  ks[is.na(ks)] <- ks_min
  ks[tmean>=thresh] <- 1
  
  # convert to stomatal resistance.
  sr <- sr/ks
  
  # ra is aerodynamic resistance, rs is bulk surface resistance
  ra <- 208/wind # 
  rs <- sr/(0.5*24*0.12) # value of 70 when sr=100
  
  # Saturation vapor pressure , 
  es <- 0.6108*exp(tmin*17.27/(tmin+237.3))/2+0.6108*exp(tmax*17.27/(tmax+237.3))/2  # saturation vapor pressure
  ea <- 0.6108*exp((dpt)*17.27/((dpt)+237.3))                                        # actual vapor pressure
  vpd <- es - ea
  vpd[vpd<0] <- 0    
  
  # delta - Slope of the saturation vapor pressure vs. air temperature curve
  delta <- (4098 * es)/(tmean + 237.3)^2  
  P <- 101.3*((293-0.0065*elev)/293)^5.26  # Barometric pressure in kPa
  lambda <- 2.501-2.361e-3*tmean # latent heat of vaporization    
  cp <- 1.013*10^-3 # specific heat of air
  gamma <- cp*P/(0.622*lambda) # Psychrometer constant (kPa C-1)
  pa <- P/(1.01*(tmean+273)*0.287) # mean air density at constant pressure
  
  # Calculate potential max solar radiation or clear sky radiation.
  GSC = 0.082      # MJ m -2 min-1 (solar constant)
  phi <- pi*lat/180 
  dr <- 1+0.033*cos(2*pi/365*doy)      
  delt <- 0.409*sin(2*pi/365*doy-1.39)     
  omegas <- acos(-tan(phi)*tan(delt))     
  Ra <- 24*60/pi*GSC*dr*(omegas*sin(phi)*sin(delt) +cos(phi)*cos(delt)*sin(omegas)) # daily extraterrestrial radiation
  Rso <- Ra*(0.75+2e-5*elev)     #For a cloudless day, Rs is roughly 75% of extraterrestrial radiation (Ra)
  
  # radfraction is a measure of relative shortwave radiation, or of
  # possible radiation (cloudy vs. clear-sky), needs to be less than 1
  radfraction <- radiation/Rso
  radfraction[radfraction>1] <- 1
  
  # longwave radiation
  longw <- 4.903e-9*n_days*((tmax+273.15)^4+(tmin+273.15)^4)/2*(.34-.14*sqrt(ea))*(1.35*radfraction-.35)     
  
  # net radiation
  netrad <- radiation*n_days*(1-albedo)-longw     
  
  # ET from long-form P-M eqn.
  et0 <- .408*((delta*(netrad-G))+(pa*cp*vpd/ra*3600*24*n_days))/(delta+gamma*(1+rs/ra))
  return(et0)
} 


aetmod <- function(et0,input,awc,soil_prev=NULL){
  # This function computes AET given ET0, H2O input, soil water capacity, and beginning-of-month soil moisture
  # Arguments:
  # et0: vector of monthly reference evapotranspiration in mm
  # input: vector of monthly water input to soil in mm
  # awc: vector of soil water capacity in mm
  # soil_prev: vector of soil water content for the previous month (mm).  If left NULL this is assigned to be zero.
  #
  # Value:
  # returns a data frame with columns for AET, deficit, end-of-month soil moisture, and runoff.
  
  N <- length(et0)
  runoff <- def <-  aet <- soil <- rep(NA,N) # 
  if(is.null(soil_prev)) soil_prev <- rep(0,N)
  
  deltasoil <- input-et0 # positive=excess H2O, negative=H2O deficit
  
  # Case when there is a moisture surplus:
  Case <- deltasoil>=0
  if(sum(Case)>0){
    aet[Case] <- et0[Case]
    def[Case] <- 0
    soil[Case] <- pmin(soil_prev[Case]+deltasoil[Case],awc[Case])	# increment soil moisture, but not above water holding capacity
    runoff[Case] <- pmax(soil_prev[Case]+deltasoil[Case]-awc[Case],0) # when awc is exceeded, send the rest to runoff
  }
  
  # Case where there is a moisture deficit:  soil moisture is reduced
  Case <- deltasoil<0
  if(sum(Case)>0){
    soildrawdown <- soil_prev[Case]*(1-exp(-(et0-input)[Case]/awc[Case]))	# this is the net change in soil moisture (neg)
    aet[Case] <- pmin(input[Case] + soildrawdown,et0[Case])
    def[Case] <- et0[Case] - aet[Case]
    soil[Case] <- soil_prev[Case]-soildrawdown
    runoff[Case] <- 0
  }
  
  return(data.frame(aet=aet,def=def,soil=soil,runoff=runoff))
  
}


