#####This the code of the function definition of the model in Microrisk Lab#####

#Required Packages
library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(rhandsontable)
library(mvtnorm)
library(mc2d)
library(mgcv)
library(ggplot2)
library(Metrics)
library(shinyalert)
library(plotly)
library(nlstools)

#####################################  Models  #####################################
#Primary models (referenced to Tabel 1)

##Growth models
###Complete Gompertz model
Gompertz = function(t, mumax, lag, Y0, Ymax){
  Yt = Y0 + (Ymax - Y0) * exp(-exp(mumax*exp(1) * (lag - t)/(Ymax - Y0) + 1))
  return(Yt)
}
###Complete Baranyi model
Baranyi = function(t, mumax, lag, Y0, Ymax){
  At = mumax*t+log(exp(-mumax*t)+exp(-mumax*lag)-exp(-mumax*t-mumax*lag))
  Yt = Y0 + At-log(1+(exp(At)-1)/exp(Ymax-Y0))
  return(Yt)
}
###Complete Buchanan model
Buchanan = function(t, mumax, lag, y0, ymax){
  yt = y0 + (t >= lag) * (t <= (lag + (ymax - y0)*log(10)/mumax)) * (mumax/log(10)) * (t -lag) + (t >=lag) * (t >(lag + (ymax - y0)*log(10)/mumax)) * (ymax - y0)
  return(yt)
}
###Lag-Logistic model
Laglogistic = function(t, mumax, lag, Y0, Ymax){
  Yt = (t <= lag) * Y0 + (t >lag) * (Ymax - log(1 + (exp(Ymax - Y0) - 1) * exp(-mumax * (t-lag))))
  return(Yt)
}
###Complete Huang model
Huang = function(t, mumax, lag, Y0, Ymax){
  Yt = Y0 +Ymax-log(exp(Y0)+(exp(Ymax)-exp(Y0))*exp(-mumax *(t+log((1+exp(-4*(t-lag)))/(1+exp(4*lag)))/4)))
  return(Yt)
}
###Logistic model
Logistic = function(t, mumax, Y0, Ymax){
  Yt =Y0+Ymax-log(exp(Y0)+(exp(Ymax)-exp(Y0))*exp(-mumax*t))
  return(Yt)
}
###No lag Buchanan model
Buchanannolag = function(t, mumax, y0, ymax){
  yt =  y0 + (t <= ((ymax - y0)*log(10)/mumax)) * (mumax/log(10)) * t+ (t > ((ymax - y0)*log(10) /mumax)) * (ymax - y0)
  return(yt)
}
###Reduced Baranyi model
ReducedBaranyi = function(t, mumax, lag, Y0){
  At = (lag>=0)*(mumax*t+log(exp(-mumax*t)+exp(-mumax*lag)-exp(-mumax*t-mumax*lag)))
  Yt = Y0 + At
  return(Yt)
}
###Reduced Buchanan model
ReducedBuchanan = function(t, mumax,lag, y0){
  yt =  (t <= lag) * y0 + (t >lag) * (y0+(mumax/log(10))* (t - lag))
  return(yt)
}
###Reduced Huang model
ReducedHuang = function(t, mumax, lag, Y0){
  Yt = Y0+mumax*(t+log((1+exp(-4*(t-lag)))/(1+exp(4*lag)))/4)
  return(Yt)
}
###Linear model
Linear = function(t, mumax, y0){
  yt = y0+(mumax/log(10))*t
  return(yt)
}

##Inactivation models
###Complete Geeraerd model
Geeraerd = function(t,L,k,y0,yres){
  yt = yres + log10((10^(y0 - yres) - 1) * exp(k* L)/(exp(k* t) + (exp(k* L) - 1)) + 1)
  return(yt)
}
###Three-phase model
Three = function(t,L,k,y0,yres){
  yt = y0 - (t >= L) * (t <= (L + (y0 - yres) * log(10)/k)) * k * (t - L)/log(10) + (t >= L) * (t > (L + (y0 - yres) * log(10)/k)) * (yres - y0)
  return(yt)
}
###Weibull-tail model
Weibulltail = function(t,delta,p,y0,yres){
  yt=yres + log10((10^(y0 - yres) - 1) * 10^(-(t/delta)^p) + 1)
  return(yt)
}
###No shoulder Geeraerd model
Geeraerdnos = function(t,k,y0,yres){
  yt = yres + log10(1 + (10^(y0 - yres) - 1) * exp(-k * t))
  return(yt)
}
###No shoulder two-phase model
Twonos = function(t,k,y0,yres){
  yt = y0 - (t <= ((y0 - yres) * log(10)/k)) * k * t/log(10) + (t > ((y0 - yres) * log(10)/k)) * (yres - y0)
  return(yt)
}
###No tail Geeraerd model
Geeraerdnot = function(t,L,k,y0){
  yt = y0 - k * t/log(10) + log10(exp(k * L)/(1 + (exp(k * L) - 1) * exp(-k * t)))
  return(yt)
}
###No tail two-phase model
Twonot = function(t,L,k,y0){
  yt = (t <= L) * y0 + (t > L) * (y0 - k/log(10) * (t - L))
  return(yt)
}
###Weibull model
Weibull = function(t,delta,p,y0){
  yt=y0 - (t/delta)^p
  return(yt)
}
###Bigelow model
Bigelow = function(t,D,y0){
  yt = y0-(t/D)
  return(yt)
}


#Secondary model (referenced to Tabel 2)

##Temperature model
###Suboptimal square-root model
sSQR = function(T, a,Tmin){
  rate=(a*(T-Tmin))^2
  return(rate)
}
###Full square-root model
fSQR = function(T, a,b,Tmin,Tmax){
  rate=((a*(T-Tmin))*(1-exp(b*(T-Tmax))))^2
  return(rate)
}
###Suboptimal Huang square-root model
sHSQR = function(T, a,Tmin){
  rate=((a*(T-Tmin)^0.75))^2
  return(rate)
}
###Full Huang square-root model 
fHSQR = function(T, a,b,Tmin,Tmax){
  rate=((a*(T-Tmin)^0.75)*(1-exp(b*(T-Tmax))))^2
  return(rate)
}
###Cardinal parameter model
CPMT = function(T, muopt,Topt,Tmin,Tmax){
  rate=(muopt*(T-Tmax)*((T-Tmin)^2))/((Topt-Tmin)*((Topt-Tmin)*(T-Topt)-(Topt-Tmax)*(Topt+Tmin-2*T)))
  return(rate)
}

##pH model
###Cardinal 3-parameter model
CPMPH1 = function(pH, muopt,pHopt,pHmin){
  rate=((pH >= pHmin) & (pH <= (2 * pHopt - pHmin))) * muopt * (pH - pHmin) * (pH - (2 * pHopt - pHmin))/((pH - pHmin) * (pH - (2 * pHopt - pHmin)) - (pH - pHopt)^2)
  return(rate)
}
###Cardinal 4-parameter model
CPMPH2 = function(pH, muopt,pHopt,pHmin,pHmax){
  rate=((pH >= pHmin) & (pH <= pHmax)) * muopt * (pH - pHmin) * (pH - pHmax)/((pH - pHmin) * (pH - pHmax) - (pH - pHopt)^2)
  return(rate)
}
####Quasi-mechanistic model
QMPH = function(pH, muopt,pHopt,pHmin){
  rate= muopt * (1-(10^(pHmin-pH)))
  return(rate)
}

#Aw model
#Cardinal 2-parameter model
CPMAW1 = function(aw, muopt,awmin){
  rate=(aw >= awmin) * muopt * (aw - awmin)^2/(1 - awmin)^2
  return(rate)
}
#Cardinal 3-parameter model 
CPMAW2 = function(aw, muopt,awopt,awmin){
  rate=(aw >= awmin) * muopt * (aw - 1) * (aw - awmin)^2/((awopt - awmin) * ((awopt - awmin) * (aw - awopt) - (awopt - 1) * (awopt + awmin - 2 * aw)))
  return(rate)
}


#Complex models (referenced to Tabel 3)

##Two flora competition growth model
###Jameson - No lag Buchanan model
JamBunolag=function(t, Flora,mumax1, y01, tmax,mumax2, y02){
  yt=(Flora == 1) * ((t < tmax) * (y01 + mumax1/log(10) *t) + (t >= tmax) * (y01 + mumax1/log(10) * tmax)) + 
    (Flora == 2) * ((t < tmax) * (y02 + mumax2/log(10) * t) + (t >= tmax) * (y02 + mumax2/log(10) * tmax))
  
  return(yt)
}
###Jameson - Buchanan model
JamBu=function(t, Flora,lag1,mumax1, y01,tmax,lag2, mumax2, y02){
  yt=(Flora == 1) * ((t <= lag1) * y01 + ((t > lag1) & (t < tmax)) * (y01 + mumax1/log(10) * (t - lag1)) + (t >= tmax) * (y01 + mumax1/log(10) * (tmax - lag1))) + 
    (Flora == 2) * ((t <= lag2) * y02 + ((t > lag2) & (t < tmax)) * (y02 + mumax2/log(10) * (t - lag2)) + (t >= tmax) * (y02 + mumax2/log(10) * (tmax - lag2)))
  return(yt)
}


##Non-isothermal growth model (referenced to Tabel 3)
###Cardinal parameter model for dynamic
dCPM = function(Temp,muopt,Topt,Tmin,Tmax){
  len = length(Temp)
  for(i in 1:len){
    if(Temp[i]>Tmin&Temp[i]<Tmax){
      rate[i] =muopt * (Temp[i] - Tmax) * (Temp[i]  - Tmin)^2/((Topt - Tmin) * ((Topt - Tmin) * (Temp[i]  - Topt) - (Topt - Tmax) * (Topt + Tmin - 2*Temp[i] )))
    }
    else{
      rate[i] = 0
    }
  }
  return(rate)
}
###Baranyi - Cardinal parameter model
dBaranyiModel = function(t, Y, mu, Q, Ymax){
  dY = Q*mu*(1-exp(Y-Ymax))/(1+Q)
  return(dY)
}
dQ = function(t,Q,mu){
  dQ = mu*Q
  return(dQ)
}
dBaranyiCPMNLS = function(t,ts,Temp, muopt,Topt,Tmin,Tmax, Q0, Y0, Ymax){
  len_t=length(t)-1
  LenY=length(ts)
  Q[1] = Q0
  mu = dCPM(Temp, muopt,Topt,Tmin,Tmax)
  Zout[1] = Y0
  for(i in 1:len_t){
    k1Q = dQ(t[i], Q[i], mu[i])
    k2Q = dQ(t[i]+0.5*dt, Q[i]+0.5*dt*k1Q, mu[i])
    k3Q = dQ(t[i]+0.5*dt, Q[i]+0.5*dt*k2Q, mu[i])
    k4Q = dQ(t[i]+dt, Q[i]+dt*k3Q, mu[i])
    Q[i+1] = Q[i] +dt*(k1Q+2*k2Q+2*k3Q+k4Q)/6
    k1 = dBaranyiModel(t[i], Zout[i], mu[i], Q[i], Ymax)
    k2 = dBaranyiModel(t[i]+0.5*dt, Zout[i]+0.5*dt*k1, mu[i], Q[i], Ymax)
    k3 = dBaranyiModel(t[i]+0.5*dt, Zout[i]+0.5*dt*k2, mu[i], Q[i], Ymax)
    k4 = dBaranyiModel(t[i]+dt, Zout[i]+dt*k3, mu[i], Q[i], Ymax)
    Zout[i+1] = Zout[i] +dt*(k1+2*k2+2*k3+k4)/6
  }
  for(i in 1:length(t)){
    for(j in 1:LenY){
      if(t[i]==ts[j]){
        Yout[j]=Zout[i]
      }
    }
  }
  return(Yout)
}
dBaranyiCPM = function(t,Temp,muopt,Topt,Tmin,Tmax, Q0, Y0, Ymax){
  len_t=length(t)-1
  Q[1] = Q0
  mu = dCPM(Temp, muopt,Topt,Tmin,Tmax)
  Zout[1] = Y0
  for(i in 1:len_t){
    k1Q = dQ(t[i], Q[i], mu[i])
    k2Q = dQ(t[i]+0.5*dt, Q[i]+0.5*dt*k1Q, mu[i])
    k3Q = dQ(t[i]+0.5*dt, Q[i]+0.5*dt*k2Q, mu[i])
    k4Q = dQ(t[i]+dt, Q[i]+dt*k3Q, mu[i])
    Q[i+1] = Q[i] +dt*(k1Q+2*k2Q+2*k3Q+k4Q)/6
    k1 = dBaranyiModel(t[i], Zout[i], mu[i], Q[i], Ymax)
    k2 = dBaranyiModel(t[i]+0.5*dt, Zout[i]+0.5*dt*k1, mu[i], Q[i], Ymax)
    k3 = dBaranyiModel(t[i]+0.5*dt, Zout[i]+0.5*dt*k2, mu[i], Q[i], Ymax)
    k4 = dBaranyiModel(t[i]+dt, Zout[i]+dt*k3, mu[i], Q[i], Ymax)
    Zout[i+1] = Zout[i] +dt*(k1+2*k2+2*k3+k4)/6
  }
  return(Zout)
}
###Huang - Cardinal parameter model
dHuangModel = function(t, Y, mu,A,m,Ymax){
  lag=exp(A)/((mu)^m)
  dY = mu*(1-exp(Y-Ymax))/(1+exp(-4*(t-lag)))
  return(dY)
}
dHuangCPMNLS = function(t,ts,Temp,muopt,Topt,Tmin,Tmax,Y0,A,m,Ymax){
  len_t=length(t)-1
  LenY=length(ts)
  Zout[1] = Y0
  mu = dCPM(Temp, muopt,Topt,Tmin,Tmax)
  for(i in 1:len_t){
    k1 = dHuangModel(t[i], Zout[i], mu[i], A,m,Ymax)
    k2 = dHuangModel(t[i]+0.5*dt, Zout[i]+0.5*dt*k1, mu[i],A,m,Ymax)
    k3 = dHuangModel(t[i]+0.5*dt, Zout[i]+0.5*dt*k2, mu[i],A,m,Ymax)
    k4 = dHuangModel(t[i]+dt, Zout[i]+dt*k3, mu[i],A,m,Ymax)
    Zout[i+1] = Zout[i] +dt*(k1+2*k2+2*k3+k4)/6
  }
  for(i in 1:length(t)){
    for(j in 1:LenY){
      if(t[i]==ts[j]){
        Yout[j]=Zout[i]
      }
    }
  }
  return(Yout)
}
dHuangCPM = function(t,Temp,muopt,Topt,Tmin,Tmax,Y0,A,m,Ymax){
  len_t=length(t)-1
  Zout[1] = Y0
  mu = dCPM(Temp, muopt,Topt,Tmin,Tmax)
  for(i in 1:len_t){
    k1 = dHuangModel(t[i], Zout[i], mu[i], A,m,Ymax)
    k2 = dHuangModel(t[i]+0.5*dt, Zout[i]+0.5*dt*k1, mu[i],A,m,Ymax)
    k3 = dHuangModel(t[i]+0.5*dt, Zout[i]+0.5*dt*k2, mu[i],A,m,Ymax)
    k4 = dHuangModel(t[i]+dt, Zout[i]+dt*k3, mu[i],A,m,Ymax)
    Zout[i+1] = Zout[i] +dt*(k1+2*k2+2*k3+k4)/6
  }
  return(Zout)
}
###Logistic - Cardinal parameter model
dLogisticModel = function(t, Y, mu, Ymax){
  dY = mu*(1-exp(Y-Ymax))
  return(dY)
}
dLogisticCPMNLS = function(t,ts,Temp,muopt,Topt,Tmin,Tmax, Y0, Ymax){
  len_t=length(t)-1
  LenY=length(ts)
  mu = dCPM(Temp, muopt,Topt,Tmin,Tmax)
  Zout[1] = Y0
  for(i in 1:len_t){
    k1 = dNoLagModel(t[i], Zout[i], mu[i],  Ymax)
    k2 = dNoLagModel(t[i]+0.5*dt, Zout[i]+0.5*dt*k1, mu[i], Ymax)
    k3 = dNoLagModel(t[i]+0.5*dt, Zout[i]+0.5*dt*k2, mu[i], Ymax)
    k4 = dNoLagModel(t[i]+dt, Zout[i]+dt*k3, mu[i], Ymax)
    Zout[i+1] = Zout[i] +dt*(k1+2*k2+2*k3+k4)/6
  }
  for(i in 1:length(t)){
    for(j in 1:LenY){
      if(t[i]==ts[j]){
        Yout[j]=Zout[i]
      }
    }
  }
  return(Yout)
}
dLogisticCPM = function(t,Temp,muopt,Topt,Tmin,Tmax, Y0, Ymax){
  len_t=length(t)-1
  mu = dCPM(Temp, muopt,Topt,Tmin,Tmax)
  Zout[1] = Y0
  for(i in 1:len_t){
    k1 = dNoLagModel(t[i], Zout[i], mu[i],  Ymax)
    k2 = dNoLagModel(t[i]+0.5*dt, Zout[i]+0.5*dt*k1, mu[i], Ymax)
    k3 = dNoLagModel(t[i]+0.5*dt, Zout[i]+0.5*dt*k2, mu[i], Ymax)
    k4 = dNoLagModel(t[i]+dt, Zout[i]+dt*k3, mu[i], Ymax)
    Zout[i+1] = Zout[i] +dt*(k1+2*k2+2*k3+k4)/6
  }
  return(Zout)
}


##Non-isothermal inactivation model (referenced to Tabel 3)
###Dynamic Bigelow model 
dBigelowModel = function(t,Temp,Y,Tref,Z,Dref){
  dY=-1/(Dref*(10^(-(Temp-Tref)/Z)))
  return(dY)
}
dBigelowNLS = function(t,Temp,ts,Tref,Z,Dref,Y0){
  len_t=length(t)-1
  Zout[1] = Y0
  leny = length(ts)
  for(i in 1:len_t){
    k1 = dBigelowModel(t[i], Temp[i],Zout[i],Tref,Z,Dref)
    k2 = dBigelowModel(t[i]+0.5*dt, Temp[i], Zout[i]+0.5*dt*k1,Tref,Z,Dref)
    k3 = dBigelowModel(t[i]+0.5*dt, Temp[i], Zout[i]+0.5*dt*k2,Tref,Z,Dref)
    k4 = dBigelowModel(t[i]+dt, Temp[i], Zout[i]+dt*k3,Tref,Z,Dref)
    Zout[i+1] = Zout[i] +dt*(k1+2*k2+2*k3+k4)/6
  }
  for(i in 1:length(t)){
    for(j in 1:leny){
      if(t[i]==ts[j]){
        Yout[j]=Zout[i]
      }
    }
  }
  return(Yout)
}
dBigelow = function(t,Temp,Tref,Z,Dref,Y0){
  len_t=length(t)-1
  Zout[1] = Y0
  for(i in 1:len_t){
    k1 = dBigelowModel(t[i], Temp[i],Zout[i],Tref,Z,Dref)
    k2 = dBigelowModel(t[i]+0.5*dt, Temp[i], Zout[i]+0.5*dt*k1,Tref,Z,Dref)
    k3 = dBigelowModel(t[i]+0.5*dt, Temp[i], Zout[i]+0.5*dt*k2,Tref,Z,Dref)
    k4 = dBigelowModel(t[i]+dt, Temp[i], Zout[i]+dt*k3,Tref,Z,Dref)
    Zout[i+1] = Zout[i] +dt*(k1+2*k2+2*k3+k4)/6
  }
  return(Zout)
}