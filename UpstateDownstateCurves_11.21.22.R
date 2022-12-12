# Explore the differences in temperature-trait R0 based on data from two
# Culex pipiens mosquito populations.

# Set working directory
setwd("~/OneDrive/Documents/Arbo Lab/Life history/Model")

# Load R0 functions
source("R0_functions2_fay_20211006.R")

# load package for minor tick
library("Hmisc")

# Read in trait measurements
in.csv = "TraitInput_2pop_9.14.22.csv"
trait.data2 = read.csv(in.csv)

#
head(trait.data2)
#
View(trait.data2)

# Generate 'typical' Cx. pipiens trait profiles based on Shocket's fit
Temps = seq(0,40, 0.1)
R0s = sapply(Temps, calc.R0.v2, 1, 1)

LS = sapply(Temps, LifeSpan)
BR = sapply(Temps, BitingRate)
EFGC = sapply(Temps, calc_EFGC)
MDR = sapply(Temps, calc_MDR)
pLAs = sapply(Temps, calc_pLA)

#UNY = Upstate New York
#DNY = Downstate New York

#plot pO
pop.pO = trait.data2[trait.data2$Trait == "ProportionOvipositioning", ]
pop.pO$COLOR = "#CC79A7"
#'purple3'

pop.pO$COLOR[pop.pO$Population == "Albany"] = "#56B4E9"
#'Green3'
pop.pO$Temperature = as.numeric(pop.pO$Temperature)
pop.pO$Value = as.numeric(pop.pO$Value)
plot(pop.pO$Temperature, pop.pO$Value, col = pop.pO$COLOR, main = "Percent Ovipositing (pO)" ,xlim = c(0,40), ylim = c(0,0.6), xlab = "Temperature (°C)", ylab = "Percent ovipositing", pch = 19, cex.lab=1.5, cex.main=1.5)
axis(1, xlim = c(0,40), at = seq(0,40,5))

#DNY.pOs = sapply(Temps, ) quadratic
SNY.pOs = sapply(Temps, Quadratic, 0.046, 21.5, 28.5)
par(new = TRUE)
plot(Temps, SNY.pOs, xlim = c(0,40), ylim = c(0,0.6), xlab = "Temperature (°C)", ylab = "Percent ovipositing", col = "#CC79A7", cex.lab= 1.5, cex.main=1.5)
legend("topleft", 
       legend = c("Downstate", "Upstate"), 
       col = c("#CC79A7", "#56B4E9"), 
       pch = c(19))

#UNY pOs quadratic 
ENY.pOs = sapply(Temps, Quadratic, 0.010000, 17, 28.5)
par(new = TRUE)
plot(Temps, ENY.pOs, xlim = c(0,40), ylim = c(0,0.6), xlab = "Temperature (°C)", ylab = "Percent ovipositing", col = "#56B4E9", xaxp = c(0, 40, 8),cex.lab=1.5, cex.main=1.5)
legend("topleft", 
       legend = c("Downstate", "Upstate"), 
       col = c("#CC79A7", "#56B4E9"), 
       pch = c(19))
minor.tick(nx=10, ny=2, tick.ratio=0.5, x.args = list(), y.args = list())

# plot lf 
# From Table S1 in Mordecai et al. 2019
Reverse.Briere = function(Temp, q, Tmin, Tmax){
  
  #stop("This implementation does not actually capture a Reverse Briere Shape")
  
  
  # Suppress warnings for sqrt of negative numbers
  if (Tmax < Temp){
    value = 0
  }else{
    equiv.Temp = (Tmax - Temp) + Tmin # Shift to be the same deviation from Tmin
    value = q*equiv.Temp *(equiv.Temp - Tmin)*sqrt(Tmax - equiv.Temp)
    
    #value = q*Temp *-sqrt(Temp - Tmin)*(Temp - Tmax) # This didn't work
  }
  
  # Constrain value to not be NA or negative
  if(is.na(value)){ value = 0 }
  if(value < 0){    value = 0 }
  # Constrain function to return 0's for negative inputs
  if(Temp  < 0){    value = 0 }
  
  return(value)
}

#plot lf
#max.ls = 230
pop.lifespan = trait.data2[trait.data2$Trait == "LifeSpan", ]
pop.lifespan$COLOR = "#CC79A7"
#par(new = TRUE)


# Data taken from Ciota et al. 2011 https://doi.org/10.2987/8756-971X-27.1.21
# Add a point for known low-temperature limits based on over-wintering
tunnel.1.temp = mean(c(9.41, 11.38))
# 18% drop from Dec - Jan. Do we assume this is a monthly mortality?
# Nov 28 - Jan 7 and Dec. 5 - Jan 8
# 40 days for tunnel 1, 34 for tunnel 2
# So daily mortality should be
#tunnel.1.mortality = 0.18 / 40 # Not correct, proportion changes are multiplicative
#tunnel.1.mortality.v2 = 0.18**(1/40) # NOpe, that's not right - need to apply the exponent to survival
tunnel.1.survival = 1 - 0.18
tunnel.1.daily.survival = tunnel.1.survival**(1/40)
tunnel.1.mortality = 1 - tunnel.1.daily.survival
# Life span is inverse of daily mortality rate, correct?
tunnel.1.lifespan = 1 / tunnel.1.mortality
# Populations mean differences are assumed be the same in the absence of direct measurements.
Albany.mean.val = sum(as.numeric(pop.lifespan$Value[pop.lifespan$Population == "Albany"])) / length(pop.lifespan$Value[pop.lifespan$Population == "Albany"])
suffolk.mean.val = sum(as.numeric(pop.lifespan$Value[pop.lifespan$Population == "Suffolk"])) / length(pop.lifespan$Value[pop.lifespan$Population == "Suffolk"])
mean.pop.diff = Albany.mean.val - suffolk.mean.val
# Note that this population is from Cohoes, NY (I thought the study was at Fort Totten? near Queens)
added.record.1 = c(tunnel.1.temp, "LifeSpan", tunnel.1.lifespan, "Suffolk", NA, NA, "#CC79A7")
added.record.2 = c(tunnel.1.temp, "LifeSpan", tunnel.1.lifespan- mean.pop.diff, "Albany", NA, NA, "#56B4E9")
pop.lifespan = rbind(pop.lifespan, added.record.1)
pop.lifespan = rbind(pop.lifespan, added.record.2)

pop.lifespan = trait.data2[trait.data2$Trait == "LifeSpan", ]
pop.lifespan$COLOR = "#CC79A7"
pop.lifespan$COLOR[pop.lifespan$Population == "Albany"] = "#56B4E9"

#par(new = TRUE)

# DNY LS
plot(pop.lifespan$Temperature, pop.lifespan$Value, xlim = c(0,40), ylim = c(0, 250), col = pop.lifespan$COLOR,main = "Adult Lifespan (lf)", xlab = "Temperature (°C)", ylab = "Time (days)", pch=19, cex.lab=1.5, cex.main=1.5)
axis(1, xlim = c(0,40), at = seq(0,40,5))
# Fit to middle point and low temperature point
q.Rev = 0.15
Tmin.Rev = 0
Tmax.Rev = 31
test.ls = sapply(Temps, Reverse.Briere, q.Rev, Tmin.Rev, Tmax.Rev)
par(new = TRUE)
plot(Temps, test.ls, xlim = c(0,40), ylim = c(0,250), col = "#CC79A7", xlab = "Temperature (°C)", ylab = "Time (days)", cex.lab=1.5, cex.main=1.5)
legend("topright", 
       legend = c("Downstate","Upstate"), 
       col = c("#CC79A7", "#56B4E9"), 
       pch = c(19))

# UNY LS
q.Rev.2 = 0.17
Tmin.Rev.2 = 0
Tmax.Rev.2 = 30
test.ls.2 =  sapply(Temps, Reverse.Briere, q.Rev.2, Tmin.Rev.2, Tmax.Rev.2)
par(new = TRUE)
plot(Temps, test.ls.2, xlim = c(0,40), ylim = c(0,250), col = "#56B4E9", xlab = "Temperature (°C)", ylab = "Time (days)", xaxp = c(0, 40, 8), cex.lab=1.5, cex.main=1.5)
legend("topright", 
       legend = c("Downstate", "Upstate"), 
       col = c("#CC79A7", "#56B4E9"), 
       pch = c(19))
minor.tick(nx=10, ny=2, tick.ratio=0.5, x.args = list(), y.args = list())

# plot a
BR = sapply(Temps, BitingRate)

#plot(Temps, BR, xlim = c(0,40), ylim = c(0, max(BR)), main = "Biting Rate (a)", xlab = "Temperature (°C)", ylab = "Rate (day-1)")
pop.bitingrate = trait.data2[trait.data2$Trait == "BitingRate", ]
pop.bitingrate$COLOR = "#CC79A7"
pop.bitingrate$COLOR[pop.bitingrate$Population == "Albany"] = "#56B4E9"
#par(new = TRUE)

plot(pop.bitingrate$Temperature, pop.bitingrate$Value, xlim = c(0,40), ylim = c(0, 0.37), col = pop.bitingrate$COLOR, main = "Biting Rate (a)", xlab = "Temperature (°C)", ylab = "Rate (day-1)", pch=19, xaxt="n", cex.lab=1.5, cex.main=1.5)
axis(1, xlim = c(0,40), at = seq(0,40,5))

# DNY change fit Briere FITS BEST
SNY.BR= sapply(Temps, Briere, 0.00085, 20.5, 31.5)
par(new = TRUE)
plot(Temps, SNY.BR, xlim = c(0,40), ylim = c(0, 0.37), col = "#CC79A7", xlab = "", ylab = "")

#UNY 
ENY.BR= sapply(Temps, Briere, 0.00048, 15, 28.2)
par(new = TRUE)
plot(Temps, ENY.BR, xlim = c(0,40), ylim = c(0, 0.37), col = "#56B4E9", xlab = "", ylab = "", xaxp = c(0, 40, 8), cex.lab=1.5, cex.main=1.5)
legend("topleft", 
       legend = c("Downstate", "Upstate"), 
       col = c("#CC79A7", "#56B4E9"), 
       pch = c(19))
minor.tick(nx=10, ny=2, tick.ratio=0.5, x.args = list(), y.args = list())


#plot MDR
MDR = sapply(Temps, calc_MDR)

pop.MDR = trait.data2[trait.data2$Trait == "MDR", ]
pop.MDR$COLOR = "#CC79A7"
pop.MDR$COLOR[pop.MDR$Population == "Albany"] = "#56B4E9"

plot(pop.MDR$Temperature, pop.MDR$Value, xlim = c(0,40), ylim = c(0, 0.14), col = pop.MDR$COLOR,
     main = "Mosquito Development Rate (MDR)",  xlab = "Temperature (°C)", ylab = "Rate (day-1)",pch=19, cex.lab=1.5, cex.main=1.5)
axis(1, xlim = c(0,40), at = seq(0,40,5))

# DNY change fit Briere FITS  Best
SNY.MDR= sapply(Temps, Briere, 0.00021, 17, 31)
par(new = TRUE)
plot(Temps, SNY.MDR, xlim = c(0,40), ylim = c(0, 0.14), col = "#CC79A7", xlab = "", ylab = "", cex.lab=1.5, cex.main=1.5)
legend("topleft", 
       legend = c("Downstate","Upstate"), 
       col = c("#CC79A7", "#56B4E9"), 
       pch = c(19))

# UNY change fit Briere FITS  Best
ENY.MDR= sapply(Temps, Briere, 0.00024, 16, 30 )
par(new = TRUE)
plot(Temps, ENY.MDR, xlim = c(0,40), ylim = c(0, 0.14), col = "#56B4E9", xlab = "", ylab = "", cex.lab=1.5, cex.main=1.5,  
     xaxp = c(0, 40, 8))
legend("topleft", 
       legend = c("Downstate", "Upstate"), 
       col = c("#CC79A7", "#56B4E9"), 
       pch = c(19))
minor.tick(nx=10, ny=2, tick.ratio=0.5, x.args = list(), y.args = list())

#plot pLa
pLAs = sapply(Temps, calc_pLA)

pop.pLa = trait.data2[trait.data2$Trait == "pLa", ]
pop.pLa$COLOR = "#CC79A7"
pop.pLa$COLOR[pop.pLa$Population == "Albany"] = "#56B4E9"

plot(pop.pLa$Temperature, pop.pLa$Value, xlim = c(0,40), ylim = c(0, 1), col = pop.pLa$COLOR,
     main = "Larval-to-Adult Survival (pLA)", xlab = "Temperature (°C)", ylab = "Survival Probability",pch=19, cex.lab=1.5, cex.main=1.5)
axis(1, xlim = c(0,40), at = seq(0,40,5))

#DNY.pLa = sapply(Temps, ) quadratic BEST
SNY.pLa = sapply(Temps, Quadratic, 0.010, 15, 35)
par(new = TRUE)
plot(Temps, SNY.pLa, xlim = c(0,40), ylim = c(0,1), xlab = "Temperature (°C)", ylab = "Survival Probability", col = "#CC79A7", cex.lab=1.5, cex.main=1.5)
legend("topleft", 
       legend = c("Downstate", "Upstate"), 
       col = c("#CC79A7", "#56B4E9"), 
       pch = c(19))

#UNY pLa quadratic 
ENY.pLa = sapply(Temps, Quadratic, 0.014, 15, 31)
par(new = TRUE)
plot(Temps, ENY.pLa, xlim = c(0,40), ylim = c(0,1), xlab = "Temperature (°C)", ylab = "Survival Probability", col = "#56B4E9", xaxp = c(0, 40, 8), cex.lab=1.5, cex.main=1.5)
legend("topleft", 
       legend = c("Downstate", "Upstate"), 
       col = c("#CC79A7", "#56B4E9"), 
       pch = c(19))
minor.tick(nx=10, ny=2, tick.ratio=0.5, x.args = list(), y.args = list())

#plot ER
# define ER
ER.q = 6.36 * 10**-1 
ER.Tmin = 5.0        
ER.Tmax = 37.7

ER = sapply(Temps, Quadratic,ER.q, ER.Tmin, ER.Tmax)

#plot(Temps, ER, xlim = c(0,40), ylim = c(0, max(ER)), xlab = "Temperature (°C)", ylab = "Eggs per raft", pch=19)
pop.EggsPerRaft = trait.data2[trait.data2$Trait == "EggsPerRaft", ]
pop.EggsPerRaft$COLOR = "#CC79A7"
pop.EggsPerRaft$COLOR[pop.EggsPerRaft$Population == "Albany"] = "#56B4E9"
#par(new = TRUE)

plot(pop.EggsPerRaft$Temperature, pop.EggsPerRaft$Value, xlim = c(0,40), ylim = c(0, 200), col = pop.EggsPerRaft$COLOR,main = "Fecundity (ER)",  xlab = "Temperature (°C)", ylab = "Eggs per raft", pch=19, cex.lab=1.5, cex.main=1.5)


# DNY fit quadratic 
Suffolk.EggsPerRaft= sapply(Temps, Briere, 0.35, 16, 29)
par(new = TRUE)

plot(Temps, Suffolk.EggsPerRaft, xlim = c(0,40), ylim = c(0, 200), col = "#CC79A7", xlab = "", ylab = "", cex.lab=1.5, cex.main=1.5)


#UNY fit quadratic 
Albany.EggsPerRaft= sapply(Temps, Briere, 0.36, 15, 29)
par(new = TRUE)
plot(Temps, Albany.EggsPerRaft, xlim = c(0,40), ylim = c(0, 200), col = "#56B4E9", xlab = "", ylab = "",cex.lab=1.5, cex.main=1.5,
     xaxp = c(0, 40, 8))
legend("topleft", 
       legend = c("Downstate", "Upstate"), 
       col = c("#CC79A7", "#56B4E9"), 
       pch = c(19))
minor.tick(nx=10, ny=2, tick.ratio=0.5, x.args = list(), y.args = list())

#RO calculation

# New code for species setup for UNY  
source("R0_functions2_fay_20211006.R")

test1 = species.setup("Culex pipiens")

Temps = seq(0,40,0.1)

##UNY R0 all set up
# Run them each separtely and then use par(new = TRUE)
new.all = species.setup("Culex pipiens")
new.all[1] = 0 
new.all[34] = 0.36
new.all[35] = 15
new.all[36] = 29
new.all[30] = 1
new.all[31] = 0.010000
new.all[32] = 17
new.all[33] = 28.5
new.all[2] = 0.00048
new.all[3] = 15
new.all[4] = 28.2
new.all[23] = 0.00024
new.all[24] = 16
new.all[25] = 30
new.all[37] = 0.17
new.all[38] = 0
new.all[39] = 30
new.all[16] = 0.014
new.all[17] = 15
new.all[18] = 31

#par(new = TRUE)
R0s.all = sapply(Temps, calc.R0.v2, 1,1, setup = new.all)
plot(Temps, R0s.all, xlim = c(0,40), ylim = c(0, 18), xlab = "Temperature (°C)", ylab = "Relative R0", col = "#56B4E9", xaxp = c(0, 40, 8), cex.lab=1.5, cex.main=1.5)
legend("topleft", 
       legend = c("Downstate", "Upstate"), 
       col = c("#CC79A7", "#56B4E9"), 
       pch = c(19))
minor.tick(nx=10, ny=2, tick.ratio=0.5, x.args = list(), y.args = list())

# DNY all set up
# Run them each separtely and then use par(new = TRUE)
new.all1 = species.setup("Culex pipiens")
new.all1[1] = 0 
new.all1[34] = 0.35
new.all1[35] = 16
new.all1[36] = 29
new.all1[30] = 1
new.all1[31] = 0.046
new.all1[32] = 21.5
new.all1[33] = 28.5
new.all1[2] = 0.00085
new.all1[3] = 20.5
new.all1[4] = 31.5
new.all1[23] = 0.00021
new.all1[24] = 17
new.all1[25] = 31
new.all1[37] = 0.15
new.all1[38] = 0
new.all1[39] = 31
new.all1[16] = 0.010
new.all1[17] = 15
new.all1[18] = 35

par(new = TRUE)
R0s.all = sapply(Temps, calc.R0.v2, 1,1, setup = new.all1)
plot(Temps, R0s.all, xlim = c(0,40), ylim = c(0, 18), xlab = "Temperature (°C)", ylab = "Relative R0", col = "#CC79A7", xaxp = c(0, 40, 8), cex.lab=1.5, cex.main=1.5)
legend("topleft", 
       legend = c("Downstate", "Upstate"), 
       col = c("#CC79A7", "#56B4E9"), 
       pch = c(19))
minor.tick(nx=10, ny=2, tick.ratio=0.5, x.args = list(), y.args = list())
