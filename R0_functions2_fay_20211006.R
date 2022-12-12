# SUPPLEMENT S3 to "West Nile virus is predicted to be more geographically
# widespread in New York State and Connecticut under future climate change"

# Provide R functions to access the mean fitted relationships from Shocket et al. 2020
# Author: A.C. Keyel <akeyel@albany.edu>
# Equations and parameters taken from Mordecai et al. 2019 and Shocket et al. 2020.

# DISCLAIMER
# This script is provided as is, with no warranty of fitness or suitability for any purpose

#' Implement a Briere functional relationship
#' 
#' NOTE: Switch to use K instead of C if temps below 0 C are desired as inputs with non-zero outputs
#' 
#' @param Temp The input temperature in C. Should be non-negative, as negative inputs will create an output of 0.
#' @param q scaling coefficient
#' @param Tmin Lower thermal limit
#' @param Tmax Upper thermal limit
#'
#' From Table S1 in Mordecai et al. 2019
#'  
Briere = function(Temp, q, Tmin, Tmax){
  
  # Suppress warnings for sqrt of negative numbers
  if (Tmax < Temp){
    value = 0
  }else{
    value = q*Temp *(Temp - Tmin)*sqrt(Tmax - Temp)
  }
  
  # Constrain value to not be NA or negative
  if(is.na(value)){ value = 0 }
  if(value < 0){    value = 0 }
  # Constrain function to return 0's for negative inputs
  if(Temp  < 0){    value = 0 }
    
  return(value)
}

#' Implement a Reverse Briere functional relationship
#' 
#' NOTE: Switch to use K instead of C if temps below 0 C are desired as inputs with non-zero outputs
#' 
#' @param Temp The input temperature in C. Should be non-negative, as negative inputs will create an output of 0.
#' @param q scaling coefficient
#' @param Tmin Lower thermal limit
#' @param Tmax Upper thermal limit
#'
#' From Table S1 in Mordecai et al. 2019
#'  
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



#' Implement a quadratic function
#' 
#' @param Temp The input temperature in C. Should be non-negative, as negative inputs will create an output of 0.
#' @param q scaling coefficient
#' @param Tmin Lower thermal limit
#' @param Tmax Upper thermal limit
#'
Quadratic = function(Temp, q, Tmin, Tmax){
  
  value = -q*(Temp - Tmin)*(Temp - Tmax)
  
  # Constrain value to not be NA or negative
  if(is.na(value)){ value = 0 }
  if(value < 0){    value = 0 }
  if(Temp  < 0){    value = 0 }
  
  return(value)
}

#https://www.thoughtco.com/normal-distribution-bell-curve-formula-3126278
#**# Should probably look for a more authoritative source, but want to see if it works first
#' Return the shape of a normal distribution
#' 
#' @param x Equivalent to Temp, the temperature at which to evaluate the function
#' @param in.mean The mean for the normal distribution
#' @param in.sd The standard deviation for the normal distribution
#' @param scale.factor A scaling factor to convert from a probability density to the units of interest
#'
normal.curve = function(x, in.mean, in.sd, scale.factor){
  term1 = 1 / (in.sd * sqrt(2*pi))
  term2 = exp(-((x - in.mean)**2/(2*in.sd**2)))
  y = term1 * term2
  y = y * scale.factor
  return(y)
}



#' Implement a linear relationship
#' 
#' @param Temp The input temperature in C. Should be non-negative, as negative inputs will create an output of 0.
#' @param m Slope
#' @param z Intercept
#' 
Linear = function(Temp, m, z){
  value = -m*Temp + z
  
  # Constrain value to be non-negative
  if (value <= 0){ value = 0  }
  if (Temp <= 0) { value = 0  }
  
  return(value)
}


#' Biting Rate a (default= Cx. pipiens & WNV)
#' 
BitingRate = function(Temp, q = 1.70*10**-4, Tmin = 9.4, Tmax = 39.6){
  a = Briere(Temp, q, Tmin, Tmax)
  return(a)
}


#ER.q = 6.36 * 10**-1 # DATA from Cx. quinquefasciatus
#ER.Tmin = 5.0        # DATA from Cx. quinquefasciatus
#ER.Tmax = 37.7 

#' Vector competence (b*c) (default = Cx.pipiens & WNV)
#'
VectorCompetence = function(Temp, q = 3.05*10**-3, Tmin = 16.8, Tmax = 38.9){
  vc = Quadratic(Temp, q, Tmin, Tmax) 
  return(vc)  
}

#' Life span (lf) (default = Cx. pipiens & WNV)
#' 
LifeSpan = function(Temp, m = 4.86, z = 169.8){
  lf = Linear(Temp, m, z)
  return(lf)
}

#' Parasite Development Rate (default = Cx. pipiens & WNV)
#' 
ParasiteDevelopmentRate = function(Temp, q = 7.38*10**-5, Tmin = 11.4, Tmax = 45.2){
  PDR = Briere(Temp, q, Tmin, Tmax) 
  return(PDR)
}

#' Eggs Per Female per Gonotrophic Cycle
#' 
calc_EFGC = function(Temp, q = 0.598, Tmin = 5.3, Tmax= 38.9){
  EFGC = Quadratic(Temp, q, Tmin, Tmax) # 
  return(EFGC)
}

# For Cx. pipiens, Immature survival is given as pLA, rather than pEA.
#' Calculate pLA
#'
calc_pLA = function(Temp, q = 3.60*10**-3, Tmin = 7.8, Tmax = 38.4){
  pLA = Quadratic(Temp, q, Tmin, Tmax) 
  return(pLA)
}

#' Egg Viability
#'
#' EV.type: 1 for Quadratic, 2 for Briere
calc_EV = function(Temp, q = 2.11*10**-3, Tmin = 3.2, Tmax = 42.6, EV.type = 1){
  if (EV.type == 1){ EV = Quadratic(Temp, q, Tmin, Tmax) }
  if (EV.type == 2){ EV = Briere(Temp, q, Tmin, Tmax)    }
  return(EV)  
}

# Egg-to-adult development rate
calc_MDR = function(Temp, q = 3.76*10**-5, Tmin = 0.1, Tmax = 38.5){
  MDR = Briere(Temp, q, Tmin, Tmax) # 
  return(MDR)  
}

# new code form sasha 9.22
#' Calculate R0 for selected species of mosquito and West Nile virus
#' 
#' Note that r and N are not important if a relative R0 is being calculated,
#' as they will cancel out during the scaling process.
#' 
#' @param Temp Temperature in Celsius
#' @param r Host Recovery rate
#' @param N Host density
#' @param species Species of mosquito. Only used if setup is not specified
#' @param R0.type VA: Virus and Abundance components,
#'                V:Virus portion only,
#'                A: Abundance portion only
#' @param setup A vector mirroring the species.setup function output.
#'              Set species to "setup" to use this option
#' 
#' @return R0
calc.R0.v2 = function(Temp, r, N, species = "Culex pipiens", R0.type = "VA",
                      setup = NA){
  
  supported.types = c("VA", "V", "A")
  if (!R0.type %in% supported.types){
    stop(sprintf("%s R0.type is not supported. Supported options are %s",
                 R0.type, paste(supported.types, collapse = ", ")))
  }
  
  #return(c(use.EFGC, a.q, a.Tmin, a.Tmax, vc.q, vc.Tmin, vc.Tmax, lf.m, lf.z, PDR.q,
  #         PDR.Tmin, PDR.Tmax, EFGC.q, EFGC.Tmin, EFGC.Tmax, pLA.q, pLA.Tmin, pLA.Tmax,
  #         EV.type, EV.q, EV.Tmin, EV.Tmax, MDR.q, MDR.Tmin, MDR.Tmax, lf.type, lf.mean,
  #         lf.sd, lf.scale, pO.type, pO.q, pO.Tmin, pO.Tmax, ER.q, ER.Tmin, ER.Tmax))
  
  # If setup is not specified
  if (length(setup) == 1){ setup = species.setup(species) }
  
  use.EFGC =  setup[1]
  a.q =       setup[2]
  a.Tmin =    setup[3]
  a.Tmax =    setup[4]
  vc.q =      setup[5]
  vc.Tmin =   setup[6]
  vc.Tmax =   setup[7]
  lf.m =      setup[8]
  lf.z =      setup[9]
  PDR.q =     setup[10]
  PDR.Tmin =  setup[11]
  PDR.Tmax =  setup[12]
  EFGC.q =    setup[13]
  EFGC.Tmin = setup[14]
  EFGC.Tmax = setup[15]
  pLA.q =     setup[16]
  pLA.Tmin =  setup[17]
  pLA.Tmax =  setup[18]
  EV.type =   setup[19]
  EV.q =      setup[20]
  EV.Tmin =   setup[21]
  EV.Tmax =   setup[22]
  MDR.q =     setup[23]
  MDR.Tmin =  setup[24]
  MDR.Tmax =  setup[25]
  lf.type =   setup[26]
  lf.mean =   setup[27]
  lf.sd =     setup[28]
  lf.scale =  setup[29]
  pO.type =   setup[30]
  pO.q =      setup[31]
  pO.Tmin =   setup[32]
  pO.Tmax =   setup[33]
  ER.q =      setup[34]
  ER.Tmin =   setup[35]
  ER.Tmax =   setup[36]
  lf.q =      setup[37]
  lf.Tmin =   setup[38]
  lf.Tmax =   setup[39]

  # Biting Rate
  a = BitingRate(Temp, a.q, a.Tmin, a.Tmax) 
  # Data given for lower and upper 95% CI's.
  #a.95L = Briere(Temp, 1.18*10**-4, 2.8, 37.9) #Topt 31.3
  #a.95H = Briere(Temp, 2.29*10**-4, 13.4, 40.6) #Topt 33.6
  #**# NOT CURRENTLY USED
    
  # b and c were not used for Cx. pipiens - used vector capacity
  # Proportion of mosquitoes that become infectious
  #b.mean.tarsalis = Quadratic(Temp, 2.94*10**-3, 11.3, 41.9) # 26.6 T.optimum
  
  # Proportion of exposed mosquitoes that become infected
  #c = Quadratic(Temp, 2.56*10**-3, 15.6, 52.2) #33.9 T.optimum
  
  # Vector competence (b*c)
  vc = VectorCompetence(Temp, vc.q, vc.Tmin, vc.Tmax)
  
  if (is.na(lf.type) | lf.type == 1){
    # Life Span #NOTE: in Shocket et al. 2020, this is constrained to be constant below 14 C for Cx. pipiens
    lf = LifeSpan(Temp, lf.m, lf.z) 
 
  }
  if(!is.na(lf.type)) {
    if (lf.type == 2){
      lf = normal.curve(Temp, lf.mean, lf.sd, lf.scale)
    }
    if (lf.type == 3){
      lf = Reverse.Briere(Temp, lf.q, lf.Tmin, lf.Tmax)
    } 
    
  }
  
  # Throw an informative error if lf is not defined
  if (!exists("lf")){
    stop(sprintf("Life Span not correctly processed. Code entered %s. Supported
                 types are 'linear', 'normal' and 'reversebriere'. Please check
                 the species' source setup code for species %s and fix the error",
                 lf.type, species))
  }
  
  # mosquito adult mortality rate (inverse of lifespan)
  u = 1/lf
  
  # Parasite development rate (1/(extrinsic incubation period))
  PDR = ParasiteDevelopmentRate(Temp, PDR.q, PDR.Tmin, PDR.Tmax)
  
  if (use.EFGC == 1){
    # Gives fecundity for pipiens as (EFGC) Eggs per female per gonotrophic cycle
    EFGC = calc_EFGC(Temp, EFGC.q, EFGC.Tmin, EFGC.Tmax)
  
    # Convert to Eggs per female per day
    EFD = EFGC * a # a is the inverse of gonotrophic cycle
  
  # Otherwise use eggs per raft (for Cx. quinquefasciatus)
  }else{
    # Gives proportion ovipositing.
    if (pO.type == 1){ p.ovipositing = Quadratic(Temp, pO.q, pO.Tmin, pO.Tmax) }
    if (pO.type == 2){ p.ovipositing = Briere(Temp, pO.q, pO.Tmin, pO.Tmax)    }
  
    ER = Quadratic(Temp, ER.q, ER.Tmin, ER.Tmax)
    
    
  
    # Eggs per female per day
    EFD = p.ovipositing * ER * a    
  }
  
  # Larval surival probability
  pLA = calc_pLA(Temp, pLA.q, pLA.Tmin, pLA.Tmax)
  
  # Egg Viability
  EV = calc_EV(Temp, EV.q, EV.Tmin, EV.Tmax, EV.type)
  
  # pEA = survival probability (eggs to adult)
  pEA = EV * pLA
  
  # Egg-to-adult development rate
  MDR = calc_MDR(Temp, MDR.q, MDR.Tmin, MDR.Tmax)
  
  
  if (R0.type == "VA"){
    term1 = a**2 * vc * exp(-u/PDR) * EFD * pEA * MDR
    term2 = N * r * u**3
  }
  
  if (R0.type == "V"){
    term1 = a**2 * vc * exp(-u/PDR)
    term2 = u
  }
  if (R0.type == "A"){
    term1 = EFD * pEA * MDR
    term2 = u**2
  }
  
  R0 = sqrt(term1 / term2)
  
  # Ensure that R0 is not NaN (if temp is 0, R0 will return NaN due to exp(-u/PDR) returning inf)
  if (is.na(R0)){ R0 = 0 }
  
  return(R0)  
}

#' Helper function to move species-specific info out of the R0 function
#' 
#' @param species The species of interest.
#' @param lf.type Linear fit for life span: NA or 1. Normal fit for life span: 2. Reverse Briere fit: 3
#' 
species.setup = function(species){
  supported.species = c("Culex pipiens", "Culex restuans", "Culex univittatus",
                        "Culex tarsalis", "Culex quinquefasciatus", "Culiseta melanura", 
                        "Culex pipiens Syracuse all", "Culex pipiens Suffolk BR", "Culex pipiens Syracuse ls","Culex pipiens Syracuse br ls mdr",
                          "Culex pipiens Suffolk ls","Culex pipiens Suffolk br ls mdr pO ER", "Culex pipiens Syracuse all traits",
                         "Culex pipiens Suffolk all traits" )
  # "Culex tarsalis": Not all traits fit, decided to do pipiens with univittatus vector competence (but not PDR)
  # Not all traits fit, old world species, similar pattern for Cx. tarsalis, so I filled in that species instead.
  if (!species %in% supported.species){
    stop(sprintf("Species %s is not supported. Supported species are %s", 
                 species, paste(supported.species, collapse = ", ")))
  }
  
  # Set added parameters to NA for all other species
  lf.type = NA
  lf.mean = NA
  lf.sd = NA
  lf.scale = NA
  EV.type = 1 # Set default to quadratic (1)
  pO.type = NA
  pO.q = NA
  pO.Tmin = NA
  pO.Tmax = NA
  ER.q = NA
  ER.Tmin = NA
  ER.Tmax = NA
  EFGC.q = NA
  EFGC.Tmin = NA
  EFGC.Tmax = NA
  lf.q = NA
  lf.Tmin = NA
  lf.Tmax = NA
  
  # Syracuse with all 
  if (species == "Culex pipiens Syracuse all traits"){
    lf.type = 3 #"reversebriere"
    use.EFGC = 1
    a.q = 0.00155
    a.Tmin = 20.5
    a.Tmax = 28.5
    vc.q = 3.05*10**-3
    vc.Tmin = 16.8
    vc.Tmax = 38.9
    lf.m = 4.86
    lf.z = 169.8
    lf.q = 0.16
    lf.Tmin = 0
    lf.Tmax = 30
    PDR.q = 7.38*10**-5
    PDR.Tmin = 11.4
    PDR.Tmax = 45.2
    EFGC.q = 0.598
    EFGC.Tmin = 5.3
    EFGC.Tmax = 38.9
    pLA.q = 3.60*10**-3
    pLA.Tmin = 7.8
    pLA.Tmax = 38.4
    EV.q = 2.11*10**-3
    EV.Tmin = 3.2
    EV.Tmax = 42.6
    MDR.q = 0.011
    MDR.Tmin = 5
    MDR.Tmax = 29
    ER.q = 0.43
    ER.Tmin = 16
    ER.Tmax = 29
    pO.q = 0.023
    pO.Tmin = 21.5
    pO.Tmax = 28.5
  } 
 
 # Culex pipiens SYR LS = "reversebriere" with biting rate
  if (species == "Culex pipiens Syracuse BR"){
    use.EFGC = 1
    a.q = 0.00155
    a.Tmin = 20.5
    a.Tmax = 28.5
    vc.q = 3.05*10**-3
    vc.Tmin = 16.8
    vc.Tmax = 38.9
    lf.m = 4.86
    lf.z = 169.8
    PDR.q = 7.38*10**-5
    PDR.Tmin = 11.4
    PDR.Tmax = 45.2
    EFGC.q = 0.598
    EFGC.Tmin = 5.3
    EFGC.Tmax = 38.9
    pLA.q = 3.60*10**-3
    pLA.Tmin = 7.8
    pLA.Tmax = 38.4
    EV.q = 2.11*10**-3
    EV.Tmin = 3.2
    EV.Tmax = 42.6
    MDR.q = 3.76*10**-5
    MDR.Tmin = 0.1
    MDR.Tmax = 38.5
    # ER =  
  }
  
  # Culex pipiens SYR LS = "reversebriere" with biting rate
  if (species == "Culex pipiens Syracuse ls br"){
    lf.type = 3 # "reversebriere"
    use.EFGC = 1
    a.q = 0.0017
    a.Tmin = 20.5
    a.Tmax = 28.1
    vc.q = 3.05*10**-3
    vc.Tmin = 16.8
    vc.Tmax = 38.9
    lf.m = NA
    lf.z = NA
    lf.q = 0.16
    lf.Tmin = 0
    lf.Tmax = 30
    PDR.q = 7.38*10**-5
    PDR.Tmin = 11.4
    PDR.Tmax = 45.2
    EFGC.q = 0.598
    EFGC.Tmin = 5.3
    EFGC.Tmax = 38.9
    pLA.q = 3.60*10**-3
    pLA.Tmin = 7.8
    pLA.Tmax = 38.4
    EV.q = 2.11*10**-3
    EV.Tmin = 3.2
    EV.Tmax = 42.6
    MDR.q = 3.76*10**-5
    MDR.Tmin = 0.1
    MDR.Tmax = 38.5
  }
  
  # Culex pipiens SYR LS = "reversebriere" NO biting rate
  if (species == "Culex pipiens Syracuse ls"){
    lf.type = 3 # "reversebriere"
    use.EFGC = 1
    a.q = 1.70*10**-4
    a.Tmin = 9.4
    a.Tmax = 39.6
    vc.q = 3.05*10**-3
    vc.Tmin = 16.8
    vc.Tmax = 38.9
    lf.m = NA
    lf.z = NA
    lf.q = 0.16
    lf.Tmin = 0
    lf.Tmax = 30
    PDR.q = 7.38*10**-5
    PDR.Tmin = 11.4
    PDR.Tmax = 45.2
    EFGC.q = 0.598
    EFGC.Tmin = 5.3
    EFGC.Tmax = 38.9
    pLA.q = 3.60*10**-3
    pLA.Tmin = 7.8
    pLA.Tmax = 38.4
    EV.q = 2.11*10**-3
    EV.Tmin = 3.2
    EV.Tmax = 42.6
    MDR.q = 3.76*10**-5
    MDR.Tmin = 0.1
    MDR.Tmax = 38.5
  }
  #0.00085, 20.5, 31.5 BR
  if (species == "Culex pipiens Suffolk BR"){
    use.EFGC = 1
    a.q = 0.00085
    a.Tmin = 20.5
    a.Tmax = 31.5
    vc.q = 3.05*10**-3
    vc.Tmin = 16.8
    vc.Tmax = 38.9
    lf.m = 4.86
    lf.z = 169.8
    PDR.q = 7.38*10**-5
    PDR.Tmin = 11.4
    PDR.Tmax = 45.2
    EFGC.q = 0.598
    EFGC.Tmin = 5.3
    EFGC.Tmax = 38.9
    pLA.q = 3.60*10**-3
    pLA.Tmin = 7.8
    pLA.Tmax = 38.4
    EV.q = 2.11*10**-3
    EV.Tmin = 3.2
    EV.Tmax = 42.6
    MDR.q = 3.76*10**-5
    MDR.Tmin = 0.1
    MDR.Tmax = 38.5
  }
  #q.Rev = 0.15 Tmin.Rev = 0 Tmax.Rev = 31
  if (species == "Culex pipiens Suffolk ls"){
    lf.type = 3 # "reversebriere"
    use.EFGC = 1
    a.q = 1.70*10**-4
    a.Tmin = 9.4
    a.Tmax = 39.6
    vc.q = 3.05*10**-3
    vc.Tmin = 16.8
    vc.Tmax = 38.9
    lf.m = NA
    lf.z = NA
    lf.q = 0.15
    lf.Tmin = 0
    lf.Tmax = 31
    PDR.q = 7.38*10**-5
    PDR.Tmin = 11.4
    PDR.Tmax = 45.2
    EFGC.q = 0.598
    EFGC.Tmin = 5.3
    EFGC.Tmax = 38.9
    pLA.q = 3.60*10**-3
    pLA.Tmin = 7.8
    pLA.Tmax = 38.4
    EV.q = 2.11*10**-3
    EV.Tmin = 3.2
    EV.Tmax = 42.6
    MDR.q = 3.76*10**-5
    MDR.Tmin = 0.1
    MDR.Tmax = 38.5
  }
  # suffolk br and ls ALL 
  if (species == "Culex pipiens Suffolk br ls mdr pO ER"){
    lf.type = 3 # "reversebriere"
    use.EFGC = 1
    a.q = 0.00085
    a.Tmin = 20.5
    a.Tmax = 31.5
    vc.q = 3.05*10**-3
    vc.Tmin = 16.8
    vc.Tmax = 38.9
    lf.m = NA
    lf.z = NA
    lf.q = 0.15
    lf.Tmin = 0
    lf.Tmax = 31
    PDR.q = 7.38*10**-5
    PDR.Tmin = 11.4
    PDR.Tmax = 45.2
    EFGC.q = 0.598
    EFGC.Tmin = 5.3
    EFGC.Tmax = 38.9
    pLA.q = 3.60*10**-3
    pLA.Tmin = 7.8
    pLA.Tmax = 38.4
    EV.q = 2.11*10**-3
    EV.Tmin = 3.2
    EV.Tmax = 42.6
    MDR.q = 0.013
    MDR.Tmin = 5
    MDR.Tmax = 29
    ER.q = 0.35
    ER.Tmin = 16
    ER.Tmax = 29
    pO.q = 0.046
    pO.Tmin = 21.5
    pO.Tmax = 28.5
  }
  if (species == "Culex pipiens"){
    use.EFGC = 1
    a.q = 1.70*10**-4
    a.Tmin = 9.4
    a.Tmax = 39.6
    vc.q = 3.05*10**-3
    vc.Tmin = 16.8
    vc.Tmax = 38.9
    lf.m = 4.86
    lf.z = 169.8
    PDR.q = 7.38*10**-5
    PDR.Tmin = 11.4
    PDR.Tmax = 45.2
    EFGC.q = 0.598
    EFGC.Tmin = 5.3
    EFGC.Tmax = 38.9
    pLA.q = 3.60*10**-3
    pLA.Tmin = 7.8
    pLA.Tmax = 38.4
    EV.q = 2.11*10**-3
    EV.Tmin = 3.2
    EV.Tmax = 42.6
    MDR.q = 3.76*10**-5
    MDR.Tmin = 0.1
    MDR.Tmax = 38.5
  }
  
  if (species == "Culex restuans"){
    # Traits with no data used values for Cx. pipiens
    # Except vector competence, which used values for Cx. univittatus
    use.EFGC = 1
    a.q = 1.70*10**-4
    a.Tmin = 9.4
    a.Tmax = 39.6
    vc.q = 2.32*10**-3 # Cx. univittatus
    vc.Tmin = 4.2      # Cx. univittatus
    vc.Tmax = 45.2     # Cx. univittatus
    lf.m = NA          # Use a normal distribution
    lf.z = NA          # Use a normal distribution
    PDR.q = 7.54*10**-5# Cx. univittatus
    PDR.Tmin = 10.2    # Cx. univittatus
    PDR.Tmax = 34.4    # Cx. univittatus
    EFGC.q = 0.598
    EFGC.Tmin = 5.3
    EFGC.Tmax = 38.9
    pLA.q = 3.60*10**-3 # Manual fit to Ciota et al. 2014 - good enough fit, rather than 'best' fit
    pLA.Tmin = 4.5      # Manual fit to Ciota et al. 2014 - good enough fit, rather than 'best' fit
    pLA.Tmax = 35.2     # Manual fit to Ciota et al. 2014 - good enough fit, rather than 'best' fit
    EV.q = 2.11*10**-3
    EV.Tmin = 3.2
    EV.Tmax = 42.6
    MDR.q = 3.76*10**-5 # As bad for Cx. restuans as it was for Cx. pipiens based on Ciota et al. 2014
    MDR.Tmin = 0.1      # As bad for Cx. restuans as it was for Cx. pipiens based on Ciota et al. 2014
    MDR.Tmax = 38.5     # As bad for Cx. restuans as it was for Cx. pipiens based on Ciota et al. 2014
    lf.type = 1      # Fit manually to Ciota et al. 2014 - good enough fit, rather than 'best' fit
    lf.mean = 6.5      # Fit manually to Ciota et al. 2014 - good enough fit, rather than 'best' fit
    lf.sd = 9          # Fit manually to Ciota et al. 2014 - good enough fit, rather than 'best' fit
    lf.scale = 3000    # Fit manually to Ciota et al. 2014 - good enough fit, rather than 'best' fit
    
  }
  
  if (species == "Culex univittatus"){
    # ALL TRAITS FROM Cx pipiens EXCEPT vc and PDR
    #**# NOT CERTAIN THIS WAS Mordecai/Shocket's approach, but see 
    # Appendix 1 Fig. 14 from Shocket. But elsewhere, they say they fit
    # traits with trait values from all species except the focal species as
    # priors, and then re-fit with the species-specific values.
    use.EFGC = 1
    a.q = 1.70*10**-4
    a.Tmin = 9.4
    a.Tmax = 39.6
    vc.q = 2.32*10**-3
    vc.Tmin = 4.2
    vc.Tmax = 45.2
    lf.m = 4.86
    lf.z = 169.8
    # PDR 
    PDR.q = 7.54*10**-5
    PDR.Tmin = 10.2
    PDR.Tmax = 34.4
    EFGC.q = 0.598
    EFGC.Tmin = 5.3
    EFGC.Tmax = 38.9
    pLA.q = 3.60*10**-3
    pLA.Tmin = 7.8
    pLA.Tmax = 38.4
    EV.q = 2.11*10**-3
    EV.Tmin = 3.2
    EV.Tmax = 42.6
    MDR.q = 3.76*10**-5
    MDR.Tmin = 0.1
    MDR.Tmax = 38.5
    
  }
  
  if (species == "Culiseta melanura"){
    # When data were missing for a vector–virus pair, we substituted the most
    # conservative (i.e. least restrictive of transmission) trait thermal response from a vector
    # that occurs within the geographic range of disease transmission.
    # Commented lines are from Cx. pipiens
    use.EFGC = 0
    a.q = 1.87*10**-4
    a.Tmin = 7.8
    a.Tmax = 31.8
    vc.q = 1.51*10**-3 # DATA for EEE from Ae. triseriatus
    vc.Tmin = 7.0      # DATA for EEE from Ae. triseriatus
    vc.Tmax = 50.3     # DATA for EEE from Ae. triseriatus
    lf.m = 4.86        # DATA FROM Cx. pipiens
    lf.z = 169.8       # Data from Cx. pipiens
    PDR.q = 7.05*10**-5 # DATA for EEE from Ae. triseriatus
    PDR.Tmin = 11.6     # DATA for EEE from Ae. triseriatus
    PDR.Tmax = 44.8     # DATA for EEE from Ae. triseriatus
    EFGC.q = NA
    EFGC.Tmin = NA
    EFGC.Tmax = NA
    pLA.q = 3.03*10**-3
    pLA.Tmin = 10.1
    pLA.Tmax = 36.2
    EV.q = 2.11*10**-3 # DATA FROM Cx. pipiens
    EV.Tmin = 3.2      # DATA FROM Cx. pipiens
    EV.Tmax = 42.6     # DATA from Cx. pipiens
    MDR.q = 2.74*10**-5
    MDR.Tmin = 8.6
    MDR.Tmax = 37.6
    #**# These require eggs per raft to be used, and those exist only for Cx. quinquefasciatus
    pO.type = 1 # 1 for Quadratic relationship, 2 for Briere
    pO.q = 6.31*10**-3 
    pO.Tmin = 8.7
    pO.Tmax = 33.6
    ER.q = 6.36 * 10**-1 # DATA from Cx. quinquefasciatus
    ER.Tmin = 5.0        # DATA from Cx. quinquefasciatus
    ER.Tmax = 37.7       # DATA from Cx. quinquefasciatus
  }
  
  
  # Uses same shape equations as Cx. pipiens
  if (species == "Culex tarsalis"){
    a.q = 0.000167
    a.Tmin = 2.3
    a.Tmax = 32.0
    #Used vc from Cx. pipiens, as only b was available for Cx. tarsalis
    vc.q = 3.05*10**-3
    vc.Tmin = 16.8
    vc.Tmax = 38.9
    lf.m = 1.69
    lf.z = 69.6
    PDR.q = 0.0000657
    PDR.Tmin = 11.2
    PDR.Tmax = 44.7
    # EFGC used pipiens fit (see appendix 1 Fig 13)
    use.EFGC = 1 
    EFGC.q = 0.598
    EFGC.Tmin = 5.3
    EFGC.Tmax = 38.9
    pLA.q = 0.00212
    pLA.Tmin = 5.9
    pLA.Tmax = 43.1
    # EV used pipiens fit (see appendix 1 Fig. 13)
    EV.q = 2.11*10**-3
    EV.Tmin = 3.2
    EV.Tmax = 42.6
    MDR.q = 0.0000412
    MDR.Tmin = 4.3
    MDR.Tmax = 39.9
  }
  

  if (species == "Culex quinquefasciatus"){
    use.EFGC = 0
    a.q = 0.0000728
    a.Tmin = 3.1
    a.Tmax = 39.3
    # Used vc from Cx. pipiens, not available for Cx. quinquefasciatus
    vc.q = 3.05*10**-3
    vc.Tmin = 16.8
    vc.Tmax = 38.9
    lf.m = 3.80
    lf.z = 136.3
    PDR.q = 0.0000712
    PDR.Tmin = 19.0
    PDR.Tmax = 44.1
    # Use pipiens fit (see appendix 1 Fig 13)
    ER.q = 6.36 * 10**-1 
    ER.Tmin = 5.0        
    ER.Tmax = 37.7 
    pO.type = 2 # 2 for Briere "B"
    pO.q = 0.000667
    pO.Tmin = 1.7
    pO.Tmax = 31.8
    pLA.q = 0.00426
    pLA.Tmin = 8.9
    pLA.Tmax = 37.7
    # Use pipiens fit (see appendix 1 Fig. 13)
    EV.type = 2 #"B" 1 for Quadratic, 2 for Briere
    EV.q = 0.0047
    EV.Tmin = 13.6
    EV.Tmax = 38.0
    MDR.q = 0.0000414
    MDR.Tmin = 0.1
    MDR.Tmax = 38.6
    
  }
    
  return(c(use.EFGC, a.q, a.Tmin, a.Tmax, vc.q, vc.Tmin, vc.Tmax, lf.m, lf.z, PDR.q,
         PDR.Tmin, PDR.Tmax, EFGC.q, EFGC.Tmin, EFGC.Tmax, pLA.q, pLA.Tmin, pLA.Tmax,
         EV.type, EV.q, EV.Tmin, EV.Tmax, MDR.q, MDR.Tmin, MDR.Tmax, lf.type, lf.mean,
         lf.sd, lf.scale, pO.type, pO.q, pO.Tmin, pO.Tmax, ER.q, ER.Tmin, ER.Tmax,
         lf.q, lf.Tmin, lf.Tmax))
}


 # END OF FILE