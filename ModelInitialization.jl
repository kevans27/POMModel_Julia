#= Carbon dynamics Model per PARTICLE base
Model input parameters

Modified from the model by Trang Nguyen
Modified: May 2, 2024

-------------------------------------------------------
 LOAD IN USER DEFINED PARAMETERS
-------------------------------------------------------=#

# Input parameters
global r = 0.05
global mybeta = 150

#you can either set total cells per particle or cell density at formation
global cell_per_particle = 1.5E4 #initial cells particle-1 
global cell_density_init = cell_per_particle/(4*pi*r*r*1E-4) #cells.m-2.particle-1

global F0_max = 5e10#5e10; #cells/m3  concentration of free living bact community 

global umax = 3.7  #day-1

global runName = "noAttach"

#=-------------------------------------------------------
           MODEL PARAMS                               
-------------------------------------------------------

Constants =#
day_to_sec = 24*60*60;
D = 6.7E-6*3600*24;#cm2.day-1 #Diffusion of monomer (D of glucose)
viscosity = 1E-2;#cm2.s-1 #Kioebe 2001 
sc = viscosity/(D/3600/24); #Schmidt number

#Depth and Particle Radius
depth = 100; #m  Depth of POM formation
rcell = 6.2e-5; #cm  for Vol_cell = 1um3

# POM density and sinking rate calculation
Af = 1.49; #unitless
D_Omand = 1.85;#1.85; #aggregate;  = 3 for small solid particle
D_s = 3; #small particle
c_avg_pa_ref = 1.23; #g/cm3, McDonell 2010, Sicko-Goad 1980, Cram 2018 - Measured pom density for small solid particles
fluid_density_ref = 1.028;#g/cm3, Omand 2020, at the base of euphotic zone, assumed to be constant in the density equation
g = 9.81; #m/s2
mu = 1e-3; #kg/m/s
rs_3_minus_d = 0.0012*18*mu/(g*Af*(c_avg_pa_ref-fluid_density_ref)*1000*(1e-3)^(D_Omand-1));#rs is in meter; in which the reference particle is 1e-3 m with 0.0012 m/s sinking rate
cut_off = rs_3_minus_d^(1/(3-D_Omand))*100;#cm, aka radius of the sub-unit
volu = 4/3*pi*r^3;

#Were originally treated as arrays, but current use is restricted to single values
f = Af*rs_3_minus_d*(r/100).^(D_Omand-3);#unitless, Omand 2020, Fraction of POM in an aggregate
f = ifelse(f >= 1, 1, f)
f = ifelse(r <= cut_off, 1, f)


c_avg_pa_omand = f*(c_avg_pa_ref-fluid_density_ref)+fluid_density_ref; #g/cm3, Omand 2020, POM Note4 Trang
c_avg_pa_omand = ifelse(r <= cut_off, c_avg_pa_ref, c_avg_pa_omand)
c_avg_pa_omand = ifelse(c_avg_pa_omand >= c_avg_pa_ref, c_avg_pa_ref, c_avg_pa_omand)#Assuming all particles > constant density to have constant density
pseudo_pom = f.*c_avg_pa_ref.*volu/12*1000;#mmolC/particle, Omand 2020, POM Note4 Trang
convert_Alldredge_Omand = 0.0494;#hard code the B ~ regression between Omand C content and Alldredge C content to account for non-org C fraction in a particle
c_avg_pa = pseudo_pom*convert_Alldredge_Omand;#save this B: 0.0494

w_aggregate = g*Af*(c_avg_pa_ref-fluid_density_ref)*1000*(r/100).^(D_Omand-1)*rs_3_minus_d/(18*mu)*100;#cm/s, Omand 2020, POM Note4 Trang, Sinking rate based on density at certain depth and based on the change of fluid density at depth
w_smallpa = g*Af*(c_avg_pa_ref-fluid_density_ref)*1000*(r/100).^(D_s-1)/(18*mu)*100;#cm/s, Sinking rate for solid particle fraction dimension D=3 aka a_s^(D-3)=1
if r <= cut_off
    Ws_cms = w_smallpa; #cm/s
else
    Ws_cms = w_aggregate;#cm/s
end
Ws = Ws_cms*864; #m/day
Ws1 = Ws;

rey = Ws_cms * r / viscosity; #unitless, Reynold number #Reynold Number (Moradi 2018, describe how viscous the particle is compared to the environment, depends on the diffusion and advection on the particle #Re=sinking*radius/water viscosity (high Re=High sinking)
sc = viscosity/(D/3600/24); #Schmidt number, unitless, Bianche 2018 Nat Geo (Supp)
sh = 1+0.619*rey^0.412*sc^(1/3);#Sherwood number, unitless, Bianche 2018 Nat Geo (Supp)
diffloss = 4*pi*D*sh*r*1e-6;#m3/day/particle

#= Bacterial initialization

Initial population =#
nbft = 2;#number of bacterial groups so that initial colonizers and 'recruits' tracked seperately
cell_c_content = 20E-15/12*1E3; #mmolC.cell-1 #Lee and Fuhrman 1987, unit fg into mmol C
initba = cell_per_particle*cell_c_content;#mmolC.particle-1

#Encounter rate
newAttachedCells = F0_max*cell_c_content;#mmolC of cell/m3
#newAttachedCells = 0;#mmolC of cell/m3


#= Bacterial dynamics and behavior - Global variables

Motility (For Encounter Rate calculation)=#
motility = (1e-9 + 1e-10)/2*60; #m2/minute

#Mortality rates
mlin = 2.5;#linear mortality for b

#Detachment (Leaving) rate
L0 = 0.25;

#Metabolism for Aerobic Heterotrophs
ymom = 0.2; #mmolC_cells/mmolC_MOM - aerobic bacteria yield
vmom_max = umax/ymom;#day-1
kmom = 7.2;#mmolC_MOM.m-3; Half saturation for MOM
 
# Initialize following recruitment experiment
betas = mybeta * ones(1, nbft);
bo_pas = [initba, cell_c_content]#early colonizers, recruited cells ~ late comers
#=
-------------------------------------------------------
           PARAMS REGARDING TO SOLVING EQUATIONS       
-------------------------------------------------------
 TimeSteps=#
ndays = 5;#day
dt = 1/24/6;#1/24/6/1; #day/step ~ apx every 10 mins
timestep = round(ndays/dt);

# Temperature Function through depth
varTemp = true;   #1 means variable temp, 0 means constant temp
TempAeArr = -4E3;  
TemprefArr = 293.15 ;      
Tkel = 273.15 ;  
TempCoeffArr = 0.8;

Temp = 12*exp(-depth/150)+12*exp(-depth/500)+2;
TempFun = TempCoeffArr *exp(TempAeArr *(1 /(Temp+Tkel)-1 /TemprefArr));

if varTemp == 0
    TempFun = 1;
end


#=
-------------------------------------------------------
           INITIALIZE VALUES FOR PATCH                
-------------------------------------------------------=#

tsInt = Integer(timestep)

PATCH = Dict(
    "pom" => zeros(tsInt+1),              # Preallocate a vector for `pom` with `timesteps` elements
    "mom" => zeros(tsInt+1),              # Preallocate a vector for `mom` with `timesteps` elements
    "z" => zeros(tsInt+1),                # Preallocate a vector for `z` (depth) with `timesteps` elements
    "bo_pa" => zeros(length(bo_pas), tsInt+1),  # Preallocate a matrix for `bo_pa` with `length(bo_pas)` rows and `tsInt+1` columns
    "celldensity" => fill(cell_density_init / nbft, nbft),
    "beta" => copy(betas),  # or betas[:]
    "pool" => nbft,
    "depth" => zeros(tsInt+1),            # Preallocate a vector for `depth` with `timesteps` elements
    "time" => zeros(tsInt+1),             # Preallocate a vector for `time` with `timesteps` elements
    "Sh" => zeros(tsInt),             # Preallocate a vector for `time` with `timesteps` elements
    "r" => zeros(tsInt+1),                # Preallocate a vector for `r` with `timesteps` elements
    "newcells" => newAttachedCells,         # Preallocate a vector for `newcells` with `timesteps` elements
    "varTemperature" => varTemp,   # Preallocate a vector for `varTemperature` with `timesteps` elements
    "vmax" => vmom_max * ymom - mlin,             # Preallocate a vector for `vmax` with `timesteps` elements
    "km" => kmom,               # Preallocate a vector for `km` with 1 element
    "mlinear" => mlin           # Preallocate a vector for `mlinear` with `timesteps` elements
)

PATCH["pom"][1] = c_avg_pa #mmolC.particle-1
PATCH["mom"][1] = 0 #mmolC.m-3
PATCH["bo_pa"][:, 1] = bo_pas #mmolC.m-3 # Assign `bo_pas` to the first column of the matrix
PATCH["z"][1] = depth #m
PATCH["depth"][1] = depth #m
PATCH["time"][1] = 0 #day
PATCH["r"][1] = r #cm
#PATCH["newcells"][1] = newAttachedCells #mmolC of cell/m3
#PATCH["varTemperature"][1] = TempFun #unitless
#PATCH["vmax"][1] = vmom_max * ymom - mlin #mmolC.m-3.day-1
#PATCH["km"][1] = kmom #mmolC.m-3
#PATCH["mlinear"][1] = mlin


#Global variables
module ModelInitialization
export nbft, initba, c_avg_pa, rcell
export mlin, D, kmom, ymom, vmom_max, varTemp
export depth, L0, newAttachedCells, cell_c_content, cell_per_particle
export day_to_sec, motility
export Af, D_Omand, D_s, c_avg_pa_ref, fluid_density_ref, g, mu, rs_3_minus_d, cut_off, volu, convert_Alldredge_Omand #For Omand density and sinking equation
export viscosity, sc
export betas, Ws1
export ndays, timestep, dt, TempFun
export TempAeArr, TemprefArr, Tkel, TempCoeffArr
export PATCH

end
