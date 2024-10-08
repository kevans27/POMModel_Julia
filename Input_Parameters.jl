# Input parameters

global r = 0.05
global mybeta = 150

#you can either set total cells per particle or cell density at formation
global cell_per_particle = 1.5E4 #initial cells particle-1 
global cell_density_init = cell_per_particle/(4*pi*r*r*1E-4) #cells.m-2.particle-1

global F0_max = 5e10#5e10; #cells/m3  concentration of free living bact community 

global umax = 3.7  #day-1

global runName = "test"