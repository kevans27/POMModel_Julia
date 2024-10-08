#=-----------------------------------------------------%
% Modified createStructure for Trang simplified model %
%                Not working yet                      %
%-----------------------------------------------------=#

function createStructure(NPATCH, dCdt, t)
global nbft

NPATCH["pom"][t]=dCdt[1][1]   #total pool
NPATCH["mom"][t]=dCdt[1][2]   #total pool
NPATCH["z"][t]=dCdt[1][3]    #depth over time, inside myRK
#NPATCH.pomz[t]=dCdt[4] #POM over depth over time, inside myRK

for ibft in 1:nbft 
    NPATCH["bo_pa"][ibft,t]=dCdt[1][3+ibft]   #dCdt(4+ibft)
end

return NPATCH
end