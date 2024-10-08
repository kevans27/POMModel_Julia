#=--------------------------------------------------%
% Modified createVector for Trang simplified model %
%                Not working yet                   %
%--------------------------------------------------=#

#include("ModelInitialization.jl") 


function createVector(NPATCH, j)

    dCdt = zeros(Float64, 3 + nbft) # Initialize the array with zeros

    dCdt[1]= NPATCH["pom"][j]; #total pool
    dCdt[2]= NPATCH["mom"][j]; #total pool
    dCdt[3]= NPATCH["z"][j]; #depth over time inside myRK
    #dCdt(4)= NPATCH.pomz(i); #POM over depth over time inside myRK

    for ibft=1:nbft 
        dCdt[3+ibft]= NPATCH["bo_pa"][ibft,j] #dCdt(4+ibft)= NPATCH.bo_pa(ibft,end);
    end

    return dCdt
end