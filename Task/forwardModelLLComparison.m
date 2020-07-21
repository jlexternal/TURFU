forwardInfNoiseLL = 0;
postdicInfNoiseLL = 0;
forwardSelNoiseLL = 0;

for i = 1
    infNoiseGlazeForward_sim;
    forwardInfNoiseLL = forwardInfNoiseLL + totalLL;
    
    infNoiseGlaze_sim;
    postdicInfNoiseLL = postdicInfNoiseLL + totalLL;
    
    selNoiseGlazeForward_sim;
    forwardSelNoiseLL = forwardSelNoiseLL + totalLL;
    
    if mod(i,100) == 0
        disp(['Simulation number: ' num2str(i)]);
    end
end

forwardInfNoiseLL = forwardInfNoiseLL.*-1;
postdicInfNoiseLL = postdicInfNoiseLL.*-1;
forwardSelNoiseLL = forwardSelNoiseLL.*-1;

exp(forwardInfNoiseLL./postdicInfNoiseLL)

bar([1:3], [forwardInfNoiseLL postdicInfNoiseLL forwardSelNoiseLL]);