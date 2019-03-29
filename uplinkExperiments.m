clear all
%this version is paired down and optimized. use old one for more free
%params
tic

%% Save file name
%it's always a good idea to save your data

saveFileName = 'myExperiment';

%% Waveform parameters
%you can change these params

f = 42.37e9; %all calculations are done wrt. to center frequency
            %this is, of course, an approximation.
lambda = 3e8/f;
bandwidth = 2e9; %passband

txPowerUE = 10;

polarizationAngle = pi/4; %pi/4 = circular, only used for rainfall 
                          %attenuation estimates
  
%% Antenna system parameters
%you can change these params, be sure to update the unique string!

patchElement= rectangularPatch(f,9.8,.025*lambda); %alumina substrate
patchPatternSq = @(phi,theta) patchElement.patternSq(phi,theta);

%isotropicPatternSq = @(phi,theta) 1;
%dipolePatternSq =...
...@(phi,theta) (cos(pi*cos(theta)/2)./sin(theta)).^2; %z oriented

patchUniqueString = ...
sprintf('%f_Hz,%f_perm,%f*lambda_thickness_CAB_optimalDesign_',f,9.8,.025);

numAPAntY = 4; %access point subarray dimensions
numAPAntZ = 4;

numUEAntY = 2; %UE subarray dimensions
numUEAntZ = 2;

%% Noise Temperature of Access Point Recievers
%you can change these params

rxNoiseTemp = 1200;

%% Environmental parameters 
%you can change these params

temperature = 295; %ambient temperature in Kelvin

surfaceTemperature = 295; %surface temperature in Kelvin

dryPressure = 1013.15;%atmospheric pressure due to dry air in hPa
                      %approximately equal to total barometric pressure
                      %but if you want something more exact or use
                      %radiosonde data see the conversion scripts.

RH = 50; %relative humidity in percent.

rainRate = 0;%rate of rainfall in mmph
fogDensity = 0;%gpm^3 


%% Monte Carlo and stochastic geometry setup
%you can change these params

numTrials = 100;

maxNumUsers = 25; %max number of users in the network
coneMaxAngle = pi/6; %angular of spherical sector served
maxRadius = 10;  %max radius of spherical sector served
minRadius = .001; %min radius of spherical sector served

maxTilt = pi/3; %tilt of UE array planes
                %more documentation on this feature forthcoming.
                %we used maxTilt = 0 in our paper.
                %max tilt should be less than pi/2 
                %no guarantee this works otherwise... 



%% Derived parameters for attenuation and noise.

%We overestimate the noise from atmospheric effects via assuming that the
%antenna noise temperature is equal to the atmospheric mean radiating
%temperature. This is accurate for beams steered along the horizontal or at
%frequencies close to 60 GHz. We estimate the mean radiating temperature
%from the surface temperature via an ITU standard. Also, the same standard
%indicates that the noise from terrestrial sources is comparable. 

%Reference: "Radio noise", ITU Std. ITU-R P.372-12, 2015.

meanRadiatingTemp = surfaceTemperature*.81+37.34; %Estimates noise from 
                                                  %atmospheric radiators
                                                 
rxNoisePower = (bandwidth*(meanRadiatingTemp+rxNoiseTemp)*1.38064852e-23);

txSNR = txPowerUE/rxNoisePower;

[atmosphericSA,~] = atmosphericAttenuation(f,dryPressure,relativeHumidityDryPressureToVaporPressure(RH, dryPressure,temperature), temperature);
fogSA = fogAttenuation( f, temperature, fogDensity);


%% Just some initialization
SINRAnalogOnly = zeros(numTrials,maxNumUsers);
perUserCapacityAnalogOnly = zeros(numTrials,maxNumUsers);
sumRateCapacityAnalogOnly = zeros(numTrials,maxNumUsers);

SINRZF = zeros(numTrials,maxNumUsers);
perUserCapacityZF = zeros(numTrials,maxNumUsers);
sumRateCapacityZF = zeros(numTrials,maxNumUsers);

SINRMMSE = zeros(numTrials,maxNumUsers);
perUserCapacityMMSE = zeros(numTrials,maxNumUsers);
sumRateCapacityMMSE = zeros(numTrials,maxNumUsers);

%% Monte Carlo trials 

for monteCarloTrial = 1:numTrials
    if(mod(monteCarloTrial,50)==0)
        fprintf('Trial: %d/%d \n',monteCarloTrial,numTrials)
    end
    
    
    %We begin each trial by specifying an aerial access point with the
    %specified antenna system. There are initially no users associated with
    %the access point.
    ap = accessPoint(f, numAPAntY, ...
        numAPAntZ,lambda/2,patchPatternSq,patchUniqueString);
    
    
    %For each user, this loop generates a position according to a uniform
    %distribution over the spherical sector specified above. It also
    %generates a uniformly random tilt of the array plane, and a uniform
    %rotation of the array system within the array plane itself.
    %It then instantiates an uplinkUser object and associates it with the
    %access point.
    
    for idx = 1:maxNumUsers

        % generate random position with positive x.    thetaPrime = acos(1-rand*(1-cos(maxConeAngle))); %inverse transform sampling.
        thetaPrime = acos(1-rand*(1-cos(coneMaxAngle))); %inverse transform sampling.
        rPrime = nthroot((maxRadius^3-minRadius^3)*rand+minRadius^3,3);
        phiPrime = 2*pi*rand;
    
        %we redefined coordinates  for this section only just to make
        %generating the points easier.
        x = rPrime*cos(thetaPrime);
        y = rPrime*sin(thetaPrime)*cos(phiPrime);
        z = rPrime*sin(thetaPrime)*sin(phiPrime);
        displacementVector = [x;y;z];    
      
        tiltAngle = pi-rand*maxTilt; 
        rotationAngle = rand*2*pi;
        
        xhat = [cos(tiltAngle);sin(tiltAngle)*cos(rotationAngle);sin(tiltAngle)*sin(rotationAngle)];
        yhat = [0;1;0] - (xhat'*[0;1;0])*xhat;
        yhat = yhat/sqrt(yhat'*yhat);
        zhat = cross(xhat,yhat);
        cartesianToLocal = [xhat, yhat, zhat];

        inPlaneRotation = rand*2*pi;
        R = [1 0 0; 0 cos(inPlaneRotation) -sin(inPlaneRotation); 0 sin(inPlaneRotation) cos(inPlaneRotation) ];
        cartesianToLocal = cartesianToLocal*R;
        
        %set up uplink user and associate with the access point
        uplinkUser(ap, displacementVector, numUEAntY, numUEAntZ, lambda/2, cartesianToLocal, patchPatternSq,patchUniqueString);

    end

    
    %This script generates the channel according to the channel model in
    %the paper. The analog beamforming is built into the uplinkUser and
    %accessPoint objects. Right now, only the heuristic beamforming
    %strategy from our paper is available. This will hopefully be updated
    %in the future, or feel free to mess with the code yourself. 
    
    H = zeros(maxNumUsers,maxNumUsers);
    
    %this outer loop iterates through AP subarrays
    for steeredAt = 1:1:maxNumUsers
        
        ap.steerAtUEidx(steeredAt);
        
        %this inner loop iterates through UEs
        for lossIdx = 1:1:maxNumUsers
            
            uex= ap.ueList(lossIdx);
            FSPL = fspldB( lambda, uex.range ); %this includes path loss
                                                %as well as the "shrinking
                                                %aperature" term from the
                                                %Friis formula.         
                                                
            %Rain specific attenuation calculation depends on the geometry 
            %and is in general different for different users.
            rainSA = ...
                rainAttenuation( f, ...
                rainRate,uex.slantAngle, polarizationAngle );
            
            %total atmospheric loss
            atmosphericLoss = (fogSA + rainSA  + atmosphericSA)*uex.range;  
            
            %subarray gains
            gains = 10*log10(uex.steeredGain) +...
                10*log10(ap.rxGains(lossIdx));
            %this is just the Friis formula
            signalLoss =  gains + FSPL + atmosphericLoss;   
            
            %applies our digital channel model, assumes each path from user
            %i to subarray j has a phase that is iid uniformly random.
            H(steeredAt,lossIdx) = ...
                sqrt(10^(signalLoss/10))*...
                exp(1i*(2*pi/lambda)*uex.range)*exp(1i*2*pi*rand);
            
        end
        
    end
    
    % zero forcing calculation
    for numUsers = 1:maxNumUsers
       Hsub = H(1:numUsers,1:numUsers);
       diagonalComponents = diag(inv(Hsub'*Hsub));
       noisePowerZF = diagonalComponents*rxNoisePower;
       SINRZF(monteCarloTrial,numUsers) = mean(txPowerUE./noisePowerZF);    
       sumRateCapacityZF(monteCarloTrial,numUsers) = bandwidth*sum(log2(1+txPowerUE./noisePowerZF));
       perUserCapacityZF(monteCarloTrial,numUsers) = sumRateCapacityZF(monteCarloTrial,numUsers)/numUsers;     
    end
    
    %This block computes capacities for various digital beamforming
    %strategies.
    
    % MMSE calculation
    for numUsers = 1:maxNumUsers
       Hsub = H(1:numUsers,1:numUsers); 
       den = diag(inv(eye(numUsers)+txSNR*(Hsub'*Hsub))); %cite http://www.utdallas.edu/~aria/papers/MehanaNosratinia2012a.pdf, this matched your derivation.
       SINRMMSE = 1./den-1;
       sumRateCapacityMMSE(monteCarloTrial,numUsers) = ...
           sum(bandwidth*log2(1+SINRMMSE)); 
       perUserCapacityMMSE(monteCarloTrial,numUsers) = ...
           mean(bandwidth*log2(1+SINRMMSE));
    end
    
    
    % Analog beamforming only calculation ("spatial division multiple ...
    ...acccess")
    
    for numUsers = 1:maxNumUsers
        
        Hsub = H(1:numUsers,1:numUsers);
        
        for spatTrial = 1:numUsers
            powerVec =circshift( Hsub(spatTrial,:).*conj(Hsub(spatTrial,:)), -(spatTrial-1) );
            rxPower= txPowerUE*powerVec(1);
            interferencePower = sum(powerVec(2:numUsers))*txPowerUE;
            SINRAnalogOnly(monteCarloTrial,numUsers) = SINRAnalogOnly(monteCarloTrial,numUsers) + rxPower/((interferencePower+rxNoisePower)*numUsers);
            sumRateCapacityAnalogOnly(monteCarloTrial,numUsers) = sumRateCapacityAnalogOnly(monteCarloTrial,numUsers)+ bandwidth*log2(1+rxPower/(interferencePower+rxNoisePower)); 
        end
        
       perUserCapacityAnalogOnly(monteCarloTrial,numUsers)=sumRateCapacityAnalogOnly(monteCarloTrial,numUsers)/numUsers;
        
    end
    
    
    
end


avgSumCapacityZF= mean(sumRateCapacityZF);
avgPerUserCapacityZF= mean(perUserCapacityZF);
avgSumCapacityMMSE= mean(sumRateCapacityMMSE);
avgPerUserCapacityMMSE= mean(perUserCapacityMMSE);
avgPerUserCapacityD = mean(perUserCapacityAnalogOnly);
avgSumCapacityD = mean(sumRateCapacityAnalogOnly);
save(strcat('results/',saveFileName));
toc

semilogy(avgSumCapacityMMSE,'.-','LineWidth',5,'MarkerFaceColor','k','MarkerEdgeColor','k', 'MarkerSize',32);
hold on
semilogy(avgSumCapacityZF,'.-','LineWidth',5,'MarkerFaceColor','k','MarkerEdgeColor','k', 'MarkerSize',32);
semilogy(avgSumCapacityD,'.-','LineWidth',5,'MarkerFaceColor','k','MarkerEdgeColor','k', 'MarkerSize',32);
xlabel('Number of Users = Number of RF chains');
ylabel('Achievable Rate (bits/second)');
tstring = sprintf('at %.2f GHz',f/1e9);
title({'Achievable rates for several digital beamforming strategies',tstring});
set(gca,'FontSize', 18);
legend('MMSE','ZF','Analog Only (SDMA)')
hold off

