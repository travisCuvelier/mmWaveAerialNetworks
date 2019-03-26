%This function computes (an estimate) of the attenuation per kilometer 
%from absorbtion by suspended liquid water (clouds, fog).
%Valid at temperature from 265.1500 K to 298.1500 K (-8 to 25 centigrade)
%and from 30 - 100 GHz [1, 2].

%At frequencies above 100 GHz, attenuation from water vapor dominates, 
%and at lower frequencies the attenuation is not appreciable [1]. 
%The likelihood of suspended liquid water at temperatures outside the 
%range is also small [1].

%References: 
%[1] E.A. Altschuler, "A Simple Expression for Estimating Attenuation by 
%    Fog at Millimeter Wavelengths," IEEE Trans. Antennas and Propagation, 
%    Vol. AP-32, No. 7
%
%[2] E.A. Altschuler, ?Addendum to ?A simple expression for estimating 
%    attenuation by fog at millimeter wavelengths?,? IEEE Trans. Antennas 
%    Propagation, vol. 34, no. 8, pp. 1067?1067, August 1986.

%Inputs:
%f = frequency (Hz)
%temperatureK = temperature (K)
%densitygperm3 = density of liquid water (grams per cubic millimeter)

%Outputs:
%specificAttenuationdBperKm = specific attenuation due to liquid water
%absorbtion (dB per kilometer) 

function specificAttenuationdBperKm = fogAttenuation( f, temperatureK, densitygperm3)

temperatureC = temperatureK-273.15;

c = 3e8;

lambda = c./f;

lambdaMM = lambda*10^(3);

specificAttenuationdBperKm = -(-1.347+.0372*lambdaMM+18.0./lambdaMM-.022*temperatureC)*densitygperm3;


end

