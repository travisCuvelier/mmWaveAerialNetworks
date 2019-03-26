%This function converts relative humidity (%) to water vapor partial
%pressure (hPascals). This function is parameterized by the total
%atmospheric pressure.

%Note that:
%Total Atmospheric Pressure = (Partial Pressure of Dry Air) + (Partial
%Pressure of Water Vapor)

%Inputs:
%RH = relative humidity (%)
%totalAtmosphericPressurehPa = total atmospheric pressure (hPascals)
%tempK = temperature (Kelvin)

%Outputs:
%RH = relative humidity (%)

%Reference:  "The radio refractive index: its formula and refractivity data", ITU Std. ITU-R P.453-12, 2016.

function waterVaporPressurehPa = relativeHumidityToVaporPressure(RH, totalAtmosphericPressurehPa,tempK)

    tempC = tempK-273.15;

    if(tempC <-40 || tempC >50)
        error('This approximation isn''t valid in this temperature range');
    end


    EF = 1+10^(-4)*(7.2+totalAtmosphericPressurehPa.*(.0320+5.9*10^(-6)*tempC.^2));
    a = 6.1121;
    b = 18.678;
    c = 257.14;
    d = 234.5;
    
    waterVaporPressurehPa = RH.*EF*a*exp((b-tempC/d).*tempC./(tempC+c))/100;


end