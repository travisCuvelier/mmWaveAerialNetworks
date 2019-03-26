%This function converts water vapor (partial) pressures (hPascals) to
%relative humidity (%). 

%Inputs:
%waterVaporPressurehPa = partial pressure of water vapor (hPascals)
%totalAtmosphericPressurehPa = total atmospheric pressure (hPascals)
%tempK = temperature (Kelvin)

%Output:
%RH = relative humidity (%)

%Reference:  "The radio refractive index: its formula and refractivity data", ITU Std. ITU-R P.453-12, 2016.

function RH = VaporPressureToRelativeHumidity(waterVaporPressurehPa, totalAtmosphericPressurehPa,tempK)

    tempC = tempK-273.15;

    if(tempC <-40 || tempC >50)
        error('This approximation isn''t valid in this temperature range');
    end


    EF = 1+10^(-4)*(7.2+totalAtmosphericPressurehPa.*(.0320+5.9*10^(-6)*tempC.^2));
    a = 6.1121;
    b = 18.678;
    c = 257.14;
    d = 234.5;
    
    RH = waterVaporPressurehPa*100./(EF*a*exp((b-tempC/d).*tempC./(tempC+c))); 


end