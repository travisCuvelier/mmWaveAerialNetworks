%This function converts relative humidity (%) to water vapor partial
%pressure (hPascals). This function is parameterized by the partial
%pressure of "dry air" (hPascals). 

%Note that:
%Total Atmospheric Pressure = (Partial Pressure of Dry Air) + (Partial
%Pressure of Water Vapor)

%Inputs:
%RH = relative humidity (%)
%dryAtmosphericPressurehPa = partial pressure of dry air (hPascals)
%tempK = temperature (Kelvin)

%Outputs:
%waterVaporPressurehPa = Partial Pressure of Water Vapor (hPascals)

%Reference:  "The radio refractive index: its formula and refractivity data", ITU Std. ITU-R P.453-12, 2016.

function waterVaporPressurehPa = relativeHumidityDryPressureToVaporPressure(RH, dryAtmosphericPressurehPa,tempK)

    tempC = tempK-273.15;


    fprintf('Note that this is only valid for centigrade temperatures between -40 and 50\n\n\n');

    a = 6.1121;
    b = 18.678;
    c = 257.14;
    d = 234.5;
    
    
    waterVaporPressurehPa = (RH.*(1+10^(-4)*(7.2+(dryAtmosphericPressurehPa).*(.0320+5.9*10^(-6)*tempC.^2)))*a.*exp((b-tempC/d).*tempC./(tempC+c))/100)./(1-10^-4*(.0320+5.9*10^(-6)*tempC.^2).*(RH/100).*a.*exp((b-tempC/d).*tempC./(tempC+c)));


end