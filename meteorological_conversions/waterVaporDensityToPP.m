%This function converts water vapor density (grams per cubic meter) to a
%partial pressure (in hPascals). 

%Inputs:
%densitygpm3 = density of water vapor (grams per cubic meter)
%tempK = temperature (Kelvin)

%Output:
%pp = partial pressure of water vapor (hPascals)


%Reference:  "The radio refractive index: its formula and refractivity data", ITU Std. ITU-R P.453-12, 2016.
%pp = partial pressure of water vapor (hPascals)

function pp = waterVaporDensityToPP( densitygpm3, tempK )
    pp = densitygpm3.*tempK/216.7;
end

