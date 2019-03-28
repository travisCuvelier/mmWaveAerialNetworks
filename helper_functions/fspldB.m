%does not include gains. (factors that take rx gain to effective area by
%reciprocity are included.
function lossdB = fspldB( lambdam, Rkm  )

lossdB = (-20*log10(4*pi*Rkm*1000)+20*log10(lambdam));


end

