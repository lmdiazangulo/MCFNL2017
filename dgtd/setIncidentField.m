function [ e ] = setIncidentField(e, eT, eS, incField)

if (eT > eS)
	e(eT).dE(1) = e(eT).dE(1) - incField;
	e(eT).dH(1) = e(eT).dH(1) - incField / e(eT).Z;
	if (eS > 0)
    	e(eS).dE(2) = e(eS).dE(2) + incField;
    	e(eS).dH(2) = e(eS).dH(2) + incField / e(eS).Z;
	end
else
	e(eT).dE(2) = e(eT).dE(2) - incField;
	e(eT).dH(2) = e(eT).dH(2) - incField / e(eT).Z;
   	e(eS).dE(1) = e(eS).dE(1) + incField;
   	e(eS).dH(1) = e(eS).dH(1) + incField / e(eS).Z;
end

