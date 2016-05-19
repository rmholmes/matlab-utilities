
function srflux = srflux_function(cff1,cff2,srfl,Hangle,lon)

if (abs(cff1)>abs(cff2))
    if (cff1*cff2>0)
        cff = cff1*2*pi;
        srflux = max([0 srfl/cff*(cff1+cff2*cos(Hangle-lon*pi/ ...
                                                180))]);
    else
        srflux = 0;
    end
else
    cff = (cff1*acos(-cff1/cff2)+sqrt(cff2*cff2-cff1*cff1))/pi;
    srflux = max([0 srfl/cff*(cff1+cff2*cos(Hangle-lon*pi/ ...
                                                180))]);
end
