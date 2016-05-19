function waarde = get(x,n)
%geeft de nde-waarde uit het object X weer: get(X,1) geeft X in a,
%get(X,2) geeft X' in a, get(X,3) geeft X'' in a en get(X,4) geeft X''' in
%a.

if nargin == 1
    error('functie get heeft 2 inputargumenten nodig: get(X,n)')
else
    switch n
        case 1
            waarde = x.tc(1);
        
        case 2
            waarde = x.tc(2);
        
        case 3
            waarde = x.tc(3);
    
        otherwise
            error('De gevraagde waarde is geen dataveld van het object X')
        
    end
end