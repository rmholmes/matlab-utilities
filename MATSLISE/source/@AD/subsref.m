function c = subsref(a,s)
% Subscripted reference for AD objects.
switch s.type
    case '()'     
      ind=s.subs{:};
      switch ind
      case 1, c = a.tc(1);
      case 2, c = a.tc(2);
      case 3, c = a.tc(3);    
      otherwise 
        error('invalid index reference for AD')
      end
end
      