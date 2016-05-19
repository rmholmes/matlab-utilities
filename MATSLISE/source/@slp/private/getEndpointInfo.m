function [s]=getEndpointInfo(c)
%returns classification information
n=sprintf('\n');
switch c.category
    case 1
        s='The spectrum is simple, purely discrete and bounded below.';
        s=[s n 'There is no continuous spectrum.'];
        s= [s n n 'At endpoint A'];
        if c.regular(1)
           s= [s n '   The problem is regular.'];
        else
           s= [s n '   The problem is singular.'];
           s=[s n '   It is nonoscillatory for all eigenvalues.'];
           if c.limit_circle(1)>-1 
             if c.limit_circle(1)
              s= [s n '   It is limit circle.'];
             else
              s= [s n '   It is limit point.']; 
             end
           end
           %friedrichs c.bcs
           if c.bcs(2)>= 0
              s =[s n '   The Friedrichs boundary conditions are ' num2str(c.bcs(1)) 'y(a)+ ' num2str(c.bcs(2)) '(py'')(a)=0'];
           else
              if abs(c.bcs(1)) > 1e-10
                  scal=abs(1/c.bcs(1));
                  c.bcs(1)=scal*c.bcs(1);
                  c.bcs(2)=scal*c.bcs(2);
              end
              s =[s n '   The Friedrichs boundary conditions are ' num2str(c.bcs(1)) 'y(a) - ' num2str(abs(c.bcs(2))) '(py'')(a)=0'];
           end
        end
        s= [s n 'At endpoint B'];
        if c.regular(2)
           s= [s n '   The problem is regular.'];
        else
           s= [s n '   The problem is singular.'];
           s=[s n '   It is nonoscillatory for all eigenvalues.'];
           if c.limit_circle(2)>-1 
             if c.limit_circle(2)
              s= [s n '   It is limit circle.'];
             else
               s= [s n '   It is limit point.']; 
             end
           end
           if c.bcs(4)>=0
              s =[s n '   The Friedrichs boundary conditions are ' num2str(c.bcs(3)) 'y(b)+ ' num2str(c.bcs(4)) '(py'')(b)=0'];
           else
              s =[s n '   The Friedrichs boundary conditions are ' num2str(c.bcs(3)) 'y(b) - ' num2str(abs(c.bcs(4))) '(py'')(b)=0']; 
           end
        end
           %friedrichs c.bcs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 2
        s='The spectrum is simple and bounded below.';
        s=[s n 'There is a continuous spectrum in [' num2str(c.cutoff) ', infinity].'];
        if isinf(c.lastev)
            s=[s n 'There appear to be infinitely many eigenvalues below the start of the continuous spectrum.'];
        elseif c.lastev == 0
            s =[s n 'There appear to be no eigenvalues below the start of the continuous spectrum.'];
        elseif c.lastev == 1
            s=[s n 'There appears to be 1 eigenvalue below the start of the continuous spectrum.'];
        else
            s=[s n 'There appear to be ' num2str(c.lastev) ' eigenvalues below the start of the continuous spectrum.'];
        end
        s= [s n n  'At endpoint A'];
        if c.regular(1)
           s= [s n '   The problem is regular.'];
        else
           s= [s n '   The problem is singular.'];
           if c.cont_spectrum(1)
               s=[s n '   It is nonoscillatory for eigenvalues < ' num2str(c.cev(1)) ' and oscillatory otherwise.'];
           else
               s=[s n '   It is nonoscillatory for all eigenvalues.'];
               if c.bcs(2)>=0 
                  s =[s n '   The Friedrichs boundary conditions are ' num2str(c.bcs(1)) 'y(a)+ ' num2str(c.bcs(2)) '(py'')(a)=0'];
               else
                  s =[s n '   The Friedrichs boundary conditions are ' num2str(c.bcs(1)) 'y(a) - ' num2str(abs(c.bcs(2))) '(py'')(a)=0']; 
               end
           end
           if c.limit_circle(1)>-1 
            if c.limit_circle(1)
              s= [s n '   It is limit circle.'];
            else
              s= [s n '   It is limit point.']; 
            end
           end
        end
        s= [s n 'At endpoint B'];
        if c.regular(2)
           s= [s n '   The problem is regular.'];
        else
           s= [s n '   The problem is singular.'];
           if c.cont_spectrum(2)
               s=[s n '   It is nonoscillatory for eigenvalues < ' num2str(c.cev(2)) ' and oscillatory otherwise.'];
           else
               s=[s n '   It is nonoscillatory for all eigenvalues.'];
               s =[s n '   The Friedrichs boundary conditions are ' num2str(c.bcs(3)) 'y(b)+ ' num2str(c.bcs(4)) '(py'')(b)=0'];
           end
           if c.limit_circle(2)>-1 
            if c.limit_circle(2)
              s= [s n '   It is limit circle.'];
            else
              s= [s n '   It is limit point.']; 
            end
           end
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 3
        s='The spectrum is simple';
        s= [s n 'There are infinitely many negative and infinitely many positive eigenvalues.'];
        s=[s n 'There is no continuous spectrum'];
        s= [s n n  'At endpoint A'];
        if c.regular(1)
           s= [s n '   The problem is regular.'];
        else
           s= [s n '   The problem is singular.'];
           if c.oscillatory(1)
               s=[s n '   It is oscillatory for all eigenvalues.'];
           else
               s=[s n '   It is nonoscillatory for all eigenvalues.'];
           end
           if c.limit_circle(1)>-1 
            if c.limit_circle(1)
              s= [s n '   It is limit circle.'];
            else
              s= [s n '   It is limit point.']; 
            end
           end
        end
        s= [s n 'At endpoint B'];
        if c.regular(2)
           s= [s  n '   The problem is regular.'];
        else
           s= [s n '   The problem is singular.'];
           if c.oscillatory(2)
               s=[s n '   It is oscillatory for all eigenvalues.'];
           else
               s=[s n '   It is nonoscillatory for all eigenvalues.'];
           end
           if c.limit_circle(2)>-1 
            if c.limit_circle(2)
              s= [s n '   It is limit circle.'];
            else
              s= [s n '   It is limit point.']; 
            end
           end
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    case 4
        s='The spectrum is simple.';
        s=[s n 'There is a continuous spectrum consisting of the entire real line.'];
        s= [s n n  'At endpoint A'];
        if c.regular(1)
           s= [s  n '   The problem is regular.'];
        else
           s= [s n '   The problem is singular.'];
           if c.oscillatory(1)
               s=[s n '   It is oscillatory for all eigenvalues.'];
           else
               s=[s n '   It is nonoscillatory for all eigenvalues.'];
           end
           if c.limit_circle(1)>-1 
            if c.limit_circle(1)
              s= [s n '   It is limit circle.'];
            else
              s= [s n '   It is limit point.']; 
            end
           end
        end
        s= [s n 'At endpoint B'];
        if c.regular(2)
           s= [s  n '   The problem is regular.'];
        else
           s= [s n '   The problem is singular.'];
           if c.oscillatory(2)
               s=[s n '   It is oscillatory for all eigenvalues.'];
           else
               s=[s n '   It is nonoscillatory for all eigenvalues.'];
           end
           if c.limit_circle(2)>-1 
            if c.limit_circle(2)
              s= [s n '   It is limit circle.'];
            else
              s= [s n '   It is limit point.']; 
            end
           end
        end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 5
        s='The spectrum is simple.';
        s=[s n 'There are infinitely many negative and infinitely many positive eigenvalues.'];
        s=[s n 'There is no continuous spectrum'];
        s= [s n n  'At endpoint A'];
        s= [s n '   The problem is singular.'];
        s=[s n '   It is oscillatory for all eigenvalues.'];
        s= [s n '   It is limit circle.'];
        s= [s n 'At endpoint B'];
        s= [s n '   The problem is singular.'];
        s=[s n '   It is oscillatory for all eigenvalues.'];
        s= [s n '   It is limit circle.'];
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 6
        s='The nature of the spectrum is unknown. There is likely to be a continuous spectrum.';
        s= [s n n  'At endpoint A'];
        s= [s n '   The problem is singular.'];
        s=[s n '   It is oscillatory for all eigenvalues.'];
        s= [s n '   It is limit point.'];
        s= [s n 'At endpoint B'];
        s= [s n '   The problem is singular.'];
        s=[s n '   It is oscillatory for all eigenvalues.'];
        s= [s n '   It is limit point.'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 7
        s='There is a continuous spectrum consisting of the entire real line';
        s= [s n n  'At endpoint A'];
        s= [s n '   The problem is singular.'];
        s=[s n '   It is oscillatory for all eigenvalues.'];
        if c.limit_circle(1)>-1 
            if c.limit_circle(1)
              s= [s n '   It is limit circle.'];
            else
              s= [s n '   It is limit point.']; 
            end
        end
        s= [s n 'At endpoint B'];
        s= [s n '   The problem is singular.'];
        s=[s n '   It is oscillatory for all eigenvalues.'];
        if c.limit_circle(2)>-1 
            if c.limit_circle(2)
              s= [s n '   It is limit circle.'];
            else
              s= [s n '   It is limit point.']; 
            end
        end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 8
        s='The spectrum is simple.';
        s=[s n 'There is a continuous spectrum in [' num2str(c.cutoff) ', infinity].'];
        s=[s n 'At endpoint A'];
        s=[s n '   The problem is singular.'];
        if c.cont_spectrum(1)
            s=[s n '   It is nonoscillatory for eigenvalues < ' num2str(c.cev(1)) ' and oscillatory otherwise.'];
        else
            s=[s n '   It is oscillatory for all eigenvalues.'];
        end
        if c.limit_circle(1)>-1 
            if c.limit_circle(1)
              s= [s n '   It is limit circle.'];
            else
              s= [s n '   It is limit point.']; 
            end
        end
        s=[s n 'At endpoint B'];
        s=[s n '   The problem is singular.'];
        if c.cont_spectrum(2)
            s=[s n '   It is nonoscillatory for eigenvalues < ' num2str(c.cev(2)) ' and oscillatory otherwise.'];
        else
            s=[s n '   It is oscillatory for all eigenvalues.'];
        end
        if c.limit_circle(2)>-1 
            if c.limit_circle(2)
              s= [s n '   It is limit circle.'];
            else
              s= [s n '   It is limit point.']; 
            end
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 9
        s='This problem may have non-simple spectrum.';
        s=[s n 'At endpoint A'];
        s=[s n '   The problem is singular.'];
        if c.cont_spectrum(1)
            s=[s n '   It is nonoscillatory for eigenvalues < ' num2str(c.cev(1)) ' and oscillatory otherwise.'];
        else
            s=[s n '   It is oscillatory for all eigenvalues.'];
        end
        if c.limit_circle(1)>-1 
            if c.limit_circle(1)
              s= [s n '   It is limit circle.'];
            else
              s= [s n '   It is limit point.']; 
            end
        end
        s=[s n 'At endpoint B'];
        s=[s n '   The problem is singular.'];
        if c.cont_spectrum(2)
            s=[s n '   It is nonoscillatory for eigenvalues < ' num2str(c.cev(2)) ' and oscillatory otherwise.'];
        else
            s=[s n '   It is oscillatory for all eigenvalues.'];
        end
        if c.limit_circle(2)>-1 
            if c.limit_circle(2)
              s= [s n '   It is limit circle.'];
            else
              s= [s n '   It is limit point.']; 
            end
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 10
        s='This problem may have non-simple spectrum';
        s=[s n 'There is a continuous spectrum in [' num2str(c.cutoff) ', infinity].'];
        if c.lastev == -5
            s=[s n 'There appear to be infinitely many eigenvalues below the start of the continuous spectrum.'];
        elseif c.lastev == 0
            s =[s n 'There appear to be no eigenvalues below the start of the continuous spectrum.'];
        elseif c.lastev == 1
            s=[s n 'There appears to be 1 eigenvalue below the start of the continuous spectrum.'];
        else
            s=[s n 'There appear to be ' num2str(c.lastev) ' eigenvalues below the start of the continuous spectrum.'];
        end
        s=[s n 'At endpoint A'];
        s=[s n '   The problem is singular.'];
        s=[s n '   It is nonoscillatory for eigenvalues < ' num2str(c.cev(1)) ' and oscillatory otherwise.'];
        if c.limit_circle(1)>-1 
            if c.limit_circle(1)
              s= [s n '   It is limit circle.'];
            else
              s= [s n '   It is limit point.']; 
            end
        end
        s=[s n 'At endpoint B'];
        s=[s n '   The problem is singular.'];
        s=[s n '   It is nonoscillatory for eigenvalues < ' num2str(c.cev(2)) ' and oscillatory otherwise.'];
        if c.limit_circle(2)>-1 
            if c.limit_circle(2)
              s= [s n '   It is limit circle.'];
            else
              s= [s n '   It is limit point.']; 
            end
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
s=[s n n];