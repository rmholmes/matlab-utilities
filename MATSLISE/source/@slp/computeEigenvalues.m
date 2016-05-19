function [E,meshData] = computeEigenvalues(slpObject,pmin,pmax,tol,indices,varargin)
%COMPUTE_EIGENVALUES
% [E,mesh_data] = computeEigenvalues(slp,pmin,pmax,tol,indices)
% [E,mesh_data] = computeEigenvalues(slp,pmin,pmax,tol,indices,mesh_data)
% compute eigenvalues of a SLP problem, when the optional input argument (mesh_data) is
% not specified also the mesh is constructed
% 
%  Input arguments:
%                   slpObject
%                   pmin,pmax: determines which eigenvalues need to be
%                              computed (see indices)
%                   tol: input tolerance
%                   indices: boolean value
%                        indices = true: 
%                           pmin = lower index,
%                           pmax = upper index 
%                           for the eigenvalues to be searched for
%                        indices = false: 
%                           pmin = lower limit,
%                           pmax = upper limit
%                           of the interval to be scanned for eigenvalues
%                   mesh_data: if this optional argument is present, the
%                   function computeEigenvalues does not (re)compute the mesh
%                   but uses mesh_data instead
%  Output arguments:
%                   E: structure
%                      E.success:  
%                             = true if the CP=method ran successfully
%                             = false: the CP-method wasn't able to obtain the data
%                              due to an error. E.msg contains the error message.
%                       E.eigenvalues:
%                             a vector with the eigenvalues in the ascending
%                             order. 
%                       E.errors:
%                             a vector with the errors in the eigenvalues 
%                             calculated on the partition.
%                       E.indices:     
%                             a vector which contains the index n of each eigenvalue.
%                       E.status:
%                             a vector of status-flags, E.status(i) furnishes
%                             information on the quality of the eigenvalue
%                             i:
%                           E.status(i) = 0       
%                                   - correct calculation
%                           E.status(i) > 0 means that difficulties were
%                           detected, the results may be inaccurate or wrong
%
%                           E.status(i) = 1       
%                                    - the eigenvalue or its index may be inaccurate. 
%                                    This may suggest that double precision is not sufficient 
%                                    for that case.
%                           E.status(i) = 2
%                                    - the indices of the basic eigenvalue 
%                                    and the reference one are different.
%                           E.status(i) = 3        
%                                    - the index of the eigenvalue is uncertain.
%                           E.status(i) = 4        
%                                    - = E.status(i) = 2 + E.status(i) = 3
%                           E.status(i) = 5        
%                                    - the reference eigenvalue seems not to
%                                     agree with the basic one
%                           E.status(i) = 6        
%                                    - difficulties in the Prufer
%                                    computations: try a lower input
%                                    tolerance
%                           E.status(i) = 7        
%                                    - the eigenvalues are smaller than the requested 
%                                       accuracy: try a lower input tolerance
%
%                   mesh_Data: structure containing information on the
%                   constructed mesh. This structure can be passed as input
%                   argument in a next call to computeEigenvalues in order
%                   to avoid recomputation of the mesh.

if any(slpObject.classification.oscillatory)
    exception = MException('Classification:checkProblem','The problem has an oscillatory endpoint. Sorry, this problem cannot be handled by Matslise.');
    %error('Classification: the problem has an oscillatory endpoint. Sorry, this problem cannot be handled by Matslise.')
    throw(exception);
end

if indices
    if pmin < 0 || pmax < 0
        exception = MException('VerifyInput:indices','The eigenvalue indices must be positive.');
        throw(exception);
    elseif pmin > pmax 
        exception = MException('VerifyInput:indices','The lower eigenvalue index must be smaller than or equal to the upper index.');
        throw(exception);
    elseif floor(pmin) ~= pmin || floor(pmax) ~= pmax
        exception = MException('VerifyInput:indices','The eigenvalue indices must be integer values.');
        throw(exception);
    end
else
    if pmin >= pmax 
        exception = MException('VerifyInput:indices','The lower energy-value must be smaller than the upper energy-value.');
        throw(exception);
    end
end

%use classification information to check if the requested eigenvalues can
%be computed:
if slpObject.liouvilleTransformed
   pmax=check_evrange(slpObject.classificationSLP,pmin,pmax,indices);
else
   pmax=check_evrange(slpObject.classification,pmin,pmax,indices);  
end

tol=max(tol,1e-15);
iflag=0;


if nargin<6
    %mesh construction (an initial mesh is returned when the problem is singular)
    if indices
          meshData = computeMesh(slpObject,tol,pmax);
    else
          meshData = computeMesh(slpObject,tol,25);   
    end
else
    meshData=varargin{1}; %a mesh was given as input argument
end




if meshData.halfRangeReduction && slpObject.xmin==-slpObject.xmax
    %the problem is symmetric and half-range reduction has been applied
    if indices
        %compute even eigenvalues:
        slpObject2=slp(slpObject.p,slpObject.q,slpObject.w,0,slpObject.xmax,0,1,slpObject.a1,slpObject.b1,slpObject.jumps(slpObject.jumps>0));
        if ceil((pmax-1)/2)>=ceil(pmin/2)
            [E,meshData]=computeEigenvalues(slpObject2,ceil(pmin/2),ceil((pmax-1)/2),tol,indices,meshData);
        else
            E.eigenvalues=[]; E.errors=[]; E.status=[]; E.success=true; 
        end
        %compute odd eigenvalues:
        % use now min(slpObject.xmax,meshData.Infb) instead of slpObject.xmax !:
        % avoids recomputing an appropiate truncation point for infinite
        % endpoint.
        slpObject2=slp(slpObject.p,slpObject.q,slpObject.w,0,min(slpObject.xmax,meshData.Infb),1,0,slpObject.a1,slpObject.b1,slpObject.jumps(slpObject.jumps>0));
        if floor((pmax-1)/2)>=floor(pmin/2)
            [E2,meshData]=computeEigenvalues(slpObject2,floor(pmin/2),floor((pmax-1)/2),tol,indices,meshData);
        else
            E2.eigenvalues=[]; E2.errors=[]; E2.status=[]; E2.success=true; 
        end
        E.indices=pmin:pmax;
        [E.eigenvalues,ix]=sort([E.eigenvalues E2.eigenvalues]);
    else
        %half-range problem 1
        slpObject2=slp(slpObject.p,slpObject.q,slpObject.w,0,slpObject.xmax,0,1,slpObject.a1,slpObject.b1,slpObject.jumps(slpObject.jumps>0));
        [E,meshData]=computeEigenvalues(slpObject2,pmin,pmax,tol,indices,meshData);
        %half-range problem 2
        slpObject2=slp(slpObject.p,slpObject.q,slpObject.w,0,min(slpObject.xmax,meshData.Infb),1,0,slpObject.a1,slpObject.b1,slpObject.jumps(slpObject.jumps>0));
        E2=computeEigenvalues(slpObject2,pmin,pmax,tol,indices,meshData);
        %combine:
        E.indices=min(E.indices(1)*2,E2.indices(1)*2+1):max(E.indices(end)*2,E2.indices(end)*2+1);
        [E.eigenvalues,ix]=sort([E.eigenvalues E2.eigenvalues]);
    end  
    if E.success && E2.success
        E.eigenvalues=E.eigenvalues(1:length(E.indices));
        errs=[E.errors E2.errors];
        status=[E.status E2.status];
        E.errors=errs(ix(1:length(E.indices)));
        E.status=status(ix(1:length(E.indices)));
    else
        E.success=false;
    end
    return;
end




% get_intervals returns the energy range [emin,emax] and the index range 
% [nmin,nmax] of the eigenvalues to be obtained. 
[emin,emax,nmin,nmax,err] = getIntervals(slpObject,pmin,pmax,indices,meshData,tol);


if err.fail
   E.eigenvalues=[];
   E.indices=[];
   E.errors=[];
   E.status=[]; 
   E.success = false;
   E.msg = err.msg ;
elseif nmax - nmin < 0 &&  ~(isinf(slpObject.xmin) || isinf(slpObject.xmax))
   E.eigenvalues=[];
   E.indices=[];
   E.errors=[];
   E.status=[]; 
   E.success = false;
   E.msg ='no eigenvalue is found in the given interval' ;
else
   % The effective search of the eigenvalues with index n from nmin to nmax is
   % now done.

    if slpObject.classification.cont_spectrum(1) || slpObject.classification.cont_spectrum(2)
       %there is a continuous spectrum
       if indices
            nmax=pmax;
       end
       c=slpObject.classification.cutoff; %start of continuous spectrum
       % enlarge the truncated interval (at the infinite endpoint)
       % until the index of c computed on
       % this truncated interval is larger than or equal to nmax
       if isinf(slpObject.xmin) || isinf(slpObject.xmax)
           n=locateEigenvalue(slpObject,c,meshData,false);
           while n<nmax
             meshData=enlargeMesh(slpObject,meshData,tol,c,nmax-n);
             n=locateEigenvalue(slpObject,c,meshData,false);
           end
           meshData=enlargeMesh(slpObject,meshData,tol,c,10);
       end
    end

   E=[]; 
   if meshData.trunca || meshData.truncb
       % finite singularity and not a radial schrodinger problem
       % -> the problem needs to be truncated
       % choice of truncation point depends on which eigenvalue needs to be
       % computed
       [E,meshData]=truncateSingularity(slpObject,meshData,emin,emax,nmax,tol);
       if ~(isinf(slpObject.xmin) || isinf(slpObject.xmax)) % no infinite problem
           if nmin==nmax %there needs to be computed only one eigenvalue
              if abs(E.errors(1))>0.1
                 [emin,emax,nmin,nmax,err] = getIntervals(slpObject,pmin,pmax,indices,meshData,tol);
                 E=calculateEigenvalues(slpObject,emin,emax,nmax,nmax,tol,meshData);
                 if E.status(1)>0 || abs(E.errors(1)>0.1)
                     exception = MException('ComputeEigenvalues:Prufer','Difficulties detected in the Prufer process: try to choose a lower input tolerance.');
                     throw(exception);
                 end
              end 
              E.success=true;
              E.status=iflag;
              return;
           end
          [emin,emax,nmin,nmax,~] = getIntervals(slpObject,pmin,pmax,indices,meshData,tol);
       end
   end
   if isinf(slpObject.xmin) || isinf(slpObject.xmax) %infinite problem
       [E,meshData]=truncateInfiniteEndpoint(slpObject,meshData,indices,E,emax,nmax,tol);
       if nmin==nmax && indices
           %only one eigenvalue was requested
           E.success=true;
           E.status=iflag;
           if abs(E.eigenvalues(end))<tol
                E.success=false;
                E.status=7;
                E.msg='The absolute value of the eigenvalue is smaller than the requested accuracy. Try a smaller tolerance. ';
           end
           return;
       end
       [emin,emax,nmin,nmax,~] = getIntervals(slpObject,pmin,pmax,indices,meshData,tol); %als weg dan problemen bij harm osc indices between 1000 en 1010: status>0
   end
   % a regular problem remains:
   if nmax - nmin < 0 
       E.eigenvalues=[];
       E.indices=[];
       E.errors=[];
       E.status=[]; 
       E.success = false;
       E.msg ='no eigenvalue is found in the given interval' ;
       return;
   end
   
   %compute eigenvalues of regular problem:
   [E,meshData]=calculateEigenvalues(slpObject,emin,emax,nmin,nmax,tol,meshData);
   
   for i=1:length(E.status)
       E.status(i)=max(E.status(i),iflag);
   end
   E.success=true;
   if abs(E.eigenvalues(end))<tol
       E.success=false;
       inds=abs(E.eigenvalues)<tol;
       E.status(inds)=7;
       E.msg='The absolute value of the eigenvalue is smaller than the requested accuracy. Try a smaller tolerance. ';
   end
end


function [E,meshData]=truncateInfiniteEndpoint(slpObject,meshData,indices,E,emax,nmax,tol)
%enlarges the mesh until successive eigenvalue approximations remain
%similar
   c=slpObject.classification.cutoff; %start of continuous spectrum 
   if ~indices
      nmax1=locateEigenvalue(slpObject,emax,meshData,false);
      nmax=nmax1+10;
      while round(nmax)~=round(nmax1)
        nmax=nmax1;
        meshData=enlargeMesh(slpObject,meshData,tol);
        nmax1=locateEigenvalue(slpObject,emax,meshData,false);
      end
   end   
   if  meshData.trunca || meshData.truncb
       meshData=enlargeMesh(slpObject,meshData,tol);
       E1=E;
       E=calculateRefinedEigenvalue(slpObject,E,tol,meshData,nmax);
   else 
       [emint,emax,~,nmax,~] = getIntervals(slpObject,nmax,nmax,true,meshData,tol);
       E1=calculateEigenvalues(slpObject,emint,emax,nmax,nmax,tol,meshData);
       meshData=enlargeMesh(slpObject,meshData,tol);
       E=calculateEigenvalues(slpObject,emint,emax,nmax,nmax,tol,meshData);
       E.eigenvalues=min(c,E.eigenvalues);
   end
   it=1;
   while abs(E.eigenvalues-E1.eigenvalues)>tol*max(1,abs(E.eigenvalues)) || (meshData.LNF && abs(E.errors)>tol*max(5,abs(E.eigenvalues)))
       E1=E;
       meshData=enlargeMesh(slpObject,meshData,tol,E1.eigenvalues(1));
       E=calculateRefinedEigenvalue(slpObject,E1,tol,meshData,nmax);
       it=it+1;
       if it>10 && abs(E.eigenvalues-E1.eigenvalues)<tol*max(1,abs(E.eigenvalues))
           break;
       end
       if ~indices && E.eigenvalues>c
           break;
       end
   end

function [E,meshData,iflag]=truncateSingularity(slpObject,meshData,emin,emax,nmax,tol)
%find good truncated value for singular endpoint.
%mesh in meshData is repeatedely refined near the singularity until
%some successive approximations for E remain similar
    iflag=0;
    c=slpObject.classification.cutoff; %start of continuous spectrum
    if slpObject.classification.cont_spectrum(1) || slpObject.classification.cont_spectrum(2)
           %there is a continuous spectrum
           %move truncated endpoint closer to the (finite) singular endpoint:
           %until the index of c remains the same
           meshData=refineSingularMesh(slpObject,meshData,tol);
           n=locateEigenvalue(slpObject,c,meshData,false);
           n1=n+10;
           while round(n)~=round(n1)
               n1=n;
               meshData=refineSingularMesh(slpObject,meshData,tol);
               n=locateEigenvalue(slpObject,c,meshData,false);
           end
    end
    %compute eigenvalue with index nmax:
    E1=calculateEigenvalues(slpObject,emin,emax,nmax,nmax,tol,meshData);
    meshData=refineSingularMesh(slpObject,meshData,tol);
%     if E1.status==5 %something wrong
%            %it may be necessary to recompute the requested energy-range or
%            %index-range over the enlarged interval
%             [emin,emax,nmin,nmax,err] = getIntervals(slpObject,pmin,pmax,indices,meshData,tol);
%             E1=calculateEigenvalues(slpObject,emin,emax,nmax,nmax,tol,meshData);
%             meshData=refineSingularMesh(slpObject,meshData,tol);
%     end
    E=calculateRefinedEigenvalue(slpObject,E1,tol,meshData,nmax);
    while abs(E.eigenvalues-E1.eigenvalues)>tol*max(1,abs(E.eigenvalues)) || (meshData.LNF && abs(E.errors)>tol*max(5,abs(E.eigenvalues))) %|| abs(E.eigenvalues-E1.eigenvalues)>100*tol*max(1,abs(E.eigenvalues))))
               E1=E;
               meshData=refineSingularMesh(slpObject,meshData,tol);
               E=calculateRefinedEigenvalue(slpObject,E1,tol,meshData,nmax);
               if meshData.trunca && meshData.h(1)<1e-13 || meshData.truncb && meshData.h(end)<1e-13
                   iflag=7;
                   break;
               end
               if E.eigenvalues==c && E1.eigenvalues==c
                   break;
               end
    end
   
    

function [emin,emax,nmin,nmax,err]=getIntervals(slp,pmin,pmax,indices,meshData,tol)
% GET_INTERVALS
% returns in which energy-range and which index-range the required eigenvalues are
% situated
% 
if meshData.LNF %problem in Liouville normal form
    vmin=meshData.V0(meshData.imatch);
    if meshData.radial
        vmin=min(slp.q(meshData.r0/2),vmin);
    end
    % we would like accepte= energy value < lowest eigenvalue
    accepte = -2*abs(vmin);
    if accepte<-1e5 && (isinf(slp.xmin) || isinf(slp.xmax))
        accepte=vmin;
    end
else    
    accepte = min(meshData.Q0./meshData.W0);
end
err.fail = false;


if indices
    % the energy range [emin,emax] consistent with the input index range is
    % determined 
    if ~meshData.LNF && meshData.Q0(1)/meshData.W0(1)==accepte
          tmp=meshData.Q0./meshData.W0;
          tmp2=(tmp(2:end)-tmp(1:length(tmp)-1))./meshData.h(2:end);
          tmp3=find(abs(tmp2)<1e4);
          ind=tmp3(1);
          accepte=tmp(ind);
          acceptn = locateEigenvalue(slp,accepte,meshData,false);
          while ((acceptn - pmin -1) > eps) && ind>1
              ind=ind-1;
              accepte = tmp(ind);
              acceptn = locateEigenvalue(slp,accepte,meshData,false);
          end
    else
        acceptn = locateEigenvalue(slp,accepte,meshData,false);
        if ((acceptn - pmin -1) > eps)
              % first eigenvalue may be smaller than minimum of potential
               if abs(accepte)<1
                   accepte=accepte-2;
               else
                   accepte=accepte-abs(accepte)/2;
               end
               acceptn = locateEigenvalue(slp,accepte,meshData,false);
               it=1;
               while ((acceptn - pmin -1) > eps) && it<15
                   if abs(accepte)<1
                       accepte=accepte-2;
                   else
                       accepte=accepte-abs(accepte)/2;
                   end
                   acceptn = locateEigenvalue(slp,accepte,meshData,false);
                   it=it+1; 
               end
               if it>15
                   err.fail = true; 
                   err.msg = sprintf(' The minimal eigenvalue index to be accepted is %d.',floor(acceptn));
                   emin = accepte;
                   emax = min(0,slp.classification.cutoff);
                   nmin = 0;
                   nmax = 0;
                   return;
              end
        end
    end
    e = accepte;
    scal = min([abs(e)/(pmax+1-pmin),2000,abs(e)/10]);
    emin = convert_n2e(slp,e,scal,meshData,pmin);
    nmin = pmin;
    c=slp.classification;
    if ~isempty(c) && (c.cont_spectrum(1) || c.cont_spectrum(2))
        scal = min(abs(c.cutoff-emin)/(pmax+1-pmin),2000);
    end
    emax = convert_n2e(slp,emin,scal,meshData,pmax+1);
    if emax == 0
        emax = -eps;
    end
    nmax = pmax;
else
    % produce the index of the lowest and of the highest eigenvalue in the
    % given energy interval (= nmin, nmax resp. )
    emin = pmin;
    emax = pmax;
    nmin0 =locateEigenvalue(slp,emin,meshData,false);
    nmin = floor(nmin0);
    if(abs(nmin-nmin0)<tol)
     warning('getIntervals:warning','it may be safer to decrease tol')    
    end
    nmax = locateEigenvalue(slp,emax,meshData,false);
    nmax = floor(nmax - 1);
end



function [e] = convert_n2e(slp,e,he,part,n)
%returns an e-value corresponding to an input index n
if he == 0
    he = 1;
end
c=slp.classification;
cs=any(slp.classification.cont_spectrum);
rn = locateEigenvalue(slp,e,part,false);
p1 = n+1-rn-eps;
while p1 <= 0
    e = e - he;
    rn = locateEigenvalue(slp,e,part,false);	         	
    p1 = n+1-rn-eps;
end
it=0;
while true && it<40
  it=it+1;
  rn = locateEigenvalue(slp,e,part,false);      	
  p1 = n+1-rn-eps;
  if p1 > 0 
	  e0 = e;
	  p0 = p1;
	  e = e0 + he;
      if cs && e<c.cutoff
         he=(c.cutoff-e)/5;  
      else
     	 he = 2*he;
      end
  else
	  e1 = e;
	  break;
  end
end
if it==40
  e1=e;   
end
et=e0+(e1-e0)*0.95;
rn =  locateEigenvalue(slp,et,part,false);
pt = n+1-rn-eps;
if pt>0
    e0=et;
    p0=pt;
else
    e1=et;
end
while true
   if p0 + 1e-14 > 1 
	  et = (e0 + e1)/2;
	  scalev = max(1,abs(et));	
	  if abs(e1-e0) < (eps*scalev)
          e = e0;
          break;
      else
	      rn =  locateEigenvalue(slp,et,part,false);
	      pt = n+1-rn-eps;
	      if p0*pt < 0 
	        e1 = et;
	      else
	        e0 = et;
	        p0 = pt;
	      end
	  end
   else
	  e = e0;
	  break;
   end
end 





%---------------------------------------------
function [emax]=check_evrange(class,emin,emax,indices)
%check if the requested eigenvalue range is possible
 if class.category>=3 && class.category<=7 && isempty(strfind(class.msg,'uncertain'))
     exception = MException('VerifyInput:indices','The eigenvalue indices must be positive.');
     throw(exception);
     error('Classification: the eigenvalues of this problem can not be labeled by the number of zeros in (a, b).')
 end
if indices
  if class.category==2 || class.category==10
   lastev=class.lastev;
   if lastev>=0 && emax>=lastev
        if emin < lastev
            emax=lastev-1;
            if lastev == 0
               warning('There appear to be no eigenvalues below the start of the continuous spectrum.');
            else
               warning(['There are only ' num2str(lastev) ' eigenvalues in the discrete spectrum. Only the eigenvalues with index < ' num2str(lastev) ' are calculated.']);
            end
        else
            if lastev == 0
              exception = MException('VerifyInput:classification','There appear to be no eigenvalues below the start of the continuous spectrum..');
               throw(exception);
            else
               s ='There are only ' + num2str(lastev) + ' eigenvalues in the discrete spectrum.';
               exception = MException('VerifyInput:classification', s);
               throw(exception);
            end
        end
   end
  end
else
    cutoff=class.cutoff;
    if emax > cutoff
        if emin < cutoff
            emax = cutoff;
        else
           exception = MException('VerifyInput:classification', 'The requested energy range is situated in the continuous spectrum.');
           throw(exception);
        end
    end
    if class.category==2 || class.category==10
        lastev=class.lastev;
        if lastev == 0
            exception = MException('VerifyInput:classification', 'There appear to be no eigenvalues below the start of the continuous spectrum.');
            throw(exception);
        end
    end
end




function [E,mesh]=calculateEigenvalues(slpObject,emin,emax,nmin,nmax,tol,mesh)
% search of the eigenvalues with index n from nmin to nmax 
e_low = emin;
n_low = locateEigenvalue(slpObject,emin,mesh,false); %niet noodzakelijk gelijk aan nmin
number_eigv = 0; 
h_en = (emax-emin)/(5*(nmax-nmin+1));
h_mem = h_en;

%preallocation
eigv = zeros(1,nmax-nmin+1);
index = zeros(1,nmax-nmin+1);
indexr = zeros(1,nmax-nmin+1);
err = zeros(1,nmax-nmin+1);
status = zeros(1,nmax-nmin+1);

refined=false;


for n=nmin:nmax  %for each eigenvalue index
	number_eigv = number_eigv + 1;
    status(number_eigv)=0;
    h_en = abs(min(h_mem,h_en));
    % an energy interval [e_low,e_up] is searched such that any e in
    % this interval represents a good initial guess for the Newton
    % iteration.
    while true
      %h_en = abs(min(h_mem,h_en));  
	  e_up = e_low + h_en;
      n_up = locateEigenvalue(slpObject,e_up,mesh,false); %return index of e_up
      if h_en > 0 && n_up < n_low && abs(n_up-n_low)>tol
          % something wrong
          if ~refined
              mesh=refineMesh(slpObject,mesh);
              refined=true;
          end
          e_up = e_low + 2*h_en;
          n_up = locateEigenvalue(slpObject,e_up,mesh,false);
      end
      rz2 = n_up-n-1;
      if rz2 < 0 
	     e_low = e_up;
         n_low = n_up;
         h_en = h_en*2;
      else
	     break; % stop while	
      end
    end %while
     % calculate the eigenvalue corresponding to the basic
     % method:
     e_lows=e_low;
     n_lows=n_low;
     [ev,e_low,n_low,e_up,n_up,evpruf0,fail] = ...
        calculateEigenvalueBP(slpObject,e_low,n_low,e_up,n_up,n,mesh,tol);
     if fail && (e_low >= e_up || round(n_low -1) ~= n) 
         eigv(number_eigv) = ev;
         index(number_eigv) = n;
         indexr(number_eigv) = -1;
         continue;
     end
      % check whether ev is actually inside [e_low,e_up]:
      if (ev-e_low)*(ev-e_up) > 0 || round(evpruf0 -1) ~= n || abs(evpruf0-1-n)>1e-3
          if n_lows>n+1 && n_lows>n_low
             e_lows=e_low; n_lows=n_low;
          end
          [ev,e_low,n_low,e_up,n_up,evpruf0]=...
              calculate_eigenvalue_lin_interpolation(slpObject,e_lows,n_lows,e_up,n_up,n,mesh,false);
      end
	 % compute reference eigenvalue: the basic one is used as
     % initial guess
      [evx,evpruf]= calculateEigenvalueRP(slpObject,ev,mesh,tol/10,true);
      if round(evpruf -1) ~= n || abs(evpruf-1-n)>1e-3
          evxs=evx;
          [evx,e_low,n1,e_up,n_up,evpruf]=...
              calculate_eigenvalue_lin_interpolation(slpObject,e_low,n_low,e_up,n_up,n,mesh,true);
          if abs((evx-evxs)/evx)>tol*100
              status(number_eigv)=5; %something wrong
              %als SLP en geen jumps dan Liouville transf toepassen
              if ~mesh.LNF && isempty(slpObject.jumps)
                  schrod=slp(slpObject.p,slpObject.q,slpObject.w,mesh.h(1),slpObject.xmin+sum(mesh.h),slpObject.a0,slpObject.b0,slpObject.a1,slpObject.b1);
                  if exist('meshSchrod','var')
                    [E,meshSchrod] = computeEigenvalues(schrod,n,n,tol,true,meshSchrod);
                  else
                    [E,meshSchrod] = computeEigenvalues(schrod,n,n,tol,true);
                  end
                  ev=E.eigenvalues(1); evpruf0=E.indices(1)+1; evpruf=evpruf0; evx=ev+E.errors(1); status(number_eigv)=E.status(1);
              end
          end
      end
      eigv(number_eigv) = ev;
      index(number_eigv) = abs(round(evpruf0 - 1));
      indexr(number_eigv) = abs(round(evpruf - 1));	
	  err(number_eigv) = (ev-evx); %error estimation
      if abs(err(number_eigv))<eps
          err(number_eigv) = eps;
      end
      h_en = max(min(e_up-e_low,h_en),tol);
	  e_low = e_up;
	  n_low = n_up;
end %for

%   The results of the previous calculation are now scanned once again
%   and corrected if necessary. 
for i=1:number_eigv
    if index(i) ~= indexr(i) || index(i)  ~= nmin+i-1
	  n = nmin+i-1;
      index(i) = n;
	  e1 = eigv(i);
	  h_en = 2*tol; 
	  while true
	    e1 = e1 - h_en;
        pruf1 = locateEigenvalue(slpObject,e1,mesh,false);       	
	    h_en = 2*h_en;
        if (pruf1-n-1) <= 0
            break;
        end
	  end %while
      while true
	    e2 = e1 + h_en/2;
        pruf2 = locateEigenvalue(slpObject,e2,mesh,false);
        if (pruf2-n-1) >= 0 
          [ev0,t1,t2,t3,t4,pruf0]= calculate_eigenvalue_lin_interpolation(slpObject,e1,pruf1,e2,pruf2,n,mesh,false);
          [ev,t1,pruf1,t2,t3,pruf] = calculateEigenvalueBP(slpObject,e1,pruf1,e2,pruf2,n,mesh,tol);
          if abs(ev0- ev) > tol
              if abs(round(pruf-1)-n) > abs(round(pruf0)-1)
                  ev = ev0;
                  pruf = pruf0;
              end
          end
          if indexr(i) ~= -1
              indexr(i) = round(pruf-1);
          end
          [evx,evpruf]= calculateEigenvalueRP(slpObject,ev,mesh,tol,true);
          err(i) = (ev-evx);
          if abs(err(i))<eps
            err(i) = eps;
          end
          eigv(i) = ev;
          break;    
        else
	      e1 = e2;
	      h_en = 2*h_en;
	      pruf1 = pruf2;
        end
      end %while
    end
end %for

for i=1:number_eigv
    if indexr(i) == -1
        status(i)=1;
        continue;
    elseif indexr(i) ~= index(i)
        status(i)=2;
    elseif index(i)~= nmin+i-1
        status(i)=3;
    end
    if status(i) == 2 && index(i) ~= nmin+i-1 
        status(i) = 4;
    end
    if abs(err(i))>1
        status(i) = 6;
    end  
end %for
E.eigenvalues=eigv;
E.errors=err;
E.indices=index;
E.status=status;



function E=calculateRefinedEigenvalue(slp,E,tol,mesh,n)    
%computes new eigenvalue approximation when a previous approximation is
%available
cutoff=slp.classification.cutoff;
if cutoff<1e100
     pruf = locateEigenvalue(slp,cutoff,mesh,false);
     if pruf<n
         E.eigenvalues=cutoff;
         E.errors=1;
         E.indices=n;
        return
    end
end

% compute new eigenvalue
[ev,evpruf0]= calculateEigenvalueRP(slp,E.eigenvalues(1),mesh,tol,false);
eigv = ev;
index = abs(round(evpruf0 - 1));


%   The results of the previous calculation are now scanned once again
%   and corrected if necessary. 
if index  ~= n || abs(evpruf0-1-n)>1e-3
      index = n;
	  e1 = eigv;
	  h_en = tol; 
	  while true
	    e1 = e1 - h_en;
        pruf1 = locateEigenvalue(slp,e1,mesh,false);
	    h_en = 4*h_en;
        if (pruf1-n-1) <= 0
            break;
        end
	  end %while
      while true
	    e2 = e1 + h_en/4;
        pruf2 = locateEigenvalue(slp,e2,mesh,false);
        if (pruf2-n-1) >= 0 
          [ev,~,~,~,~,pruf0]= calculate_eigenvalue_lin_interpolation(slp,e1,pruf1,e2,pruf2,n,mesh,false);
          eigv = ev;
          break;    
        else
	      e1 = e2;
	      h_en = 2*h_en;
	      pruf1 = pruf2;
        end
      end %while
end
% compute reference eigenvalue: the basic one is used as
% initial guess
[evx,~]= calculateEigenvalueRP(slp,ev,mesh,tol,true);
err = (ev-evx); %error estimation
if abs(err)<eps
   err=eps;
end
E.eigenvalues=eigv;
E.errors=err;
E.indices=index;




%--------------------------------------------------------------------------
function [e,e1,n1,e2,n2,pruf,small]= calculateEigenvalueBP(slp,e1,n1,e2,n2,n,part,tol0)
%   This calculates the eigenvalue (e_n) by the Newton iteration procedure.
%     e1 / e2     - in input: lower and upper limits for the eigenvalue e_n 
%                 - in output: sharper lower and upper limits
%     n1 / n2     - these are the pruefer theta/pi values corresponding to e1 and e2.
%     part        - data concerning the partition (mesh) of the equation domain.
%     tol         - the requested tolerance.
%
%      pruf       - the Pruefer theta/pi corresponding to the eigenvalue e.
%      small      - small = false indicates that there was no difficulty 
%                          to localise ev. 
%                   small = true says that the Pruefer phase varies
%                          with e so violently that the original [e1,e2] had
%                          to be diminished up to a width consistent with 
%                          eps.

% init_ev : generates a value ev in [e1,e2] to be used as the initial guess
 % for the Newton procedure 
[e,e1,n1,e2,n2,small] = init_ev(slp,e1,n1,e2,n2,n,part,tol0);
[e,pruf]= calculateEigenvalueRP(slp,e,part,tol0,false);

%--------------------------------------------------------------------------
function [e,pruf]= calculateEigenvalueRP(slp,e,part,tol0,reference)
%computes eigenvalue using a newton iteration process
%input e is a good approximation to start the newton iteration process with
newton=1e300;
%initial conditions:
a0 = slp.a0;b0 = slp.b0;
a1 = slp.a1;b1 = slp.b1;
sq1 = sqrt(a0^2+b0^2);
sq2 = sqrt(a1^2+b1^2);
if part.radial
    [yi,yei]=initstep(slp,e,part);
else
    yi = [-b0 / sq1;a0 / sq1]; %initial values
    yei=[0; 0];
end
yf = [-b1 / sq2;a1 / sq2];
yef=[0; 0];
it = 0;
while true
   %Newton iteration
    it = it+1;
    if it >= 70
        break;
    end
    scalev = max(1,abs(e));
    tol = scalev*tol0;
    %forward:
    [yl,yel] = propagateSolution(e,yi,yei,part,true,reference);
    %backward:
    [yr,yer] = propagateSolution(e,yf,yef,part,false,reference);
    %the data obtained at each side of xmatch are combined:
    yl1=yl(1); yr1=yr(1);
    yl2=yl(2); yr2=yr(2);
    yel1=yel(1);  yer1=yer(1);
    yel2=yel(2);  yer2=yer(2);
    phi = yl1*yr2-yr1*yl2; %miss-match 
    dphi = yel1*yr2+yl1*yer2-yer1*yl2-yr1*yel2;
    aux = yr1*yl1;
    auxp = yr2*yl2;
    daux = yer1*yl1+yr1*yel1;
    dauxp = yer2*yl2+yr2*yel2;
    newton0=newton;
    if abs(aux) >= abs(auxp)
        newton = -phi*aux/(dphi*aux-daux*phi);
    else %(abs(aux) < abs(auxp))
        newton = -phi*auxp/(dphi*auxp-dauxp*phi);
    end
    if abs(dphi)>abs(phi)
        testnewton=-phi/dphi;
        if newton*testnewton<0
            newton=testnewton;
        end
    end
    if isnan(newton)
        newton= abs(e)/100000;
    end
    if abs(newton)>abs(newton0) && abs(newton/newton0)>10
        break;
    end
    e = e + newton;
    if abs(newton)/min(abs(e),1) <= tol
        break;
    end
    if part.radial
        [yi,yei]=initstep(slp,e,part);
    end
end %while  
pruf = locateEigenvalue(slp,e,part,reference);



%--------------------------------------------------------------------------
function [ev,e1,pruf1,e2,pruf2,small]= init_ev(slp,e1,pruf1,e2,pruf2,n,part,tol)
%@CPM\INIT_EV
% [ev,e1,pruf1,e2,pruf2,small]= init_ev(e1,pruf1,e2,pruf2,n,part)
%  For given input e1 and e2 for which it is known that the eigenvalue
%  e_n is between e1 and e2, function init_ev furnishes a value ev to be 
%  used as the initial value in the Newton procedure of
%  calculate_eigenvalue for locating e_n. 
%  ev is calculated by alternative use of the linear interpolation and 
%  of halving.
%
%   INPUT:
%      n                 - the value of the eigenvalue index
%    part                - the set of data corresponding to the 
%                          partition of the domain of the equation.
%     e1,e2              - lower and upper limits for the 
%                          eigenvalue e_n;
%     pruf1,pruf2       -  these are the pruefer 
%                          theta/pi values corresponding to e1 and e2.
%	                       Conditions pruf1-n-1 < 0 and pruf2-n-1 > 0 are 
%                          assumed to be satisfied in the input.
%
%   OUTPUT:
%      ev                - the approximate value of the eigenvalue e_n.
%      small             - small = false indicates that there was no difficulty 
%                          to localise ev. 
%                          small = true says that the Pruefer phase varies
%                          with e so violently that the original [e1,e2] had
%                          to be diminished up to a width consistent with 
%                          eps.
%      e1,e2             - in output, sharper lower and upper limits.
%      pruf1,pruf2       - the pruefer theta/pi values corresponding to e1
%      and e2.
% 

if part.LNF
   vmin = part.V0(part.imatch);
else
   vmin=min(part.Q0./part.W0);
   if part.trunca || part.truncb
     vmin=max(e1,e2);
   end
end
small = false;
phi1 = pruf1-(n+1);
phi2 = pruf2-(n+1);
it = 0;
while true
   ec = (phi2*e1-phi1*e2)/(phi2-phi1); % linear interpolation
   if (abs(phi1) <= 0.1 && abs(phi2) <= 0.1)  %#1
       ev = ec;
	   return;		    	
   else
       scalev = max(1,max(abs(vmin-ec),abs(ec)));
       if scalev*eps  < (e2-e1) %#2
            it = it + 1;
            if rem(it,2) == 0  %even 
                ec = (e1+e2)/2;
            end
            prufc = locateEigenvalue(slp,ec,part,false);
            if (abs(phi1)+abs(phi2) < 0.25 && (ec-e1)*(ec-e2) < 0) 	%#3
               [~,ind]=min([prufc-(n+1),phi1,phi2]);
               tmp=[ec,e1,e2];
               ev=tmp(ind);
	           %ev = ec;
	           return;
	        else
	           phic = prufc-(n+1);
               if phic < (phi1-max((n+1)*eps,tol)) || phic > (phi2+max((n+1)*eps,tol))
         	        %warning(sprintf('INIT_EV: Inaccurate computation of the Pruefer phase is detected at e = %s ',num2str(ec)))
               end
               if phi1*phic<0 && abs(phi1)<1e-12 && abs(phic)<1e-12 
                  %uncertain, therefore try another ec:
                   ec = e1+(e2-e1)*3/4;
                   prufc = locateEigenvalue(slp,ec,part,false);
                   phic = prufc-(n+1);
                   phic=phic-1e-14;
               end
               if phic>0 %phi1*phic <= 0 problemen bij harm osc indices between 1000 and 1010
	              e2 = ec;
	              phi2 = phic;
	              pruf2 = prufc;
	           else
		          e1 = ec;
		          phi1 = phic;	
		          pruf1 = prufc;
               end
             end	%#3
             if abs(e1 - e2)/min(max(abs(e1),abs(e2)),1) < tol || (abs(e1-e2)<1e-3*min(1,abs(ec)) && abs(phi1)<1 && abs(phi2)<1) || abs(e1)<1e-14 && abs(e2)<1e-14
                 ev = ec;
                 return;
             end
	    else %#2
	         small = true;		
	         ev = ec;
	         return	      
        end%#2
    end %#1
 end %while


 %-------------------------------------------------------------------------
 function [e,e_low,n_low,e_up,n_up,pruf]=calculate_eigenvalue_lin_interpolation(slp,e_low,n_low,e_up,n_up,n,part,reference)
%   For given index n and lower and upper limits e_low and e_up, this 
%   function calculates the eigenvalue e_n (argument e) by using 
%   the Pruefer phase function. Since this function uses linear 
%   interpolation, it is more robust but slower than calculate_eigenvalue(). 
%   For these reasons calculate_eigenvalue_lin_interpolation() is called only 
%   when there are some hints that the results produced by 
%   calculate_eigenvalue() may be inaccurate. 
if part.LNF
    vmin=part.V0(part.imatch);
else
    [vmin,~]=max(((e_low+e_up)/2*part.W0-part.Q0).*part.P0);
end
tol0=eps*5;
p_low = n_low-(n+1);
p_up = n_up-(n+1);
it=0;	
while true
    it = it + 1;
	e = (p_up*e_low-p_low*e_up)/(p_up-p_low); %linear interpolation
    if rem(it,2) == 1
        e = (e_low+e_up)/2;
    elseif rem(it,11) == 1 
        et = ((p_up+1)*e_low-(p_low+1)*e_up)/(p_up-p_low);
        if e_low<et && et<e_up
            e=et;
        end
    end
    pruf = locateEigenvalue(slp,e,part,reference);
	pc = pruf-(n+1);
    if it>500 && rem(it,2) == 1 %added to avoid problems when 2 eigenvalues are very close
        return;
    end
%   if pc < p_low-(n+1)*eps  || pc > p_up+(n+1)*eps
 %      warning(['calculate_eigenvalue_lin_interpolation: Inaccurate computation of the Pruefer phase is detected at e = ' num2str(e) ])
 %   end
    if abs(pc)<eps*max(abs(e),1)
        return;
    end
	scalev = max(1,max(abs(vmin-e),abs(e)));
	tol = scalev*tol0;
    if abs(e_up-e_low) > tol %#1
        if p_low*pc < 0 %#2
         p_up = pc;
	     e_up = e;
         n_up = pruf;
	    else %#2
          if p_low == 0
            e = e_low;
            pruf = n_low;
            return;
          end
          if abs(pc) <eps
            return;
          end
	      e_low = e;
          n_low = pruf;
	      p_low = pc;
        end %#2
	else %#1
	  break;
    end %#1
end %while


