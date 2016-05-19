classdef slp
    %SLP Sturm-Liouville problem : class definition file
    %   A class that implements a data type for Sturm-Liouville problems
    %
    % slpObject = slp(p,q,w,xmin,xmax,a0,b0,a1,b1,jumps)
    
    properties (Hidden)
        liouvilleTransformed
        % data associated to Liouville's transformation: 
        x2r
        r2x
        SLPp
        SLPw
    end
    
    properties (SetAccess = private, GetAccess = public)
       % coefficient functions (function_handles)
        p
        q 
        w  
       % endpoints integration interval
        xmin  
        xmax
       % boundary conditions
        a0 = 1;
        b0 = 0;
        a1 = 1;
        b1 = 0;
       % classification information
        classification
        classificationSLP = [];
        % list of x-values where jumps are
        % located in one of the coefficient functions
        jumps=[]
    end
    
    methods
        function obj = slp(varargin)
            % Construct a slp object
            % The slp constructor method can be called with the following
            % arguments:
            % obj=slp(p,q,w,xmin,xmax)
            % obj=slp(p,q,w,xmin,xmax,a0,b0,a1,b1)
            if isa(varargin{1},'slp') % copy-constructor
                obj.p = varargin{1}.p;  obj.q = varargin{1}.q;  obj.w = varargin{1}.w;
                obj.xmin = varargin{1}.xmin;  obj.xmax = varargin{1}.xmax;
                obj.a0 = varargin{1}.a0; obj.b0 = varargin{1}.b0;
                obj.a1 = varargin{1}.a1; obj.b1 = varargin{1}.b1;
                obj.classification = varargin{1}.classification;
                obj.liouvilleTransformed = varargin{1}.liouvilleTransformed;
            else
                if nargin~=5 && nargin~=10 && nargin~=9
                    error('MATLAB:slp:WrongArgument','slp-constructor: wrong number of input arguments')
                end
                if isa(varargin{1},'function_handle')  
                   obj.p = varargin{1};
                elseif isa(varargin{1},'char')
                   obj.p = eval(['@(x)' varargin{1}]); %inline(varargin{1}); %werkt niet met arrayfun
                else
                   error('MATLAB:slp:WrongArgument','slp-constructor: arg1: You need to pass the coefficient function p as string or function handle')
                end
                if isa(varargin{2},'function_handle')  
                   obj.q = varargin{2};
                elseif isa(varargin{2},'char')
                   obj.q = eval(['@(x)' varargin{2}]);
                else
                   error('MATLAB:slp:WrongArgument','slp-constructor: arg2: You need to pass the coefficient function q as string or function handle')
                end 
                if isa(varargin{3},'function_handle')  
                   obj.w = varargin{3};
                elseif isa(varargin{3},'char')
                   obj.w = eval(['@(x)' varargin{3}]);
                else
                   error('MATLAB:slp:WrongArgument','slp-constructor: arg3: You need to pass the coefficient function w as string or function handle')
                end 
                if isa(varargin{4},'double')
                   obj.xmin = varargin{4};
                else
                   error('MATLAB:slp:WrongArgument','slp-constructor: arg4: You need to pass the endpoints of the integration interval as double') 
                end
                if isa(varargin{5},'double')
                   obj.xmax = varargin{5};
                else
                   error('MATLAB:slp:WrongArgument','slp-constructor: arg5: You need to pass the endpoints of the integration interval as double') 
                end
                if obj.xmax <= obj.xmin
                    error('MATLAB:slp:WrongArgument',...
                        'slp-constructor: arg4-arg5: The endpoint xmax of the integration interval needs to be larger than the endpoint xmin of the integration interval') 
                end
                if nargin>5
                  if isa(varargin{6},'double')
                     obj.a0 = varargin{6};
                  else
                     error('MATLAB:slp:WrongArgument','slp-constructor: arg6: You need to pass the boundary condition coefficients as double') 
                  end
                  if isa(varargin{7},'double')
                     obj.b0 = varargin{7};
                  else
                     error('MATLAB:slp:WrongArgument','slp-constructor: arg7: You need to pass the boundary condition coefficients as double') 
                  end
                  if isa(varargin{8},'double')
                     obj.a1 = varargin{8};
                  else
                     error('MATLAB:slp:WrongArgument','slp-constructor: arg8: You need to pass the boundary condition coefficients as double') 
                  end
                  if isa(varargin{9},'double')
                     obj.b1 = varargin{9};
                  else
                     error('MATLAB:slp:WrongArgument','slp-constructor: arg9: You need to pass the boundary condition coefficients as double') 
                  end
                  if obj.b0==0 && obj.a0==0
                     error('MATLAB:slp:WrongArgument','slp-constructor: arg6-7: invalid boundary condition') 
                  end
                  if obj.b1==0 && obj.a1==0
                     error('MATLAB:slp:WrongArgument','slp-constructor: arg8-9: invalid boundary condition') 
                  end
                  if nargin>9
                      obj.jumps= varargin{10};
                  else
                      obj.jumps = [];
                  end
                end
                [class,obj] = classify(obj);
                obj.a0=class.bcs(1);
                obj.b0=class.bcs(2);
                obj.a1=class.bcs(3);
                obj.b1=class.bcs(4);
                obj.classification=class;
                
                
                obj.liouvilleTransformed=false;
%                pryce24 is een voorbeeld van een singulier slp die beter
%                niet met liouville transfo wordt gedaan
                if isempty(obj.jumps) && ~class.LNF 
                    if all(class.regular) && ~(isinf(obj.w(obj.xmin)/obj.p(obj.xmin)) || isinf(obj.w(obj.xmax)/obj.p(obj.xmax)))
                        %regular problem without jumps:
                        %apply Liouville transformation
                        %leads to more accurate results
                        he=eps;
                        while isinf(obj.w(obj.xmin)/obj.p(obj.xmin)) || isnan(obj.w(obj.xmin)/obj.p(obj.xmin))
                            obj.xmin=obj.xmin+he;
                            he=he*2;
                        end  
                        he=eps;
                        while isinf(obj.w(obj.xmax)/obj.p(obj.xmax)) || isnan(obj.w(obj.xmax)/obj.p(obj.xmax))
                            obj.xmax=obj.xmax-he;
                            he=he*2;
                        end  
                        obj=liouvilletransformation(obj);
                        obj.liouvilleTransformed=true;
                        obj.classificationSLP = class; %original SLP classification
                    elseif  ~isinf(obj.xmin) && ~isinf(obj.xmax) && ~(isinf(obj.p(obj.xmin)) || isinf(obj.p(obj.xmax)) || isinf(obj.w(obj.xmin)) || ...
                        isinf(obj.w(obj.xmax)) || (obj.p(obj.xmin)==0 && obj.w(obj.xmin)~=0)  || (obj.p(obj.xmax)==0 && obj.w(obj.xmax)~=0)) 
                
                        if ~class.regular(1) 
                            obj.xmin=obj.xmin+1e-13;
                        end
                        if ~class.regular(2) 
                            obj.xmax=obj.xmax-1e-13;
                        end
                        obj=liouvilletransformation(obj);
                        obj.classificationSLP = class; %original SLP classification
                        obj.liouvilleTransformed=true;
                    end
                end
             end
        end
        
        function s=classificationInfo(slp)
            if slp.liouvilleTransformed
                classif=slp.classificationSLP;
            else
                classif=slp.classification;
            end
             s='';
             if length(classif.msg)>1
                    n=sprintf('\n');
                    s='Warning: Problems were detected in the endpoint classification algorithm. ';
                    s=[s 'There may be insufficient information to make a reliable classification. '];
                   % s=[s n classif.msg];
                    s=[s n n]; 
             end
             s=[s getEndpointInfo(classif)];
        end
        
        function obj = set.liouvilleTransformed(obj,value)
            obj.liouvilleTransformed=value;
        end
            
        function disp(obj)
            %DISP Display slp-object
            disp(['  slp-object: (p(x)*y''(x))''+q(x)*y(x)=E*w(x)*y(x) :']);
            disp(['     p(x) = ']) 
            disp(obj.p)
            disp(['     q(x) = ']) 
            disp(obj.q)
            disp(['     w(x) = ']) 
            disp(obj.w)
            
            disp(['     xmin = ' num2str(obj.xmin)])
            disp(['     xmax = ' num2str(obj.xmax)])
            disp(' ')
            disp(['     Boundary conditions:'])
            disp(['        a0*y(xmin)+b0*p(xmin)*y''(xmin)=0 :  '])
            disp(['          a0 = ' num2str(obj.a0) ' , b0 = ' num2str(obj.b0)])
            disp(['        a1*y(xmax)+b1*p(xmax)*y''(xmax)=0 :   '])
            disp(['          a1 = ' num2str(obj.a0) ' , b1 = ' num2str(obj.b0)])
            disp(' ');
        end
        
        function bool=liouvilleNormalForm(obj)
            %returns true if the slp-object is in the Schrodinger form
            %(i.e. p(x)=1=w(x))
            bool= obj.classification.LNF;
        end
        
            
        
    end
    
    
    methods (Access = private, Hidden = true )
        function obj = addLiouvilleTransInfo(obj,fcn,fcn2,SLPp,SLPw)
            obj.x2r=fcn;
            obj.r2x=fcn2;
            obj.SLPp=SLPp;
            obj.SLPw=SLPw;
        end
        end
    
end

