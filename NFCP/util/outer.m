function zmat  = outer( vec, varargin ) 
    %{ 
    OUTER 
    
    A function to apply arbitrary function to "cross-product" of inputs
    Concept based on the  R-language 'outer' function.
    copyright (C) 2012 by Carl Witthoft. Code may be be freely distributed
    and incorporated into other functions, so long as it
    is not part of any product which is sold.

    USAGE: outer(vector1 [,vector2]  [, function] )

    The optional 'function' can be either a string, e.g., 'atan2,' or a
    function handle.
    The function must accept two arguments and produce a vector of
    results the length of the first input. There are a few cases, e.g.
    f(x,y) = numel([ x y ]) which may succeed but these are not
    guaranteed.
    Examples: 
        
        outer(([1 2 3 4],5,'plus') 
        outer({ 'a' 'b' 'c' 'd'} , 'c','strcmp')
        foo = @atan2; outer([1:5],foo)
    
    The default function is "*" , i.e. outer product of two numeric vectors.
    If 'vec2' is not supplied, "function(vec,vec)" is assumed.
    Inputs can be vectors of numbers, vectors of characters, or cell vectors
    of characters/char strings. 
    Note that "*" is illegal for characters, so if either vec or
    vec2 is non-numeric, a function must be supplied.
    %}

    strop = '*'; %default value of function
    vec2 = vec; % will overwrite later if vec2 is provided
    if nargin==2
        % if "vec" not number then it's string, so next arg MUST be function name 
         if ~isnumeric(vec)
              strop = (varargin{1});
              if (~iscell(vec))
                   bb = textscan(num2str(vec),'%s');
                  vec=bb{1} ;
                  vec2=vec;
              end;
          else % vec is number; next arg is either number or MUST be function name
              if ~isnumeric(varargin{1}) 
                  strop = (varargin{1});
              else 
                 vec2 = (varargin{1}); %both are numeric; use default operator 
              end
          end
    end
    if nargin ==3 %arg3 is always function name
          vec2 = (varargin{1});
          strop = (varargin{2});
        if (~(isnumeric(vec) && (isnumeric(vec2))))
       %not both numeric, so convert both to cell array
            if(~iscell(vec)) 
       %  This will NOT split a char string input (except at spaces)
                bb = textscan(num2str(vec),'%s');
                  vec=bb{1} ;
            end
            if(~iscell(vec2)) 
                bb = textscan(num2str(vec2),'%s');
                  vec2=bb{1} ;
            end
        end
    end
    %  MatLab handles A'*B as a matrix op. All else must be longhand.
      if (strcmp(strop , '*' ) ) 
        zmat = vec(:)*vec2(:)';
      else 
          %check for handleism
          if (~isa(strop,'function_handle')) strop  = str2func(strop);end;
      % two cases:  vec, vec2 are either numeric or cells
          if isnumeric(vec) % by design vec2 isnum as well
      %check whether the  result is a string
            if isnumeric(strop(vec(1),vec2))          
                zmat = zeros(length(vec), length(vec2));        
                for ii = 1:length(vec)   
                 zmat(ii,:) = strop(vec(ii),vec2);  
                end
            else
                zmat = cell(length(vec),length(vec2));
            % if output is not cell, have to make it so
                if iscell(strop(vec(1),vec2))
                    for ii = 1:length(vec)
                        zmat(ii,:) = strop(vec(ii),vec2);              
                    end
                else
                    for ii = 1: length(vec)
                       bb = textscan(num2str(strop(vec(ii),vec2)),'%s');
                       zmat(ii,:) = bb{1};
                    end
                end
            end
          else
     % Do the same 'switching' for nonnumeric vectors
            zmat = cell(length(vec),length(vec2));
            if iscell(strop(vec{1},vec2))
                for ii = 1:length(vec)
                    zmat(ii,:) = strop(vec{ii},vec2);              
                end
            else
                for ii = 1: length(vec)
                    zmat(ii,:) = num2cell(strop(vec{ii},vec2));
                end
            end
          end 
       end 
    end

