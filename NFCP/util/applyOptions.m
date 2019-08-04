function options = applyOptions(defaults, options, verbose)
    %{
    Convert default arguments into a struct, to create a single convention
    for handing optional arguments and `varargin`
    
    Default arguments can be passed as sequence of key,value pairs as 
    in Matlab convention (like plot), or as a struct, or a cell array 
    with (key,value) pairs. 
    
    If the optional `verbose` flag is enabled, 
    this routine will print the values assigned to each argument
    as well as whether said value is default or user-supplied 
    
    There are many conventions for passing optional arguments in Matlab.
    The one used by Mathworks is to pass optional arguments via `varagin`
    as a 1 x N cell array, as a sequence of key,value entries. 
    Other reasonable conventions would be to pass options as a list of
    key,value pairs in a cell array, or as a struct whose fields are the
    defaults.

    
    Parameters
    ----------
    defaults: `cell`
        Ndefaults x 3 cell array containing a list of key, value, message
        tuples.
    options: `struct` or `cell`
        Optional arguments passed to a function. Most likely the contents
        of `varargin`. Can be a struct or cell array. If a cell array, it
        should either be a 1x1 cell array containing a struct or a cell 
        array, or a cell array that is a Noptions x 2 list of key,value
        pairs or a 1 x (Noptions*2) sequence of said pairs. 
    verbose : `bool`, default `false`
        Whether to print logging and debugging options
        
    Returns
    -------
    options : `struct`
        Returns a struct with `key`→`value` mapping. Any keys not defined
        in the `options` argument will be filled in with default values
        from the `defaults` argument.
    %}

    if nargin<3, verbose=false; end
    
    if nargin<2,
        % The user has omitted the options argument
        % This is ok, this function can also be used to convert
        % a cell array of kay,value or key,value,message tuples
        % into a struct format.
        % Initialize options to an empty struct
        options = {};
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % Interpret the options argument. 
    % Possible cases:
    %   Empty cell array (no options set)
    %   Cell array containing a single struct entry (passed via varargin)
    %   Cell array containing a single cell array (passed via varargin)
    %   Cell array noptions x 2 of keys, values
    %   Flattened cell array with noptions * 2 entries
    %   Struct mapping keys to values
    if iscell(options),
        % optional arguments appear as cell array.
        % This can happen if arguments were passed as varargin,
        % in which case we expect a flattened list of key,value 
        % pairs. It is also possible that the user passed a 
        % prepared `options` struct or cell array, which would
        % appear as a single entry in a {1,1} cell array in a 
        % function expecting `varargin`. 
        noptions = numel(options);
        if noptions==0,
            options = struct;
        elseif noptions==1, 
            options = options{1,1};
        end
    end
    if iscell(options),
        % optional arguments appear in cell array
        % can be a flat sequence, or list of (key,value) pairs
        % convert it to a struct format
        structoptions = struct;
        [nrow,ncol]   = size(options);
        if ncol>2 & ~mod(ncol,2) & nrow==1,
            nrow = floor(ncol/2);
            ncol = 2;
            options = reshape(options,ncol,nrow)';
        end
        if ncol~=2,
            error('Options should be struct or list of key,value pairs');
        end
        for i=1:nrow,
            [k,v]=options{i,1:end};
            structoptions.(k)=v;
        end
        options = structoptions;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    % Interpret the defaults argument
    if iscell(defaults),
        % default arguments provided in cell array format
        % expect nargumentx x 2 or narguments x 3 cell array, where
        % each row is either (argument name, value) or 
        % (argument name, value, message)
        [ndefaults,ncols] = size(defaults);
        if ncols<2 | ncols>3,
            error('Specify defaults as list of key,val or k,v,msg tuples');
        end
        % re-pack (key,value) pairs into a cell array
        for ii=1:ndefaults,
            kv(ii,1:2) = defaults(ii,1:2);
        end
    elseif isstruct(defaults),
        % default argument values provided as struct
        % field names are interpreted as argument names
        % field values are interpreted as argument values
        ndefaults = numel(fieldnames(defaults));
        k = fieldnames(defaults);
        ncols = 2;
        % re-pack (key,value) pairs into a cell array
        for ii=1:ndefaults,
            kv(ii,1:2) = {k{ii}; defaults.(k{ii})};
        end
    else
        error('Defaults should be specified as a struct or cell array.');
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    % Apply the default values to the provided options
    % if provided options are missing keys
    for ii=1:ndefaults,
        % for every option not supplied by the user,
        % set it to the value in defaults
        [k,v] = kv{ii,1:2};
        mode = 'passed ';
        if ~isfield(options,k),
            options.(k)=v;
            mode = 'default';
        end
        if verbose && strcmp(mode,'passed '),
            m = '';
            if ncols==3, m=defaults{ii,3}; end
            valuename = options.(k);
            % Matlab is in serious need of better "toString" functionality
            toString  = @(var) evalc(['disp(var)']);
            if isnumeric(valuename) || islogical(valuename), 
                if isscalar(valuename),
                    valuename = num2str(valuename);
                else,
                    x1 = num2str(valuename(1));
                    x2 = num2str(valuename(2));
                    valuename = sprintf('%s %s …',x1,x2);
                end
            end
            if iscell(valuename),
                valuename = valuename(:);
                valuename = toString(valuename(1:max(numel(valuename),2)));
            end
            if ~isstring(valuename),
                try,
                    valuename = toString(valuename);
                catch,
                    valuename = '?';
                end
            end
            valuename = strtrim(valuename);
            if numel(valuename)>40,
                valuename = valuename(:);
                valuename = valuename(1:min(numel(valuename),40));
                valuename = sprintf('%s...',valuename);
            end
            fprintf('Using %s options.%s = %s %s\n',mode,k,valuename,m);
        end
    end

