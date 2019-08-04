function varargout=voke(f,varargin)
    %{
    cached_invoke is a way of memoizing large matlab functions
    
    All arguments should be integers or strings -- the sort of thing you'd
    want to see in a file name.     These arguments will be used to build a
    .mat file containg the output of fname called on the arguments.
    
    If fname requires large amounts of data, as is likely if you've found need
    of this functionality, suggest having fname load this data from disk.
    If you require parameterization over a finite number of values of type
    that is neither integer nor string, suggest passing an index to an array
    stored on dist.
    
    This caching function should only be used for intensive computations
    where loading arguments from disk (slow) is acceptable.
    
    This function cannot distinguish between different versions of a 
    function with the same name, which can cause invalid cache retrievals.

    The default cache location is the local directory.
    
    Parameters
    ----------
    f : function to invoke
    varargin : arguments to forward to f
    
    Returns
    result : 
        f(varargin{:})
    %}
    display(['( ' func2str(f) ' )'])
    fname = func2str(f);
    cachedir = '.';
    var      = ['__' fname '_cache'];
    for i=[1:size(varargin,2)]
        cachedir=[cachedir '/' var];
        if exist(cachedir)~=7
            mkdir(cachedir);
        end
        var = varargin{i};
        if ~isa(var,'char')
            var=int2str(var);
        end
        display(['Argument ' int2str(i) ' is ' var]);
    end
    cachename = fname;
    for i=[1:size(varargin,2)]
        var = varargin{i};
        if ~isa(var,'char')
            var=int2str(var);
        end
        cachename = [cachename '_' var];
    end
    cachename = strrep(cachename,'/','_');
    dispname  = strrep(cachename,'_',' ');
    display([': ' dispname])
    cachename = [cachedir '/' cachename '.mat'];
    display(['Cache name is ' cachename]);    
    display(['nargsout is ' int2str(nargout)])
    if ~exist(cachename,'file')
        disp(['did not find  ' cachename]);
    end
    try
        load(cachename);
        disp(['loaded        ' cachename]);
    catch err
        disp(['computing     ' cachename]);
        if nargout<=1
            varargout={f(varargin{:})};
        else
            c = cell(nargout,1);
            [c{:}]=f(varargin{:});
            varargout=c;
        end
        display(['output type: ' class(varargout)])
        %varargout
        save(cachename,'varargout');
    end



