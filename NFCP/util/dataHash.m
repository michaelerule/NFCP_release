function H = dataHash(Data)
    %{
    Jan Simon's struct hash function
    mathworks.com/matlabcentral/answers/3314-hash-function-for-matlab-struct
    
    Parameters
    ----------
    Data : struct
    
    Returns
    -------
    H : hash value of data
    %}
    Engine = java.security.MessageDigest.getInstance('MD5');
    H = CoreHash(Data, Engine, 0);
    H = sprintf('%.2x', H);   % To hex string

function H = CoreHash(Data, Engine, recursion_depth)
    %{
    Recursively walk through a `struct` and build a hash value.
    Warning: this function appears to behave erratically in Matlab
    versions older than r2017a; Have not identified the bug yet
    
    Parameters
    ----------
    Data : any variable within the structure
    Engine : Hashing engine
    recursion_dept : recusion depth counter
    
    Returns
    -------
    H : hash value of data
    %}
    if recursion_depth>100000,
        H = 0;
        return
    end
    
    % Consider the type of empty arrays:
    S = [class(Data), sprintf('%d ', size(Data))];
    Engine.update(typecast(uint16(S(:)), 'uint8'));
    H = double(typecast(Engine.digest, 'uint8'));
    if isa(Data, 'struct')
        n = numel(Data);
        if n == 1  % Scalar struct:
           F = sort(fieldnames(Data));  % ignore order of fields
           for iField = 1:length(F)
              H = bitxor(H, CoreHash(Data.(F{iField}), Engine, recursion_depth+1));
           end
        else  % Struct array:
           for iS = 1:n
              H = bitxor(H, CoreHash(Data(iS), Engine, recursion_depth+1));
           end
        end
    elseif isempty(Data)
        % No further actions needed
    elseif isnumeric(Data)
        Data = full(Data); % no sparse arrays allowed
        Engine.update(typecast(Data(:), 'uint8'));
        H = bitxor(H, double(typecast(Engine.digest, 'uint8')));
    elseif ischar(Data)  % Silly TYPECAST cannot handle CHAR
        Engine.update(typecast(uint16(Data(:)), 'uint8'));
        H = bitxor(H, double(typecast(Engine.digest, 'uint8')));
    elseif iscell(Data)
        for iS = 1:numel(Data)
           H = bitxor(H, CoreHash(Data{iS}, Engine, recursion_depth+1));
        end
    elseif islogical(Data)
        Engine.update(typecast(uint8(Data(:)), 'uint8'));
        H = bitxor(H, double(typecast(Engine.digest, 'uint8')));
    elseif isa(Data, 'function_handle')
        H = bitxor(H, CoreHash(functions(Data), Engine, recursion_depth+1));
    else
        warning(['Type of variable not considered: ', class(Data)]);
    end

