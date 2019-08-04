function ntoclear = overprint(msg,ntoclear)
    %{
    Erase current line then print message
    
    Parameters
    ----------
    msg : `string`
        Message to print
    ntoclear : `int`
        Number of characters to clear
        
    Returns
    -------
    ntoclear : `int`
        Length of `msg`; number of characters to clear on next iteration
    %}
	fprintf(repmat('\b', 1, ntoclear));
	fprintf(msg);
	ntoclear = numel(msg);
