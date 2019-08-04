function op = sparsifyUnitary(op,truncate)
    %{
    Creates a sparisfied version of op by setting values smaller than
    threshold to zero. Attempts to ensure that probability is conserved
    by re-normalizing the resuling matrix rows to sum to one. Repeats this
    procedure if re-normalization violates the sparsity constraint. 
    
    Parameters
    ----------
    op : square matrix,
        Probability-preserving operator to sparsify
    truncate : number, default `1e-4`
        Small positive number below which (relative to mean value) we
        should set entry to zero
        
    Returns
    -------
    op : matrix
        Operator with values smaller than truncate*max(abs(operator)) set
        to zero.
    %}
    
    % Round values smaller than `truncate` down to zero to create a 
    % sparse operator.
    if nargin<2, truncate = 1e-4; end

    % normalize rows so density is conserved
    op  = bsxfun(@rdivide, op', sum(op,2))';
    %remove small entries that more or less don't matter
    big = max(abs(op(:)));
    toosmall = truncate*big;
    
    for iter=1:3,
        % normalize rows so density is conserved
        op = bsxfun(@rdivide, op', sum(op,2))';
        remove = abs(op)<toosmall;
        remove = remove.*remove';
        op(remove>0) = 0.0;
        % sparsify
        op = sparse(op);
    end
end
