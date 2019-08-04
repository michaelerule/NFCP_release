function x = min0rows(x)
    %{
    Shift rows so that minimum value is at zero
    %}
    x = bsxfun(@minus,x,min(x,[],2));
