function x = normrows(x)
    %{
    Ensure all rows sum to one (differs from matlab normr)
    %}
    x = bsxfun(@rdivide,x,sum(x,2));


