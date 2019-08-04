function [x] = nearestPositivePoint(x0,P)
    %{
    Function nearestPositivePoint finds a point `x` based on the 
    multivariate Gaussian with mean `x0` and precision `P` such that
    `x` is the most probably point such that all entries in `x` are 
    positive.
    
    This is done using quadratic programming. 
    
    First, we transform into a coordinate system where Euclidian distance
    reflects distance in probability. 
    
    Then, we re-code the positivity constraint as a collection of 
    constraints in this tranformed space, building an open simplex.
    
    x > 0
    %}
    
    % Function to minimize is Euclidian norm in transformed space
    % Get as close to mean as possible
    % minimize (x-x0)' P (x-x0)
    % x' P x - x0' P x - x' P x0 + x0' P x0
    % 0.5 x' P x - x0' P x
    % 
    % Constaint is just positivity
    % x >= 0
    H = P;
    f = -x0'*P;
    m = size(P,1);
    A = -eye(m);
    b = zeros(m,1);
    options =  optimoptions('quadprog','Display','off');
    %x = quadprog(H,f,A,b,Aeq,beq,lb,ub,x0,options);
    x = quadprog(H,f,A,b,[],[],[],[],[],options);
        
    
