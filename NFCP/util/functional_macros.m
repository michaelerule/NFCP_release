%{
Functional programming macros for Matlab. `macros.m` should be called as
a script to load these macros as lambda expressions into the current 
workspace. 

These macros make it easier to define subroutines from within a script, 
without creating a new dedicated file. This allows for more compact and
maintainable scripts.

The main use case is in defining objective functions to be passed to 
optimization routines. These functions often need to close over scop
(e.g. large datasets), and are often very short and awkward to maintain
in a separate file. 

Matlab's restrictions on syntax lead to some limitaions, and these macros
might be more difficult to read in some use cases. 
%}

% call function with returned results f(g())
apply = @(f,args) f(args{:});

% lazy evaulation
evaluate = @(f)f();

% index cell arrays within expressions
cellitem = @(x,i) x{i};       

% index vectors within expressions
matitem = @(x,i) x(i);

% ternary if (non-lazy)
select = @(b,t,f) cellitem({t,f},2-(~~b));

% ternary if (lazy)
ifelse = @(b,t,f) evaluate(select(b,t,f));     

% safe cell2mat
any2mat = @(x) ifelse(iscell(x),@()cell2mat(x),@()x); 

% inline argmax
argmax = @(x) find(x==max(x));
