%


fprintf('Inferring states...\n');
tic()
model.safety =0 ;
[llsum,infstate,margvar,infe] = stateInfer(ini,model,xydata,false,{...
    'doplot'      ,false 
    'peakactivity',false
    'showprog'    ,true
    });
toc()
fprintf('Log-likelihood is %f\n',llsum);


