%{
Run state inference while plotting at a fine time resolution
%}

fprintf('Inferring states...\n');
tic()
model.safety =0 ;
[llsum,infstate,margvar,infe] = stateInfer(ini,model,xydata,false,{...
    'doplot'      ,true 
    'upscale'     ,8
    'skipinf'     ,5
    'showduration',300/dt
    'showmaxy'    ,1
    'ratescale'   ,A_color_scale
    'peakactivity',false
    'points'      ,false
    'showprog'    ,true
    'softmax'     ,true
    'MarkerSize'  ,1
    });
toc()

fprintf('Log-likelihood is %f\n',llsum);


