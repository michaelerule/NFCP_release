%

%rE = 1;
%rA = 0.5;
%rR = 0.1;

fprintf('======================================================================\n');
fprintf('Loading file \n%s\n',fn);
fprintf('----------------------------------------------------------------------\n');
load(filepath);

bisa = bias(:);
inhomogeneous_gain = inhomogeneous_gain(:);

fprintf('Original sample rate : %0.1f Hz\n', Fs_Hz_original);
fprintf('Frames in recording  : %d      \n', stop_frame_original);
fprintf('Recording duration   : %d s    \n', duration_seconds);
fprintf('Spikes binned at     : %0.1f s \n', binsize_seconds);
fprintf('Frames, downsampled  : %d      \n', size(spikexy,2));
fprintf('The grid size is     : %dx%d   \n', N,N);
fprintf('Dispersion (1=Pois.) : %f      \n', alpha);
fprintf('Estimated A->Q rate  : %f 1/s  \n', rE);
fprintf('Estimated A->R rate  : %f 1/s  \n', rA);
fprintf('Estimated R    rate  : %f 1/s  \n', rR);
fprintf('Recovery time        : %f s    \n', recovery_timescale);
fprintf('Recovery states      : %f      \n', recovery_nstates);
fprintf('Average gain of      : %f      \n', gain);
fprintf('======================================================================\n');


N = double(N);

NFRAMES = size(spikexy,2);

xydata  = {};
for i=1:NFRAMES,
    xydata{i} = (spikexy{i}-1)./(64+1e-7);
    assert(all(all(xydata{i}>=0)));
    assert(all(all(xydata{i}<=1)));
end

% bias should be positive
bias(bias<0)=0;



%{
+--------------------------+----------------+---------------------+
| Variables                | Units...       | ...per units        | 
+==========================+================+=====================+
| model.dt                 | seconds      / | bin_t               |
+--------------------------+----------------+---------------------+
| model.dx                 | array_x²     / | bin_x²              |
+--------------------------+----------------+---------------------+
| model.n                  | bin_x        / | array_x             |
+--------------------------+----------------+---------------------+
| model.nn                 | bin_x²       / | array_x²            |
+--------------------------+----------------+---------------------+
| model.volume             | s∙array_x²   / | bin_x²∙bin_t        |
+--------------------------+----------------+---------------------+
| model.bias               | spikes       / | bin_x²∙bin_t        |
+--------------------------+----------------+---------------------+
| model.gain               | spikes       / | s∙array_x²          |
+--------------------------+----------------+---------------------+
| model.alpha              | spike mean   / | spike standard dev. |
+--------------------------+----------------+---------------------+
| model.inhomogeneous_gain | unitless                             |
+--------------------------+----------------+---------------------+
| model.adjusted_bias      | unitless                             |
+--------------------------+----------------+---------------------+
| model.adjusted_gain      | alpha∙spikes / | bin_x²∙bin_t        |
+--------------------------+----------------+---------------------+
| rate (inside obj. fun.)  | spikes       / | bin_x²∙bin_t        |
+--------------------------+----------------+---------------------+
%}    

%model.bias = model.bias .* rebin;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Heuristically coarsen time-binning 
% We selected 0.1 second bins because some datasets contain fast waves
% However, most datasets contain slower waves
% We should coursen the timescales to handle this
%{
rebin   = min(10,max(1,floor(recovery_timescale/dt/10)));
if rebin>1,
    fprintf(1,'Recovery time %0.1f suggests coarser timescale\n',recovery_timescale);
    fprintf(1,'... Downsampling by a factor of %d\n',rebin);
    xydata  = reformatSpikes(xydata,rebin);
    NFRAMES = numel(xydata);
    dt      = dt*rebin;
    binsize_seconds = dt;
    model.dt = dt;
    bias = bias./rebin;
    NFRAMES = size(xydata,2);
    fprintf(1,'... There are now %d frames\n',NFRAMES);
    fprintf(1,'... Time bin size is now %0.2f s\n',dt);
end
%}

fprintf('----------------------------------------------------------------------\n');

