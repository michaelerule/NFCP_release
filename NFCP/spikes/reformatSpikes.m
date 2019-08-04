function xydata = reformatSpikes(xydata,rebin,verbose)
    %{
    Convert spikes into a standardized format in order to handle slightly 
    different formats in some of the data.
    
    Parameters
    ----------
    xydata : `cell`
        Point data retrieved from one of the retinal wave datasets. 
        Points should be provided as a cell array of one (possibly empty)
        `Nspikes x 2` array of (x,y) locations. This function will 
        automatically detect transposed inputs. Spikes locations should be
        normalized so that all points like within [0,1]. If this is not
        the case, the function will automatically rescale points to fit.
    rebin : `int`, default 1
        Whether to collapse spikes into a courser timescale. Default of 1
        preserves original time binning. 
    verbose : `bool`, default `false`
        Whether to print detailed logging information
    
    Returns
    -------
    xydata : `cell`
        Corrected point data. 1D cell-array with one entry per time point.
        Each entry is either the empty matrix or an `Nspikes x 2` array of
        spikes with (x,y) locations in (0,1)²
    %}
    if nargin<2,
        rebin=1;
    end
    if nargin<3,
        verbose=false;
    end
    Nsample = size(xydata,2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % If data is not in [0,1]^2 recenter and rescale it
    for i=1:Nsample,
        transposed{i} = xydata{i}';
    end
    allxy = cell2mat(transposed);
    mins  = min(allxy'); 
    maxs  = max(allxy'); 
    if (any(mins<0) || any(maxs>1)),
        if verbose,
            fprintf(2,'Spikes outside of [0,1] interval; rescaling');
        end
        transformed = {};
        xscale = 1.0/(maxs(1)-mins(1));
        yscale = 1.0/(maxs(2)-mins(2));
        scale  = min(xscale,yscale);
        for i=1:Nsample,
            pxy = xydata{i};
            if numel(pxy)>0,
                px = (pxy(1:end,1) - mins(1)).*scale;
                py = (pxy(1:end,2) - mins(2)).*scale;
                transformed{i} = [px py];
            else,
                transformed{i} = pxy;
            end
        end
        xydata = transformed;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Optionally bin spikes at a courser time-scale
    % Can speed things up if slow-timescale is quite slow
    assert(rebin>0);
    newxydata = {};
    for i=1:Nsample,
        pxy = xydata{i};
        if numel(pxy)>0,
            j=floor((i-1)/rebin)+1;
            if j<=numel(newxydata),
                pxy2 = newxydata{j};
                pxy2 = [pxy2; pxy];
                newxydata{j} = pxy2;
            else,
                newxydata{j} = pxy;
            end
        end
    end
    xydata  = newxydata;
    Nsample = numel(xydata);
