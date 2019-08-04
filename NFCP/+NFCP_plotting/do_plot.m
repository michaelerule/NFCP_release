function [inf_means,inf_errs] = do_plot(model,i,M,C,xydata,true_states,saved_means,inf_means,inf_errs,opt)
    %{
    Plotting helper for state inference. No user-accessible functionality.
    
    Parameters
    ----------
    model : struct
        pre-initialized model structure; see `initializeModel`
    i : int
        Current filtering iteration
    M : vector  
        Current filtered mean estiamte
    C : matrix
        Current filtered covariance estimate
    xydata :
    true_states :
    saves_means :
    inf_means :
    inf_errs :
    opt : struct
        Options struct, see `stateInfer` documentation
    %}

    % Store spatially averaged means for display
    inf_means(i,1:opt.K) = opt.meanproj*M(:);
    if opt.peakactivity,
        % Special case: visualize max over spacefor the active A state
        [maxval, maxidx] = nanmax(M(model.nn+1:model.nn*2));
        inf_means(i,2) = maxval*opt.ratescale;
    end
    
    %{
    % Estimate expected activation directly from counts
    % This will be plotted as a sanity check against inferred A density
    count = pointsToHistogram(xydata{i},model.n,model.interpolate,0);
    count = model.blur(count.*model.alpha);
    if strcmp(model.link,'linear'),
        Aest = count(:)./model.premultiplied_gain-model.premultiplied_bias;
    else,
        Aest = log(count(:)./(model.volume.*model.alpha.*model.bias))./model.gain;
    end
    % remove dead regions
    Aest  = Aest(model.premultiplied_gain>0.0); 
    % use the maximum over the spatial domain to let localized events show through
    if opt.peakactivity, f=@nanmax; else, f=@nanmean; end;
    %}
    
    var = diag(opt.meanproj*C*opt.meanproj');
    inf_errs(i,1:opt.K) = sqrt(var).*1.96;
    if opt.peakactivity,
        % Special case: visualize max over spatial area for active A state
        idx = maxidx + model.nn;
        inf_errs(i,2) = sqrt(C(idx,idx)).*1.96*opt.ratescale;
    end    
    
    % If the opt.doplot flag is true, show progress as a plot
    if mod(i,opt.skipinf)==0
        NFCP_plotting.update_plot(...
            model,i,M,xydata,true_states,saved_means,inf_means,inf_errs,opt);
    end
end


