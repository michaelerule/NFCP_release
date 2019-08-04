function update_plot(model,i,M,xydata,true_states,saved_means,inf_means,inf_errs,opt)
    %{
    Plotting helper for state inference. No user-accessible functionality.
    
    Parameters
    ----------
    model : struct
        pre-initialized model structure; see `initializeModel`
    i : int
        Current filtering iteration
    M : vector  
        Current filtered mean estiamtes
    xydata 
    true_states
    saves_means
    inf_means
    inf_errs
    opt : struct
        Options struct, see `stateInfer` documentation
    
    %}
    NFCP_plotting.fix_figure(opt);
    newplot; hold off; 
    
    % x-axis limits
    a = max(0,(i-opt.showduration)*abs(model.dt));
    b = a + opt.showduration*abs(model.dt);
    
    if opt.mergeRstates, nplot = 4;
    else nplot = model.nstates+1; end
    
    if opt.have_groundtruth,
        % Show simulated (ground truth) in left subplot
        % Plot average intensities over time (true; if available)
        trueM = real(true_states{i});
        NFCP_plotting.fix_figure(opt);
        subplot(231); 
        cla;
        trueRGB = NFCP_plotting.fieldsToRGB(trueM,model.n,model.cscale,opt.upscale,opt.softmax);
        NFCP_plotting.showField(trueRGB,sprintf('True (%d out of %d)',i,model.Nsample));
        if opt.points, 
            NFCP_plotting.add_points(i,xydata,max(1,i-opt.skipinf),opt); 
        end
        % Plot ground truth 
        NFCP_plotting.fix_figure(opt); subplot(2,3,[2 3]); cla;
        plot(real(saved_means(1:i,1:nplot)),'LineWidth',1);
        ylabel('Mean concentration'); 
        if opt.mergeRstates, legend([model.names(1:2) 'R' 'Total']);
        else                 legend([model.names 'Total']); end
        grid on; xlim([a b]); ylim([0 opt.showmaxy]); set(gca,'xticklabel',{[]});
    end
    
    % Plot the spatiotemporal states
    inferredRGB = NFCP_plotting.fieldsToRGB(M,model.n,model.cscale,opt.upscale,false);
    NFCP_plotting.fix_figure(opt);
    if opt.have_groundtruth, subplot(234); % show both true and inferred
    else subplot(131); end % only show inferred
    cla; 
    NFCP_plotting.showField(inferredRGB,'Reconstructed');
    if opt.points, 
        NFCP_plotting.add_points(i,xydata,max(1,i-min(opt.poolpoints,opt.skipinf)),opt); 
    end
    
    % Plot the true values, if available (check that these lie within
    % confidence intervals), 
    % or plot the inferred means if the true means are not given
    times = (1:i)*abs(model.dt);

    % First line group: mean lines (ground truth if available)
    % these will show the average concentrations in each state, 
    % plus the total number of agents estimated    
    NFCP_plotting.fix_figure(opt);
    if opt.have_groundtruth, subplot(2,3,[5 6]); useforplot=saved_means;
    else                     subplot(1,3,[2 3]); useforplot=inf_means; end
    cla; 
    
    % There should be NSTATES + 2 (or 3+2 if refractory are merged)
    % different signals to plot. The first ones reflect the concentrations
    % The rest...
    
    if opt.mergeRstates, ntoplot = 3;
    else,                ntoplot = model.nstates; end
    for j=1:ntoplot,
        plot(times,real(useforplot(1:i,j)),'color',model.colors(j,1:end),'LineWidth',1);
        hold on; 
    end
    for j=ntoplot:(size(useforplot,2)-1),
        plot(times,real(useforplot(1:i,j)),'LineWidth',1);
        hold on; 
    end
    xlabel('Time step'); ylabel('Inferred means');
    
    % Plot observed counts averaged over arrays
    c = cell2mat(cellfun(@(x)size(x,1),xydata,'uni',false));
    % Convert counts to (heuristic) activation
    c = (c./model.nn-mean(model.bias(:)))./(model.gain*model.volume);
    c = c.*opt.ratescale;
    c = c.*model.alpha;
    plot(times,c(1:i));
    
    % Label lines
    names = model.names;
    if opt.mergeRstates, names=[names(1:2),'R']; end
    legend([names,'Total','Empirical activation']);
    
    % Shade confidence intervals
    xx    = [1:i, fliplr(1:i)].*abs(model.dt);
    % Plot quiescent state with error bar
    meanQ = inf_means (1:i,1);
    errQ  = inf_errs(1:i,1);
    fill(xx,[meanQ-errQ;flipud(meanQ+errQ)]',1, 'facecolor','blue','edgecolor','none','facealpha',.3);
    % Plot active state with error bar
    meanA = inf_means (1:i,2);
    errA  = inf_errs(1:i,2);
    fill(xx,[meanA-errA;flipud(meanA+errA)]',1, 'facecolor','red','edgecolor','none','facealpha',.3);
    % Plot refractory states with error bar
    % TODO: must handle merged refractory states properly here!
    meanR = inf_means (1:i,3);
    errR  = inf_errs(1:i,3);
    fill(xx,[meanR-errR; flipud(meanR+errR)]', 1, 'facecolor','green','edgecolor','none','facealpha', 0.3);
    %{
    % Plot total count with error bar
    meanC = inf_means (1:i,4);
    errC  = inf_errs(1:i,4);
    plot(real(meanC),'LineWidth',1);
    fill(xx,[meanC-errC; flipud(meanC+errC)]', 1,'facecolor','magenta','edgecolor','none','facealpha', 0.3);
    %}
    
    % Clean up plot   
    NFCP_plotting.fix_figure(opt); 
    grid on; xlim([a b]); ylim([0 opt.showmaxy]); 
    hold off; f=getframe(gcf); clear f; 
    
    if opt.save_figure,
        if opt.have_groundtruth, NFCP_plotting.save_figure(i,inferredRGB,trueRGB,opt);
        else                     NFCP_plotting.save_figure(i,inferredRGB,false,opt); end
    end
end




