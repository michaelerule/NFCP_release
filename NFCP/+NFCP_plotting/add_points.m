function add_points(i,xydata,last_frame,opt)
    %{
    Plotting helper for state inference. No user-accessible functionality.
    
    Parameters
    ----------
    i : int
        current frame
    xydata : cell
        point proces data
    last_frame : int
        most recent frame rendered
    opt : dict
        options
    %}
    NFCP_plotting.fix_figure(opt);
    while last_frame<=i,
        if (size(xydata{last_frame},1)>0)
            plot(xydata{last_frame}(:,1),...
                 xydata{last_frame}(:,2),'.w','MarkerSize',opt.MarkerSize);
        end
        last_frame = last_frame+1;
    end
end



