function save_figure(i,inferredRGB,trueRGB,opt)
    %{
    Plotting helper for state inference. No user-accessible functionality.
    
    Parameters
    ----------
    trueRGB
    inferredRGB: RGB image array
    trueRGB: RGB image array
    options: struct
    %}
    NFCP_plotting.fix_figure(opt);
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    saveas(gcf         ,sprintf('inference_fields_frame_%d.png',i));
    saveas(gcf         ,sprintf('inference_fields_frame_%d.svg',i));
    %if isfield(opt,'have_groundtruth') && opt.have_groundtruth && trueRGB,
    if trueRGB,
        imwrite(trueRGB,sprintf('inference_fields_RGB_true__frame_%d.png',i));
    end
    if inferredRGB,
        imwrite(inferredRGB,sprintf('inference_fields_RGB_infer_frame_%d.png',i));
        %fprintf('Saved file inference_fields_RGB_infer_frame_%d.png',i);
    end
end
