function printModel(model)
    %{
    Print a model, displaying only those parameters over which we are 
    optimizing. Used for debugging and monitoring of parameter fits
    
    Parameters
    ----------
    model : struct
        Pre-initialized model struct; see `initializeModel(model)`
    %}
    keys = strsplit('rAA sigma thr gain alpha gamma ss_rescale');
    s = ['rates=' mat2str(model.linearRates(:)') '; '];
    for k=keys,
        k = k{1,1};
        s = strcat(s,sprintf('%s=%0.3e;',k,model.(k)));
    end
    fprintf(1,[s '\n']);
end
