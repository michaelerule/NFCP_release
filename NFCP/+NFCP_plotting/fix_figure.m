function fix_figure(opt)
    %{
    Helper function for the stateInfer plotting routines, no user
    accessible functionality. 
    
    Tries to workaround Matlab figure focusing issues, and stop the
    animation from drawing to the wrong plot if user changes 
    focus. 
    
    TODO: replace this with axis handles; do it properly!
    %}
    f = gcf;
    if ~((f==opt.figure)||(f.Number==opt.figure)),
        % window focus has changed? 
        % This steals focus, annoying but what can you do?
        %figure(opt.figure);
        set(0,'CurrentFigure',opt.figure);
        fprintf(2,'hey the figure window changed?\n');
    end
end
    




