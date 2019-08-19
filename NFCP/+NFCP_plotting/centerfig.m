function f = centerfig(f,w,h)
    %{
    Centers a figure on screen

    Parameters
    ----------
    f: figure
        Figure to center
    w: int
        Desired width of figure
    h: int
        Desired height of figure
    %}
    scsize = get( groot, 'Screensize' );
    sw = scsize(3);
    sh = scsize(4);
    x = (sw-w)/2;
    y = (sh-h)/2;
    set(f,'Position',[x y w h]);
