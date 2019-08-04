function f = centerfig(f,w,h)
        scsize = get( groot, 'Screensize' );
        sw = scsize(3);
        sh = scsize(4);
        x = (sw-w)/2;
        y = (sh-h)/2;
        set(f,'Position',[x y w h]);
