function showField(imagedata,titlestring)
    %{
    Plotting helper: show square rendered image on a [0 1]x[0 1] axis
    
    Parameters
    ----------
    imagedata: n x n x 3 array
        `n x n x [R G B]` image; RGB in [0,1]
    titlestring: string
        Title of plot
    %}
    imagesc([0 1], [0 1], imagedata); 
    hold on; 
    axis([0 1 0 1]);
    axis square;
    title(titlestring);
    set(gca,'yticklabel',{[]});
    set(gca,'xticklabel',{[]});
