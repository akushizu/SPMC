clear

parpool

for index = 1:2

    switch index

        case 1

            plot_rho()

        case 2

            plot_fw()

    end
end

delete(gcp('nocreate'))



function plot_rho()

    load adata

    %%{
    v = VideoWriter( 'muscato_rho.avi', 'Uncompressed AVI' ) ;
    v.FrameRate = 10 ;
    open(v)
    %}

    figure
    set(gcf,'color','w');

    fileID = fopen( 'rho.dat', 'r' ) ;
    textscan( fileID, '%*s\t%*s', 1 ) ;

    N = 0;
    
    while ~feof(fileID)
    
        textscan( fileID, '%*s', 1 ) ;
    
        C = textscan( fileID, '# = %d of %f', 1 ) ;
        i = C{1} ;
    
        C = textscan( fileID, 't = %f', 1 ) ;
        t = C{1} ;
    
        C = textscan( fileID, '%f\t%f' ) ;
        x = C{1, 1}' ;
        rho = C{1, 2}' ;
    
        if i == 0
    
            N = sum( rho ) ;
    
        end
    
        rho = rho / N / ( x(2) - x(1) ) ;
    
        p1 = plot( x, rho ) ;
        hold on
        plot( adata(i+1).x, 0.05 * exp( -(adata(i+1).x - 0).^2 / 2 / 1^2 ) )
        p2 = plot( adata(i+1).x , abs( adata(i+1).Amp ).^2 ) ;
        hold off
        title( "t = " + t + " fs" )
        xlabel( 'x (nm)' )
        ylabel( 'density' )
        axis([-30 30 -inf inf])
        legend( [ p1 p2 ], 'Monte Carlo', 'backward difference' )
        set(gca, 'FontName', 'DejaVu Sans')
        ax = gca ;
        ax.FontSize = 20 ;
        p1.LineWidth = 2 ;
        p2.LineWidth = 2 ;
    
        writeVideo( v , getframe(gcf) )
    
    end
    
    saveas( gcf, 'rho', 'svg' )
    fclose(fileID) ;

    close(v)

end

function plot_fw()

    load adata
    
    %%{
    v = VideoWriter( 'muscato_fw.avi', 'Uncompressed AVI' ) ;
    v.FrameRate = 10 ;
    open(v)
    %}
    
    figure%( 'Position' , [ 1 1 2 1 ] .* get( 0, 'defaultfigureposition' ) )
    set(gcf,'color','w');
    
    fileID = fopen( 'Fw.dat', 'r' ) ;
    textscan( fileID, '%*s', 1 ) ;
    
    C = textscan( fileID, '%f' ) ;
    x = C{1}' ;
    
    textscan( fileID, '%*s', 1 ) ;
    
    C = textscan( fileID, '%f' ) ;
    p = flipud( C{1} ) ;
    
    textscan( fileID, '%*s', 1 ) ;
    
    N = 0;
    %load bwr bwr
    cMap = getPyPlot_cMap( 'bwr' );

    while ~feof(fileID)
    
        textscan( fileID, '%*s', 1 ) ;
    
        C = textscan( fileID, '# = %d of %f', 1 ) ;
        i = C{1} ;
    
        C = textscan( fileID, 't = %f', 1 ) ;
        t = C{1} ;
    
        C = textscan( fileID, '%f' ) ;
    
        fw = reshape( C{1}, [ length(x) length(p) ] )' ;
    
        if i == 0
    
            N = sum( sum(fw) ) ;
    
        end
    
        fw = fw / N / ( x(2) - x(1) ) / ( p(1) - p(2) ) ;
    
        [ xx, pp ] = meshgrid( x, p ) ;
    
        contourf( xx, pp, fw )
        colormap( cMap )
        colorbar
        clim([-0.3 0.3])
        title( "fw, t = " + t + " fs" )
        xlabel('x (nm)')
        ylabel('p')
        set(gca, 'FontName', 'DejaVu Sans')
        ax = gca ;
        ax.FontSize = 16 ;
    
        writeVideo( v , getframe(gcf) )

    end
    
    saveas( gcf, 'Fw', 'svg' )
    fclose(fileID) ;
    
    close(v)

end
