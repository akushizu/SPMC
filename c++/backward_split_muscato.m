clear

parpool

parfor index = 1:2

    switch index

        case 1

            backward();

        case 2

            split();

    end
end

compare();

delete(gcp('nocreate'))


function backward()

    hb = 6.62607015 / 1.602176634 / 2 / pi ;
    m = 0.32 ;
    e = 1 ;
    k0 = 0.7 ;
    s0 = 2.852 ;
    x0 = -15 ;
    xc = 0 ;
    a = 0.3 ;
    s = 1 ;
    
    %--------step size----------------------
    dx = 0.005 ;
    dt = 0.001 ;
    l = ( 0.5 * sqrt(-1) * hb / m ) * ( dt / dx^2 ) ;
    
    x = -40 : dx : 40 ;
    t = 0 : dt : 20 ;
    
    %------initial & boundaty condition------------------
    w = sqrt(sqrt( 1 / (2 * pi * s0^2) )) * exp( -(x - x0).^2 / (4 * s0^2) ) .* exp( sqrt(-1) * k0 * (x - xc) ) ;
    
    N = 1 ; %sqrt( dx * sum( abs(w).^2 ) ) ;
    
    w = w / N ;
    
    %---------potential------------------------
    Vw = e * a * exp( -(x - xc).^2 / (2 * s^2) ) ;
    
    %---------backward difference------------------------
    W = spdiags( ( 1 + 2 * l ) * ones( length(x) , 1 ) + ( sqrt(-1) * dt / hb ) * Vw' , 0 , speye( length(x) ) ) ;
    W = spdiags( -l * ones( length(x) , 1 ) , 1 , W ) ;
    W = spdiags( -l * ones( length(x) , 1 ) , -1 , W ) ;
    
    adata( ( length(t) - 1 ) / 50 + 1 ) = struct( 't' , [] , 'x' , [] , 'Amp' , [] , 'time' , [] ) ;
    adata( 1 ) = struct( 't' , t(1) , 'x' , x , 'Amp' , w , 'time' , 0 ) ;
    
    tic
    
    for i = 2:length(t)
        
        w1 = w / W ;
        
        N = 1 ; %sqrt( dx * sum( abs(w1).^2 ) ) ;
        
        w = w1 / N ;
        
        if mod( ( i - 1 ) , 50 ) == 0
            
            adata( ( i - 1 ) / 50 + 1 ) = struct( 't' , t(i) , 'x' , x , 'Amp' , w , 'time' , toc ) ;
            
        end
        
    end
    
    toc
    
    %----------save----------------------------
    save adata adata

    %------------------------------------------
    x = adata(1).x ;
    t = [ adata.t ] ;
    
    figure( 'Position' , [ 0.5 0.5 2 2 ] .* get( 0, 'defaultfigureposition' ) )
    
    col = lines ;
    
    %%{
    v = VideoWriter( 'wob.avi', 'Uncompressed AVI' ) ;
    v.FrameRate = 12 ;
    %v.Quality = 100 ;
    open(v)
    %}
    
    for i = 1:length(t)
        
        subplot( 2 , 2 , 1 ) ;
        plot( x , 0.5*Vw , '--' , 'Color' , col( 2 , : ) ) ;
        hold on
        plot( x , abs( adata(i).Amp ).^2 , 'Color' , col( 1 , : ) ) ;
        hold off
        axis( [ x(1) x(end) -inf inf ] ) ;
        xlabel('x') ;
        ax = gca ;
        ax.FontSize = 16 ;
        title( [ 'den, t = ' num2str( t(i) ) 'fs' ] ) ;
        ax = gca ;
        ax.FontSize = 16 ;
        
        subplot( 2 , 2 , 2 ) ;
        plot( x , Vw , '--' , 'Color' , col( 2 , : ) ) ;
        hold on
        plot( x , abs( adata(i).Amp ) , 'Color' , col( 1 , : ) ) ;
        hold off
        axis( [ x(1) x(end) -inf inf ] ) ;
        xlabel('x') ;
        title( 'abs' ) ;
        ax = gca ;
        ax.FontSize = 16 ;
        
        subplot( 2 , 2 , 3 ) ;
        plot( x , Vw , '--' , 'Color' , col( 2 , : ) ) ;
        hold on
        plot( x , abs( adata(i).Amp ) , 'k:' , x , -abs( adata(i).Amp ) , 'k:' ) ;
        plot( x , real( adata(i).Amp ) , 'Color' , col( 1 , : ) ) ;
        hold off
        axis( [ x(1) x(end) -inf inf ] ) ;
        xlabel('x') ;
        title( 'real' ) ;
        ax = gca ;
        ax.FontSize = 16 ;
        
        subplot( 2 , 2 , 4 ) ;
        plot( x , Vw , '--' , 'Color' , col( 2 , : ) ) ;
        hold on
        plot( x , abs( adata(i).Amp ) , 'k:' , x , -abs( adata(i).Amp ) , 'k:' ) ;
        plot( x , imag( adata(i).Amp ) , 'Color' , col( 1 , : ) ) ;
        hold off
        axis( [ x(1) x(end) -inf inf ] ) ;
        xlabel('x') ;
        title( 'imag' ) ;
        ax = gca ;
        ax.FontSize = 16 ;
        
        writeVideo( v , getframe(gcf) )
        
    end
    
    close(v)

end

function split()

    hb = 6.62607015 / 1.602176634 / 2 / pi ;
    m = 0.32 ;
    e = 1 ;
    k0 = 0.7 ;
    s0 = 2.852 ;
    x0 = -15 ;
    xc = 0 ;
    a = 0.3 ;
    s = 1 ;
    
    %--------step size----------------------
    N = 16384 ;
    dt = 0.001 ;
    
    x = linspace( -40, 40, N ) ;
    dx = x(2) - x(1) ;
    k = linspace( -pi / dx, pi / dx, N ) ;
    t = 0 : dt : 20 ;
    
    %------initial & boundaty condition------------------
    psi = sqrt(sqrt( 1 / (2 * pi * s0^2) )) * exp( -(x - x0).^2 / (4 * s0^2) ) .* exp( sqrt(-1) * k0 * (x - xc) ) ;
    V = e * a * exp( -(x - xc).^2 / (2 * s^2) ) ;
    
    Ut = exp( -sqrt(-1) * (hb^2 * k.^2 / m / 2) * dt / hb ) ;
    Uv = exp( -sqrt(-1) * V * dt / 2 / hb ) ;
    
    edata( ( length(t) - 1 ) / 50 + 1 ) = struct( 't' , [] , 'x' , [] , 'Amp' , [] , 'time' , [] ) ;
    edata( 1 ) = struct( 't' , t(1) , 'x' , x , 'Amp' , psi , 'time' , 0 ) ;
    
    tic
    
    for i = 2:length(t)
        
        psi = Uv .* psi ;
        psih = fftshift( fft(psi) ) ;
        psih = Ut .* psih ;
        psi = ifft( ifftshift(psih) ) ;
        psi = Uv .* psi ;
        
        if mod( ( i - 1 ) , 50 ) == 0
            
            edata( ( i - 1 ) / 50 + 1 ) = struct( 't' , t(i) , 'x' , x , 'Amp' , psi , 'time' , toc ) ;
            
        end
        
    end
    
    toc
    
    %----------save----------------------------
    save edata edata
    
    %------------------------------------------
    x = edata(1).x ;
    t = [ edata.t ] ;
    
    figure( 'Position' , [ 0.5 0.5 2 2 ] .* get( 0, 'defaultfigureposition' ) )
    
    col = lines ;
    
    %%{
    v = VideoWriter( 'wos.avi', 'Uncompressed AVI' ) ;
    v.FrameRate = 12 ;
    %v.Quality = 100 ;
    open(v)
    %}
    
    for i = 1:length(t)
        
        subplot( 2 , 2 , 1 ) ;
        plot( x , 0.5*V , '--' , 'Color' , col( 2 , : ) ) ;
        hold on
        plot( x , abs( edata(i).Amp ).^2 , 'Color' , col( 1 , : ) ) ;
        hold off
        axis( [ x(1) x(end) -inf inf ] ) ;
        xlabel('x') ;
        ax = gca ;
        ax.FontSize = 16 ;
        title( [ 'den, t = ' num2str( t(i) ) 'fs' ] ) ;
        ax = gca ;
        ax.FontSize = 16 ;
        
        subplot( 2 , 2 , 2 ) ;
        plot( x , V , '--' , 'Color' , col( 2 , : ) ) ;
        hold on
        plot( x , abs( edata(i).Amp ) , 'Color' , col( 1 , : ) ) ;
        hold off
        axis( [ x(1) x(end) -inf inf ] ) ;
        xlabel('x') ;
        title( 'abs' ) ;
        ax = gca ;
        ax.FontSize = 16 ;
        
        subplot( 2 , 2 , 3 ) ;
        plot( x , V , '--' , 'Color' , col( 2 , : ) ) ;
        hold on
        plot( x , abs( edata(i).Amp ) , 'k:' , x , -abs( edata(i).Amp ) , 'k:' ) ;
        plot( x , real( edata(i).Amp ) , 'Color' , col( 1 , : ) ) ;
        hold off
        axis( [ x(1) x(end) -inf inf ] ) ;
        xlabel('x') ;
        title( 'real' ) ;
        ax = gca ;
        ax.FontSize = 16 ;
        
        subplot( 2 , 2 , 4 ) ;
        plot( x , V , '--' , 'Color' , col( 2 , : ) ) ;
        hold on
        plot( x , abs( edata(i).Amp ) , 'k:' , x , -abs( edata(i).Amp ) , 'k:' ) ;
        plot( x , imag( edata(i).Amp ) , 'Color' , col( 1 , : ) ) ;
        hold off
        axis( [ x(1) x(end) -inf inf ] ) ;
        xlabel('x') ;
        title( 'imag' ) ;
        ax = gca ;
        ax.FontSize = 16 ;
        
        writeVideo( v , getframe(gcf) )
        
    end
    
    close(v)

end

function compare()

    load adata
    load edata
    
    e = 1 ;
    xc = 0 ;
    a = 0.12 ;
    s = 1 ;
    
    x = -40 : 0.005 : 40 ;
    V = e * a * exp( -(x - xc).^2 / (2 * s^2) ) ;
    
    figure
    
    col = lines ;
    
    v = VideoWriter( 'woc.avi', 'Uncompressed AVI' ) ;
    v.FrameRate = 12 ;
    %v.Quality = 100 ;
    open(v)
    
    for i = 1:length(adata)
              
        plot( x , 0.5*V , '--' , 'Color' , col( 3 , : ) ) ;
        hold on
        p1 = plot( adata(i).x , abs( adata(i).Amp ).^2 , 'Color' , col( 1 , : ) ) ;
        p2 = plot( edata(i).x , abs( edata(i).Amp ).^2 , 'Color' , col( 2 , : ) ) ;
        hold off
        axis( [ x(1) x(end) -inf inf ] ) ;
        xlabel('x') ;
        title( [ 'den, t = ' num2str( adata(i).t ) 'fs' ] ) ;
        legend( [ p1 p2 ], 'backward difference', 'split operator' )
        
        writeVideo( v , getframe(gcf) )
        
    end
    
    close(v)

end
