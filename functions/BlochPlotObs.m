function BlochPlotObs(obs)
    sig0 = eye(2);
    sig1 = [0+0i,1+0i;1+0i,0+0i];
    sig2 = [0+0i,-1i;1i,0+0i];
    sig3 = [1+0i,0+0i;0+0i,-1+0i] ;
    
    r0 = trace(obs' * sig0);
    r1 = trace(obs' * sig1);
    r2 = trace(obs' * sig2);
    r3 = trace(obs' * sig3);
    
    % taken from https://www.youtube.com/watch?v=qysCuvPdX6E
    
    figure(1)
    hold on
    
    
    [Xs, Yx, Zx] = sphere(25);
    mySphere = surf( Xs, Yx, Zx );
    axis equal
    shading interp
    mySphere.FaceAlpha = 0.25;
    
    line([-1 1], [0 0], [0 0])
    line([0 0], [-1 1], [0 0])
    line([0 0], [0 0], [-1 1])
    
    text(0, 0, 1.1, '$\left| 0 \right>$', 'Interpreter', 'latex', 'FontSize', 20, 'HorizontalAlignment', 'Center')
    text(1.1, 0, 0, '$\left| + \right>$', 'Interpreter', 'latex', 'FontSize', 20, 'HorizontalAlignment', 'Center')
    text(-1.1, 0, 0, '$\left| - \right>$', 'Interpreter', 'latex', 'FontSize', 20, 'HorizontalAlignment', 'Center')    
    text(0, 0, -1.1, '$\left| 1 \right>$', 'Interpreter', 'latex', 'FontSize', 20, 'HorizontalAlignment', 'Center')

    view([60 12])
    
    p1 = [0 0 0];
    p2 = [r1 r2 r3];
    dp = p2-p1;
    quiver3(p1(1),p1(2),p1(3),dp(1),dp(2),dp(3),'LineWidth',3,'Color','r')
    hold off
end

