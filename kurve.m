%*********************************************************************** 
%
%     File_name :   KURVE.M
%
%    Dieses Programm bestimmt die Kurvenwerte
%    n2 und lab und par muessen bekannt sein
%    lab = 1 : Ellipse
%    lab = 2 : gerundetes Rechteck
%    lab = 3 : Drachen
%    lab = 4 : Erdnuss
%    Fuer jede Figuration wird zunaechst die Standardform berechnet.
%    Hierfuer werden die Parameter par(1), par(2) benutzt. Danach kann
%    um par(5) (in Grad, mathematisch positiv) gedreht und um Vektor
%    (par(3),par(4)) verschoben werden.
%***********************************************************************

function [x,dx,ddx,dddx] = kurve(n2, object, par)
pin = 2*pi/n2;
Ix = [0:n2-1]'*pin;
rot = par(5)*pi/180;
R = [ cos(rot) -sin(rot) ; sin(rot) cos(rot) ];

switch object
    
    case 'circle'
        %******************************************************
        %     Ellipse mit Halbachsen al=par(1) un2 be=par(2)
        %******************************************************
        al = par(1);
        be = par(2);
        x(:,1)  =  al*cos(Ix);
        x(:,2)  =  be*sin(Ix);
        dx(:,1) = -al*sin(Ix);
        dx(:,2) =  be*cos(Ix);
        ddx(:,1) = -x(:,1);
        ddx(:,2) = -x(:,2);
        dddx(:,1) = -dx(:,1);
        dddx(:,2) = -dx(:,2);
        
    case 'rectangle'
        %*****************************************************************
        %     gerundetes Rechteck mit Halbachsen al=par(1) und be=par(2)
        %*****************************************************************
        al = par(1);
        be = par(2);
        nr = 10;
        cs  = cos(Ix)/al;
        sn  = sin(Ix)/be;
        h1  = cs.^nr + sn.^nr;
        f   = h1.^(1/nr);
        h2  = ( sn.^(nr-2)/(be*be) - cs.^(nr-2)/(al*al) )./h1;
        s2  = sin(2*Ix);
        df  = -f/2.*s2.*h2;
        h3  = ( sn.^(nr-4)/(be^4) + cs.^(nr-4)/(al^4) )./h1;
        ddf = -( 5.5*df.*s2 + f.*cos(2*Ix) ).*h2 - 2*f.*s2.*s2.*h3;
        cs = cs*al;
        sn = sn*be;
        x(:,1)  =  f.*cs;
        x(:,2)  =  f.*sn;
        dx(:,1) =  df.*cs - x(:,2);
        dx(:,2) =  df.*sn + x(:,1);
        ddx(:,1) = ddf.*cs - df.*sn - dx(:,2);
        ddx(:,2) = ddf.*sn + df.*cs + dx(:,1);
                
    case 'kite'
        %-------------------------------------------------------------------
        %        Berechnet die Daten eines Drachens (kite) mit Parameter
        %        al=par(1) und be=par(2) (in Colton/kress: al=1.5, be=0.65)
        %-------------------------------------------------------------------
        al = par(1);
        be = par(2);
        x(:,1)  =  cos(Ix) + be*( cos(2*Ix) - 1 );
        x(:,2)  =  al*sin(Ix);
        dx(:,1) = -sin(Ix) - 2*be*sin(2*Ix);
        dx(:,2) =  al*cos(Ix);
        ddx(:,1) = -cos(Ix) - 4*be*cos(2*Ix);
        ddx(:,2) = -x(:,2);
        dddx(:,1) = sin(Ix) + 8*be*sin(2*Ix);
        dddx(:,2) = -dx(:,2);
                
    case 'nut'
        %------------------------------------------------------------------
        %        Berechnet die Daten einer Erdnuss (peanut) mit Parameter
        %        al=par(1) und be=par(2)
        %------------------------------------------------------------------
        al = par(1);
        be = par(2);
        f0 = sqrt( al*cos(Ix).^2 + be*sin(Ix).^2 );
        x(:,1)  =  f0.*cos(Ix);
        x(:,2)  =  f0.*sin(Ix);
        f1 = .5*(be-al)*sin(2*Ix)./f0;
        dx(:,1) = f1.*cos(Ix) - x(:,2);
        dx(:,2) = f1.*sin(Ix) + x(:,1);
        f2 = ( (be-al)*cos(2*Ix) - f1.*f1 ) ./ f0;
        ddx(:,1) = ( f2-f0).*cos(Ix) - 2*f1.*sin(Ix);
        ddx(:,2) = ( f2-f0).*sin(Ix) + 2*f1.*cos(Ix);
        f3 = ( 2*f1.*(f1.*sin(2*Ix) - 2*f0.*cos(2*Ix)) ) ./ f0.^3;
        f3 = f3 - ( f2.*sin(2*Ix) + 4*f0.*sin(2*Ix) ) ./ f0.^2;
        f3 = f3 * (be-al)/2;
        dddx(:,1) = f3.*cos(Ix) - 3*f2.*sin(Ix) - 2*f1.*cos(Ix) - dx(:,1);
        dddx(:,2) = f3.*sin(Ix) + 3*f2.*cos(Ix) - 2*f1.*sin(Ix) - dx(:,2);
        
end

scaling = par(6);
x   = scaling * (R*x')' + ones(n2,1)*[par(3) par(4)];
dx  = scaling * (R*dx')';
ddx = scaling * (R*ddx')';
dddx = scaling * (R*dddx')';



