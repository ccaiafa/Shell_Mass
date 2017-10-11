function [semimajor_axis, semiminor_axis, x0, y0, phi,error,error1] = ellipse_fit_opt_V2(x, y)
D1 = [x.^2,x.*y, y.^2];
D2 = [x,y, ones(size(x))];
S1 = D1'*D1;
S2 = D1'*D2;
S3 = D2'*D2;
T = -inv(S3)*S2';
M = S1 + S2*T;
M = [M(3,:)./2; -M(2,:); M(1,:)./2];
[evec,eval] = eig(M);
cond = 4*evec(1,:).*evec(3,:) - evec(2,:).^2;
vect = evec(:,find(cond>0));

if ~isempty(vect)

    B = T*vect;
    a = vect(1);
    b = vect(2)/2;
    c = vect(3);
    d = B(1)/2;
    f = B(2)/2;
    g = B(3);

    delta = b^2-a*c;

    x0 = (c*d - b*f)/delta;
    y0 = (a*f - b*d)/delta;

    phi = 0.5 * acot((c-a)/(2*b));

    nom = 2 * (a*f^2 + c*d^2 + g*b^2 - 2*b*d*f - a*c*g);
    s = sqrt(1 + (4*b^2)/(a-c)^2);

    a_prime = sqrt(nom/(delta* ( (c-a)*s -(c+a))));

    b_prime = sqrt(nom/(delta* ( (a-c)*s -(c+a))));

    semimajor_axis = max(a_prime, b_prime);
    semiminor_axis = min(a_prime, b_prime);

    if (a_prime < b_prime)
        phi = pi/2 - phi;
    end
    % Given the quadratic form of an ellipse: 
    %  
    %       a*x^2 + 2*b*x*y + c*y^2  + 2*d*x + 2*f*y + g = 0   (1)
    % Calcular error
    %Construct M


    % Rotate points angle phi
    R = [cos(phi), sin(phi);
         -sin(phi), cos(phi)];

    for n = 1:size(x,1)
        rot = R*[(x-x0)';(y-y0)'];    
    end

    x = rot(1,:)';
    y = rot(2,:)';

    % scale
    x = x./semimajor_axis;
    y = y./semiminor_axis;

    error = x.^2 + y.^2 -1;
    error = sum(error.^2)/size(error,1);


    Npoints = 100;
    xref = zeros(Npoints,1);
    yref = zeros(Npoints,1);

    for n = 1:Npoints
        ang = n*2*pi/Npoints;
        xref(n) = cos(ang);
        yref(n) = sin(ang);
    end

    error1 = 0;
    for n = 1:size(x,1);
        d = (xref - x(n)).^2 +(yref - y(n)).^2;
        [C,I] = min(d); 

        error1 = error1 + C;
    end

    error1 = error1/size(x,1);

    % fig1 = figure
    % scatter(x,y)
    % ellipse(1,1,0,0,0,'r');
    % annotation(fig1,'textbox',...
    %     [0.726785714285714 0.253761904761905 0.0883928571428571 0.0619047619047619],...
    %     'String',{num2str(error)});



    % error = [x.^2 2*x.*y y.^2 2*x 2*y ones(size(x))]*[a;b;c;d;f;g];
    % error = sum(error.^2)/size(error,1);
    % error = error/(semimajor_axis*semiminor_axis);
else %% ERROR
    semimajor_axis = 0;
    semiminor_axis = 0;
    x0 = 0;
    y0 = 0;
    phi = 0;
    error = Inf;
    error1 = Inf;
    
end



end
