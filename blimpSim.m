%% Blimp simulation
% By Piero Alessandro Riva Riquelme
% Undergraduate student.
% Universidad de Concepción, Chile
% pieroriva@udec.cl
% November of 2020

% Based on work of Michale Frye, Stephen Gammon and Chunjiang Qian
% "The 6-DOF dynamic model and simulation of the tri-turbofan remote-controlled airship".
% Based on work of S.B.V. Gomes 
% "An investigation of the flight dynamics of airships with application to the YEZ-2A".
% "Airship dynamic modeling for autonomous operation".

%% Step 0: Cleaning
clc
clear all
%close all

%% Step 1: Parameters
  fprintf('\nStep1: Parameters...');
  % a) Physical constants
    rhoAir  = 1.205; % Density of air at NTP (20°C, 1atm)
    rhoHe   = 0.1664; % Density of helium at NTP (20°C, 1atm)
    g       = 9.80665;    % acceleration of gravity
    gradISA = -0.0065; % Temperature gradient versus height K/m 
    Temp    = 20; %°C
    Pressure = 1; %atm
    DegToRad = pi/180;
    RadToDeg = DegToRad^-1;
    msTokmh = 60^2/1000;
    % b) Vehicle geometry and parameters (Ellipsoid) [m,Kg]
    a = 0.9;        % a length
    b = 0.45;       % b length
    c = 0.45;       % c length
    Vol = 4*pi*a*b*c/3; % Volume
    Area = Vol^(2/3);   % Area
    mBlimp = Vol*(rhoAir-rhoHe);   % Blimp mass
    mHe = rhoHe*Vol;        % Mass due to Helium
    mTotal = mBlimp+mHe;    % Total mass
    dx = a/8;   % X axis distance from CV (Center of Volume) to propellers
    dy = b/2;   % Y axis distance from CV to propellers
    dz = c;     % Z axis distance from CV to propellers
    ax = 0;     % X axis distance from center of gravity CG and center of volume CV
    ay = 0;     % Z axis distance from center of gravity CG and center of volume CV
    az = -(0.2*mBlimp)*b/mTotal;  % Z axis distance from center of gravity CG and center of volume CV
    displacedAirMass = Vol*rhoAir;

    % c) Masses and inertias
    Ix = mTotal*(b^2+c^2)/5; % Inertia of an ellipsoid along the x axis
    Iy = mTotal*(c^2+a^2)/5; % Inertia of an ellipsoid along the y axis
    Iz = mTotal*(a^2+b^2)/5; % Inertia of an ellipsoid along the z axis
    
    % c.1) Lamb's coefficients of inertia for a prolate ellipsoid
    eLamb      = sqrt(1-c^2/a^2);
    gammaLamb  = 2*(1-eLamb^2)*(0.5*log((1+eLamb)/(1-eLamb))-eLamb)/eLamb^3;
    alphaLamb  = (1/eLamb^2)-(1-eLamb^2)/(2*eLamb^3)*log((1+eLamb)/(1-eLamb));

    k1 = gammaLamb/(2-gammaLamb);
    k2 = alphaLamb/(2-alphaLamb);
    k_ = (eLamb^4*(alphaLamb-gammaLamb))/((2-eLamb^2)*(2*eLamb^2-(2-eLamb^2)*(alphaLamb-gammaLamb)));
    
    % c.2) Tuckerman for a prolate ellipsoid
    eTuckerman = eLamb;
    alphaTuckerman = (1-eTuckerman^2)*(log((1+eTuckerman)/(1-eTuckerman))-2*eTuckerman)/eTuckerman^3;
    betaTuckerman = (1-eTuckerman^2)*((eTuckerman/(1-eTuckerman^2))-0.5*log((1+eTuckerman)/(1-eTuckerman)))/eTuckerman^3;
    gammaTuckerman = betaTuckerman;
    
    K1 = Vol*(alphaTuckerman/(2-alphaTuckerman));
    K2 = Vol*(betaTuckerman/(2-betaTuckerman));
    K3 = Vol*(gammaTuckerman/(2-gammaTuckerman));
    K1_ = Vol*Ix*(((b^2-c^2)/(b^2+c^2))^2*((gammaTuckerman-betaTuckerman)/(2*((b^2-c^2)/(b^2+c^2))-(gammaTuckerman-betaTuckerman))));
    K2_ = Vol*Iy*(((c^2-a^2)/(c^2+a^2))^2*((alphaTuckerman-gammaTuckerman)/(2*((c^2-a^2)/(c^2+a^2))-(alphaTuckerman-gammaTuckerman))));
    K3_ = Vol*Iz*(((a^2-b^2)/(a^2+b^2))^2*((betaTuckerman-alphaTuckerman)/(2*((a^2-b^2)/(a^2+b^2))-(betaTuckerman-alphaTuckerman))));
    
    % c.3) Virtual masses and inertias
        % Lamb
            %Xu = -k1*displacedAirMass; % Lambs
            %Yv = -k2*displacedAirMass;
            %Zw = -k2*displacedAirMass;
            %Lp = 0;
            %Mq = -k_*displacedAirMass*Iy; % Munk
            %Nr = -k_*displacedAirMass*Iz;
        % Tuckerman
            Xu = -K1*rhoAir;
            Yv = -K2*rhoAir;
            Zw = -K3*rhoAir;
            Lp = 0;
            Mq = -K2_*rhoAir;
            Nr = -K3_*rhoAir;
        % Others (Gomes)
            Mu = 0;
            Lv = 0;
            Nv = 0;
            Mw = 0;
            Yp = 0;
            Xq = 0;
            Zq = 0;
            Yr = 0;
        % Group
            mx = mTotal-Xu;
            my = mTotal-Yv;
            mz = mTotal-Zw;
            Jx = Ix-Lp;
            Jy = Iy-Mq;
            Jz = Iz-Nr;
            Jxz = 0;

    % d) M matrix
        M = [mx             0                   0               0               mTotal*az-Xq    0;...
             0              my                  0               -mTotal*az-Yp   0               mTotal*ax-Yr;...
             0              0                   mz              0               -mTotal*ax-Zq   0;...
             0              -mTotal*az-Lv       0               Ix-Lp           0               -Jxz;...
             mTotal*az-Mu   0                   -mTotal*ax-Mw   0               Iy-Mq           0;...
             0              mTotal*ax-Nv         0               -Jxz            0               Iz-Nr];
     
    % d.1) M inverse
        invM = inv(M);
    
%% Step 2: Flight simulation config
  fprintf('\nStep2: Initializing the tail configuration...');    
    % a) Tools
    h = @(t) heaviside(t);      % Heaviside
    r = @(t) heaviside(t)*t;    % Ramp
    
    % b) Time definition
    ti = 10.1;         % Initial sim time
    step = 0.1;        % Step size
    tf = 560;          % Final sim time
    tt = ti:step:tf;
    nData = length(tt);
    
    % c) Process configuration and declaration
    nInputs = 3;
    nStates = 6;
    tol = 1e-6;
    
    u = zeros(nInputs,nData);
    X = zeros(nStates,nData);   % Zero initial conditions
    Y = zeros(nStates,nData);   % Zero initial conditions
    Vearth = zeros(3,nData);    % Zero initial conditions
    Xearth = zeros(3,nData);    % Zero initial conditions
    DCM = zeros(3);             % DCM
    lambda = zeros(3);          % Lambda's
    dX = zeros(nStates,nData);  % Zero initial conditions
    Vtot = zeros(1,nData);      % Rectiliniar speed
    wTot = zeros(1,nData);      % Angular speed
    G  = zeros(nStates,nData);  % Gravitational force
    A  = zeros(nStates,nData);  % Aerodynamic force
    Fd = zeros(nStates,nData);  % Dynamic forces
    P  = zeros(nStates,nData);  % Propulsion forces
    
    % d) Input generation
    Amp  = 0.05;            % Motor maximum forces [N]
    Amp2 = 90*DegToRad;     % Motor maximum angle rotation mu
    u1 = @(t) h(t)*-Amp*0+h(t-100)*Amp*0.9-h(t-200)*Amp*0.9;
    u2 = @(t) h(t)*-Amp*0+h(t-100)*Amp*0.9-h(t-200)*Amp*0.9;
    u3 = @(t) h(t)*Amp2*0; 
    for i=1:nData
        u(1,i) = (u1(tt(i)))+(randn-0.5)*Amp*0.025*0;
        u(2,i) = (u2(tt(i)))+(randn-0.5)*Amp*0.025*0;
        u(3,i) = (u3(tt(i)))+(randn-0.5)*Amp2*0.025*0;
        if u(1,i)<-Amp
            u(1,i)=-Amp;
        end
        if u(1,i)>Amp
            u(1,i)=Amp;
        end
        if u(2,i)<-Amp
            u(2,i)=-Amp;
        end
        if u(2,i)>Amp
            u(2,i)=Amp;
        end
        if u(3,i)<-Amp2
            u(3,i)=-Amp2;
        end
        if u(3,i)>Amp2
            u(3,i)=Amp2;
        end
    end
    

%% Step 3: Simulation start

  for n=2:nData
    % a) Dynamics vector, Fd
    f1 = -mz*X(3,n-1)*X(5,n-1)+my*X(6,n-1)*X(2,n-1)+mTotal*(ax*(X(5,n-1)^2+X(6,n-1)^2)-az*X(6,n-1)*X(4,n-1));
    f2 = -mx*X(1,n-1)*X(6,n-1)+mz*X(4,n-1)*X(3,n-1)+mTotal*(-ax*X(4,n-1)*X(5,n-1)-az*X(6,n-1)*X(5,n-1));
    f3 = -my*X(2,n-1)*X(4,n-1)+mx*X(5,n-1)*X(1,n-1)+mTotal*(-ax*X(6,n-1)*X(4,n-1)+az*(X(5,n-1)^2+X(4,n-1)^2));
    f4 = -(Jz-Jy)*X(6,n-1)*X(5,n-1)+Jxz*X(4,n-1)*X(5,n-1)+mTotal*az*(X(1,n-1)*X(6,n-1)-X(4,n-1)*X(3,n-1));
    f5 = -(Jx-Jz)*X(4,n-1)*X(6,n-1)+Jxz*(X(6,n-1)^2-X(4,n-1)^2)+mTotal*(ax*(X(2,n-1)*X(4,n-1)-X(5,n-1)*X(1,n-1))-az*(X(3,n-1)*X(5,n-1)-X(6,n-1)*X(2,n-1)));
    f6 = -(Jy-Jx)*X(5,n-1)*X(4,n-1)-Jxz*X(5,n-1)*X(6,n-1)+mTotal*(-ax*(X(1,n-1)*X(6,n-1)-X(4,n-1)*X(3,n-1)));
    Fd(:,n-1) = [f1 f2 f3 f4 f5 f6]';
    
    % b) Propulsion vector, P
    P1 = (u(1,n-1)+u(2,n-1))*cos(u(3,n-1));
    P2 = 0;
    P3 = -(u(1,n-1)+u(2,n-1))*sin(u(3,n-1));
    P4 = (u(2,n-1)-u(1,n-1))*sin(u(3,n-1))*dy;
    P5 = (u(1,n-1)+u(2,n-1))*(dz*cos(u(3,n-1))-dx*sin(u(3,n-1)));
    P6 = (u(2,n-1)-u(1,n-1))*cos(u(3,n-1))*dy; 
    P(:,n-1) =[P1 P2 P3 P4 P5 P6]';
    
    
    % c) Aerodynamic force vector, A
    
    % c.2) Calculate longitudinal and lateral incidences: alpha and beta (not in use)
    Vtot(n-1) = sqrt(X(1,n-1)^2+X(2,n-1)^2+X(3,n-1)^2);
    wTot(n-1) = sqrt(X(4,n-1)^2+X(5,n-1)^2+X(6,n-1)^2);
    % Easy way
      CD = 0.9;
      CY = 0.9;
      CL = 0.9;
      Cl = 0.9;
      Cm = 0.9;
      Cn = 0.9; 
    
    % c.4) Aerodynamic forces A
    coefA1 = 0.5*rhoAir*Vtot(n-1)^2*Area;
    coefA2 = 0.5*rhoAir*Vtot(n-1)^2*Vol;
    
    coefB1 = 0.5*rhoAir*X(1,n-1)^2*sign(X(1,n-1))*Area;
    coefB2 = 0.5*rhoAir*X(2,n-1)^2*sign(X(2,n-1))*Area;
    coefB3 = 0.5*rhoAir*X(3,n-1)^2*sign(X(3,n-1))*Area;
    coefB4 = 0.5*rhoAir*X(4,n-1)^2*sign(X(4,n-1))*Vol;
    coefB5 = 0.5*rhoAir*X(5,n-1)^2*sign(X(5,n-1))*Vol;
    coefB6 = 0.5*rhoAir*X(6,n-1)^2*sign(X(6,n-1))*Vol;
    
    A1 = -CD*coefB1;
    A2 = -CY*coefB2;
    A3 = -CL*coefB3;
    A4 = -Cl*coefB4;
    A5 = -Cm*coefB5;
    A6 = -Cn*coefB6;
    A(:,n-1) = [A1 A2 A3 A4 A5 A6]';
    
    % d) Gravitational force vector, G 
    % d.1) Calculate the direction cosines matrix (DCM)
    lambda(1,1) = cos(Y(5,n-1))*cos(Y(6,n-1));
    lambda(1,2) = cos(Y(5,n-1))*sin(Y(6,n-1));
    lambda(1,3) = sin(Y(5,n-1));
    lambda(2,1) = (-cos(Y(4,n-1))*sin(Y(6,n-1))+sin(Y(4,n-1))*sin(Y(5,n-1))*cos(Y(6,n-1)));
    lambda(2,2) = (cos(Y(4,n-1))*cos(Y(6,n-1))+sin(Y(4,n-1))*sin(Y(5,n-1))*sin(Y(6,n-1)));
    lambda(2,3) = sin(Y(4,n-1))*cos(Y(5,n-1));
    lambda(3,1) = (sin(Y(4,n-1))*sin(Y(6,n-1))+cos(Y(4,n-1))*sin(Y(5,n-1))*cos(Y(6,n-1)));
    lambda(3,2) = (-sin(Y(4,n-1))*cos(Y(6,n-1))+cos(Y(4,n-1))*sin(Y(5,n-1))*sin(Y(6,n-1)));
    lambda(3,3) = cos(Y(4,n-1))*cos(Y(5,n-1));
    DCM = lambda;        
    
    % d.2) Calculate gravitational forces & moments. 
    B = rhoAir*g*Vol; % constant
    W = mTotal*g;     % constant
    
	% d.3) Using DCM elements convert these to body axes to obtain gravity vector G
        G1 = lambda(3,1)*(W-B);
        G2 = lambda(3,2)*(W-B);
        G3 = lambda(3,3)*(W-B);
        G4 = -lambda(3,2)*az*W;
        G5 = (lambda(3,1)*az-lambda(3,3)*ax)*W;
        G6 = lambda(3,2)*ax*W;
        G(:,n-1) = [G1 G2 G3 G4 G5 G6]';
        
    % e) Calculate linear and angular accelerations du, dv, dw, dp, dq, dr from [acc] =M-1[[P]+[Fd]+[A]+[G]]
    dX(:,n-1) = invM*(P(:,n-1)+Fd(:,n-1)+A(:,n-1)+G(:,n-1));
    
    % f) Integrate accelerations and obtain linear and angular body axes velocities u, v, w, p, q, r
    for i=1:nStates
        X(i,n) = sum(dX(i,1:n-1))*step+step*(dX(i,n-1)-dX(i,1))/2;
        if abs(X(i,n))<tol
            %X(i,n) = tol;
        end
    end
    
    % g) Transform linear velocities to earth axes to obtain Vnorth, Veast, Vup
    Vearth(:,n) = DCM*X(1:3,n);
    for i=1:3
        Xearth(i,n) = sum(Vearth(i,1:n-1))*step+step*(Vearth(i,n-1)-Vearth(i,1))/2;
    end
    
    % h) Calculate vehicle position in terms of displacements in the north, east and vertical directions
    for i=1:nStates
        Y(i,n) = sum(X(i,1:n-1))*step+step*(X(i,n-1)-X(i,1))/2; % (XYZ) aligned with (NEU)
    end
    
    % End
    fprintf('\nn = %d',n);
  end

%% Step 4: Plot Results
fprintf('\nPlotting results...');
lineWidth = 1.5;

% a) Figure(1): Inputs
    figure(1)
    subplot(5,1,1);
    cla
    grid minor
    hold on;
    plot(tt,u(1,:),'-','LineWidth',lineWidth);
    legend('u1');
    title('Starboard motor force Tds'); %xlabel('Time (s)'); 
    ylabel('(N)');

    subplot(5,1,2);
    cla
    grid minor
    hold on;
    plot(tt,u(2,:),'-','LineWidth',lineWidth);
    legend('u2');
    title('Port motor force Tsp'); %xlabel('Time (s)'); 
    ylabel('(N)');

    subplot(5,1,3);
    cla
    grid minor
    hold on;
    plot(tt,u(3,:).*RadToDeg,'-','LineWidth',lineWidth);
    legend('u3');
    title('Propeller angle \mu'); %xlabel('Time (s)');  
    ylabel('(°)');

    subplot(5,1,4);
    cla
    grid minor
    hold on;
    plot(tt,Vtot,'-','LineWidth',lineWidth);
    %legend('');
    title('Net speed'); xlabel('Time (s)');  ylabel('(m/s)');

    subplot(5,1,5);
    cla
    grid minor
    hold on;
    plot(tt,wTot.*RadToDeg,'-','LineWidth',lineWidth);
    %legend('');
    title('Net angular speed'); xlabel('Time (s)');  ylabel('(°/s)');

% b) Figure(2): Velocities
    figure(2)
    subplot(3,2,1);
    cla
    grid minor
    hold on;
    plot(tt,X(1,:),'-','LineWidth',lineWidth);
    title('u = X axis velocity'); xlabel('Time (s)'); ylabel('Speed (m/s)');
    
    subplot(3,2,3);
    cla
    grid minor
    hold on;
    plot(tt,X(2,:),'-','LineWidth',lineWidth);
    title('v = Y axis velocity'); xlabel('Time (s)'); ylabel('Speed (m/s)');
    
    subplot(3,2,5);
    cla
    grid minor
    hold on;
    plot(tt,X(3,:),'-','LineWidth',lineWidth);
    title('w = Z axis velocity'); xlabel('Time (s)'); ylabel('Speed (m/s)');
    
    subplot(3,2,2);
    cla
    grid minor
    hold on;
    plot(tt,X(4,:).*RadToDeg,'-','LineWidth',lineWidth);
    title('p = X axis angular velocity'); xlabel('Time (s)'); ylabel('Angular speed (°/s)');
    
    subplot(3,2,4);
    cla
    grid minor
    hold on;
    plot(tt,X(5,:).*RadToDeg,'-','LineWidth',lineWidth);
    title('q = Y axis angular velocity'); xlabel('Time (s)'); ylabel('Angular speed (°/s)');
    
    subplot(3,2,6);
    cla
    grid minor
    hold on;
    plot(tt,X(6,:).*RadToDeg,'-','LineWidth',lineWidth);
    title('r = Z axis angular velocity'); xlabel('Time (s)'); ylabel('Angular speed (°/s)');
    
% c) Figure(3): Position and orientation
    figure(3)
    subplot(3,2,1);
    cla
    grid minor
    hold on;
    plot(tt,Y(1,:),'-','LineWidth',lineWidth);
    title('x = X axis Position'); xlabel('Time (s)'); ylabel('Position (m)');
    
    subplot(3,2,3);
    cla
    grid minor
    hold on;
    plot(tt,Y(2,:),'-','LineWidth',lineWidth);
    title('y = Y axis Position'); xlabel('Time (s)');  ylabel('Position (m)');
    
    subplot(3,2,5);
    cla
    grid minor
    hold on;
    plot(tt,Y(3,:),'-','LineWidth',lineWidth);
    title('z = Z axis Position'); xlabel('Time (s)');  ylabel('Position (m)');
    
    subplot(3,2,2);
    cla
    grid minor
    hold on;
    plot(tt,Y(4,:).*RadToDeg,'-','LineWidth',lineWidth);
    title('\phi = X axis angle'); xlabel('Time (s)'); ylabel('Angle (°)');
    
    subplot(3,2,4);
    cla
    grid minor
    hold on;
    plot(tt,Y(5,:).*RadToDeg,'-','LineWidth',lineWidth);
    title('\theta = Y axis angle'); xlabel('Time (s)'); ylabel('Angle (°)');
    
    subplot(3,2,6);
    cla
    grid minor
    hold on;
    plot(tt,Y(6,:).*RadToDeg,'-','LineWidth',lineWidth);
    title('\psi = Z axis angle'); xlabel('Time (s)'); ylabel('Angle (°)');

% d) Figure(4): Forces
    figure(4)
    subplot(2,2,1);
    cla
    grid minor
    hold on;
    plot(tt,Fd(1,:),'-','LineWidth',lineWidth);
    plot(tt,Fd(2,:),'-','LineWidth',lineWidth);
    plot(tt,Fd(3,:),'-','LineWidth',lineWidth);
    plot(tt,Fd(4,:),'-','LineWidth',lineWidth);
    plot(tt,Fd(5,:),'-','LineWidth',lineWidth);
    plot(tt,Fd(6,:),'-','LineWidth',lineWidth);
    title('Fd(X(t))'); xlabel('Time (s)'); ylabel('Dynamic forces (N)');
    legend('x','y','z','roll','pitch','yaw');
    
    subplot(2,2,2);
    cla
    grid minor
    hold on;
    plot(tt,G(1,:),'-','LineWidth',lineWidth);
    plot(tt,G(2,:),'-','LineWidth',lineWidth);
    plot(tt,G(3,:),'-','LineWidth',lineWidth);
    plot(tt,G(4,:),'-','LineWidth',lineWidth);
    plot(tt,G(5,:),'-','LineWidth',lineWidth);
    plot(tt,G(6,:),'-','LineWidth',lineWidth);
    title('G(X(t))'); xlabel('Time (s)'); ylabel('Gravity force (N)');
    legend('x','y','z','roll','pitch','yaw');
    
    subplot(2,2,3);
    cla
    grid minor
    hold on;
    plot(tt,A(1,:),'-','LineWidth',lineWidth);
    plot(tt,A(2,:),'-','LineWidth',lineWidth);
    plot(tt,A(3,:),'-','LineWidth',lineWidth);
    plot(tt,A(4,:),'-','LineWidth',lineWidth);
    plot(tt,A(5,:),'-','LineWidth',lineWidth);
    plot(tt,A(6,:),'-','LineWidth',lineWidth);
    title('A(X(t))'); xlabel('Time (s)'); ylabel('Aerodynamic forces (N)');
    legend('x','y','z','roll','pitch','yaw');
    
    subplot(2,2,4);
    cla
    grid minor
    hold on;
    plot(tt,P(1,:),'-','LineWidth',lineWidth);
    plot(tt,P(2,:),'-','LineWidth',lineWidth);
    plot(tt,P(3,:),'-','LineWidth',lineWidth);
    plot(tt,P(4,:),'-','LineWidth',lineWidth);
    plot(tt,P(5,:),'-','LineWidth',lineWidth);
    plot(tt,P(6,:),'-','LineWidth',lineWidth);
    title('P(u(t))'); xlabel('Time (s)'); ylabel('Propulsion forces (N)');
    legend('x','y','z','roll','pitch','yaw');
      
%% END
  save("blimpSim.mat",'tt','u','dX','X','Y');
  fprintf('\n\nEND\n');
    
    
    
    
    
    