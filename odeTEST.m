function dYdt = odeTEST(t,Y)  


    global sys 

    Kij = bdGetValue(sys.pardef,'Kij');
    % aee = bdGetValue(sys.pardef,'aee');
    aee = 0.36;
    aei = bdGetValue(sys.pardef,'aei');
    aie = bdGetValue(sys.pardef,'aie');
    ane = bdGetValue(sys.pardef,'ane');
    ani = bdGetValue(sys.pardef,'ani');
    b = bdGetValue(sys.pardef,'b');
    C = bdGetValue(sys.pardef,'C');
    r = bdGetValue(sys.pardef,'r');
    phi = bdGetValue(sys.pardef,'phi');
    Gion = bdGetValue(sys.pardef,'Gion');
    Vion = bdGetValue(sys.pardef,'Vion');
    thrsh = bdGetValue(sys.pardef,'thrsh');
    delta = bdGetValue(sys.pardef,'delta');
    VT = bdGetValue(sys.pardef,'VT');
    ZT = bdGetValue(sys.pardef,'ZT');
    deltaV = bdGetValue(sys.pardef,'deltaV');
    deltaZ = bdGetValue(sys.pardef,'deltaZ');
    I = bdGetValue(sys.pardef,'I');

% V1 = -0.01; V2 = 0.15; V3 = 0; V4 = 0.3; V5 = 0; V7 = 0; V9 = 0.3; V8 = 0.15;
% gCa = 1; gK = 2.0; gL = 0.5; gNa = 6.7;
% VK = -0.7; VL = -0.5; I = 0.3; b = 0.1; phi = 0.7; VNa = 0.53; VCa = 1;
% ani = 0.4; vs = 1; aei = 2; aie = 2; aee = 0.36; ane = 1; rnmda = 0.25;


    % Extract incoming values from Y
    % keyboard
    Y = reshape(Y,[],3);        % reshape Y to 3 columns
    V = Y(:,1);                 % 1st column of Y contains vector V
    W = Y(:,2);                 % 2nd column of Y contains vector W    
    Z = Y(:,3);                 % 3rd column of Y contains vector Z

    % Extract conductance parameters
    % keyboard
    gCa = Gion(1);
    gK  = Gion(2);
    gNa = Gion(3);
    gL  = Gion(4);
    
    % Extract Nerst potentials
    VCa = Vion(1);
    VK  = Vion(2);
    VNa = Vion(3);
    VL  = Vion(4);
    
    % Extract Gain threshold parameters
    TCa = thrsh(1);
    TK  = thrsh(2);
    TNa = thrsh(3);
    
    % Extract Gain slope parameters
    deltaCa = delta(1);
    deltaK  = delta(2);
    deltaNa = delta(3);
    % keyboard
    % Compute Firing-rate functions
    Qv = gain(V, VT, deltaV);       % (nx1) vector
    % keyboard
    Qz = gain(Z, ZT, deltaZ);       % (nx1) vector

    % Compute fraction of open channels
    mCa = gain(V, TCa, deltaCa);    % (nx1) vector
    mK  = gain(V, TK,  deltaK );    % (nx1) vector
    mNa = gain(V, TNa, deltaNa);    % (nx1) vector
    
    % Compute Mean firing rates
    k = sum(Kij,2);                 % (nx1) vector
    % QvMean = (Kij*Qv)./k;           % (nx1) vector
    QvMean = ((Qv'*Kij)')./k;           % (nx1) vector
    QvMean(isnan(QvMean)) = 0;    
    
    % Excitatory cell dynamics
    dV = -(gCa + (1-C).*r.*aee.*Qv + C.*r.*aee.*QvMean).*mCa.*(V-VCa) ...
         - gK.*W.*(V-VK) ...
         - gL.*(V-VL) ... 
         - (gNa.*mNa + (1-C).*aee.*Qv + C.*aee.*QvMean).*(V-VNa) ...
         + ane.*I ...
         - aie.*Qz.*Z;
     
    % K cell dynamics
    dW = phi.*(mK-W);
    
    % Inhibitory cell dynamics
    dZ = b.*(ani.*I + aei.*Qv.*V);

    % Return a column vector
    dYdt = [dV; dW; dZ]; 
end

% Non-linear gain function
function f = gain(VAR,C1,C2)
    f = 0.5*(1+tanh((VAR-C1)./C2));
end