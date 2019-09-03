

function dx = DNA_damage_de(t,x,P) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DNA damage differential equations
%
% t: time
% x: vector or variables
% P: parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%
% Parameter values
%%%%%%%%%%%%%%%%%%%%%%

Vo = 2*10^(-3);
Vi = 7*10^(-6);
Kd = 0.0143;

k1t = 0.0440;    % Gammaoi
k2t = 29.086;    % Gammaio

k1 = 0.1270;     % Damage + N7 -> C7
k2 = 0.0973;     % Damage + N7 <- C7
k5 = 0.0130;     % C7 -> A7

k3 = 0.5003;     % Damage + O6 -> C6
k4 = 0.0457;     % Damage + O6 <- C6
k6 = 0.0159;     % C6 -> A6




%%%%%%%%%%%%%%%%%%%%%%
% Equations
%%%%%%%%%%%%%%%%%%%%%%

% TMZ Extracellular
F = -k1t*x(1)+k2t*(Vi/Vo)*x(2)-Kd*x(1);

% TMZ Intracellular
GG = k1t*(Vo/Vi)*x(1) - k2t*x(2)-Kd*x(2)-k1*x(2)*x(3)-k3*x(2)*x(5)+k2*x(4)+k4*x(6);

% N7 sites on DNA
G = -k1*x(2)*x(3)+k2*x(4);

% DNA adducts at N7 sites (C7)
H = k1*x(2)*x(3)-k2*x(4)-k5*x(4);

% O6 sites on DNA
I = -k3*x(2)*x(5)+k4*x(6);

% DNA adducts at O6 sites (C6)
J = k3*x(2)*x(5)-k4*x(6)-k6*x(6);

% Repaired N7
M = k5*x(4);

% Repaired O6
L = k6*x(6);


dx = [F; GG; G; H; I; J; M; L];

end



function F = DNA_damage_soln(P,TMZ,tdata,x0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DNA damage solution of differential equations
%
% P: parameters
% TMZ: dose of TMZ
% tdata: time
% x0 : initial conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%
% Initial Conditions
%%%%%%%%%%%%%%%%%%%%%%

x0 = [TMZ 0 13.95*TMZ 0 7.3196*TMZ 0 0 0]; 

%%%%%%%%%%%%%%%%%%%%%%
% Solve the equations
%%%%%%%%%%%%%%%%%%%%%%

[t,F1] = ode23s(@DNA_damage_de,tdata,x0,[],P); 


%%%%%%%%%%%%%%%%%%%%%%
% Desired outputs
%%%%%%%%%%%%%%%%%%%%%%

% Total DNA
N = F1(:,3)+F1(:,4)+F1(:,5)+F1(:,6)+F1(:,7)+F1(:,8);

F = [F1(:,7) F1(:,8) N];

end



function dx = DNA_repair_de(t,x,P) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DNA repair differential equations
%
% t: time
% x: vector or variables
% P: parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%
% Parameter values
%%%%%%%%%%%%%%%%%%%%%%

kfa = 0.025-0.001;     % D7 + APNG -> C7
kra = 0.001;           % D7 + APNG <- C7
kpa = 0.35;            % C7 -> DNAr

kfm = 0.0916-0.001;    % D6 + MGMT -> C6
krm = 0.001;           % D6 + MGMT <- C6
kpm = 0.1717;          % C6 -> DNAr

Sa = P(1);             % Production of APNG
Sm = P(2);             % Production of MGMT

lambdaA = P(3);      % Degradation rate of APNG
lambdaM = P(4);      % Degradation rate of MGMT

%%%%%%%%%%%%%%%%%%%%%%
% Equations
%%%%%%%%%%%%%%%%%%%%%%


% N7-meG adducts
M = -kfa*x(3)*x(8)+kra*x(5);

% O6-meG adducts
L = -kfm*x(4)*x(7)+krm*x(6);

% N7-TMZ complex
N = kfa*x(3)*x(8)-kra*x(5)-kpa*x(5);

% O6-TMZ complex
O = kfm*x(4)*x(7)-krm*x(6)-kpm*x(6);

% Repaired N7 sites on DNA
G = kpa*x(5);

% Repaired O6 sites on DNA
I = kpm*x(6);

% MGMT
P = Sm - lambdaM*x(7)-kfm*x(4)*x(7)+krm*x(6);

% APNG 
Q = Sa - lambdaA*x(8)-kfa*x(3)*x(8)+kra*x(5);



%%%%%%%%%%%%%%%%%%%%%%
% Output
%%%%%%%%%%%%%%%%%%%%%%

dx = [G; I; M; L; N; O; P; Q];

end




function F = DNA_repair_soln(P,TMZ,tdata,x0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DNA repair solution of differential equations
%
% P: parameters
% TMZ: dose of TMZ
% tdata: time
% x0 : initial conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%
% Solve the equations
%%%%%%%%%%%%%%%%%%%%%%

[t,F1] = ode23s(@DNA_repair_de,tdata,x0,[],P); 


%%%%%%%%%%%%%%%%%%%%%%
% Output
%%%%%%%%%%%%%%%%%%%%%%

% Total DNA
N = F1(:,3)+F1(:,4)+F1(:,5)+F1(:,6)+F1(:,1)+F1(:,2);

F = [F1(:,1) F1(:,2) N F1(:,7) F1(:,8) F1(:,3) F1(:,4)];

end



function Y= Population_Model(Na,Mt,t,con,alpha,K,APNg,MGMt,lambdaA,lambdaM)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Age-structure Tumor Population
%
% Na: Number of age points
% Mt: Number of time points
% t: time vector
% con: dose of TMZ
% alpha: Proliferation rate
% K: carrying capacity
% APNg: Production of APNG
% MGMt: Production of MGMT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Time steps for got variables
ht = 1/Mt;
ha = 1/Na;

% Create age vector
for m=1:Na
    age(m) = (m-1)*ha;
end

%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%%%


mu1 = 174.2872;    % Rate of cell arrest
rho1 = 2.9288;     % Rate of cell repair

b1 = 2.6411;       % Rate of cell death caused by N7-meG damage
b2 = 7.2196;       % Rate of cell death caused by O6-meG damage

c1 = 24.383;       % Rate of cell death caused by long arrest
c2 = 55.0068;      % Effect of MGMT on cell death
c3 = 0.0112;       % Effect of N7-meG adducts on cell death

a0 = age(Na/2);    % Characteristic waiting time of arrested cells



%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initiation fo Variables
%%%%%%%%%%%%%%%%%%%%%%%%%%

A = zeros(Na,Mt);     % Arrested Cells

% Only Time Dependent Variables

P = zeros(1,Mt);      % Proliferative cells
N = zeros(1,Mt);      % Total number of cells


% Initial Conditions
P(1) = 1e+05;
N(1) = 1e+05;



if con ~=0 % If treatment is given
    
F2 = DNA_damage_soln(1,con,t,x0_damage);

mu = mu1*((F2(:,1)+F2(:,2))./F2(:,3));   % Cell arrest

% Explicit recursion in time
for m=1:(Mt-1)
    clear F;
    totalDNA = (F2(m,1)+F2(m,2));
    
    % DNA repair
    Ninit = F2(m,1);
    Oinit = F2(m,2);
    MGMT_init = F2(m,7);
    APNG_init = F2(m,8);
    
    % Comment when doing local sensitivity analysis
    lambdaA = 0.0057;      % Degradation rate of APNG
    lambdaM = 0.0043;      % Degradation rate of MGMT


    x0_repair = [0 0 Ninit  Oinit 0 0 MGMT_init APNG_init];
    
    F = DNA_repair_soln([APNg MGMt lambdaA lambdaM],con,age,x0_repair);
    
    if totalDNA==0 
        rho=0*F(:,1);  % If there is no damage DNA, repair is 0
    else
        rho = rho1*(((F(:,1)+F(:,2))./totalDNA));  % Cell repair
    end

    % Time iteration in Proliferating cells
    P(m+1) = P(m)+ht*(alpha*P(m)*(1-(P(m)/(K)))-mu(m)*P(m)+sum(rho.*A(:,m))*ha);
    
    % Boundary condition in Arrested cells
    A(1,m+1) = mu(m)*P(m);

   % Explicit recursion in age 
   for l=2:Na-1
      
       delta6 = b1*(1./(1+exp(-c1*(age(l)-a0))))*exp(-c2*F(l,4));          % Death by O6-meG adducts
       delta7 = b2*(1./(1+exp(-c1*(age(l)-a0))))*(1./(1+exp(-c3*F(l,6)))); % Death by N7-meG adducts

       A(l,m+1)= A(l,m)+ht*(-(1/(2*ha))*(A(l+1,m)-A(l-1,m))-delta6*A(l,m)-delta7*A(l,m)-rho(l)*A(l,m));

   end 
   
   % Boundary condition in Arrested cells
   A(end,m+1)=A(end-1,m+1);

   % Total tumor cells
    N(m+1) = P(m+1)+sum(A(:,m+1))*ha;
end

else % If there is no treatment
    for m=2:Mt
    P(m)=P(m-1)+ht*(alpha*P(m-1)*(1-(P(m-1)/K)));
    N(m)=P(m);
    end
end


Y = [P' N'];

end


