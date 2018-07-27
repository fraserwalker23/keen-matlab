function out = keen2()
%%
close all
clc
warning off

%% Parameters attribution

nu = 3; %capital to output ratio
alpha = 0.025; %productivity growth rate
beta = 0.02; %population growth rate
delta = 0.01; %depreciation rate
r = 0.03; %real interest rate
%r_d = 0.02; % interest rate on deposits
phi0 = 0.04/(1-0.04^2); %Phillips curve parameter
phi1 = 0.04^3/(1-0.04^2); %Phillips curve parameter
kappa_eq = nu*(alpha+beta+delta); 
%kappa_prime_eq = 5;
kappa0 = -0.0065 ;
kappa1 = exp(-5);
kappa2 = 20;
l=[];

%% Other definitions
TOL = 1E-7;
options = odeset('RelTol', TOL);
txt_format = '%3.3g';
marker = 2+2;

% random numebr generation
ran = makedist('Rayleigh','B',0.75);

%% Functions Definition
fun_kappa = @(x) kappa0+kappa1.*exp(kappa2.*x);
fun_kappa_prime = @(x) kappa1.*kappa2.*exp(kappa2.*x);
fun_kappa_inv = @(x) log((x-kappa0)./kappa1)./kappa2;

fun_Phi = @(x) phi1./(1-x).^2- phi0;
fun_Phi_prime = @(x) 2*phi1./(1-x).^3;
fun_Phi_inv = @(x) 1 - (phi1./(x+phi0)).^(1/2);
%% Sample paths in keen2

% Additional parameters

k_r = 0.08; % target capital adequacy ratio, 0 < k_r < 1
eta_p = 4; % speed, eta_p > 0
% pi_eq = fun_kappa_inv(kappa_eq);
% omega_eq = 1 - pi_eq - r*((kappa_eq - pi_eq)/(alpha + beta));
% m = 1 + (1-omega_eq)/(omega_eq);
m = 1.2; % markup, m >= 1
gamma = 0.8; % paramter for money illusion
epsilon = 0.005; %numerical precision parameter

% Banking network parameters
N = 3; % size of network
a = [0.4, 0.35, 0.25]; % propotion of loans distributed to each bank. must sum to 1
b = [0.5, 0.3, 0.2]; % proportion of deposits distributed to each bank. must sum to 1

% Initial conditions
omega01 = 0.9;
lambda01 = 0.9;
l01 = 0.1;
Y01 = 100;
p01 = 1.0;
d01 = (1 - k_r)*l01;
T = 100;
run = 1; % number of loops
q = 0;

omega0 = omega01;
lambda0 = lambda01;
l0 = l01;
Y0 = Y01;
p0 = p01;
d0 = d01;

% array to store all information from model
output = zeros(3000*run, 11);

% set seed for reproducability
rng(123)

for i = 1:run
    
    % Dynamical System
    [t,y] = ode15s(@(t,y) keen(t,y), [0 T], [convert([omega0, lambda0, l0]), p0, d0], options);
    y = [retrieve(y(:,1:3)),y(:,4),y(:,5)]; % convert back to original vars
    s = size(y,1);
    ind = (q+1 : q+s); % index currently working over  
    output(ind,1) = t + (i-1)*T; % t 
    output(ind,2:6) = y(:,1:5); % omega (1), lambda (2), l (3), p (4), d (5)
    output(ind,7) = Y0*y(:,2)/lambda0.*exp((alpha+beta)*t); % Y_output 
    output(ind,8) = y(:,3).*output(ind,7).*y(:,4); % \Lambda loans
    output(ind,9) = 1-(y(:,5)./y(:,3)); % k_eff
    output(ind,10) = y(:,5).*output(ind,7).*y(:,4);% \Delta deposits
    output(ind,11) = output(ind,8) - output(ind,6); % Equity X_b
    %output(ind,11) = k_r.*y(:,4).*output(ind,7).*(fun_kappa(1-y(:,1)-r.*y(:,3))-1+y(:,1)+r.*y(:,3)); % S_b
    %output(ind,12) = r.*output(ind,8) - r_d.*output(ind,6) - output(ind,11); % Dividends
  
    % Distribute total Loand and Deposits among nodes
    a = randfixedsum(3,1,1,0,1); % propotion of loans distributed to each bank. must sum to 1
    b = a; % proportion of deposits distributed to each bank. must sum to 1
    externalAssets = (output(ind(end),8)*a)';
    externalLiabilities = (output(ind(end),10)*b)';
    dif  = externalAssets - externalLiabilities;    
    if sum(dif) < N
        disp('error: bank with negative equity being attributed with assets')
    end
    % Define the interbank Liabilities Matrix    
    interbankLiabilitiesMatrix = abs([0 (0.9 + 0.2*rand)*dif(1) 0; ...
                              0 0 (0.9 + 0.2*rand)*dif(2);... 
                              (0.9 + 0.2*rand)*dif(3) 0 0]);     
    
    % Initial shock to the system
    uniformShock = 0.1;
    shock = uniformShock*externalAssets;
    externalAssetsShocked  =  externalAssets - shock;
    
    % Eisenberg & Noe clearing payment function
    [equityLoss equity equityZero paymentVector error] = eisenbergnoe(interbankLiabilitiesMatrix,externalAssetsShocked,externalLiabilities);

%     % initial conditions for next run
%     omega0 = output(ind(end),2);
%     lambda0 = output(ind(end),3);
%     l0 = f_loans.*output(ind(end),4);
%     p0 = output(ind(end),5);
%     d0 = f_deposits.*output(ind(end),6);
%     Y0 = output(ind(end),7);
    
    q = q + s;
    q
end

figure(1)
plot(output(1:q,1),output(1:q,9)), hold on;
xlabel('time')
ylabel('Equity Ratio k (k_r = 0.08)')
title(['\omega_0 = ', num2str(omega01, txt_format), ...
    ', \lambda_0 = ', num2str(lambda01, txt_format), ...
    ', l_0 = ', num2str(l01, txt_format), ...
    ', p_0 = ', num2str(p01, txt_format), ...
    ', d_0 = ', num2str(d01, txt_format), ...
    ', Y_0 = ', num2str(Y01, txt_format)])

figure(2)
yyaxis right
plot(output(1:q,1),output(1:q,4),'r', output(1:q,1),output(1:q,6),'r-.')
ylabel('loan ratio l and deposit ratio d')

yyaxis left
plot(output(1:q,1),output(1:q,2),'b', output(1:q,1),output(1:q,3),'b-.')
ylabel('\omega, \lambda')

legend('\omega','\lambda', 'l', 'd')
title(['\omega_0 = ', num2str(omega01, txt_format), ...
    ', \lambda_0 = ', num2str(lambda01, txt_format), ...
    ', l_0 = ', num2str(l01, txt_format), ...
    ', p_0 = ', num2str(p01, txt_format), ...
    ', d_0 = ', num2str(d01, txt_format), ...
    ', Y_0 = ', num2str(Y01, txt_format)])


% figure(3)
% plot(output(1:q,1),output(1:q,2)), hold on;
% plot(output(1:q,1),output(1:q,3))
% xlabel('time t')
% ylabel('\omega, \lambda')
% legend('\omega', '\lambda')
% title(['\omega_0 = ', num2str(omega01, txt_format), ...
%     ', \lambda_0 = ', num2str(lambda01, txt_format), ...
%     ', l_0 = ', num2str(l01, txt_format), ...
%     ', p_0 = ', num2str(p01, txt_format), ...
%     ', d_0 = ', num2str(d01, txt_format), ...
%     ', Y_0 = ', num2str(Y01, txt_format)])

figure(4)
plot(output(1:q,1),eta_p*(m*output(1:q,2) - 1))
xlabel('time t')
ylabel('Inflation Rate')

% output(ind(end),2)
% output(ind(end),3)
% output(ind(end),4)
% output(ind(end),6)
% eta_p*(m*output(ind(end),2) - 1)
%% Auxiliary functions
% Instead of integrating the system in terms of omega, lambda and l, we
% decided to use:
%
% log_omega = log(omega), 
% tan_lambda = tan((lambda-0.5)*pi), and 
% pi_n=1-omega-r*l
% log_p = log(p)
%
% Experiments show that the numerical methods work better with
% these variables.

function new = convert(old)
    % This function converts
    % [omega, lambda, l, p]
    % to
    % [log_omega, tan_lambda, pi_n, log_p]
    % It also work with the first 2 or 3 elements alone.
    % each of these variables might also be a time-dependent vector
    
    n = size(old,2); %number of variables
    new = zeros(size(old));
    
    new(:,1) = log(old(:,1));
    new(:,2) = tan((old(:,2)-0.5)*pi);
    if n>2
        new(:,3) = 1-old(:,1)-r*old(:,3);
        if n==4
            new(:,4) = log(old(:,4));    
        end
    end
end

function old = retrieve(new)
    % This function retrieves [omega, lambda, l, p]
    % from
    % [log_omega, tan_lambda, pi_n, log_p]
    %  
    %  It also work with the first 2 or 3 elements alone.
    % each of these variables might also be a time-dependent vector
    
    n = size(new,2); %number of variables    
    old = zeros(size(new));
    
    old(:,1) = exp(new(:,1));
    old(:,2) = atan(new(:,2))/pi+0.5;
    if n>2
        old(:,3) = (1-old(:,1)-new(:,3))/r;
        if n==4
            old(:,4) = exp(new(:,4));
        end
    end
end

function f = keen(t,y)
   f = zeros(5,1);
   
   % variables
   log_omega = y(1); % omega
   tan_lambda = y(2); % lambda
   pi_n = y(3); % pi (profit) - l can be calculated from this
   p = y(4); % p
   d = y(5);
   
   % transformed vars
   lambda = atan(tan_lambda)/pi+0.5;
   omega = exp(log_omega);
   i_omega = eta_p*(m*omega - 1); % function i(omega) eqn (15) in Lipton

   indicator = 1 - r*d0/l0;
   k_n = 1 - (r*d)/(1-omega-pi_n);
   if ge(indicator,k_r - epsilon)
       fun_R = @(x) 1;
       g_Y = fun_kappa(pi_n)/nu-delta; % g(omega,l) from eqn (31) from Costa Lima
   else if (indicator < (k_r - epsilon)) && (ge(k_n,0))
           fun_R = @(x) x/k_r;
           g_Y = (k_n*(fun_kappa(pi_n)-pi_n)+pi_n)/(k_r*nu)-delta;
       else if(indicator < 0)
               fun_R = @(x) 0;
               g_Y = -delta;
           end
       end
   end
   
   % ODEs system
   f(1) = fun_Phi(lambda)-alpha-(1-gamma)*i_omega; %d(log_omega)/dt
   f(2) = (1+tan_lambda^2)*pi*lambda*(g_Y-alpha-beta); %d(tan_lambda)/dt
   f(3) = -omega*f(1) - r*fun_R(k_n)*(fun_kappa(pi_n)-pi_n) + (1-omega-pi_n)*(i_omega+g_Y); %d(pi_n)/dt
   f(4) = p*i_omega; % d(p)/dt
   f(5) = -d*(i_omega + g_Y) + (1-k_r)*fun_R(k_n)*(fun_kappa(pi_n) - pi_n); %d(d)/dt
end

end
