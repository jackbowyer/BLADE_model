%                          
% Weinberg et al. 2017. BLADE platform.

% Experimental data and adapted angular metric scores
Data = [67 76 60 2 6 13 77 0    2.04;
        63 73 1  7 4 88 0  92   2.22].';

% Circuits to be simulated
GFPmCherry = 4;
mCherry = 3;
GFP = 2;
termination = 1;

Vin = [GFPmCherry  GFP         termination GFP       ;
       GFPmCherry  termination mCherry     mCherry].';


% Ideal truth table for each circuit
A = [0 0 1 0 0 1 1 1].';
IdealTT = zeros(8,2);       
for h = 1:2
IdealTT(:,h) = [A((Vin(1,h)*2)-1) A(Vin(1,h)*2)...
                A((Vin(2,h)*2)-1) A(Vin(2,h)*2)...
                A((Vin(3,h)*2)-1) A(Vin(3,h)*2)...
                A((Vin(4,h)*2)-1) A(Vin(4,h)*2)].';
end


%%%%%%%%%%%%%%%%%%%%%% GILLESPIE ALGORITHM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for m = 1:4 
   
c1 = [0 1 0 1].'; 
c2 = [0 0 1 1].';

M = 17;                     % Number of reaction pathways
N = 10;                     % Number of molecular species considered

% STEP 0
nmols = 1000;
% c = [a1/B a2/A a3/B].';   % In general, a vector of length M of stochastic reaction constants.
Y0 = zeros(N,1);            % Initial number of molecules of species X. In general, a vector of length N of initial population numbers for each species (X1,X2,...,XN).
Y0(3) = nmols;
Yplot = (Y0./nmols)*100;    % For plotting.
t = 0;                      % Time variable.
tplot = zeros(1);           % For plotting.
n = 0;                      % Reaction counter.
n_max = 10000;
% rng(2)                      % specify a random number generator for reproducible results

% Optimal parameter set from the Genetic Algorithm
OptPars = exp([-0.463464805119107  -6.380196305530780 -10.078440710821784 -10.542821382522401  -9.415239471967475...
-5.711702227971008  -8.703783350567580]).';

alpha_Cre = c1(m)*OptPars(1);
alpha_Flp = c2(m)*OptPars(1);
beta_p = OptPars(2);
k1c = OptPars(3);
km1c = OptPars(5);
k2c = OptPars(3);
km2c = OptPars(5);
k1f = OptPars(4);
km1f = OptPars(5);
k2f = OptPars(4);
km2f = OptPars(5);
delta_X = OptPars(6);
delta_D = OptPars(7);


% STEP 1
while n < n_max
%     h = [].';                % Number of molecular reactant combinations available in current state. In general, h is a vector of length M consisting of combinatorial functions of molecular population numbers X (which is of length N).
    a = [alpha_Cre beta_p*Y0(1) alpha_Flp beta_p*Y0(2) k1c*Y0(1)*Y0(3) km1c*Y0(4)*Y0(5) k1f*Y0(2)*Y0(3) km1f*Y0(6)*Y0(7) k2c*Y0(1)*Y0(6) km2c*Y0(8)*Y0(10) k2f*Y0(2)*Y0(4) km2f*Y0(8)*Y0(9) delta_D*Y0(3) delta_X*Y0(5) delta_X*Y0(7) delta_X*Y0(9) delta_X*Y0(10)].';       % a is the propensity of the reaction pathway in current state. In general, a vector of length M   kX*Y0(3) kX*Y0(4) kX*Y0(5) kX*Y0(6)
   
    a0 = sum(a);               % a0 is total propensity that anything happens. This number emerges more out of mathematical necessity than physical intuition.
    
    % STEP 2
    r = rand(2,1);
    tau = -log(r(1))/a0;
    csa = cumsum(a);
    check = sum([[0; csa(1:M-1)] < r(2)*a0 r(2)*a0 <= csa].');  
    mu = find(check==max(check));
    
    % STEP 3
    t = t + tau;
    % Adjust population levels based on reaction formula(e)
    switch round(mean(mu))
        case 1
            Y0(1) = Y0(1) + 1;
            Y0(2) = Y0(2);
            Y0(3) = Y0(3);
            Y0(4) = Y0(4);
            Y0(5) = Y0(5);
            Y0(6) = Y0(6);
            Y0(7) = Y0(7);
            Y0(8) = Y0(8);
            Y0(9) = Y0(9);
            Y0(10) = Y0(10);
        case 2
            Y0(1) = Y0(1) - 1;
            Y0(2) = Y0(2);
            Y0(3) = Y0(3);
            Y0(4) = Y0(4);
            Y0(5) = Y0(5);
            Y0(6) = Y0(6);
            Y0(7) = Y0(7);
            Y0(8) = Y0(8);
            Y0(9) = Y0(9);
            Y0(10) = Y0(10);
        case 3
            Y0(1) = Y0(1);
            Y0(2) = Y0(2) + 1;
            Y0(3) = Y0(3);
            Y0(4) = Y0(4);
            Y0(5) = Y0(5);
            Y0(6) = Y0(6);
            Y0(7) = Y0(7);
            Y0(8) = Y0(8);
            Y0(9) = Y0(9);
            Y0(10) = Y0(10);
        case 4
            Y0(1) = Y0(1);
            Y0(2) = Y0(2) - 1;
            Y0(3) = Y0(3);
            Y0(4) = Y0(4);
            Y0(5) = Y0(5);
            Y0(6) = Y0(6);
            Y0(7) = Y0(7);
            Y0(8) = Y0(8);
            Y0(9) = Y0(9);
            Y0(10) = Y0(10);
        case 5
            Y0(1) = Y0(1) - 2;
            Y0(2) = Y0(2);
            Y0(3) = Y0(3) - 1;
            Y0(4) = Y0(4) + 0.75;
            Y0(5) = Y0(5) + 0.25;
            Y0(6) = Y0(6);
            Y0(7) = Y0(7);
            Y0(8) = Y0(8);
            Y0(9) = Y0(9);
            Y0(10) = Y0(10);
        case 6
            Y0(1) = Y0(1) + 2;
            Y0(2) = Y0(2);
            Y0(3) = Y0(3) + 1;
            Y0(4) = Y0(4) - 0.75;
            Y0(5) = Y0(5) - 0.25;
            Y0(6) = Y0(6);
            Y0(7) = Y0(7);
            Y0(8) = Y0(8);
            Y0(9) = Y0(9);
            Y0(10) = Y0(10);
        case 7
            Y0(1) = Y0(1);
            Y0(2) = Y0(2) - 1;
            Y0(3) = Y0(3) - 1;
            Y0(4) = Y0(4);
            Y0(5) = Y0(5);
            Y0(6) = Y0(6) + 0.75;
            Y0(7) = Y0(7) + 0.25;
            Y0(8) = Y0(8);
            Y0(9) = Y0(9);
            Y0(10) = Y0(10);
        case 8
            Y0(1) = Y0(1);
            Y0(2) = Y0(2) + 1;
            Y0(3) = Y0(3) + 1;
            Y0(4) = Y0(4);
            Y0(5) = Y0(5);
            Y0(6) = Y0(6) - 0.75;
            Y0(7) = Y0(7) - 0.25;
            Y0(8) = Y0(8);
            Y0(9) = Y0(9);
            Y0(10) = Y0(10);
        case 9
            Y0(1) = Y0(1) - 1;
            Y0(2) = Y0(2);
            Y0(3) = Y0(3);
            Y0(4) = Y0(4);
            Y0(5) = Y0(5);
            Y0(6) = Y0(6) - 1;
            Y0(7) = Y0(7);
            Y0(8) = Y0(8) + 0.75;
            Y0(9) = Y0(9);
            Y0(10) = Y0(10) + 0.25;
        case 10
            Y0(1) = Y0(1) + 1;
            Y0(2) = Y0(2);
            Y0(3) = Y0(3);
            Y0(4) = Y0(4);
            Y0(5) = Y0(5);
            Y0(6) = Y0(6) + 1;
            Y0(7) = Y0(7);
            Y0(8) = Y0(8) - 0.75;
            Y0(9) = Y0(9);
            Y0(10) = Y0(10) - 0.25;
        case 11
            Y0(1) = Y0(1);
            Y0(2) = Y0(2) - 1;
            Y0(3) = Y0(3);
            Y0(4) = Y0(4) - 1;
            Y0(5) = Y0(5);
            Y0(6) = Y0(6);
            Y0(7) = Y0(7);
            Y0(8) = Y0(8) + 0.75;
            Y0(9) = Y0(9) + 0.25;
            Y0(10) = Y0(10);
        case 12
            Y0(1) = Y0(1);
            Y0(2) = Y0(2) + 1;
            Y0(3) = Y0(3);
            Y0(4) = Y0(4) + 1;
            Y0(5) = Y0(5);
            Y0(6) = Y0(6);
            Y0(7) = Y0(7);
            Y0(8) = Y0(8) - 0.75;
            Y0(9) = Y0(9) - 0.25;
            Y0(10) = Y0(10);
        case 13
            Y0(1) = Y0(1);
            Y0(2) = Y0(2);
            Y0(3) = Y0(3) - 1;
            Y0(4) = Y0(4);
            Y0(5) = Y0(5);
            Y0(6) = Y0(6);
            Y0(7) = Y0(7);
            Y0(8) = Y0(8);
            Y0(9) = Y0(9);
            Y0(10) = Y0(10);
        case 14
            Y0(1) = Y0(1);
            Y0(2) = Y0(2);
            Y0(3) = Y0(3);
            Y0(4) = Y0(4);
            Y0(5) = Y0(5) - 1;
            Y0(6) = Y0(6);
            Y0(7) = Y0(7);
            Y0(8) = Y0(8);
            Y0(9) = Y0(9);
            Y0(10) = Y0(10);
        case 15
            Y0(1) = Y0(1);
            Y0(2) = Y0(2);
            Y0(3) = Y0(3);
            Y0(4) = Y0(4);
            Y0(5) = Y0(5);
            Y0(6) = Y0(6);
            Y0(7) = Y0(7) - 1;
            Y0(8) = Y0(8);
            Y0(9) = Y0(9);
            Y0(10) = Y0(10);
        case 16
            Y0(1) = Y0(1);
            Y0(2) = Y0(2);
            Y0(3) = Y0(3);
            Y0(4) = Y0(4);
            Y0(5) = Y0(5);
            Y0(6) = Y0(6);
            Y0(7) = Y0(7);
            Y0(8) = Y0(8);
            Y0(9) = Y0(9) - 1;
            Y0(10) = Y0(10);
        case 17
            Y0(1) = Y0(1);
            Y0(2) = Y0(2);
            Y0(3) = Y0(3);
            Y0(4) = Y0(4);
            Y0(5) = Y0(5);
            Y0(6) = Y0(6);
            Y0(7) = Y0(7);
            Y0(8) = Y0(8);
            Y0(9) = Y0(9);
            Y0(10) = Y0(10) - 1;
    end
    n = n + 1;
    % At this point, all the physics has been simulated, it only remains to
    % record the new values of t and X to vectors to we can plot it later.   
    Yplot(:,n+1) = (Y0./nmols)*100;
    tplot(n+1) = t;  
    
end


G12 = Yplot(3,:)*IdealTT(1,1)+Yplot(4,:)*IdealTT(3,1)+Yplot(6,:)*IdealTT(5,1)+Yplot(8,:)*IdealTT(7,1);        % GFP
R12 = Yplot(3,:)*IdealTT(2,1)+Yplot(4,:)*IdealTT(4,1)+Yplot(6,:)*IdealTT(6,1)+Yplot(8,:)*IdealTT(8,1);        % mCherry
G14 = Yplot(3,:)*IdealTT(1,2)+Yplot(4,:)*IdealTT(3,2)+Yplot(6,:)*IdealTT(5,2)+Yplot(8,:)*IdealTT(7,2);        % GFP
R14 = Yplot(3,:)*IdealTT(2,2)+Yplot(4,:)*IdealTT(4,2)+Yplot(6,:)*IdealTT(6,2)+Yplot(8,:)*IdealTT(8,2);        % mCherry

if m == 1
    Sim1 = G12;
    Sim2 = R12;
    Sim9 = G14;
    Sim10 = R14;
    time1 = tplot;
elseif m == 2
    Sim3 = G12;
    Sim4 = R12;
    Sim11 = G14;
    Sim12 = R14;
    time2 = tplot;
elseif m == 3
    Sim5 = G12;
    Sim6 = R12;
    Sim13 = G14;
    Sim14 = R14;
    time3 = tplot;
else
    Sim7 = G12;
    Sim8 = R12;
    Sim15 = G14;
    Sim16 = R14;
    time4 = tplot;
end



% figure(10)
% subplot(2,2,m)
% stairs(tplot,Yplot(1,:),'m--')
% hold on
% stairs(tplot,Yplot(2,:),'m')
% stairs(tplot,Yplot(3,:),'b')
% stairs(tplot,Yplot(4,:),'r')
% stairs(tplot,Yplot(5,:),'c--')
% stairs(tplot,Yplot(6,:),'g')
% stairs(tplot,Yplot(7,:),'c--')
% stairs(tplot,Yplot(8,:),'r')
% stairs(tplot,Yplot(8,:),'g--')
% stairs(tplot,Yplot(9,:),'c--')
% stairs(tplot,Yplot(10,:),'c--')

figure(20)
subplot(2,2,m)
stairs(tplot,G12,'g')
hold on
stairs(tplot,R12,'r--')
% stairs(tplot,Yplot(1,:),'m')
% stairs(tplot,Yplot(2,:),'m')

figure(21)
subplot(2,2,m)
stairs(tplot,G14,'g')
hold on
stairs(tplot,R14,'r--')
% stairs(tplot,Yplot(1,:),'m')
% stairs(tplot,Yplot(2,:),'m')

% figure(30)
% subplot(4,4,m+8)
% stairs(tplot,G12,'g')
% hold on
% stairs(tplot,R12,'r--')
% 
% figure(31)
% subplot(4,4,m+12)
% stairs(tplot,G14,'g')
% hold on
% stairs(tplot,R14,'r--')

end


T = min([time1(end) time2(end) time3(end) time4(end)]);
T1 = find(abs(time1-T)==min(abs(time1-T)));
T2 = find(abs(time2-T)==min(abs(time2-T)));
T3 = find(abs(time3-T)==min(abs(time3-T)));
T4 = find(abs(time4-T)==min(abs(time4-T)));

SimTT = [Sim1(T1) Sim2(T1) Sim3(T2) Sim4(T2) Sim5(T3) Sim6(T3) Sim7(T4) Sim8(T4)...
         Sim9(T1) Sim10(T1) Sim11(T2) Sim12(T2) Sim13(T3) Sim14(T3) Sim15(T4) Sim16(T4)];

   
theta = [acosd((dot(IdealTT(:,1),abs(SimTT(1:8))))/(sqrt(sum(IdealTT(:,1).^2)*sum(SimTT(1:8).^2))))/nnz(IdealTT(:,1))...
         acosd((dot(IdealTT(:,2),abs(SimTT(9:16))))/(sqrt(sum(IdealTT(:,2).^2)*sum(SimTT(9:16).^2))))/nnz(IdealTT(:,2))];
disp([theta, sum(theta)])
disp(abs(Data(end,1)-theta(1)) + abs(Data(end,2)-theta(2)))
