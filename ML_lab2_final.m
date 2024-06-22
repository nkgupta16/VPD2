clear all

voltages = [10,15,20,25,30,35,40,45,50];

L = 0.0047; % inductance
i = 48; % gear ratio of the motor
m = 0.016;  % rotor mass
r = 1.15/100;  % rotor radius
Jed = (m*r^2)/2;  % moment of inertia
J = i^2*Jed % reduced moment of inertia

files_n = ["data-10", "data-15","data-20", "data-25", "data-30", "data-35", "data-40", "data-45", "data-50"];
files_P = ["data10", "data15","data20", "data25", "data30", "data35", "data40", "data45", "data50"];

voltage_n = [-0.39,-0.58,-0.68,-0.82,-0.95,-1.09,-1.24,-1.40,-1.53];
voltage_p = [0.38,0.59,0.75,0.85,0.94,1.11,1.24,1.41,1.53];
current_n = [-0.06,-0.1,-0.13,-0.17,-0.21,-0.24,-0.28,-0.31,-0.34];
current_p = [0.06,0.1,0.14,0.17,0.20,0.24,0.27,0.31,0.34];

% Graph current v/s voltage (POSITIVE)
figure(1)
fun = @(par, current_p)par(1)*current_p;
par0 = 0;
par  = lsqcurvefit(fun, par0, current_p,voltage_p);
resistance_p = par(1);
voltage_p_aprx = current_p*resistance_p;
set(gcf, "Position", [100,100,800,600]);
plot(current_p,voltage_p,current_p,voltage_p_aprx,"LineWidth",2);
xlabel("I, A");
ylabel("U, B");
grid on;
set(gca,"GridAlpha",0.7);
fontsize(gcf,20,"points");
title("U(I) Positive Voltage");
legend("Measurement","Approximation","Location","best");

% Graph current v/s voltage (NEGATIVE)
figure(2)
fun = @(par, current_n)par(1)*current_n;
par0 = 0;
par  = lsqcurvefit(fun, par0, current_n,voltage_n);
resistance_n = par(1);
voltage_n_aprx = current_n*resistance_n;
set(gcf, "Position", [100,100,800,600]);
plot(current_n,voltage_n,current_n,voltage_n_aprx,"LineWidth",2);
xlabel("I, A");
ylabel("U, B");
grid on;
set(gca,"GridAlpha",0.7);
fontsize(gcf,20,"points");
title("U(I) Negative Voltage");
legend("Measurement","Approximation","Location","best");

R = (resistance_n + resistance_p)/2;

% T_m calculation
T_m_p = []

for i=1:9
    data = readmatrix(files_P{i}+".csv");
    time =  data(:,1);
    angle = data(:,2)*pi/180;
    U_pr = voltages(i);
    par0 = [0.1;0.000157];
    fun = @(par,time)U_pr*par(1)*(time-par(2)*(1-exp(-time/par(2))));
    par = lsqcurvefit(fun,par0,time,angle);
    k = par(1);
    T_m_p(i) = par(2)
    ks(i,:) = k;
    w_nls_P = U_pr*k;
    w_all_positive(i,:) = w_nls_P;
end

T_m_n = []

for i=1:9
    data = readmatrix(files_n{i}+".csv");
    time =  data(:,1);
    angle = data(:,2)*pi/180;
    U_pr = voltages(i);
    par0 = [0.1;0.000157];
    fun = @(par,time)U_pr*par(1)*(time-par(2)*(1-exp(-time/par(2))));
    par = lsqcurvefit(fun,par0,time,angle);
    k = par(1);
    T_m_n(i) = par(2)
    ks(i,:) = k;
    w_nls_P = U_pr*k;
    w_all_positive(i,:) = w_nls_P;
end



% Angular Speed Calculation
omega_p_sr = [];
speed = cell(1, 9);
time = cell(1, 9);

for i=1:9
    data = readmatrix(files_P(i)+".csv");
    omega = data(:,3)*pi/180;
    delt_omega = omega(end-19:end);
    omega_p_sr(i) = sum(delt_omega)/length(delt_omega);
    speed{i} = omega;
    time{i} = data(:,1);
end

for i=1:9
    data = readmatrix(files_n(i)+".csv");
    omega = data(:,3)*pi/180;
    delt_omega = omega(end-19:end);
    omega_n_sr(i) = sum(delt_omega)/length(delt_omega);
    speed{i+9} = omega;
    time{i+9} = data(:,1);
end

% Calculation of k_e
fun = @(par,omega_p_sr)par(1)*omega_p_sr;
par = lsqcurvefit(fun, par0, omega_p_sr,voltage_p);
k_e_p = par(1);
figure(3)
set(gcf,"Position",[100,100,800,600])
plot(omega_p_sr,voltage_p,omega_p_sr,omega_p_sr.*k_e_p,"LineWidth",2)
xlabel("Angular Speed, rad/sec")
ylabel("U, Volts")
grid on;
set(gca, "GridAlpha",0.7)
set(gca, "LineWidth", 1.1)
fontsize(gcf, 20, "points")
title("U(w) Positive Voltage")
legend("Measurement","Approximation","Location","best")

fun = @(par,omega_n_sr)par(1)*omega_n_sr;
par = lsqcurvefit(fun, par0, omega_n_sr,voltage_n);
k_e_n = par(1);
figure(4)
set(gcf,"Position",[100,100,800,600])
plot(omega_n_sr,voltage_n,omega_n_sr,omega_n_sr.*k_e_n,"LineWidth",2)
xlabel("Angular Speed, rad/sec")
ylabel("U, Volts")
grid on;
set(gca, "GridAlpha",0.7)
set(gca, "LineWidth", 1.1)
fontsize(gcf, 20, "points")
title("U(w) Negative Voltage")
legend("Measurement","Approximation","Location","best")

k_e = (k_e_n + k_e_p)/2;
k_m = k_e;

% Positive data simulation
figure(5)
set(gcf,"Position",[100,100,1000,600])
hold on; grid on;
set(gca,"GridAlpha",0.7);
set(gca, "LineWidth", 1.1);
fontsize(gcf,20,"points")
ylabel("Current (I), A")
xlabel("Time (t), s")
title("Simulation of Voltage")

merged_figure = figure(6);
set(merged_figure, "Position", [100,100,1200,600])
hold on; grid on;
set(gca, "GridAlpha", 0.7);
set(gca, "LineWidth", 1.1);
fontsize(merged_figure, 20, 'points')
ylabel("Angular Speed (w), rad/s")
xlabel("Time (t), s")
title("Measurement and Simulation of Voltages")

positive_colors = lines(9); 
negative_colors = parula(9); 

for i = 1:9
    U = voltage_p(i);
    simulin_data = sim("simulinkLAB2.slx");

    figure(5)
    plot(simulin_data.I.Time, simulin_data.I.Data, 'LineWidth', 2, "DisplayName", "Simulation "+ voltages(i) + "%","Color", positive_colors(i,:))

    figure(6)
    plot(simulin_data.w.Time, simulin_data.w.Data, "LineWidth", 2, "DisplayName","Simulation "+ voltages(i) + "%", "Color", positive_colors(i,:))
    plot(time{i}, speed{i}, "LineWidth", 2, "DisplayName","Measurement "+ voltages(i) + "%", "Color", positive_colors(i,:))
end

for i = 1:9
    U = voltage_n(i);
    simulin_data = sim("simulinkLAB2.slx");

    figure(5);
    plot(simulin_data.I.Time, simulin_data.I.Data, "LineWidth", 2, "DisplayName", "Simulation "+ -voltages(i) + "%","Color", negative_colors(i,:));
    
    figure(6)
    plot(simulin_data.w.Time, simulin_data.w.Data, "lineWidth", 2, "DisplayName", "Simulation "+ -voltages(i) + "%", "Color", negative_colors(i,:))
    plot(time{i+9}, speed{i+9}, "LineWidth", 2, "DisplayName", "Measurement "+ -voltages(i) + "%", "Color", negative_colors(i,:));
end

figure(5)
title("Simulation of Voltage")
legend("Location","eastoutside");
xlim([0, 4]);
xticks(0:0.5:4);

figure(6)
title("Measurement and Simulation of Voltages")
legend("NumColumns", 2, "Location", "eastoutside");
xlim([0, 4]);
xticks(0:0.5:4);
