clear
close all;
clc

%Importing data from Excel(avoiding first 2 row)
A=readmatrix('Data.csv','NumHeaderLines',2); 
soc=A(:,1); % SOC data From Battery C1202
v=A(:,2); % OCV data of Battery C1202
eps = 0.175; 
zs = soc*(1-2*eps)+eps; % Scaled value of SOC in terms of zs

% Equation and graph for Unnewehr Universal model
p_u=[ones(length(soc),1) soc];
u=p_u';
k_u=inv(u*p_u)*(u*v); % Value of K
ocv_u= p_u*k_u; % calculating OCV for Unnewehr model

figure(1)
hold on
plot(soc, v,'b',soc, ocv_u,'g','linewidth',3 ); %plotting the graph
xlabel("SOC (%)");
ylabel("OCV (Volts)");
title("C1202")
legend("SOC-OCV data","Unnewehr universal model");

%Equation of Shepherd Model.
p_s=[ones(length(soc),1) 1./zs];
u=p_s';
k=inv(u*p_s)*(u*v);
ocv_s= p_s*k; % OCV value for Shepherd model.

plot(soc, v,'b',soc, ocv_s,'m','linewidth',3 );
xlabel("SOC (%)");
ylabel("OCV (Volts)");
title("C1202")
legend("SOC-OCV data","Shepherd model");

%Equation of Nernst Model
p_n=[ones(length(soc),1) log(zs) log(1-zs)];
u=p_n';
k=inv(u*p_n)*(u*v);
ocv_n= p_n*k; %OCV value of Nernst Model

plot(soc, v,'b',soc, ocv_n,'r','linewidth',3 );
xlabel("SOC (%)");
ylabel("OCV (Volts)");
title("C1202")
legend("SOC-OCV data","Nernst model");

%Equation of Combined Model
p_c=[ones(length(soc),1) 1./zs zs log(zs) log(1-zs) ];
u=p_c';
k=inv(u*p_c)*(u*v);
ocv_c= p_c*k; %OCV value of combined model
plot(soc, v,'b',soc, ocv_c,'g','linewidth',3 );
xlabel("SOC (%)");
ylabel("OCV (Volts)");
title("C1202")
legend("SOC-OCV data","Combined model");

%Equation of combined+3 model.
p_c3=[ones(length(soc),1) 1./zs 1./zs.^2 1./zs.^3 1./zs.^4 zs log(zs) log(1 -zs)];
u=p_c3';
k=inv(u*p_c3)*(u*v);
ocv_c3= p_c3*k; %OCV value of combines+3 Model
plot(soc, v,'b',soc, ocv_c3,'c','linewidth',3 );
xlabel("SOC (%)");
ylabel("OCV (Volts)");
title("C1202")
legend("SOC-OCV data","Combined+3 model");

%Equation for Polynomial value
p_poly=[ones(length(soc),1) zs zs.^2 zs.^3 zs.^-1 zs.^-2];
u=p_poly';
k=inv(u*p_poly)*(u*v); %value of K for polynomial model
ocv_poly= p_poly*k; %OCV value of polynomial model
plot(soc, v,'b',soc, ocv_poly,'y','linewidth',3 );
xlabel("SOC (%)");
ylabel("OCV (Volts)");
title("C1202")
legend("SOC-OCV data","Polynomial model");

%equation for exponential value
p_expo=[ones(length(soc),1) exp(zs) exp(zs.^2) exp(-zs)];
u=p_expo';
k=inv(u*p_expo)*(u*v); %value of K for exponential value
ocv_expo= p_expo*k; %value of OCV for exponential model
plot(soc, v,'b',soc, ocv_expo,'k','linewidth',3 );
xlabel("SOC (%)");
ylabel("OCV (Volts)");
title("C1202")
legend("SOC-OCV data","Exponential model");


title("7 models of C1205")
legend("SOC-OCV data","Unnewehr model","Shepherd Model","Nernst Model","Combined Model","Combined+3 Model","Polynomial Model","Exponential model");
hold off

% calculation for Error Metrics for Unnewehr Model
n_u = numel(ocv_u); 
Vbar_u = (1/n_u)*(norm(v));
BF_u = (1-((norm(ocv_u-v))/(norm(v-Vbar_u))))*100 %bestfit value
R_squared_u = (1 - ((norm(ocv_u-v)).^2/(norm(v-Vbar_u).^2)))*100 %R square value
Merc_u= v-ocv_u;
Maxerror_u = max(Merc_u) %Max error value
error_u = immse(v,ocv_u); 
rms_u = sqrt(error_u)

%Model evalution metrics
z_u = sum(Merc_u.^2);
M_u = numel(k_u);
AIC_u = n_u*log(z_u/n_u) + 2*(M_u + 1)

h=figure(2);hold on 
plot(soc,Merc_u,'g','linewidth',2);
xlabel("SOC (%)");
ylabel("Error (Volts)");
title("Error Graph for Unnewehr Universal Model")
legend("Error Graph");

% calculation for Error Metrics for Shepherd Model
n_s = numel(ocv_s);
Vbar_s = (1/n_s)*(norm(v));  % Value of Vbar
BF_s = (1-((norm(ocv_s-v))/(norm(v-Vbar_s))))*100 %bestfit value
R_squared_s = (1 - ((norm(ocv_s-v)).^2/(norm(v-Vbar_s).^2)))*100
Merc_s= v-ocv_s;
Maxerror_s = max(Merc_s) % Max error Value
error_s = immse(v,ocv_s);
rms_s = sqrt(error_s)

%Model evalution metrics
z_s = sum(Merc_s.^2);
M_s = numel(k);
AIC_s = n_s*log(z_s/n_s) + 2*(M_s + 1)

plot(soc,Merc_s,'m','linewidth',2);
xlabel("SOC (%)");
ylabel("Error (Volts)");
title("Error Graph for Shepherd Model")
legend("Error Graph");

% calculation for Error Metrics for Nernst Model
n_n = numel(ocv_n);
Vbar_n = (1/n_n)*(norm(v)); % value of Vbar
BF_n = (1-((norm(ocv_n-v))/(norm(v-Vbar_n))))*100 %best fit value
R_squared_n = (1 - ((norm(ocv_n-v)).^2/(norm(v-Vbar_n).^2)))*100
Merc_n= v-ocv_n;
Maxerror_n = max(Merc_n) % Maximum Error
error_n = immse(v,ocv_n);
rms_n = sqrt(error_n)

%Model evalution metrics
z_n = sum(Merc_n.^2);
M_n = numel(k);
AIC_n = n_n*log(z_n/n_n) + 2*(M_n + 1)

plot(soc,Merc_n,'r','linewidth',2);
xlabel("SOC (%)");
ylabel("Error (Volts)");
title("Error Graph for Nernst Model")
legend("Error Graph");

% calculation for Error Metrics for Combined Model
n_c = numel(ocv_c);
Vbar_c = (1/n_c)*(norm(v));  % Value of V bar
BF_c = (1-((norm(ocv_c-v))/(norm(v-Vbar_c))))*100 % Best fit value
R_squared_c = (1 - ((norm(ocv_c-v)).^2/(norm(v-Vbar_c).^2)))*100
Merc_c= v-ocv_c;
Maxerr_c = max(Merc_c) % Maximum error
error_c = immse(v,ocv_c);
rms_c = sqrt(error_c)

% Model Evaluation metrics
z_c = sum(Merc_c.^2);
M_c = numel(k);
AIC_c = n_c*log(z_c/n_c) + 2*(M_c + 1)

plot(soc,Merc_c,'c','linewidth',2)
xlabel("SOC (%)");
ylabel("Error (Volts)");
title("Error Graph for Combined Model")
legend("Error Graph");

% calculation for Error Metrics for Combined+3 Model
n_c3 = numel(ocv_c3);
Vbar_c3 = (1/n_c3)*(norm(v)); % Value of V bar
BF_c3 = (1-((norm(ocv_c3-v))/(norm(v-Vbar_c3))))*100 % Best fit value
R_squared_c3 = (1 - ((norm(ocv_c3-v)).^2/(norm(v-Vbar_c3).^2)))*100
Merc_c3= v-ocv_c3;
Maxerror_c3 = max(Merc_c3) % Max error
error_c3 = immse(v,ocv_c3);
rms_c3 = sqrt(error_c3)

% Model Evaluation metrics
z_c3 = sum(Merc_c3.^2);
M_c3 = numel(k);
AIC_c3 = n_c3*log(z_c3/n_c3) + 2*(M_c3 + 1)

plot(soc,Merc_c3,'k','linewidth',2);
xlabel("SOC (%)");
ylabel("Error (Volts)");
title("Error Graph for Combined+3 Model")
legend("Error Graph");

% calculation for Error Metrics for Polynomial Model
n_poly = numel(ocv_poly);
Vbar_poly = (1/n_poly)*(norm(v)); % V bar value
BF_poly = (1-((norm(ocv_poly-v))/(norm(v-Vbar_poly))))*100 % best fit value 
R_squared_poly = (1 - ((norm(ocv_poly-v)).^2/(norm(v-Vbar_poly).^2)))*100
Merc_poly= v-ocv_poly;
Maxerror_poly = max(Merc_poly)  % Maximum error
error_poly = immse(v,ocv_poly);
rms_poly = sqrt(error_poly)

% Model Evaluation metrics
z_poly = sum(Merc_poly.^2);
M_poly = numel(k);
AIC_poly = n_poly*log(z_poly/n_poly) + 2*(M_poly + 1)

plot(soc,Merc_poly,'b','linewidth',2);
xlabel("SOC (%)");
ylabel("Error (Volts)");
title("Error Graph for polynomial Model")
legend("Error Graph");

% calculation for Error Metrics for Exponential Model
n_expo = numel(ocv_expo);
Vbar_expo = (1/n_expo)*(norm(v));  % V bar value
BF_expo = (1-((norm(ocv_expo-v))/(norm(v-Vbar_expo))))*100 % best fit value
R_squared_expo = (1 - ((norm(ocv_expo-v)).^2/(norm(v-Vbar_expo).^2)))*100
Merc_expo= v-ocv_expo;
Maxerror_expo = max(Merc_expo) % Maximum Error
error_expo = immse(v,ocv_expo);
rms_expo = sqrt(error_expo)

% Model Evaluation metrics
z_expo = sum(Merc_expo.^2);
M_expo = numel(k);
AIC_expo = n_expo*log(z_expo/n_expo) + 2*(M_expo + 1)

plot(soc,Merc_expo,'y','linewidth',2);
xlabel("SOC (%)");
ylabel("Error (Volts)");
title("Error Graph for Exponential Model")
legend("Error Graph");

title("Error Graph")
xlabel("State of Charge (%)");
ylabel("OCV (Volts)");
legend("Unnweher Model","Shepherd Model","Nernst Model","Combined Model","Combined+3 Model","Polynomial Model","Exponential Model")
hold off