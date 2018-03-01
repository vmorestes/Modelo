// Dynare code for the model with complete information
// Under complete information the model collpse to a
// representative agents model. This is a preliminary version

// Things to do:
// 1. Model the exogenous process for foreign variables as a SVAR(here each one is an AR(1))
// 2. Implement exchange rate intervation on future market(this may imply on modifications in the model)
// 3. Include on the Taylor Rule the output gap


var 
R //Nominal Gross Return
c //Final Consumption
p //Consumer Price
n //labor supply
w //nominal wage
PI //  inflation domestic prices
PIstar // Diference between optmal and current price
pstar //optimal price setting 
g1 // auxiliar 1
g2 // auxiliar 2
mc // marginal cost
a // productivity
pfstar //foreign price
delta // price dispersion measure
s // nominal exchange rate
yd; // domestic production

varexo emp ea  epstar;

parameters 
beta //  discount factor
sigma //  parameter of utility function: consumption
phi //  (inverse) Frisch elasticity of labor supply 
neta //  elasticity of substitution between domestic and imported consumption goods
gamma //  1- gamma is the home bias parameter
lambdad //  elasticity of substitution between intermediate exported goods
alpha // capital share in production
theta // probability of not being able to change price(Calvo parameter)
lambdastar //  parameter of foreign demand for exported goods(its like an elasticity of substitution)
taux //hedge percentage of exporters 
lambdax // elasticity of substitution for exported intermediate goods
lambdaf // elasticity of substitution for imported intermediate goods
tauf // hedge percentage of importers
phir // Taylor Rule Parameter: interest rate smoothing 
Rbar // Steady state interest rate
PIbar // (gross) inflation target
phipi // Taylor Rule Parameter: response to inflation
phoa // AR(1) parameter: technology shock
phoy // AR(1) parameter: foreign output shock
phop // AR(1) parameter: foreign product shock
phob // AR(1) parameter: exchange rate intervention shock
phof
fxbar 
phof2
phor; //  AR(1) parameter: foreign (gross) interest rate shock


beta = 0.99;
sigma = 1.5;
phi = 3;
neta = 1.5;
gamma = 0.4;
lambdad = 6;
alpha = 1/3;
theta = 0.75;
lambdastar = 2;
taux  = 0;
lambdax  = 6;
lambdaf = 6; 
tauf  = 0;
phir = 0;
Rbar = 1/beta;
PIbar = 1;
phipi = 1.5;
phoa = 0.6;
phoy = 0.9;
phop = 0.6;
phor = 0.8;
phob  = 0.5;
phof = 0.8;
fxbar = -0.05;
phof2 = 0.8;
model;
// CONSUMER
//. 1:EULER EQUATION
1/exp(R) = beta*((exp(c(+1))/exp(c))^(-sigma) * (exp(p)/exp(p(+1))));
//.  2:LABOR SUPPLY
exp(n)^phi * exp(c) ^sigma = (exp(w)/exp(p)) ;


// Domestic goods production
//. 5-6: Definition of (gross) inflation
exp(PI) = (exp(p)/exp(p(-1)));
exp(PIstar) = exp(pstar)/exp(p);
//. 7-9: Calvo Pricing
exp(PIstar)^(1 + (alpha*lambdad)/(1 - alpha)) *exp(g2) = lambdad/(lambdad - 1) * exp(g1);
exp(g1) = exp(c)^(-sigma) * exp(mc) * exp(yd) + beta*theta* exp(PI(+1))^((lambdad)/(1 - alpha))* exp(g1(+1));
exp(g2) = exp(c)^(-sigma) * exp(yd) + beta*theta*exp(PI(+1))^(lambdad - 1)* exp(g2(+1));
//. 10: Law of Motion for domestic Prices
exp(PI)^(1- lambdad) = theta + (1-theta)*(exp(pstar)/exp(p(-1)))^(1- lambdad);
//. 11: marginal cost
exp(mc) = exp(w)/((1 - alpha)* exp(a) * exp(n)^(-alpha)*exp(delta)^alpha);
//. 12: Price dispersion
exp(delta) = (1-theta)*exp(PIstar)^(-lambdad/(1 - alpha)) + theta* exp(PI)^(lambdad/(1 - alpha))* exp(delta(-1));


//PPP
exp(p) = exp(s)*exp(pfstar);

//Central Bank 
//Taylor Rule
exp(R)/Rbar = (exp(R(-1))/Rbar)^phir*(exp(PI)/PIbar)^phipi*exp(emp);


// Equilibrium Conditions
//. 20 Domestic goods demand
exp(yd) = exp(c) ;
//. 21 Domestic goods supply
exp(yd) = exp(n)^(1-alpha)*exp(a)*exp(delta)^(alpha - 1);
//. 23 BP = 0



//exogenous processes
//. 24 Domestic Productivity
a=phoa*a(-1)+ea;
//. 26 International Price
pfstar=phop*pfstar(-1)+epstar;


end;
initval;
R = ln(1/beta);
end;


shocks; 
var  emp; stderr 0.01;
var ea; stderr 0.01;
var epstar; stderr 0.01;
end;

steady;

stoch_simul(order=1,irf=10) yd R s PI;

