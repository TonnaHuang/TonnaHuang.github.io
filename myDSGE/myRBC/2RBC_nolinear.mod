%RBC_nolinear.mod

var y c k inv w l r_k r u_a u_g;

varexo es_a es_g;

parameters beta alpha sigma chi varphi delta CY G rho_a rho_g sigma_a sigma_c;

%----------------------------------------------------------------
%\\parameters value
%----------------------------------------------------------------

beta=0.99;
alpha=0.5;
sigma=2;
chi=1;
varphi=2;
delta=0.035;
CY=0.45;
G=1.0574;

%%%%%%%%%%%%%%%%%%

%Bayes estimation
rho_a=0.8;
rho_g=0.8;

sigma_a=0.01;
sigma_g=0.01;

%%%%%%%%%%%%%%%%%%%

model;

%----------------------------------------------------------------
%\\household
%----------------------------------------------------------------

exp(c)^(-sigma)=beta*exp(r)*(exp(c(+1))^(-sigma));
exp(w)*(exp(c)^(-sigma))=chi*exp(l)^varphi;
exp(r)=exp(r_k)+1-delta;

%----------------------------------------------------------------
%\\good producer
%----------------------------------------------------------------

exp(y)=exp(u_a)*(exp(k(-1))^alpha)*(exp(l)^(1-alpha));
exp(w)=(1-alpha)*exp(y)/exp(l);
exp(r_k)=alpha*exp(y)/exp(k(-1));

%----------------------------------------------------------------
%\\market clearing
%----------------------------------------------------------------

exp(inv)=exp(k)-(1-delta)*exp(k(-1));
exp(y)=exp(c)+exp(inv)+G*exp(u_g);

%----------------------------------------------------------------
%\\shock
%----------------------------------------------------------------

u_a=rho_a*u_a(-1) + es_a;
u_g=rho_g*u_g(-1) + es_g;

end;

%initval;

steady_state_model;

r=log(1/beta);
r_k=log(exp(r)+delta-1);
KY=alpha/exp(r_k);
LY=(KY)^(alpha/(alpha-1));
w=log((1-alpha)*(LY)^(-1));
y=log(((CY)^(-sigma)*(LY)^(-varphi)*exp(w)*(chi)^(-1))^(1/(sigma+varphi)));
c=log(CY*exp(y));
l=log(LY*exp(y));
k=log(KY*exp(y));
inv=log(delta*exp(k));
G=exp(y)-exp(c)-exp(inv);
u_a = 0;
u_g = 0;
end;

steady;

check;

shocks;

var es_a = sigma_a^2;

end;

stoch_simul(order=2)y c k inv l w r_k r;










