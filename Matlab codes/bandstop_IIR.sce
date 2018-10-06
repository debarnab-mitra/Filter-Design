
m = 47;
q_m = 4; 
r_m = 7;
B_l = 2 + 0.6*q_m + 1.5*r_m;
B_h = B_l + 3;

transition_width = 1;
f_s = 90;


un_ws1 = B_l;
un_ws2 = B_h;
un_wp1  = un_ws1 - transition_width;
un_wp2  = un_ws2 + transition_width;

w_p1 = 2*%pi*un_wp1/f_s;
w_p2 = 2*%pi*un_wp2/f_s;
w_s1 = 2*%pi*un_ws1/f_s;
w_s2 = 2*%pi*un_ws2/f_s;

//analog filter specifications
W_p1 = tan(w_p1/2);
W_p2 = tan(w_p2/2);
W_s1 = tan(w_s1/2);
W_s2 = tan(w_s2/2);


W_o = sqrt(W_p1*W_p2);
B = W_p2 - W_p1

//low pass
Wl_p1 = -1;
Wl_p2 = 1;
Wl_s1 = (B*W_s1)/(-1*W_s1^2 + W_o^2);
Wl_s2 = (B*W_s2)/(-1*W_s2^2 + W_o^2);

Wl_s = min(abs(Wl_s1),abs(Wl_s2));
Wl_p = 1;


//chebyschev
delta = 0.1 
D1 = (1/((1-delta)^2)) - 1;
D2 = (1/delta^2) - 1;

//N>= Ns
Ns = (acosh(sqrt(D2/D1)))/(acosh((Wl_s/Wl_p)));
epsilon = sqrt(D1);
N = ceil(Ns);

// hanalog(s)
k_o = 1/sqrt(1+ epsilon^2);

B_k = (asinh(1/epsilon))/N;

// we get N = 4
// finding poles
k = linspace(N,(2*N-1),N);
A_k = ((2.*k +1)*%pi)/(2*N);
S_7_k =  Wl_p*(sin(A_k).*sinh(B_k)+ complex(0,1).*cos(A_k).*cosh(B_k));

// defining the polynominal
s= poly(0,'s');
a = k_o*((-1)^N)*(real(prod(S_7_k)));
b = real(poly(S_7_k,'s'));
HanalogLPF = syslin('c', a/b); 


// converting back to analog bandpass
HanalogBPF = syslin('c',horner(HanalogLPF,((B*s)/(s^2+ W_o^2))));

// converting to digital domain
z= poly(0, 'z');
p = (z-1)/(z+1);
q = (B*p)/(p^2+ W_o^2);
Hdiscrete = syslin('c',horner(HanalogLPF,q));


k = linspace(0,1023,1024);
k1 = linspace(0,1,2);
w = (%pi*k)/1024;
y1 = horner(Hdiscrete,complex(cos(w),sin(w)));
y2 = 1*ones(1024,1);
y3 = 0.9*ones(1024,1);
y4 = 0.1*ones(1024,1);
y_wp1 = w_p1*ones(2,1);
y_wp2 = w_p2*ones(2,1);
y_ws1 = w_s1*ones(2,1);
y_ws2 = w_s2*ones(2,1);
plot(w,abs(y1),'g');
plot(w,y2);
plot(w,y3);
plot(w,y4);
plot(y_wp1,k1);
plot(y_wp2,k1);
plot(y_ws1,k1);
plot(y_ws2,k1);










