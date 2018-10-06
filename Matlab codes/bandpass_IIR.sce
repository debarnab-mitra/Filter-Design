
B_l = 18.8;
B_m = 23.8;
transition_width = 1;
f_s = 140;

un_wp1 = B_l;
un_wp2 = B_m;
un_ws1  = un_wp1 - transition_width;
un_ws2  = un_wp2 + transition_width;
 

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
Wl_s1 = (W_s1^2 - W_o^2)/(B*W_s1);
Wl_s2 = (W_s2^2 - W_o^2)/(B*W_s2);

Wl_s = min(abs(Wl_s1),abs(Wl_s2));
Wl_p = 1;

//butterworth IIR
delta = 0.1 
D1 = (1/((1-delta)^2)) - 1;
D2 = (1/delta^2) - 1;

//N>= Ns
Ns = (log(D2/D1))/(2*log(Wl_s/Wl_p));
// W_P<W_c < W_S
W_P = Wl_p*(D1^(-1/(2*ceil(Ns))));
W_S = Wl_s*(D2^(-1/(2*ceil(Ns))));
N = ceil(Ns);

W_c = (W_S+W_P)/2;

s= poly(0,'s');
b= (1+ (((s/W_c))^(2*N)));  
Sp = roots(b);
Deno_HanalogSl = poly([- 1.0689279 + 0.1693015*%i,  - 1.0689279 - 0.1693015*%i,  - 0.9642938 + 0.4913322*%i,  - 0.9642938 - 0.4913322*%i,  - 0.7652679 + 0.7652679*%i,  - 0.7652679 - 0.7652679*%i,  - 0.4913322 + 0.9642938*%i , - 0.4913322 - 0.9642938*%i , - 0.1693015 + 1.0689279*%i, -0.1693015 - 1.0689279*%i],'s');
//Deno_hanalogS1 = poly(Sp,'s'); 
a = 1;
Hbutterworth = syslin('c', a/b);

HanalogLPF = syslin('c', ((W_c)^(N))/real(Deno_HanalogSl));
//plzr(a/b);
//plzr(a/real(Deno_HanalogSl));
HanalogBPF = syslin('c',horner(HanalogLPF,(s^2+ W_o^2)/(B*s)));
//bode(HanalogBPF,0.01,100);

z= poly(0, 'z');
p = (z-1)/(z+1);
q = (p^2+ W_o^2)/(B*p);
Hdiscrete = syslin('c',horner(HanalogLPF,q));
//plzr(Hdiscrete);


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


