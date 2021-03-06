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

delta_w = 2*%pi*1/140;
constant_k = 12/(2.285*delta_w);


N = 63;

wc1 = (w_s1+w_p1)/2;  //cutoff frequency of Lowpass filter 1
wc2 = (w_s2+w_p2)/2; //cutoff frequency of Lowpass filter 2
wn = zeros(1,(2*N+1)); // the impulse response 
for  i= -N:1:N
    if(i==0)
        wn(i+N+1) = (wc2-wc1)/%pi;
    else
        wn(i+N+1) = (sin(wc2*i)-sin(wc1*i))/(i*%pi);
    end
end

n = linspace(-N, N,2*N+1);
z = poly (0 , 'z','r');
pow_z = (1/z)^n;
Hdiscrete = sum(wn.*pow_z); // H(z)

k = linspace(0,1023,1024);
k1 = linspace(0,1,2);
w = (%pi*k)/1024;
y1 = horner(Hdiscrete,complex(cos(w),sin(w)));
y2 = 1.1*ones(1024,1);
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
