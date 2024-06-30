function TestSimJumpDiff();
S0  = 100;
N   = 250;
T   = 1;
mu  = 0.2;
sig = 0.2;
etau=0.05;
ldau=1;
etad=0.05;
ldad=4;


[tv,Jvu,Jvd,Xt,St] = SimJumpDiff(S0, N, T, mu,sig,etau,ldau,etad,ldad);
plot(tv,Jvu,tv,-Jvd,tv,Xt)
%plot(tv,St)