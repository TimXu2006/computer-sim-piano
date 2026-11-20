function [S, tsave, dt] = simulate(tnote, dnote, anote, inote, use_hammer, use_stiffness, use_damping, fs, refinement, duration)
% The last two inputs are optional; the original piano.m call is unchanged.
if nargin < 9, refinement = 1; end
if nargin < 10, duration = []; end
assert(isscalar(refinement) && isfinite(refinement) && refinement >= 1 && ...
       refinement == floor(refinement), 'refinement must be a positive integer');

nmax = length(tnote);

L = 1;
J = 81; dx = L/(J-1);
flow = 110;
nstrings = 25;
f = zeros(nstrings,1); tau = f; M = f; T = f; R = f; G = f; B = f; EI = f; dtmax = f;
tauhi = 0.4;
fhi   = 2000;
for i = 1:nstrings
  f(i)   = flow*2^((i-1)/12);
  M(i)   = 1;
  if use_stiffness
    B(i) = 3.5e-4*(f(i)/261.63);
  else
    B(i) = 0;
  end
  T(i)   = M(i)*(2*L*f(i))^2/(1+B(i));
  EI(i)  = B(i)*T(i)*L^2/pi^2;
  if use_damping
    tau(i) = 1.5*sqrt(440/f(i));
    xi1  = (pi/L)^2;
    xihi = (pi*(fhi/f(i))/L)^2;
    b2 = max(0, (1/tauhi - 1/tau(i))/(xihi - xi1));
    b1 = 1/tau(i) - b2*xi1;
  else
    tau(i) = 1.2*(440/f(i));
    b1 = 0;
    b2 = L^2/(tau(i)*pi^2);
  end
  G(i) = 2*M(i)*b1;
  R(i) = 2*M(i)*b2;
  a = 4*T(i)/(M(i)*dx^2) + 16*EI(i)/(M(i)*dx^4);
  d = G(i)/M(i) + 4*R(i)/(M(i)*dx^2);
  dtmax(i) = (-d + sqrt(d^2 + 4*a))/a;
end
dtmaxmin = min(dtmax);
nskip = refinement*ceil(1/(fs*dtmaxmin));
dt    = 1/(fs*nskip);

tmax = max(tnote + dnote) + 0.5;
if ~isempty(duration), tmax = duration; end
clockmax = ceil(tmax/dt);

if use_hammer
  xh = L/8;
  jstrike = round(1+xh/dx) + (-1:1);
  p    = 2.5;
  vref = 1;
  mh   = 0.4*M*L;
  tc = 1e-3*sqrt(261.63./f);
  Ip = beta(1/(p+1), 1/2)/(p+1);
  dmax = tc*vref/(2*Ip);
  K = (p+1)*mh*vref^2./(2*dmax.^(p+1));
else
  xh1 = 0.25*L; xh2 = 0.35*L;
  jstrike = ceil(1+xh1/dx):floor(1+xh2/dx);
end
nc = length(jstrike);

cT = (dt/dx^2)*(T./M);
cB = (dt/dx^4)*(EI./M);
cR = (dt/dx^2)*(R./M);
cG = dt*(G./M);
cF = dt./(M*nc*dx);

tstop = zeros(nstrings,1);
H = zeros(nstrings,J);
V = zeros(nstrings,J);
yh = zeros(nstrings,1);
vh = zeros(nstrings,1);
hon = false(nstrings,1);
j = 2:(J-1);

count = 0;
S     = zeros(1,ceil(clockmax/nskip));
tsave = zeros(1,ceil(clockmax/nskip));
n = 1;
tic
for clock = 1:clockmax
  t = clock*dt;

  while (n <= nmax) && (tnote(n) <= t)
    i = inote(n);
    if use_hammer
      yh(i)  = mean(H(i,jstrike));
      vh(i)  = anote(n)*vref;
      hon(i) = true;
    else
      V(i,jstrike) = anote(n);
    end
    tstop(i) = t + dnote(n);
    n = n + 1;
  end

  off = (t > tstop);
  H(off,:) = 0; V(off,:) = 0; hon(off) = false;

  i = find(~off);
  if ~isempty(i)
    F = zeros(length(i),1);
    k = find(hon(i));
    if ~isempty(k)
      ik = i(k);
      comp = yh(ik) - mean(H(ik,jstrike),2);
      leaving = (comp < 0);
      hon(ik(leaving)) = false;
      comp = max(comp,0);
      F(k) = K(ik).*comp.^p;
      vh(ik) = vh(ik) - dt*F(k)./mh(ik);
      yh(ik) = yh(ik) + dt*vh(ik);
    end

    V(i,j) = V(i,j) ...
      + cT(i).*(H(i,j+1) - 2*H(i,j) + H(i,j-1)) ...
      + cR(i).*(V(i,j+1) - 2*V(i,j) + V(i,j-1)) ...
      - cG(i).*V(i,j);
    if use_stiffness
      Hg = [-H(i,2), H(i,:), -H(i,J-1)];
      Hxxxx = Hg(:,j+3) - 4*Hg(:,j+2) + 6*Hg(:,j+1) - 4*Hg(:,j) + Hg(:,j-1);
      V(i,j) = V(i,j) - cB(i).*Hxxxx;
    end
    if use_hammer
      V(i,jstrike) = V(i,jstrike) + cF(i).*F;
    end
    H(i,j) = H(i,j) + dt*V(i,j);
  end

  if mod(clock,nskip) == 0
    count = count + 1;
    S(count)     = sum(H(:,2));
    tsave(count) = t;
  end
end
fprintf('Done in %.1f s of time.\n', toc);
S = S(1:count); tsave = tsave(1:count);
end
