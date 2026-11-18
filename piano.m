clear all;

beat = 0.4;
melody = [ ...
  22 1; 22 1; 23 1; 25 1;  25 1; 23 1; 22 1; 20 1;
  18 1; 18 1; 20 1; 22 1;  22 1.5; 20 0.5;   20 2;
  22 1; 22 1; 23 1; 25 1;  25 1; 23 1; 22 1; 20 1;
  18 1; 18 1; 20 1; 22 1;  20 1.5; 18 0.5;   18 2];
nmel  = size(melody,1);
tnote = ([0; cumsum(melody(1:end-1,2))]*beat)';
dnote = 0.9*melody(:,2)'*beat;
dnote(end) = 3*beat;
anote = ones(1,nmel);
inote = melody(:,1)';

chords = 'IIIVIIIVIIIVIIVI';
voicing.I = [6 10 13];
voicing.V = [1 13 17];
amp       = [0.8 0.45 0.45];
for c = 1:length(chords)
  notes = voicing.(chords(c));
  hold_time = 1.9*beat;
  if c == length(chords), hold_time = 3*beat; end
  tnote = [tnote, (c-1)*2*beat*ones(1,3)];
  dnote = [dnote, hold_time*ones(1,3)];
  anote = [anote, amp];
  inote = [inote, notes];
end

[tnote, order] = sort(tnote);
dnote = dnote(order); anote = anote(order); inote = inote(order);

fs = 8192;
models = {
  'piano_base.wav', false, false, false;
  'piano_hammer.wav', true, false, false;
  'piano_everything.wav', true, true, true };

outdir = 'results';
if ~exist(outdir, 'dir'), mkdir(outdir); end

figure(1); clf
for m = 1:size(models,1)
  [S, tsave] = simulate(tnote, dnote, anote, inote, ...
                        models{m,2}, models{m,3}, models{m,4}, fs);
  outfile = fullfile(outdir, models{m,1});
  audiowrite(outfile, S/max(abs(S)), fs);
  fprintf('Wrote %s\n\n', outfile);

  subplot(size(models,1),1,m)
  plot(tsave, S)
  title(strrep(models{m,1}, '_', '\_')); ylabel('sound level')
end
xlabel('time (s)')

for m = 1:size(models,1)
  outfile = fullfile(outdir, models{m,1});
  [S, fsread] = audioread(outfile);
  sound(S, fsread);
  pause(length(S)/fsread + 0.5);
end


function [S, tsave] = simulate(tnote, dnote, anote, inote, use_hammer, use_stiffness, use_damping, fs)

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
nskip = ceil(1/(fs*dtmaxmin));
dt    = 1/(fs*nskip);

tmax = max(tnote + dnote) + 0.5;
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