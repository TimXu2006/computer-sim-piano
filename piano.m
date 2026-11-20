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
