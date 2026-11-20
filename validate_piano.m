function validate_piano
% Run the shared piano.m solver for one full-model A4 strike at fixed J=81.
root = fileparts(mfilename('fullpath'));
outdir = fullfile(root, 'results', 'validation');
if ~exist(outdir, 'dir'), mkdir(outdir); end
fs = 8192;
tol = 0.01;                  % Chosen practical raw-waveform L2 tolerance.
max_refinement = 128;        % Resource guard, not a convergence claim.
refinements = [1 2 4];
signals = cell(1, numel(refinements));
times = cell(1, numel(refinements));
dts = zeros(1, numel(refinements));
adjacent = NaN(1, numel(refinements));
for k = 1:numel(refinements)
    [signals{k}, times{k}, dts(k)] = run_one(refinements(k), fs);
    check_signal(signals{k}, times{k}, times{1});
end
% The first two differences explicitly share norm(S_dt4).
denom4 = norm(signals{3});
assert(isfinite(denom4) && denom4 > 0, 'dt/4 signal has zero or invalid norm.');
D1 = norm(signals{1} - signals{2}) / denom4;
D2 = norm(signals{2} - signals{3}) / denom4;
adjacent(2:3) = [D1 D2];

% For later pairs, each comparison uses the finer raw signal as denominator.
% Require two consecutive differences below 1%, decreasing, before using
% the finest signal as the numerical reference.
reference_index = 0;
if D1 < tol && D2 < tol && D2 < D1
    reference_index = 3;
end
while reference_index == 0 && refinements(end) < max_refinement
    k = numel(refinements) + 1;
    refinements(k) = 2 * refinements(end);
    [signals{k}, times{k}, dts(k)] = run_one(refinements(k), fs);
    check_signal(signals{k}, times{k}, times{1});
    adjacent(k) = norm(signals{k-1} - signals{k}) / norm(signals{k});
    % The dt/2-to-dt/4 value uses the same finer denominator here.
    previous = adjacent(k-1);
    if previous < tol && adjacent(k) < tol && adjacent(k) < previous
        reference_index = k;
    end
end

relative_to_reference = NaN(1, numel(refinements));
recommended_index = 0;
if reference_index > 0
    reference = signals{reference_index};
    for k = 1:numel(refinements)
        relative_to_reference(k) = norm(signals{k} - reference) / norm(reference);
    end
    eligible = find(relative_to_reference(1:reference_index-1) < tol, 1, 'first');
    if ~isempty(eligible), recommended_index = eligible; end
end

csvfile = fullfile(outdir, 'errors.csv');
fid = fopen(csvfile, 'w');
assert(fid >= 0, 'Cannot write errors.csv.');
fprintf(fid, 'refinement,dt_seconds,adjacent_relative_L2,relative_to_reference_L2,is_reference\n');
for k = 1:numel(refinements)
    fprintf(fid, '%d,%.16g,%.16g,%.16g,%d\n', refinements(k), dts(k), ...
        adjacent(k), relative_to_reference(k), k == reference_index);
end
fclose(fid);

fig = figure('Visible', 'off', 'Color', 'w');
plot(times{1}, signals{1}, '-', times{2}, signals{2}, '--', ...
    times{3}, signals{3}, ':', 'LineWidth', 1);
hold on;
labels = {'dt', 'dt/2', 'dt/4'};
if reference_index > 3
    plot(times{reference_index}, signals{reference_index}, '-.');
    labels{end+1} = sprintf('dt/%d reference', refinements(reference_index));
end
xlabel('Time (s)'); ylabel('Raw displacement S(t)');
title('Full model, one A4 strike, J=81, 0.1 s');
legend(labels, 'Location', 'best'); grid on;
print(fig, fullfile(outdir, 'timestep_comparison.png'), '-dpng', '-r150');
close(fig);

fprintf('D1 = %.8g (%.4f%%); D2 = %.8g (%.4f%%); D2/D1 = %.4f\n', ...
    D1, 100*D1, D2, 100*D2, D2/D1);
fprintf('Convergence evidence (D2 < D1): %d\n', D2 < D1);
if reference_index == 0
    fprintf('No qualified reference through dt/%d; accuracy undecided.\n', refinements(end));
else
    fprintf('Reference: dt/%d; original dt difference = %.4f%%.\n', ...
        refinements(reference_index), 100*relative_to_reference(1));
    if recommended_index > 0
        fprintf('Coarsest tested step below 1%%: dt/%d (%.4f%%).\n', ...
            refinements(recommended_index), 100*relative_to_reference(recommended_index));
        if recommended_index > 1
            fprintf('Minimal change: pass refinement=%d as the ninth simulate input in piano.m.\n', ...
                refinements(recommended_index));
        end
    else
        fprintf('No tested step coarser than reference is below 1%%.\n');
    end
end
fprintf('Saved %s and timestep_comparison.png\n', csvfile);

reportfile = fullfile(outdir, 'validation_report.txt');
fid = fopen(reportfile, 'w');
assert(fid >= 0, 'Cannot write validation_report.txt.');
fprintf(fid, 'Validation. The full piano model was run for one A4 strike of amplitude 1 at J=81 for 0.1 s, with the key held for 1 s. ');
fprintf(fid, 'The original time step was %.9g s. Raw outputs were compared at identical physical sample times without per-run normalization. ', dts(1));
fprintf(fid, 'Using norm(S_dt4) as the common denominator, D1 = %.4f%% and D2 = %.4f%% (D2/D1 = %.3f). ', 100*D1, 100*D2, D2/D1);
if D2 < D1
    fprintf(fid, 'The decreasing differences support temporal convergence in this test. ');
else
    fprintf(fid, 'The initial differences do not support temporal convergence. ');
end
if reference_index == 0
    fprintf(fid, 'A numerical reference meeting the two-successive-refinement criterion was not obtained through dt/%d; adequacy of the original step is undetermined. ', refinements(end));
else
    fprintf(fid, 'The dt/%d result met the chosen reference criterion; the original-step difference from it was %.4f%%. ', refinements(reference_index), 100*relative_to_reference(1));
    if relative_to_reference(1) < tol
        fprintf(fid, 'The original step met the selected 1%% practical tolerance for this test. ');
    else
        fprintf(fid, 'The original step exceeded the selected 1%% practical tolerance for this test. ');
        if recommended_index > 0
            fprintf(fid, 'The coarsest tested step meeting it was dt/%d (%.4f%% from the reference). ', refinements(recommended_index), 100*relative_to_reference(recommended_index));
            fprintf(fid, 'A minimal production change would pass refinement %d as the ninth simulate input. ', refinements(recommended_index));
        else
            fprintf(fid, 'No tested step coarser than the reference met it. ');
        end
    end
end
fprintf(fid, 'The 1%% threshold is a chosen practical raw-waveform L2 criterion, neither an instructor requirement nor a perceptual standard. ');
fprintf(fid, 'This test addresses only time-step behavior for one A4 strike on a fixed grid; it does not establish correctness of the full program.\n');
fclose(fid);
end

function [S, t, dt] = run_one(refinement, fs)
% A4 is string 25: 110*2^(24/12)=440 Hz; release at 1 s > 0.1 s.
[S, t, dt] = simulate(0, 1, 1, 25, true, true, true, fs, refinement, 0.1);
end

function check_signal(S, t, reference_times)
assert(~isempty(S) && numel(S) == numel(t), 'Empty or mismatched output.');
assert(all(isfinite(S)) && all(isfinite(t)), 'Nonfinite output.');
assert(norm(S) > 0, 'Zero signal.');
assert(numel(t) == numel(reference_times), 'Sample count differs.');
assert(max(abs(t - reference_times)) < 1e-12, 'Sample times differ.');
end
