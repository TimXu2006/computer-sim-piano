# Piano time-step validation

Run in MATLAB (or Octave) from this directory:

```matlab
validate_piano
```

`piano.m` now calls the same `simulate.m` solver used by the validator. The solver's optional ninth argument refines its original time step by an integer factor; with that argument omitted, its original behavior is retained. The optional tenth argument limits simulation duration. The numerical update loop and physical parameters are unchanged.

The validator runs the full model for one amplitude-1 A4 strike on the fixed J=81 grid for 0.1 s. The note is held for 1 s, so release cannot affect the comparison. It checks finite, nonzero raw signals and matching physical sample times. It reports the dt-to-dt/2 and dt/2-to-dt/4 relative L2 differences using `norm(S_dt4)` as a shared denominator. If needed, it refines further until two consecutive adjacent differences are below the chosen 1% tolerance and decreasing. It then compares the original step with the finest qualified result and identifies the coarsest tested step that meets 1%.

After a successful run, `results/validation/` contains `timestep_comparison.png`, `errors.csv`, and `validation_report.txt` with an English paragraph based on the measured values. The 1% tolerance is a practical choice for raw waveforms, not an instructor or perceptual standard. This test covers only one strike at one fixed spatial grid.

Execution status in this environment: **not run**. MATLAB and Octave were absent. An Octave installation was attempted but stopped because available disk space was under 1 GB and the installation required many dependencies. No MATLAB/Octave results or conclusions have been fabricated.
