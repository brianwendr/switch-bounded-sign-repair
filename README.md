# Switch-Bounded Sign Repair: Replication Package

## Contents

```text
replicate_dam_frontier_tables.py      # Main replication script
expected_csv/                         # Manuscript-aligned CSV versions of Tables 1--5
trace/representative_output_trace.log # Representative console output
README.md                             # This file
```

## Requirements

- Python 3.11 or later
- Standard library only

No third-party packages are needed.

## How to run

From the repository root, run:

```bash
python replicate_dam_frontier_tables.py
```

The script regenerates the main computational outputs used for the manuscript tables, including:

- repair-complexity frontier checks,
- penalty-path calculations,
- minimum-run-length repair checks,
- scaling checks, and
- exhaustive verification for all binary sign strings up to `n <= 10` and `B <= 5`.

A successful run should end with:

```text
exhaustive_check n<=10 B<=5 True
```

## Files generated or checked

The `expected_csv/` folder contains manuscript-aligned versions of:

- `table1.csv`: positioning comparison table,
- `table2.csv`: repair-complexity frontier summary,
- `table3.csv`: representative penalty path,
- `table4.csv`: minimum-run-length repair summary,
- `table5.csv`: scaling check.

The file `trace/representative_output_trace.log` gives a representative console trace from the replication run. Runtime values are wall-clock values and may vary by machine.

## Data availability

No external empirical dataset is used. All computational values are generated from the supplied Python script.

## Notes

This package is intentionally minimal. It is designed to let readers verify the computational claims without installing additional software or downloading data.
