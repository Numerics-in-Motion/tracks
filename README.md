# Can the longer slide arrive first?

This standalone calculation compares two **frictionless, non-rolling sliding
points** released from rest at the same start and travelling to the same
endpoint: 1.000 m horizontally and 0.700 m downward. One path is straight;
the other is the endpoint-matched descending cycloid. Only path geometry
differs. It is not a rolling-ball or real-track test.

```sh
python -m pip install -r requirements.txt
python reproduce.py
```

The script solves for the cycloid parameter from the endpoint, computes path
length and first-arrival time, then checks those results by independent
quadrature and arc-length motion integration. It also checks energy and a
deliberately wrong constant-final-speed shortcut: that shortcut reverses the
winner, so it cannot support the claim. The frozen expected output and
tolerances are in `reference/result.json`; changed or inconsistent values
make the script exit nonzero.

For this ideal case the cycloid is 6.6683% longer yet arrives 13.9950%
sooner. Both reach the endpoint at the same speed because their vertical drop
is identical. The photographic video background is scenery, not evidence or
an input to this calculation.

This is a specific two-path comparison, not a claim that an experiment will
show identical timings or that every longer path is quicker. Real rolling,
sliding friction, track imperfections, and air drag are omitted.

The classical reference for the brachistochrone geometry is
[Wolfram MathWorld](https://mathworld.wolfram.com/BrachistochroneProblem.html).
