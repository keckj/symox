# Attribution of integrated R5 vector linear Atan

The same accepted R5r1 package and the linear x/(x+1) fixture used in the [paired release measurements](../integrated-r5/README.md) are profiled for 500 public log_to_atan calls after official vector import. Runtime SymPy imports are blocked; exact expected output is checked before/after. [Raw cProfile records](profile-linear-r5.json), [runner](profile-linear.py) and source identity are retained.

In this instrumented run, log_to_atan totals 81.147ms, including 73.675ms in _stock_atan_protocol and 72.232ms in _atan_stock_operators. The native _log_to_atan kernel totals 7.014ms. Forwarding-chain authentication, function-state verification and vector proof verification account for much of the guard work. Cumulative entries overlap and must not be added together.

This identifies protocol authentication as the dominant cost for this one fixture. cProfile changes timing; these values are attribution evidence, not a second benchmark or an estimate of removable overhead. Actual user callbacks, descriptor behavior and mutable class/function protocols must remain preserved by any optimization.
