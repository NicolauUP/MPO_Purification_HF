# 2D runs

Place open-square run definitions here. Begin with the documented weak-coupling
`L=4` reference before attempting larger square lattices. Keep SP2 and PM
comparisons separate from physical parameter scans when purification behaviour
is uncertain.

For each run, record horizontal and vertical bond orders separately, together
with both Fock-field components, particle number, stationarity, and bond
dimensions. The current solver has a known finite-gap square MPO-SP2
regression, so inspect `PurificationResult` rather than relying only on the
SCF return value.
