import sys
from pathlib import Path

import numpy as np

import kernelgecp
from kernelgecp import FermionicBasis, FermionicFit
from kernelgecp.integrations.green import read_green_tau

p = Path(sys.argv[1])
d = read_green_tau(p / "input.h5", p / "sim.h5", p / "bound-manifest.json")
assert "site-packages" in kernelgecp.__file__ and d.converged and d.iteration == 83
assert d.g_tau.shape[1:] == (1, 1, 2, 2)
basis = FermionicBasis(d.beta, 4.0, tol=1e-9, grid_order=8)
fit = basis.fit_tau(d.g_tau, d.tau, axis=0, rcond=1e-13)
path = Path(sys.argv[2])
fit.save(path)
restored = FermionicFit.load(path)
original = fit.evaluate_tau(d.tau)
prediction = restored.evaluate_tau(d.tau)
assert prediction.shape == d.g_tau.shape
np.testing.assert_allclose(prediction, original, atol=1e-14, rtol=1e-14)
assert np.max(np.abs(prediction - d.g_tau)) < 1e-6
print(
    "installed molecular fit and archive passed",
    float(np.max(np.abs(prediction - d.g_tau))),
)
