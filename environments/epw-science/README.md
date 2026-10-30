# Isolated EPW experiment runtime

Create the locked Python 3.12 environment with `uv sync --project environments/epw-science --frozen`. Run the EPW frequency experiment with this environment's Python, `PYTHONPATH` set to the repository root, and all five BLAS/OpenMP thread variables set to `1`. The lock is an intended runtime configuration; actual imports, calculations, and native evidence must be verified separately before claiming execution.
