# Troubleshooting

## Docs build fails with missing docstrings

Documenter can be configured to fail the build if there are docstrings in the module that are not included in the manual. If you see an error like “docstring not included”, add the symbol to a `@docs` block (for non-exported/experimental APIs) or ensure it is included by `@autodocs`.

## `PyPlot` / `matplotlib` issues

StructureSolver.jl does not depend on `PyPlot` by default.

If you install `PyPlot` yourself and see errors about `matplotlib`, the issue is typically that `PyCall` is pointed at a Python that does not have `matplotlib` installed.

Options:

- Install `matplotlib` into the Python environment used by `PyCall`, or
- Reconfigure `PyCall` to use a different Python (e.g. Conda) following the guidance printed by the error message.

## CI shows only `dev` docs, not `stable`

`/stable` docs appear after you create a version tag like `v0.1.1` (or any `v*` tag). Until then, use `/dev`.

## GitHub Pages shows only `main`

The `gh-pages` branch is created by the docs workflow after the first successful run. Trigger the workflow from the Actions tab (or push to `main`) and then configure Pages to deploy from `gh-pages`.
