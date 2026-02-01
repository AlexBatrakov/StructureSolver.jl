# Contributing

Thanks for considering contributing!

## Quick checklist

- Make sure CI is green.
- Add or update tests when changing behavior.
- Keep changes focused and documented.

## Local setup

From the repository root:

```bash
julia --project -e 'using Pkg; Pkg.instantiate()'
```

Run tests:

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

Build docs:

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

## Reporting issues

When filing a bug report, please include:

- Julia version (`versioninfo()`)
- OS and architecture
- A minimal reproducible example
- Full stacktrace/output if applicable

## Pull requests

- Prefer small PRs that do one thing.
- If you add new public API, please add docstrings and consider documenting it in the docs site.
- If you change solver behavior, add a regression test.

## Release process

See `RELEASING.md`.
