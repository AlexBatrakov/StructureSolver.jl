# Releasing

This repo uses a lightweight release process that works well for Julia packages.

## One-time setup (recommended)

- Ensure GitHub Actions are enabled.
- (Optional) If you want coverage shown publicly, create a Codecov project and add `CODECOV_TOKEN` in repo secrets.

## How to cut a release

1. Decide the next version `X.Y.Z`.
2. Update `version = "X.Y.Z"` in `Project.toml`.
3. Add an entry to `CHANGELOG.md`.
4. Tag and push:

   ```bash
   git tag vX.Y.Z
   git push origin vX.Y.Z
   ```

5. GitHub Actions:
   - `Docs` workflow runs on tags `v*` and publishes `/stable` docs for that tagged version.
   - `TagBot` creates/updates the GitHub Release page for the tag.

## What the automations do (and why)

### TagBot

- Purpose: takes your git tag (e.g. `v0.1.2`) and turns it into a GitHub Release automatically.
- Why recruiters like it: it signals a real release cadence and makes it easy to browse changes.

Trigger options:
- Push a tag `v*` (usual flow), or
- Run the workflow manually (`workflow_dispatch`), or
- Comment `@TagBot` on an issue (handy if you want to re-run it).

### CompatHelper

- Purpose: opens PRs that update `[compat]` bounds when your dependencies release new compatible versions.
- Why it matters: keeps the package installable in fresh environments and reduces “dependency rot”.
- Safety: it does not change your source code; it proposes compat-bump PRs you can review/merge.

### Coverage (Codecov)

- CI generates an `lcov.info` report (only on Julia 1.11 job) and uploads it.
- This is non-blocking (`fail_ci_if_error: false`), so it won’t break CI if Codecov isn’t configured.
- If you later decide to add a badge, we can do it once you’re happy with the percentage.
