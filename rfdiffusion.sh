#!/usr/bin/env bash
# RFdiffusion via the maintained Docker image -- comp5 / KRAS G12D pMHC wrapper.
#
# Derived from tfr1-comp4/rfdiffusion.sh. Same two invariants as every prior
# round, repeated here because getting either wrong fails before any sampling:
#
#   * Pass ONLY Hydra overrides. The image's ENTRYPOINT already invokes
#     run_inference.py, so a leading `python` or a script path makes Hydra try
#     to parse it as a config key. Do not override the container workdir either;
#     the entrypoint resolves its script path relative to /app/RFdiffusion.
#   * $PWD is mounted at /work, so every path you pass must be container-
#     absolute under /work, and you must include hydra.run.dir=/work/.hydra.
#
# Writable mounts, needed because we run as the host UID while /app is
# root-owned in the image:
#   /app/RFdiffusion/schedules - IGSO3 diffusion schedule cache, mkdir'd on
#                                first run (model_runners.py:129). Without it
#                                the run dies with PermissionError. Persisted
#                                across runs; saves several minutes of
#                                "Calculating IGSO3" every time.
#
# WEIGHTS: this follows comp4, NOT comp2. The comp2 wrapper bind-mounts
# $HOME/tnf-track3/RFdiffusion/models over /app/RFdiffusion/models, and that
# host dir holds only 3 checkpoints (Base, Complex_base, Complex_beta) against
# the image's 7 -- so the mount MASKS the image's own weights. comp5 needs
# Complex_base_ckpt.pt, which is present either way, but there is no reason to
# hide the others; we leave the image's set intact. Verified present in the
# audit at results/RUNLOG.md stage 00.
set -euo pipefail

mkdir -p "$PWD/schedules"

exec docker run --rm --gpus all \
  --user "$(id -u):$(id -g)" \
  -v "$PWD:/work" \
  -v "$PWD/schedules:/app/RFdiffusion/schedules" \
  rosettacommons/rfdiffusion:latest "$@"
