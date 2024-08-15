#!/bin/sh

NVIDIA_SMI=$(which nvidia-smi)
IMAGE_NAME="docker.io/l1drm/dreg"

export LC_ALL=C
export LANG=C

if [ -x "$NVIDIA_SMI" ]; then
    CUDA_VERSION=$(nvidia-smi --version | grep "CUDA Version" | awk '{print $NF}')
    echo "${IMAGE_NAME}:cuda${CUDA_VERSION}"
    exit 0
else
    echo "CUDA is not available" >&2
    exit 1
fi