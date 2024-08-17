#!/bin/sh

NVIDIA_SMI=$(which nvidia-smi)

export LC_ALL=C
export LANG=C

if [ -x "$NVIDIA_SMI" ]; then
    CUDA_VERSION=$($NVIDIA_SMI --version | grep "CUDA Version" | awk '{print $NF}')
    echo "$CUDA_VERSION"
    exit 0
else
    echo "CUDA is not available" >&2
    exit 1
fi