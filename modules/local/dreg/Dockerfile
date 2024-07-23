ADD file ... in /
CMD ["bash"]
ENV R_VERSION=4.0.5
ENV R_HOME=/usr/local/lib/R
ENV TZ=Etc/UTC
COPY scripts/install_R_source.sh /rocker_scripts/install_R_source.sh # buildkit
RUN /bin/sh -c /rocker_scripts/install_R_source.sh #
ENV CRAN=https://packagemanager.rstudio.com/cran/__linux__/focal/2021-05-17
ENV LANG=en_US.UTF-8
COPY scripts /rocker_scripts # buildkit
RUN /bin/sh -c /rocker_scripts/setup_R.sh #
CMD ["R"]
LABEL org.opencontainers.image.licenses=GPL-2.0-or-later org.opencontainers.image.source=https://github.com/rocker-org/rocker-versioned2 org.opencontainers.image.vendor=Rocker Project
ENV CUDA_VERSION=11.1
ENV NCCL_VERSION=2.7.8
ENV NVIDIA_VISIBLE_DEVICES=all
ENV NVIDIA_DRIVER_CAPABILITIES=compute,utility
ENV NVIDIA_REQUIRE_CUDA=cuda>=11.1 brand=tesla,driver>=418,driver<419 brand=tesla,driver>=440,driver<441 brand=tesla,driver>=450,driver<451
ENV CUDA_HOME=/usr/local/cuda
ENV LD_LIBRARY_PATH=:/usr/local/cuda/lib64:/usr/local/cuda/extras/CUPTI/lib64:/usr/local/cuda/lib64/libnvblas.so:
ENV LIBRARY_PATH=/usr/local/cuda/lib64/stubs
ENV NVBLAS_CONFIG_FILE=/etc/nvblas.conf
ENV PYTHON_CONFIGURE_OPTS=--enable-shared
ENV RETICULATE_MINICONDA_ENABLED=FALSE
ENV PATH=/usr/local/cuda/bin:/usr/local/nviida/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/local/texlive/bin/linux
RUN /bin/sh -c /rocker_scripts/install_cuda-11.1.sh
RUN /bin/sh -c /rocker_scripts/config_R_cuda.sh
RUN /bin/sh -c /rocker_scripts/install_python.sh
