FROM pytorch/pytorch:1.10.0-cuda11.3-cudnn8-runtime

# Avoid interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Install OS dependencies
RUN apt-get update && apt-get install -y \
    libboost-all-dev \
    libsm6 \
    libxrender1 \
    libxext6 \
    && rm -rf /var/lib/apt/lists/*

# Install Python packages
RUN pip install --upgrade pip setuptools wheel

# Install core libraries
RUN pip install \
    dgl-cu111 -f https://data.dgl.ai/wheels/repo.html && \
    pip install \
    rdkit-pypi \
    tqdm \
    tensorboard \
    scikit-learn==0.23.2 \
    pandas \
    numpy\
    requests

# Set default working directory
WORKDIR /app

CMD ["/bin/bash"]


# # Use the PyTorch image as the base image
# FROM pytorch/pytorch:1.10.0-cuda11.3-cudnn8-runtime

# # Avoid interactive prompts during package installation
# ENV DEBIAN_FRONTEND=noninteractive

# # Install OS dependencies (for OPERA and Python libraries)
# RUN apt-get update && apt-get install -y \
#     libboost-all-dev \
#     libsm6 \
#     libxrender1 \
#     libxext6 \
#     wget \
#     tar \
#     openjdk-11-jre \
#     && rm -rf /var/lib/apt/lists/*

# # Install Python packages for ML and cheminformatics
# RUN pip install --upgrade pip setuptools wheel
# RUN pip install \
#     dgl-cu111 -f https://data.dgl.ai/wheels/repo.html && \
#     pip install \
#     rdkit-pypi \
#     tqdm \
#     tensorboard \
#     scikit-learn==0.23.2 \
#     pandas \
#     numpy

# # Set the working directory in the container
# WORKDIR /app

# # Copy OPERA installation package into the container
# COPY OPERA2.9_CL_Par.tar.gz /app/OPERA2.9_CL_Par.tar.gz

# # Extract OPERA installation package
# RUN tar -xvzf OPERA2.9_CL_Par.tar.gz

# # # List the files to check the directory structure
# # # RUN ls -al /app | tee /dev/stdout
# # # RUN echo ls

# # # Navigate to OPERA directory and run the installation using the custom installer
# # WORKDIR /app/OPERA2.9_CL_Par
# # RUN chmod +x OPERA_P_2.9_mcr_Installer.install && \
# #     ./OPERA_P_2.9_mcr_Installer.install

# # # Set OPERA binary path to environment variable
# # ENV PATH="/app/OPERA2.9_CL_Par:$PATH"

# # # Verify OPERA installation
# # RUN OPERA --version

# # Default command to keep the container running or use bash for testing
# CMD ["/bin/bash"]