# Taken from and inspired by https://thisiskristian.substack.com/p/i-gave-claude-code-full-access-to
# For use in Podman; Claude authentication will need to be run from within the container

FROM ubuntu:24.04

# Dependencies
RUN apt-get update && \
    apt-get install -y \
    bash curl git bubblewrap socat ripgrep jq python3 python3-pip nodejs npm \
	openjdk-21-jdk podman && \
    rm -rf /var/lib/apt/lists/*

# Install Nextflow
RUN curl -fsSL https://get.nextflow.io | bash && \
    mv nextflow /usr/local/bin/

# Install Claude Code
RUN npm install -g @anthropic-ai/claude-code

# Add unpriveleged user
RUN useradd -m -s /bin/bash claude-dev

# Set working directory
USER claude-dev
WORKDIR /workspace

# Keep container running indefinitely
CMD ["tail", "-f", "/dev/null"]

# Build this Podman image using:
# podman build -t claude-dev:latest -f claude-setup.dockerfile .