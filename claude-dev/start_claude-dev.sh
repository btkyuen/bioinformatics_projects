#!/bin/bash

CONTAINER_NAME="claude-dev"

# Check if the container already exists (either running or stopped)
if podman ps -a --format '{{.Names}}' | grep -Eq "^${CONTAINER_NAME}$"; then
    echo "WARNING: A container named '${CONTAINER_NAME}' is already running or exists."
    echo "To enter the existing container, run:  podman exec -it ${CONTAINER_NAME} bash"
    echo "To clear it and mount this new folder, run:  podman rm -f ${CONTAINER_NAME}"
    exit 1
fi

echo "Starting Claude Code environment in $PWD..."

# Launch the container in the background
podman run -d \
  --name "${CONTAINER_NAME}" \
  --privileged \
  -v "$PWD":/workspace \
  -v claude-dev:/root/.claude \
  -e CLAUDE_CONFIG_DIR="/root/.claude" \
  -e CONTAINERS_STORAGE_DRIVER="vfs" \
  claude-dev:latest

# Instantly drop you into the bash shell
podman exec -it "${CONTAINER_NAME}" bash
