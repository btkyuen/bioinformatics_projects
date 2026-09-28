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
  --userns=keep-id \
  -v "$PWD":/workspace:Z \
  -v claude-dev:/home/claude-dev/.claude \
  -e CLAUDE_CONFIG_DIR="/home/claude-dev/.claude" \
  claude-env:latest >/dev/null

# Instantly drop you into the bash shell
podman exec -it "${CONTAINER_NAME}" bash
