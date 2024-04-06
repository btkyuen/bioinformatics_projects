#!/bin/bash

# Start mining operation on Solana RPC
SESSION_NAME=$(base64 /dev/urandom | head -c 10)
tmux new-session -d -s "$SESSION_NAME"
tmux send-keys -t "$SESSION_NAME" "./run_ore_mining_solana.sh" Enter

# Start mining operation on Ironforge RPC
# (taken from ORE site browser RPC)
# SESSION_NAME=$(base64 /dev/urandom | head -c 10)
# tmux new-session -d -s "$SESSION_NAME"
# tmux send-keys -t "$SESSION_NAME" "./run_ore_mining_ironforge.sh" Enter