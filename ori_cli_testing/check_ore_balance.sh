#!/bin/bash
KEY=/home/byuen/.config/solana/id.json
RPC_IRON="https://rpc.ironforge.network/mainnet?apiKey=01HTD8PPGDM1JBVQVEVJKXZ47F"
RPC_SOL="https://api.mainnet-beta.solana.com"
echo "Solana Network ORE Balance:"
ore --keypair $KEY --rpc $RPC_SOL rewards
# echo "Ironforge Network ORE Balance:"
# ore --keypair $KEY --rpc $RPC_IRON rewards
