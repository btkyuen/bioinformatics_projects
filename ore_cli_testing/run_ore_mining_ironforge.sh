#!/bin/bash
RUN_CMD="ore --rpc https://rpc.ironforge.network/mainnet?apiKey=01HTD8PPGDM1JBVQVEVJKXZ47F --keypair /home/byuen/.config/solana/id.json --priority-fee 1000 mine --threads 4"
while true; do
	echo "Starting miner."
	# $RUN_CMD
	ore \
		--rpc https://rpc.ironforge.network/mainnet?apiKey=01HTD8PPGDM1JBVQVEVJKXZ47F \
		--keypair /home/byuen/.config/solana/id.json \
		--priority-fee 100 \
		mine \
		--threads 4
		
	if [ $? -ne 0 ]; then
		echo "Mining failed; restarting."
	else
		echo "Mining successful; restarting."
	fi
	
	sleep 5
done
