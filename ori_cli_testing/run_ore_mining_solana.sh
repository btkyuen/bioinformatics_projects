#!/bin/bash
RUN_CMD="ore --rpc https://api.mainnet-beta.solana.com --keypair /home/byuen/.config/solana/id.json --priority-fee 1000 mine --threads 4"
while true; do
	echo "Starting miner."
	# $RUN_CMD
	ore \
		--rpc https://api.mainnet-beta.solana.com \
		--keypair /home/byuen/.config/solana/id.json \
		--priority-fee 1000 \
		mine \
		--threads 4

	if [ $? -ne 0 ]; then
		echo "Mining failed; restarting."
	else
		echo "Mining successful; restarting."
	fi
	
	sleep 5
done
