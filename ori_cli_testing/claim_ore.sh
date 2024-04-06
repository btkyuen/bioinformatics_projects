#!/bin/bash
while true; do
	ore \
		--keypair /home/byuen/.config/solana/id.json \
		--rpc https://mainnet.helius-rpc.com/?api-key=38087e09-fa05-45b6-9b95-c3c93f6e5836 \
		--priority-fee 100000000 \
		claim

	if [ $? -ne 0 ]; then
		echo "Claiming failed; restarting."
	else
		echo "Claiming successful; restarting."
	fi
	
	sleep 5
done