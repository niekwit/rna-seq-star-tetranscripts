#!/usr/bin/env bash

LOG=${snakemake_log[0]}
URL=${snakemake_params["url"]}
OUTPUT=${snakemake_output[0]}

curl -fL "$URL" -o "$OUTPUT.gz" 2> "$LOG" \
    || { echo "Download failed for $URL" >> "$LOG"; exit 1; }

[ -s "$OUTPUT.gz" ] \
    || { echo "Downloaded file is empty: $OUTPUT.gz" >> "$LOG"; exit 1; }

pigz -df "$OUTPUT.gz" 2>> "$LOG"
