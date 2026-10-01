#!/usr/bin/env bash
# Run queued 2403 commands (one per line of QUEUE) while fewer than MAXJOBS 'gtb solve' drivers are alive.
# Usage: setsid nohup scripts/gtb/queue.sh QUEUE_FILE [MAXJOBS] > queue_runner.out 2>&1 &
# Each line is a full shell command; its driver log goes to <workdir>.driver.log when the line has --workdir.
set -u
Q=$1; MAXJOBS=${2:-4}; LOG=${Q%.txt}.log
echo $$ > ${Q%.txt}.pid
while true; do
  live=$(ps -eo args | grep -c "^python3 -m smatrix_bootstrap.run gtb solve")
  if [ "$live" -lt "$MAXJOBS" ] && [ -s "$Q" ]; then
    cmd=$(head -n1 "$Q"); tail -n +2 "$Q" > "$Q.tmp" && mv "$Q.tmp" "$Q"
    wd=$(echo "$cmd" | grep -oE -- "--workdir [^ ]+" | head -n1 | awk '{print $2}')
    out=${wd:-/playpen1/shiqiu/gtb2403_runs/queue_job_$(date +%s)}.driver.log
    echo "$(date -u +%FT%TZ) launch (live=$live): $cmd" >> "$LOG"
    setsid nohup bash -c "cd /home/shiqiu/s-matrix-bootstrap && $cmd" > "$out" 2>&1 < /dev/null &
    sleep 30
  else
    [ -s "$Q" ] || { echo "$(date -u +%FT%TZ) queue empty, exiting" >> "$LOG"; exit 0; }
    sleep 60
  fi
done
