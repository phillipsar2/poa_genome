#!/bin/bash
count=1
while true; do
  jobs=$(squeue -o "%8i" -u arnstrm | grep -v "JOBID" | wc -l)
  if [  ${jobs} -gt 0 ]; then
    echo -e "$(date +"%D  %r")\tJobs still running. Number of jobs in queue: $jobs"
    sleep 600
  else
    echo  -e "$(date +"%D  %r")\trunning falcon now round ${count}"
    echo "command: fc_run falcon.cfg &> stdout-${count}.txt"
    fc_run falcon.cfg &> stdout-${count}.txt
    if [ $? -eq 0 ]; then
      echo  -e "$(date +"%D  %r")\tfalcon run complete"
      exit 0
    fi
    echo -e "$(date +"%D  %r")\tdone"
    ((count++))
  fi
done
