count=$(node daemon.js $@)
echo $count
while true; do
	for job in $(node daemon.js $@ $count); do
		count=$((count+1))
		echo $count
		mail -s $(hostname) JackyLeeHongJian@Gmail.com <<< $job
	done
	sleep 1m
done
