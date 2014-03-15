count=$(node find.js $@)
echo $count
while true; do
	for job in $(node find.js $@ $count); do
		count=$((count+1))
		echo $count
		mail -s $(hostname) JackyLeeHongJian@Gmail.com <<< $job
	done
	sleep 1m
done
