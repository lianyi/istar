i=0
s=0
for n in $(grep -n TORSDOF $1 | cut -d: -f1); do
	tail -n +$((s+1)) $1 | head -n $((n-s)) > $i.pdbqt
	s=$n
	i=$((i+1))
done
