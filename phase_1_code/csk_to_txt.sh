for file in all_data/norm_velocity_data/*COMMAND.txt; do
	mv "$file" "`basename $file COMMAND.txt`command.txt";
done