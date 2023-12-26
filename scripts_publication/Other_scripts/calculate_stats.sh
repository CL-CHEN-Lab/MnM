#!/bin/bash

# calculate_stats.sh
# chmod +x calculate_stats.sh
# sh ~/OneDrive\ -\ INSTITUT\ CURIE/Scripts/scCNV/calculate_stats.sh

# Check if the directory path argument is provided
if [ -z "$1" ]; then
  # Use the default path if not provided
  dir="Kronos/CNV"
else
  dir="$1"
fi

echo "Filename,Cell_count,Minimum,Maximum,Lower_Quartile_(Q1),Upper_Quartile_(Q3),Median_RPM,Median_ploidy"
for csv in $(ls "$dir"/*csv | grep -v ALL | grep -v filtered); do
	trimmedname=$(basename "${csv}")
  numbers=$(tail -n+2 $csv | awk -F',' '{print $7}')
	cellcount=$(tail -n+2 $csv | wc -l | awk '{print $1}')
  min=$(printf '%s\n' "${numbers[@]}" | sort -n | head -n 1)
  max=$(printf '%s\n' "${numbers[@]}" | sort -n | tail -n 1)
  median=$(printf '%s\n' "${numbers[@]}" | sort -n | awk 'BEGIN {count=0} {data[count]=$1; count++} END {if (count % 2 == 1) {print data[int(count/2)]} else {print (data[int(count/2)-1] + data[int(count/2)]) / 2}}')
	lower_quartile=$(printf '%s\n' "${numbers[@]}" | sort -n | awk 'BEGIN {count=0} {data[count]=$1; count++} END {idx=int(count/4); if (count % 4 == 0) {print (data[idx-1] + data[idx]) / 2} else {print data[idx]}}')
	upper_quartile=$(printf '%s\n' "${numbers[@]}" | sort -n | awk 'BEGIN {count=0} {data[count]=$1; count++} END {idx=int(3*count/4); if (count % 4 == 0) {print (data[idx-1] + data[idx]) / 2} else {print data[idx]}}')
	ploidy=$(tail -n+2 $csv | awk -F',' '{print $3}')
	median_ploidy=$(printf '%s\n' "${ploidy[@]}" | sort -n | awk 'BEGIN {count=0} {data[count]=$1; count++} END {if (count % 2 == 1) {print data[int(count/2)]} else {print (data[int(count/2)-1] + data[int(count/2)]) / 2}}')

  printf "$trimmedname,"
	printf "$cellcount,"
  printf "$min,"
  printf "$max,"
	printf "$lower_quartile,"
	printf "$upper_quartile,"
  printf "$median,"
	printf "$median_ploidy"
	echo ""
done
