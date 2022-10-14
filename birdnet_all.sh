# Iterate through all checklist subdirectories and perform BirdNet identifications
# NOTE This must be run within the BirdNET directory or the library won't work...
# Parameters: 
#    1 - Path to target directory, containing subdirectories for each checklist
#    2 - Latitude
#    3 - Longitude
#    4 - Week (of year) 
echo "Processing all files in directory $1 using Lat: $2 Lon: $3 Week: $4"
echo ""
for d in $1*/ ; do
	echo ""
	echo "Processing audio files in ${d}"
	python analyze.py --i "${d}" --lat $2 --lon $3 --week $4
done
