#! /bin/bash



#############################################################
infile=''
is_braket=false


#############################################################
while [ $# -gt 0 ]; do
	case $1 in
		-i)
			infile=$2
			shift
			;;
		--bracket)
			is_bracket=true
			;;
		*)
			infile=$1
			;;
	esac
	shift
done


if [ -z $infile ]; then
	echo "infile not given!" >&2
	exit 1
fi


#############################################################
if [ "$is_bracket" == true ]; then
	cat $infile | sed '4!d' | sed 's/.\+= //; s/: /:/g; s/\[//g; s/\]//g; s/ //g; s/&95%HPD=//g; s/{\([^,]\+\),\([^}]\+\)\+}/[\1-\2]/g'
else
	cat $infile | sed '4!d' | sed 's/.\+= //; s/: /:/g; s/\[//g; s/\]//g; s/ //g; s/&95%HPD=//g; s/{\([^,]\+\),\([^}]\+\)\+}/\1-\2/g'
fi


