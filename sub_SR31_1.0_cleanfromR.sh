#Will clean up the things from R.
list=$1;
while read file;
	do echo $file;
	tail -n +2 $file > $file.clean;
	perl -p -i -e 's/"$//g' $file.clean;
	perl -p -i -e 's/.*"//g' $file.clean;
done < $list;
