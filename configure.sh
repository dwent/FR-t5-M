#!/bin/bash

# Please edit the following paths
# ==================================================
Modeller_Path="/tjjiang/jinxiaoyang/Software/modeller9.11"

# Choose one of the sequence databases
NR_Database=""
U90_Database=""
# ==================================================

# ========== Don't Change the code below ==========

if [ ! -d $Modeller_Path ]; then
	echo "can't find Modeller directory: $Modeller_Path"
	exit 1
fi

if [ -n "$NR_Database" ]&&[ ! -d $NR_Database ]; then
	echo "Can't find nr database: $NR_Database"
	exit 1
fi

if [ -n "$U90_Database" ]&&[ ! -d $U90_Database ]; then
	echo "Can't find u90 database: $U90_Database"
	exit 1
fi

if [ -z "$NR_Database" ]&&[ -z "$U90_Database" ]; then
	echo "At least one sequence database should be specified!"
	exit 1
fi

cd ${0%/*}
mainPath=$PWD

if [ ! -d "FR-t5" ]||[ ! -d "MEFTop" ]; then
	echo "Please check your FR-t5-M package!"
	exit 1
fi

cat /dev/null > softwarePath.txt
exec 3>softwarePath.txt
echo -e "Modeller\n\t$Modeller_Path">&3
if [ -n "$NR_Database" ]; then
	echo -e "nr database\n\t$NR_Database">&3
fi
if [ -n "$U90_Database" ]; then
	echo -e "u90 database\n\t$U90_Database">&3
fi
exec 3>&-

makeprogram()
{
	make clean
	make
}

summary()
{
	if [ $2 -ne 0 ]; then
		echo "Compile Error: $1">&2
		exit 1
	fi
}	

cd $mainPath/FR-t5/src
makeprogram
summary "FR-t5" $?

cd $mainPath/MEFTop/Program/nncon
sed -i 's@^\$install_dir.*$@\$install_dir = \"'"$mainPath/MEFTop/Program/nncon"'\";@' configure.pl
if [ -n "$NR_Database" ]; then
	sed -i '0,/^\$nr_db_dir/s@^\$nr_db_dir.*$@\$nr_db_dir = \"'"$NR_Database"'/\";@' configure.pl
	sed -i '0,/^\$nr_db =/s@^\$nr_db =.*$@\$nr_db = \"nr\";@' configure.pl
else
	sed -i '0,/^\$nr_db_dir/s@^\$nr_db_dir.*$@\$nr_db_dir = \"'"$U90_Database"'/\";@' configure.pl
	sed -i '0,/^\$nr_db =/s@^\$nr_db =.*$@\$nr_db = \"u90\";@' configure.pl
fi
./configure.pl
summary "nncon" $?

cd $mainPath/MEFTop/Program/DSSP
./DsspCompileGCC
summary "DSSP" $?

cd $mainPath/MEFTop/Program/Fea37
makeprogram
summary "Fea37" $?

cd $mainPath/MEFTop/Program/HC
makeprogram
summary "HCdecoy" $?

cd $mainPath/MEFTop/Program/libsvm-weights-3.11
makeprogram
summary "libsvm" $?

echo "Finish sucessfully..."
