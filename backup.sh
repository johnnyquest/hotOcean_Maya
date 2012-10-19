#!/bin/sh
#
#		@file		backup.sh
#		@author		xy
#		@since		2011-06-16
#
#		@brief		Create a quick backup of the whole shebang.
#

ROOT_DIR_NAME=hotOcean_Maya

BASE_DIR=..
ARCHIVE_NAME=hotOcean_Maya_`date +%Y-%m-%d_%H%M_%S`_`hostname`.tar.bz2

#rm -rf mayaJournal mayaLog


OLD_DIR=`pwd`
cd $BASE_DIR

echo in dir `pwd`, archiving to file $ARCHIVE_NAME

tar cjvf $ARCHIVE_NAME $ROOT_DIR_NAME

cd $OLD_DIR
unset OLD_DIR

