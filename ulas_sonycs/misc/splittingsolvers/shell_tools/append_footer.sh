#!/bin/sh
#
# File: append_footer.sh
#
# Usage
# append_footer.sh [-debug] "name-of-footer-text-file" "input_file" "output_file"
#                   usually "FastSparseSolver_footer.txt"
# 
#

if [ "$1" = "-debug" ]
then
shift
set -x
fi

FOOTER_FILE=$1

touch $3

INPUT_FILE=$2
OUTPUT_FILE=$3

TAG="% Modified `date`"

cp $INPUT_FILE $OUTPUT_FILE # copy original message part to output destination.

# See if the message was already tagged.
grep "Part of FastSparseSolver Version:100" $INPUT_FILE >/dev/null
if [ $? -ne 0 ]
then
# add a blank line
echo "" >> $OUTPUT_FILE

# append the disclaimer
cat $FOOTER_FILE >> $OUTPUT_FILE

# append modification date
echo ${TAG} >> $OUTPUT_FILE
echo "%" >> $OUTPUT_FILE
echo "" >> $OUTPUT_FILE

fi

# end script.

