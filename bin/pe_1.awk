#!/bin/awk -f
BEGIN {
	OFS = "\t"
}

($9 == "+") {
	print $1, $2, $2 + 1, $7, $8, $9
}

($9 == "-") {
	print $1, $3 - 1, $3, $7, $8, $9
}
