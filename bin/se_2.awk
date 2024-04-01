#!/bin/awk -f
BEGIN {
	OFS = "\t"
}

($5 > 0) {
	print $0
}

($6 == "+") {
	print $1, $2, $2 + 1, $4, $5, "-"
}

($6 == "-") {
	print $1, $3 - 1, $3, $4, $5, "+"
}

