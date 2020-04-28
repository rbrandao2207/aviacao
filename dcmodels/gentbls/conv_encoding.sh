#!/bin/bash
for file in *.csv; do
	iconv -c -f ASCII -t UTF-8 "$file"
done
exit 0

## alternatively, vim *.csv followed by :argdo write ++enc=utf-8
## also needed is argdo %s/,/./ge | update
