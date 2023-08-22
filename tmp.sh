#!/bin/bash



for var in $@
do
  if [ $var == "cheese" ]
  then
	  echo "found a cheese"
  fi
done

