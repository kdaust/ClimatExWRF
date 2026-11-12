#!/bin/bash

isleap() {
   year=$1
   (( !(year % 4) && ( year % 100 || !(year % 400) ) )) &&
      echo "leap" || echo "nleap"
}

isleap 2001
