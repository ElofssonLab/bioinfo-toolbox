#!/bin/bash -x

grep -w 34 long-RBS.out | grep HIT |awk '{print $2,$3,$4,$5}' > long-RBS.dat
grep AAAGGGCTAC long-RBS.dat > narK-long-RBS.dat
grep -v AAAGGGCTAC long-RBS.dat > araH-long-RBS.dat
