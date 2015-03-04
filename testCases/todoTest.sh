#!/bin/bash

../ChromField BedToChromField testSizes.txt test1.cf test1.bed
../ChromField BedToChromField testSizes.txt test2.cf test2.bed

../ChromField ChromFieldToBed testSizes.txt test1.cf test1.back.bed
../ChromField ChromFieldToBed testSizes.txt test2.cf test2.back.bed

diff -s test1.bed test1.back.bed
diff -s test2.bed test2.back.bed

#operations
../ChromField Union testSizes.txt test1UnionTest2.cf test1.cf test2.cf
../ChromField ChromFieldToBed testSizes.txt test1UnionTest2.cf test1UnionTest2.back.bed

../ChromField Intersect testSizes.txt test1IntersectTest2.cf test1.cf test2.cf
../ChromField ChromFieldToBed testSizes.txt test1IntersectTest2.cf test1IntersectTest2.back.bed

../ChromField Subtract testSizes.txt test1SubtractTest2.cf test1.cf test2.cf
../ChromField ChromFieldToBed testSizes.txt test1SubtractTest2.cf test1SubtractTest2.back.bed