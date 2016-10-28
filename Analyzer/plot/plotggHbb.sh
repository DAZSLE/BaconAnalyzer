#!/bin/bash

selection=$1
algo=$2
jet=$3

root -b  plotggHbb.C+\(\"${selection}\"\,\"${algo}\"\,\"${jet}\",0.38,0.3\)


