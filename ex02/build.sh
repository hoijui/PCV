#!/bin/sh

SCRIPT_DIR=$(cd $(dirname $0); pwd)
TARGET_DIR="${SCRIPT_DIR}/target"

mkdir -p "${TARGET_DIR}"

cd "${TARGET_DIR}"

if [ ! -f pcv2 ]; then
	if [ ! -f CMakeCache.txt ]; then
		cmake ..
	fi
	make
fi

