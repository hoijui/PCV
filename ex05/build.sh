#!/bin/sh

SCRIPT_DIR=$(cd $(dirname $0); pwd)
TARGET_DIR="${SCRIPT_DIR}/target"

mkdir -p "${TARGET_DIR}"

cd "${TARGET_DIR}"

if [ ! -f CMakeCache.txt ]; then
	cmake -DCMAKE_BUILD_TYPE=DEBUG ..
fi
make

