#!/bin/sh

SCRIPT_DIR=$(cd $(dirname $0); pwd)
TARGET_DIR="${SCRIPT_DIR}/target"

cd "${SCRIPT_DIR}"

./build.sh

cd "${TARGET_DIR}"
./pcv2

