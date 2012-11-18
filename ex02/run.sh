#!/bin/sh

SCRIPT_DIR=$(cd $(dirname $0); pwd)
TARGET_DIR="${SCRIPT_DIR}/target"
RESOURCES_DIR="${SCRIPT_DIR}/src/main/resources"
EXECUTABLE="${TARGET_DIR}/pcv2"

cd "${SCRIPT_DIR}"

if [ ! -f "${EXECUTABLE}" ]; then
	./build.sh
fi

cd "${TARGET_DIR}"
"${EXECUTABLE}" \
	"${RESOURCES_DIR}/img/2.jpg" \
	"${RESOURCES_DIR}/img/1.jpg" \
	"${RESOURCES_DIR}/img/3.jpg"
