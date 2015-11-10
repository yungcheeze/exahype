#!/bin/bash
# back up definition files
rm -rf ${PROJECT}_bak
mkdir ${PROJECT}_bak
cp ${PROJECT}/*.def ${PROJECT}_bak
cp ${PROJECT}/PdeInfo.h ${PROJECT}_bak
cp ${PROJECT}/project.peano-specification ${PROJECT}_bak
cp -R ${PROJECT}/multiscalelinkedcell ${PROJECT}_bak

# generate glue code
java -jar $PEANO_PDT_DIR/pdt.jar --generate-gluecode \
                ${PROJECT}/project.peano-specification ${PROJECT} \
                ${PEANO_PDT_DIR}/usrtemplates:${PDT_EXT}:

# update project home dir to correct directory (is set to namespace)
sed -e 's#PROJECT_HOME = '${NAMESPACE}'#PROJECT_HOME = '${PROJECT}'#' ${PROJECT}/makefile > ${PROJECT}/makefile_1 && \
mv ${PROJECT}/makefile_1 ${PROJECT}/makefile

# grep over all generated files and replace (#include "${NAMESPACE}/..." by #include "${PROJECT}/...")
grep -rl '#include "'${NAMESPACE}'' ${PROJECT} | xargs sed -i 's,#include "'${NAMESPACE}',#include "'${PROJECT}',g'

# grep over all generated files and replace (#include "${NAMESPACE}/..." by #include "${PROJECT}/...")
grep -rl '#include "multiscalelinkedcell' ${PROJECT} | xargs sed -i 's,#include "multiscalelinkedcell,#include "'${PROJECT}'/multiscalelinkedcell,g'

# update peano home dir in makefile 
SRC="${PEANO_SRC_DIR}"
sed -e 's#PEANO_HOME   = .#PEANO_HOME = '$SRC'#' ${PROJECT}/makefile > ${PROJECT}/makefile_1 && \
mv ${PROJECT}/makefile_1 ${PROJECT}/makefile

# add current dir to include path
sed -e 's#HOME) -c#HOME) -I. -c#' ${PROJECT}/makefile > ${PROJECT}/makefile_1 && \
mv ${PROJECT}/makefile_1 ${PROJECT}/makefile
