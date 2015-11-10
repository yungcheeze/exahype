#!/bin/bash
# create project
java -jar ${PEANO_PDT_DIR}/pdt.jar --create-project ${NAMESPACE} ${PROJECT} && \

# replace YourProjectName by actual project name
sed -e 's#YourProjectName#'${PROJECT}'#' ${PROJECT}/project.peano-specification > ${PROJECT}/project.peano-specification_1 && \
mv ${PROJECT}/project.peano-specification_1 ${PROJECT}/project.peano-specification

