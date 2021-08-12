#!/bin/sh

: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
    echo Could not determine R_HOME.
    exit 1
fi

CXX17=`${R_HOME}/bin${R_ARCH_BIN}/R CMD config CXX17`
CXX17STD=`${R_HOME}/bin${R_ARCH_BIN}/R CMD config CXX17STD`
CXX17FLAGS=`${R_HOME}/bin${R_ARCH_BIN}/R CMD config CXX17FLAGS`

${CXX17} ${CXX17STD} ${CXX17FLAGS} -o /dev/null ./tuplemapr_test.cpp

