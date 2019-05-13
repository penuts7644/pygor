#!/bin/bash

if [[ $1 = "immuno-probs" ]]; then
    shift

    if [[ -z $@ ]]
    then
        exec immuno-probs -h
    else
        exec immuno-probs $@
    fi

fi

exec $@
