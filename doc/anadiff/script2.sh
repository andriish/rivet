#! /usr/bin/env bash

touch hello.md
./diffanas > hello.md
sed -i "1a ---\ntitle: 'Rivet Analysis Diffs'\n---" hello.md
