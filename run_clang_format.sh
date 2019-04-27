#!/bin/env bash

clang-format -i -style=Google kshortestpath/include/*h
clang-format -i -style=Google kshortestpath/src/*cpp
clang-format -i -style=Google include/*h
clang-format -i -style=Google src/*cpp
