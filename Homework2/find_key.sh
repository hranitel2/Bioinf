#!/bin/bash

# Ищем файл super_secret_key.txt в текущей директории
if [ -f "./super_secret_key.txt" ]; then
    echo "Found it!" > found.log
else
    echo "I did my best"
    cat ./super_secret_key.txt 2>/dev/null
fi
